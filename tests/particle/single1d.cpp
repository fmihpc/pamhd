/*
Program for propagating charged particles in 1+3d with analytic fields.

Copyright 2015, 2016, 2017 Ilja Honkonen
Copyright 2025 Finnish Meteorological Institute
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

* Neither the names of the copyright holders nor the names of their contributors
  may be used to endorse or promote products derived from this software
  without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Author(s): Ilja Honkonen
*/

#include "cmath"
#include "cstdlib"
#include "fstream"
#include "ios"
#include "iostream"
#include "random"
#include "string"

#include "boost/numeric/odeint.hpp"
#include "boost/program_options.hpp"
#include "Eigen/Core"
#include "Eigen/Geometry"
#include "mpi.h"

#include "particle/common.hpp"

constexpr double
	electron_mass = 9.10938291e-31,
	electron_charge = -1.602176565e-19,
	electron_c2m = electron_charge / electron_mass,
	vacuum_permittivity = 8.854187817e-12,
	light_speed = 299792458,
	light_speed2 = light_speed*light_speed,
	max_particle_speed = 0.999999999 * light_speed;


using namespace std;
using namespace boost::numeric::odeint;


/*
Returns electric and magnetic fields at given time and position.

length == length of simulation box,
frequency == frequency of EM wave,
e1, b1 == amplitudes of EM wave perpendicular to sim box
b0 == component of magnetic field parallel to sim box
*/
struct Fields
{
	const double
		length,
		frequency,
		e1,
		b0,
		b1,
		inv_wavelength_mid,
		inv_wavelength_ends;

	Fields(
		const double given_length,
		const double given_frequency,
		const double given_e1,
		const double given_b0,
		const double given_b1,
		const double given_wavelength_mid,
		const double given_wavelength_ends
	) :
		length(given_length),
		frequency(given_frequency),
		e1(given_e1),
		b0(given_b0),
		b1(given_b1),
		inv_wavelength_mid(1.0 / given_wavelength_mid),
		inv_wavelength_ends(1.0 / given_wavelength_ends)
	{}

	double get_base_arg(const double time, const double position) const
	{
		return
			2 * M_PI * position * this->inv_wavelength_ends
			- 2 * M_PI * frequency * time
			+ (this->inv_wavelength_ends - this->inv_wavelength_mid)
				* 2 * this->length
				* cos(M_PI * position / this->length);
	}

	std::array<std::array<double, 3>, 2> operator()(
		const double time,
		double position
	) const {
		const double arg
			= [&](){
				if (position >= 0) {
					return this->get_base_arg(time, position);
				} else {
					position -= this->length;
					// make sure arg approaches same value when position -> 0
					return
						this->get_base_arg(time, position)
						+ 2 * this->length * (
							M_PI / this->inv_wavelength_ends
							+ 2 * (this->inv_wavelength_ends - this->inv_wavelength_mid)
						);
				}
			}();

		const std::array<double, 3>
			E{
				0,
				this->e1 * sin(arg),
				this->e1 * sin(arg + M_PI_2)
			},
			B{
				this->b0,
				this->b1 * sin(arg - M_PI_2),
				this->b1 * sin(arg)
			};

		return {E, B};
	}
};


// [0] == position, [1] == velocity
using state_t = std::array<Eigen::Vector3d, 2>;

struct Propagator
{
	const double c2m;
	const Fields fields;

	// c2m = charge to mass ratio of particle
	Propagator(
		const double given_c2m,
		const Fields given_fields
	) :
		c2m(given_c2m),
		fields(given_fields)
	{}

	void operator()(const state_t& state, state_t& change, const double time) const
	{
		const double relativity_factor = sqrt(1 - state[1].squaredNorm() / light_speed2);

		const auto EB = this->fields(time, state[0][0]);
		change[0] = state[1];
		const Eigen::Map<Eigen::Vector3d>
			eb0((double*)EB[0].data()),
			eb1((double*)EB[1].data());
		change[1]
			= this->c2m * relativity_factor
			* (eb0 + state[1].cross(eb1));
	}
};



int main(int argc, char* argv[])
{
	using std::acos;
	using std::cos;
	using std::minmax;
	using std::pow;
	using std::sin;
	using std::sqrt;
	using std::to_string;

	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		std::cerr << "Couldn't initialize MPI." << std::endl;
		abort();
	}

	MPI_Comm comm = MPI_COMM_WORLD;

	int rank_tmp = 0, comm_size_tmp = 0;
	if (MPI_Comm_rank(comm, &rank_tmp) != MPI_SUCCESS) {
		std::cerr << "Couldn't obtain MPI rank." << std::endl;
		abort();
	}
	if (MPI_Comm_size(comm, &comm_size_tmp) != MPI_SUCCESS) {
		std::cerr << "Couldn't obtain size of MPI communicator." << std::endl;
		abort();
	}

	if (rank_tmp < 0) {
		std::cerr << "Negative MPI rank not supported: " << rank_tmp << std::endl;
		abort();
	}
	const size_t rank = size_t(rank_tmp);
	if (comm_size_tmp < 0) {
		std::cerr << "Negative MPI communicator size not supported: " << comm_size_tmp << std::endl;
		abort();
	}
	const size_t comm_size = size_t(comm_size_tmp);

	bool verbose = true;
	std::string outfile_name("out");
	double
		length = 100,
		time_start = 0,
		time_length = 1e4,
		save_dt = 0,
		save_dx = 1,
		time_step_factor = 1,
		em_frequency = 1e4,
		density_mid = 1e8,
		density_ends = 1e8,
		b0 = 1e-6,
		b1 = 1e-10,
		v_par_min = 1e7,
		v_par_max = 1e7,
		b_angle_min = 0,
		b_angle_max = 1e-3,
		e_angle_min = 0,
		e_angle_max = 180,
		ExB_sign = 0,
		particles = 1;
	size_t seed = 1;

	boost::program_options::options_description
		options("Usage: program_name [options], where options are");

	options.add_options()
		("help", "Print this help message")
		("quiet", "Don't save additional information")
		("outfile",
			boost::program_options::value<std::string>(&outfile_name)
				->default_value(outfile_name),
			"Append results to file arg_N, where N is MPI rank)")
		("time-start",
			boost::program_options::value<double>(&time_start)
				->default_value(time_start),
			"Start time of simulation (s)")
		("time-length",
			boost::program_options::value<double>(&time_length)
				->default_value(time_length),
			"Length of simulation in units of initial electron gyro period")
		("length",
			boost::program_options::value<double>(&length)
				->default_value(length),
			"Length of simulation box in units of smallest EM wavelength")
		("frequency",
			boost::program_options::value<double>(&em_frequency)
				->default_value(em_frequency),
			"Frequency of EM waves")
		("number-density-mid",
			boost::program_options::value<double>(&density_mid)
				->default_value(density_mid),
			"Electron number density in middle of box")
		("number-density-ends",
			boost::program_options::value<double>(&density_ends)
				->default_value(density_ends),
			"Electron number density at ends of box")
		("B0",
			boost::program_options::value<double>(&b0)
				->default_value(b0),
			"Magnetic field component parallel to x coordinate")
		("B1",
			boost::program_options::value<double>(&b1)
				->default_value(b1),
			"Magnetic field magnitude perpendicular to simulation box")
		("particles",
			boost::program_options::value<double>(&particles)
				->default_value(particles),
			"Approximate total number of particles to propagate (minimum is 1 / MPI rank)")
		("v-par-min",
			boost::program_options::value<double>(&v_par_min)
				->default_value(v_par_min),
			"Minimum magnitude of velocity of electrons parallel to magnetic field (> 0, < c, m/s)")
		("v-par-max",
			boost::program_options::value<double>(&v_par_max)
				->default_value(v_par_max),
			"Maximum magnitude of velocity of electrons parallel to magnetic field (> 0, < c, m/s)")
		("B-angle-min",
			boost::program_options::value<double>(&b_angle_min)
				->default_value(b_angle_min),
			"Minimum angle of electron velocity away from magnetic field direction (degrees)")
		("B-angle-max",
			boost::program_options::value<double>(&b_angle_max)
				->default_value(b_angle_max),
			"Maximum angle of electron velocity away from magnetic field direction (degrees)")
		("E-angle-min",
			boost::program_options::value<double>(&e_angle_min)
				->default_value(e_angle_min),
			"Minimum angle of electron velocity away from electric field direction (degrees)")
		("E-angle-max",
			boost::program_options::value<double>(&e_angle_max)
				->default_value(e_angle_max),
			"Maximum angle of electron velocity away from electric field direction (degrees)")
		("ExB-sign",
			boost::program_options::value<double>(&ExB_sign)
				->default_value(ExB_sign),
			"Sign of electron velocity component parallel to ExB (0 == all or < 0 or > 0)")
		("save-dt",
			boost::program_options::value<double>(&save_dt)
				->default_value(save_dt),
			"Save results at least every arg seconds, 0 saves "
			"initial and final states, < 0 doesn't save")
		("save-dx",
			boost::program_options::value<double>(&save_dx)
				->default_value(save_dx),
			"Save results at least every arg lengths in units of minimum EM wavelength")
		("time-step-factor",
			boost::program_options::value<double>(&time_step_factor)
				->default_value(time_step_factor),
			"Multiply maximum allowed time step with factor arg")
		("seed",
			boost::program_options::value<size_t>(&seed)
				->default_value(seed),
			"Use arg (+ MPI rank) as seed for random number generator");

	boost::program_options::variables_map option_variables;
	try {
		boost::program_options::store(
			boost::program_options::command_line_parser(argc, argv)
				.options(options)
				.run(),
			option_variables
		);
	} catch (std::exception& e) {
		if (rank == 0) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
				<< "Couldn't parse command line options: " << e.what()
				<< std::endl;
		}
		MPI_Finalize();
		return EXIT_FAILURE;
	}
	boost::program_options::notify(option_variables);

	if (option_variables.count("help") > 0) {
		if (rank == 0) {
			cout << options << endl;
		}
		MPI_Finalize();
		return EXIT_SUCCESS;
	}

	if (option_variables.count("quiet") > 0) {
		verbose = false;
	}

	if (b_angle_min < 0) {
		if (rank == 0) {
			std::cerr << "Minimum angle from magnetic field must be >= 0." << std::endl;
		}
		MPI_Finalize();
		return EXIT_FAILURE;
	}
	if (b_angle_max > 180) {
		if (rank == 0) {
			std::cerr << "Maximum angle from magnetic field must be <= 180." << std::endl;
		}
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	if (e_angle_min < 0) {
		if (rank == 0) {
			std::cerr << "Minimum angle from electric field must be >= 0." << std::endl;
		}
		MPI_Finalize();
		return EXIT_FAILURE;
	}
	if (e_angle_max > 180) {
		if (rank == 0) {
			std::cerr << "Maximum angle from electric field must be <= 180." << std::endl;
		}
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	std::tie(v_par_min, v_par_max) = minmax(abs(v_par_min), abs(v_par_max));
	b_angle_min *= M_PI / 180;
	b_angle_max *= M_PI / 180;
	e_angle_min *= M_PI / 180;
	e_angle_max *= M_PI / 180;

	if (e_angle_max + b_angle_max < M_PI / 2) {
		if (rank == 0) {
			std::cerr << "Sum of maximum angles away from B and E must be >= 90 degrees." << std::endl;
		}
		MPI_Finalize();
		return EXIT_FAILURE;
	}
	if (e_angle_min + b_angle_min < M_PI / 2) {
		if (rank == 0) {
			std::cout << "NOTE: Sum of minimum angles away from B and E < 90 degrees might slow things down significantly." << std::endl;
		}
	}

	if (v_par_min / abs(cos(b_angle_min)) >= light_speed) {
		if (rank == 0) {
			std::cerr << "Minimum electron velocity >= light speed." << std::endl;
		}
		MPI_Finalize();
		return EXIT_FAILURE;
	}
	if (v_par_max / abs(cos(b_angle_max)) >= light_speed) {
		if (rank == 0) {
			std::cout << "NOTE: Maximum electron velocity >= light speed might slow things down significantly." << std::endl;
		}
	}

	particles = std::round(particles);
	seed += rank;
	if (particles > comm_size) {
		particles = round(particles / comm_size);
	} else {
		particles = 1;
	}
	outfile_name += "_" + to_string(rank);

	std::ofstream outfile(outfile_name, std::ios::app | std::ios::ate);
	if (not outfile.good()) {
		cerr << "Couldn't open file " << outfile_name << " for writing." << endl;
		return EXIT_FAILURE;
	}
	outfile
		<< std::scientific
		<< std::setprecision(std::numeric_limits<double>::max_digits10);

	// derived parameters
	const double
		e_gyro_freq
			= abs(electron_charge)
			* sqrt(b0*b0 + 2*b1*b1)
			/ (2 * M_PI * electron_mass),
		e_plasma_freq_mid
			= sqrt(
				density_mid
				* pow(electron_charge, 2)
				/ (4 * M_PI * M_PI * electron_mass * vacuum_permittivity)
			),
		e_plasma_freq_ends
			= sqrt(
				density_ends
				* pow(electron_charge, 2)
				/ (4 * M_PI * M_PI * electron_mass * vacuum_permittivity)
			),
		wavelength_mid
			= light_speed
			/ em_frequency
			/ sqrt(
				1
				+ pow(e_plasma_freq_mid, 2)
					/ em_frequency
					/ (e_gyro_freq - em_frequency)
			),
		wavelength_ends
			= light_speed
			/ em_frequency
			/ sqrt(
				1
				+ pow(e_plasma_freq_ends, 2)
					/ em_frequency
					/ (e_gyro_freq - em_frequency)
			),
		e1 = light_speed
			* b1
			/ sqrt(
				1
				- pow(e_plasma_freq_mid / em_frequency, 2)
					/ (1 - e_gyro_freq / em_frequency)
			);

	length *= min(wavelength_mid, wavelength_ends);
	save_dx *= min(wavelength_mid, wavelength_ends);
	time_length /= e_gyro_freq;

	const double min_wavelength = std::min(wavelength_mid, wavelength_ends);

	if (verbose) {
		outfile
			<< "# Parameters\n# Start time (s): " << time_start
			<< "\n# Simulation length: " << time_length << " s, "
			<< time_length * e_gyro_freq << " electron gyro periods (g_e)"
			<< "\n# Simulation box length (m): " << length
			<< "\n# Number of particles: " << particles
			<< "\n# Minimum electron velocity parallel to magnetic field (m/s): " << v_par_min
			<< "\n# Maximum electron velocity parallel to magnetic field (m/s): " << v_par_max
			<< "\n# Minimum electron velocity angle away from magnetic field (degrees): " << b_angle_min * 180 / M_PI
			<< "\n# Maximum electron velocity angle away from magnetic field (degrees): " << b_angle_max * 180 / M_PI
			<< "\n# Minimum electron velocity angle away from electric field (degrees): " << e_angle_min * 180 / M_PI
			<< "\n# Maximum electron velocity angle away from electric field (degrees): " << e_angle_max * 180 / M_PI
			<< "\n# EM wave frequency (f_w): " << em_frequency
			<< " Hz\n# Electron number density in middle (1/m^3): " << density_mid
			<< "\n# Electron number density at ends (1/m^3): " << density_ends
			<< "\n# B0 (T): " << b0
			<< "\n# B1 (T): " << b1
			<< "\n# Random seed: " << seed
			<< "\n# Derived parameters\n# Electron gyro frequency (f_g): "
			<< e_gyro_freq << " Hz, " << e_gyro_freq / em_frequency << " f_w"
			<< "\n# Electron plasma frequency in middle: " << e_plasma_freq_mid
			<< " Hz, " << e_plasma_freq_mid / e_gyro_freq << " f_g"
			<< ", " << e_plasma_freq_mid / em_frequency << " f_w"
			<< "\n# Electron plasma frequency at ends: " << e_plasma_freq_ends
			<< " Hz, " << e_plasma_freq_ends / e_gyro_freq << " f_g"
			<< ", " << e_plasma_freq_ends / em_frequency << " f_w"
			<< "\n# EM wavelength in middle: " << wavelength_mid << " m, "
			<< wavelength_mid / length << " simulation boxes"
			<< "\n# EM wavelength at ends: " << wavelength_ends << " m, "
			<< wavelength_ends / length << " simulation boxes"
			<< "\n# EM phase speed in middle: " << wavelength_mid * em_frequency
			<< " m/s, " << wavelength_mid * em_frequency / light_speed << " c"
			<< "\n# EM phase speed at ends: " << wavelength_ends * em_frequency
			<< " m/s, " << wavelength_ends * em_frequency / light_speed << " c"
			<< "\n# E1 (V/m): " << e1
			<< "\n# Each particle on separate line, format: "
			   "r1x, v1x, v1y, v1z, Ey(r1), Ez(r1), Bx(r1), By(r1), Bz(r1), r2x..."
			<< std::endl;
	}

	std::mt19937_64 random_source;
	random_source.seed(seed);
	std::uniform_real_distribution<>
		v_par_gen(v_par_min, v_par_max),
		b_angle_gen(b_angle_min, b_angle_max),
		e_angle_gen(e_angle_min, e_angle_max),
		v_sign_gen(0, 1);

	const Fields fields(
		length,
		em_frequency,
		e1,
		b0,
		b1,
		wavelength_mid,
		wavelength_ends
	);

	const std::array<double, 3>
		initial_pos{0, 0, 0},
		init_E = fields(0, initial_pos[0])[0],
		init_B = fields(0, initial_pos[0])[1];

	Propagator propagator(electron_c2m, fields);
	runge_kutta_fehlberg78<state_t> stepper;

	size_t particle = 0;
	while (particle < particles) {

		const auto
			e_angle = e_angle_gen(random_source),
			b_angle = b_angle_gen(random_source);
		if (b_angle + e_angle < M_PI / 2) {
			continue;
		}

		const auto
			v_par = v_par_gen(random_source),
			v = v_par / abs(cos(b_angle));
		if (v >= max_particle_speed) {
			continue;
		}

		const std::array<double, 3> v_in_E_B_ExB{
			v * cos(e_angle),

			v * cos(b_angle),

			[&]() -> double {
				const auto arg_to_sqrt
					= pow(sin(e_angle) * sin(b_angle), 2)
					- pow(cos(e_angle) * cos(b_angle), 2);
				if (arg_to_sqrt < 0) {
					return 0;
				} else {
					const auto ret_val = v * sqrt(arg_to_sqrt);
					// select which side of EB plane velocity will be on
					if (ExB_sign < 0) {
						return -ret_val;
					}
					if (ExB_sign > 0) {
						return +ret_val;
					}

					auto v_sign = v_sign_gen(random_source);
					if (v_sign < 0.5) {
						return -ret_val;
					} else {
						return +ret_val;
					}
				}
			}()
		};
		if (pamhd::norm(v_in_E_B_ExB) >= max_particle_speed) {
			continue;
		}

		/*
		Rotate initial velocity to x, y, z coordinate system
		*/

		const std::array<double, 3>
			old_x = pamhd::mul(init_E, 1 / pamhd::norm(init_E)),
			old_y = pamhd::mul(init_B, 1 / pamhd::norm(init_B)),
			init_ExB = pamhd::cross(init_E, init_B),
			old_z = pamhd::mul(init_ExB, 1 / pamhd::norm(init_ExB)),
			new_x{1, 0, 0},
			new_y{0, 1, 0},
			new_z{0, 0, 1};

		// http://ocw.mit.edu/courses/aeronautics-and-astronautics/16-07-dynamics-fall-2009/lecture-notes/MIT16_07F09_Lec03.pdf
		Eigen::Matrix3d transformer;
		transformer <<
			pamhd::dot(new_x, old_x), pamhd::dot(new_x, old_y), pamhd::dot(new_x, old_z),
			pamhd::dot(new_y, old_x), pamhd::dot(new_y, old_y), pamhd::dot(new_y, old_z),
			pamhd::dot(new_z, old_x), pamhd::dot(new_z, old_y), pamhd::dot(new_z, old_z);

		const Eigen::Map<Eigen::Vector3d> v_map((double*)v_in_E_B_ExB.data());
		const Eigen::Vector3d temp = transformer * v_map;
		const std::array<double, 3> initial_vel{temp[0], temp[1], temp[2]};
		if (pamhd::norm(initial_vel) >= max_particle_speed) {
			continue;
		}

		const auto angle_V_E = acos(pamhd::dot(pamhd::norm(initial_vel), pamhd::norm(init_E)));
		if (angle_V_E < e_angle_min or angle_V_E > e_angle_max)  {
			continue;
		}
		const auto angle_V_B = acos(pamhd::dot(pamhd::norm(initial_vel), pamhd::norm(init_B)));
		if (angle_V_B < b_angle_min or angle_V_B > b_angle_max)  {
			continue;
		}

		Eigen::Vector3d temp1, temp2;
		temp1[0] = initial_pos[0];
		temp1[1] = initial_pos[1];
		temp1[2] = initial_pos[2];
		temp2[0] = initial_vel[0];
		temp2[1] = initial_vel[1];
		temp2[2] = initial_vel[2];
		state_t state{temp1, temp2};

		std::vector<double> data_to_save;
		if (save_dt >= 0 or save_dx >= 0) {
			data_to_save.push_back(state[0][0]);
			data_to_save.push_back(state[1][0]);
			data_to_save.push_back(state[1][1]);
			data_to_save.push_back(state[1][2]);
			data_to_save.push_back(init_E[1]);
			data_to_save.push_back(init_E[2]);
			data_to_save.push_back(init_B[0]);
			data_to_save.push_back(init_B[1]);
			data_to_save.push_back(init_B[2]);
		}

		double
			simulation_time = 0,
			dt = 0,
			next_save_t = time_start;

		if (save_dt < 0) {
			next_save_t += 2 * time_length;
		}

		if (save_dx < 0) {
			save_dx = 2 * length;
		} else if (save_dx == 0) {
			save_dx = length;
		}

		double previous_save_x = state[0][0];
		while (true) {

			const auto EB = fields(simulation_time, state[0][0]);
			if (
				simulation_time >= next_save_t
				or abs(previous_save_x - state[0][0]) >= save_dx
			) {
				data_to_save.push_back(state[0][0]);
				data_to_save.push_back(state[1][0]);
				data_to_save.push_back(state[1][1]);
				data_to_save.push_back(state[1][2]);
				data_to_save.push_back(EB[0][1]);
				data_to_save.push_back(EB[0][2]);
				data_to_save.push_back(EB[1][0]);
				data_to_save.push_back(EB[1][1]);
				data_to_save.push_back(EB[1][2]);

				if (save_dt == 0) {
					next_save_t += time_length;
				} else {
					next_save_t += save_dt;
				}
				previous_save_x = state[0][0];
			}

			double max_time_step = std::numeric_limits<double>::max();
			const std::array<double, 3> state1{
				state[1][0], state[1][1], state[1][2]};
			const auto step_size
				= pamhd::particle::get_step_size(
					1.0 / 32.0, min_wavelength / 8,
					electron_c2m, state1, EB[0], EB[1]
				);
			max_time_step = std::min(
				max_time_step,
				std::min(step_size.first, step_size.second)
			);
			max_time_step = std::min(max_time_step, 1 / em_frequency / 8);
			dt = max_time_step * time_step_factor;

			if (not std::isfinite(state[0][0])) {
				std::cerr << "before: "
					<< state[0][0] << ", "
					<< state[0][1] << ", "
					<< state[0][2] << ";  "
					<< state[1][0] << ", "
					<< state[1][1] << ", "
					<< state[1][2]
					<< std::endl;
				abort();
			}

			stepper.do_step(propagator, state, simulation_time, dt);
			if (state[0][0] < -1.001 * length or state[0][0] > 1.001 * length) {
				break;
			}

			if (not std::isfinite(state[0][0])) {
				std::cerr << "after: "
					<< state[0][0] << ", "
					<< state[0][1] << ", "
					<< state[0][2] << ";  "
					<< state[1][0] << ", "
					<< state[1][1] << ", "
					<< state[1][2]
					<< std::endl;
				abort();
			}
			simulation_time += dt;

			if (simulation_time > time_start + time_length) {
				std::cerr
					<< "Simulation time exceeded for particle with initial velocity "
					<< initial_vel[0] << " "
					<< initial_vel[1] << " "
					<< initial_vel[2] << ", current velocity relative to initial "
					<< state[1].norm() / pamhd::norm(initial_vel) << " and current position "
					<< state[0][0] << " "
					<< state[0][1] << " "
					<< state[0][2]
					<< std::endl;
				break;
			}
		}

		const auto EB = fields(simulation_time, state[0][0]);
		if (save_dt >= 0 or save_dx >= 0) {
			data_to_save.push_back(state[0][0]);
			data_to_save.push_back(state[1][0]);
			data_to_save.push_back(state[1][1]);
			data_to_save.push_back(state[1][2]);
			data_to_save.push_back(EB[0][1]);
			data_to_save.push_back(EB[0][2]);
			data_to_save.push_back(EB[1][0]);
			data_to_save.push_back(EB[1][1]);
			data_to_save.push_back(EB[1][2]);
		}

		for (const auto& item: data_to_save) {
			outfile << item << " ";
		}
		outfile << "\n";

		particle++;
	}
	if (not outfile.good()) {
		cerr << "Writing results to " << outfile_name << " probably failed." << endl;
	}
	outfile.close();

	MPI_Finalize();

	return EXIT_SUCCESS;
}
