/*
Reference test program for MHD solvers of PAMHD.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/


#include "array"
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "cmath"
#include "cstdlib"
#include "iostream"
#include "string"
#include "fstream"
#include "iomanip"
#include "iostream"
#include "limits"
#include "sstream"
#include "string"
#include "tuple"
#include "type_traits"

#include "gensimcell.hpp"

#include "mhd/common.hpp"
#include "mhd/hll_athena.hpp"
#include "mhd/hlld_athena.hpp"
#include "mhd/roe_athena.hpp"
#include "mhd/rusanov.hpp"
#include "mhd/variables.hpp"
#include "variables.hpp"


using namespace std;


constexpr double
	adiabatic_index = 5.0 / 3.0,
	vacuum_permeability = 4e-7 * M_PI,
	proton_mass = 1.672621777e-27;


// 1 boundary cell in each direction
constexpr size_t grid_length = 1000 + 2;

using Cell = gensimcell::Cell<
	gensimcell::Never_Transfer,
	pamhd::mhd::HD_State_Conservative,
	// used as face magnetic field since system is 1d with const Bx
	pamhd::Magnetic_Field,
	pamhd::Edge_Electric_Field,
	pamhd::mhd::HD_Flux_Conservative,
	pamhd::Magnetic_Field_Flux, // not real flux, used as dB
	pamhd::Bg_Magnetic_Field_Pos_X
>;
using Grid = std::array<Cell, grid_length>;

/*
Accessors to simulation variables.

Return a reference to data of corresponding variable in given simulation cell.
*/

// conservative MHD variables
const auto Mas
	= [](Cell& cell_data)->typename pamhd::mhd::Mass_Density::data_type&{
		return cell_data[pamhd::mhd::HD_State_Conservative()][pamhd::mhd::Mass_Density()];
	};

const auto Mom
	= [](Cell& cell_data)->typename pamhd::mhd::Momentum_Density::data_type&{
		return cell_data[pamhd::mhd::HD_State_Conservative()][pamhd::mhd::Momentum_Density()];
	};

const auto Nrj
	= [](Cell& cell_data)->typename pamhd::mhd::Total_Energy_Density::data_type&{
		return cell_data[pamhd::mhd::HD_State_Conservative()][pamhd::mhd::Total_Energy_Density()];
	};

const auto Mag
	= [](Cell& cell_data)->typename pamhd::Magnetic_Field::data_type&{
		return cell_data[pamhd::Magnetic_Field()];
	};

const auto Ele
	= [](Cell& cell_data)->typename pamhd::Edge_Electric_Field::data_type&{
		return cell_data[pamhd::Edge_Electric_Field()];
	};

const auto Bg_Mag
	= [](Cell& cell_data)->typename pamhd::Bg_Magnetic_Field_Pos_X::data_type&{
		return cell_data[pamhd::Bg_Magnetic_Field_Pos_X()];
	};

// fluxes of conservative MHD variables
const auto Mas_f
	= [](Cell& cell_data)->typename pamhd::mhd::Mass_Density::data_type&{
		return cell_data[pamhd::mhd::HD_Flux_Conservative()][pamhd::mhd::Mass_Density()];
	};

const auto Mom_f
	= [](Cell& cell_data)->typename pamhd::mhd::Momentum_Density::data_type&{
		return cell_data[pamhd::mhd::HD_Flux_Conservative()][pamhd::mhd::Momentum_Density()];
	};

const auto Nrj_f
	= [](Cell& cell_data)->typename pamhd::mhd::Total_Energy_Density::data_type&{
		return cell_data[pamhd::mhd::HD_Flux_Conservative()][pamhd::mhd::Total_Energy_Density()];
	};

const auto Mag_f
	= [](Cell& cell_data)->typename pamhd::Magnetic_Field_Flux::data_type&{
		return cell_data[pamhd::Magnetic_Field_Flux()];
	};


/*!
Returns start coordinate of the grid.
*/
constexpr double get_grid_start()
{
	return 0.0;
}

/*!
Returns end coordinate of the grid.
*/
constexpr double get_grid_end()
{
	return 1e5;
}


/*!
Returns the size of a grid cell.
*/
constexpr double get_cell_size()
{
	static_assert(
		std::tuple_size<Grid>::value > 0,
		"Grid must have at least one cell"
	);
	return (get_grid_end() - get_grid_start()) / std::tuple_size<Grid>::value;
}


/*!
Returns the center of a cell located at given index in given grid.

Index starts from 0.
Returns a quiet NaN in case of error.
*/
double get_cell_center(const size_t index)
{
	static_assert(
		std::tuple_size<Grid>::value > 0,
		"Grid must have at least one cell"
	);

	if (index >= std::tuple_size<Grid>::value) {
		return std::numeric_limits<double>::quiet_NaN();
	}

	return (0.5 + index) * get_cell_size();
}


//! Sets the initial state of simulation, zeroes fluxes
template <
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Magnetic_Field_Getter,
	class Edge_Electric_Field_Getter,
	class Mass_Density_Flux_Getter,
	class Momentum_Density_Flux_Getter,
	class Total_Energy_Density_Flux_Getter,
	class Magnetic_Field_Flux_Getter
> void initialize_mhd(
	const std::string& /*solver*/,
	Grid& grid,
	const double adiabatic_index,
	const double vacuum_permeability,
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Magnetic_Field_Getter Mag,
	const Edge_Electric_Field_Getter Ele,
	const Mass_Density_Flux_Getter Mas_f,
	const Momentum_Density_Flux_Getter Mom_f,
	const Total_Energy_Density_Flux_Getter Nrj_f,
	const Magnetic_Field_Flux_Getter Mag_f
) {
	static_assert(
		std::tuple_size<Grid>::value > 0,
		"Grid must have at least one cell"
	);

	for (size_t cell_i = 0; cell_i < grid.size(); cell_i++) {
		auto& cell_data = grid[cell_i];

		Mas_f(cell_data)    =
		Nrj_f(cell_data)    =
		Mom_f(cell_data)[0] =
		Mom_f(cell_data)[1] =
		Mom_f(cell_data)[2] =
		Mag_f(cell_data)[0] =
		Mag_f(cell_data)[1] =
		Mag_f(cell_data)[2] = 0;

		Mom(cell_data)[0]    =
		Mom(cell_data)[1]    =
		Mom(cell_data)[2]    =
		Mag(cell_data)[0]    =
		Mag(cell_data)[1]    =
		Mag(cell_data)[2]    =
		Ele(cell_data)[0]    =
		Ele(cell_data)[1]    =
		Ele(cell_data)[2]    =
		Bg_Mag(cell_data)[0] =
		Bg_Mag(cell_data)[1] =
		Bg_Mag(cell_data)[2] = 0;

		const double center = get_cell_center(cell_i);
		double pressure = -1;
		Mag(cell_data)[0] = 1.5e-9;
		if (center < (get_grid_end() - get_grid_start()) / 2) {
			Mas(cell_data) = proton_mass * 3e6;
			Mag(cell_data)[1] = 1e-9;
			pressure = 3e-12;
		} else {
			Mas(cell_data) = proton_mass * 1e6;
			Mag(cell_data)[1] = -1e-9;
			pressure = 1e-12;
		}
		// face in middle of grid has average B
		if (grid.size() % 2 == 0 and cell_i == grid.size() / 2 - 1) {
			Mag(cell_data)[1] = 0;
		}
		// cell in middle of grid has average mass, pressure
		if (grid.size() % 2 == 1 and cell_i == grid.size() / 2) {
			Mas(cell_data) = proton_mass * 2e6;
			pressure = 2e-12;
		}

		Nrj(cell_data)
			= pamhd::mhd::get_total_energy_density(
				Mas(cell_data),
				Mom(cell_data),
				pressure,
				Mag(cell_data),
				adiabatic_index,
				vacuum_permeability
			);
	}
}

/*!
Advances the MHD solution for one time step of length dt with given solver.

Returns the maximum allowed length of time step for the
next step.
*/
template <
	class Mass_Density_Getter,
	class Momentum_Density_Getter,
	class Total_Energy_Density_Getter,
	class Magnetic_Field_Getter,
	class Electric_Field_Getter,
	class Mass_Density_Flux_Getter,
	class Momentum_Density_Flux_Getter,
	class Total_Energy_Density_Flux_Getter,
	class Magnetic_Field_Flux_Getter
> double solve_mhd(
	const pamhd::mhd::Solver solver,
	Grid& grid,
	const double dt,
	const double adiabatic_index,
	const double vacuum_permeability,
	const Mass_Density_Getter Mas,
	const Momentum_Density_Getter Mom,
	const Total_Energy_Density_Getter Nrj,
	const Magnetic_Field_Getter Mag,
	const Electric_Field_Getter Ele,
	const Mass_Density_Flux_Getter Mas_f,
	const Momentum_Density_Flux_Getter Mom_f,
	const Total_Energy_Density_Flux_Getter Nrj_f,
	const Magnetic_Field_Flux_Getter Mag_f
) {
	static_assert(
		std::tuple_size<Grid>::value >= 3,
		"Grid must have at least three cells"
	);

	// shorthand notation for variables used internally
	const pamhd::mhd::Mass_Density mas{};
	const pamhd::mhd::Momentum_Density mom{};
	const pamhd::mhd::Total_Energy_Density nrj{};
	const pamhd::Magnetic_Field mag{};

	const double
		cell_size = get_cell_size(),
		face_area = 1.0;

	double max_dt = std::numeric_limits<double>::max();

	// internal data type used by MHD solvers
	using MHD = gensimcell::Cell<
		gensimcell::Never_Transfer,
		pamhd::mhd::Mass_Density,
		pamhd::mhd::Momentum_Density,
		pamhd::mhd::Total_Energy_Density,
		pamhd::Magnetic_Field
	>;

	// calculate fluxes
	for (size_t cell_i = 0; cell_i < std::tuple_size<Grid>::value - 1; cell_i++) {

		Cell
			&cell = grid[cell_i],
			&neighbor = grid[cell_i + 1];

		double max_vel = -1;
		MHD state_neg, state_pos, flux;
		state_neg[mas] = Mas(cell);
		state_neg[mom] = Mom(cell);
		state_neg[nrj] = Nrj(cell);
		state_neg[mag] = Mag(cell);
		state_pos[mas] = Mas(neighbor);
		state_pos[mom] = Mom(neighbor);
		state_pos[nrj] = Nrj(neighbor);
		state_pos[mag] = Mag(neighbor);
		const auto bg_magnetic_field = Bg_Mag(cell);

		#define SOLVER(name) \
			name< \
				pamhd::mhd::Mass_Density, \
				pamhd::mhd::Momentum_Density, \
				pamhd::mhd::Total_Energy_Density, \
				pamhd::Magnetic_Field \
			>( \
				state_neg, \
				state_pos, \
				bg_magnetic_field, \
				adiabatic_index, \
				vacuum_permeability \
			)
		switch (solver) {
		case pamhd::mhd::Solver::rusanov:
			std::tie(flux, max_vel) = SOLVER(pamhd::mhd::get_flux_rusanov);
			break;
		case pamhd::mhd::Solver::hll_athena:
			std::tie(flux, max_vel) = SOLVER(pamhd::mhd::athena::get_flux_hll);
			break;
		case pamhd::mhd::Solver::hlld_athena:
			std::tie(flux, max_vel) = SOLVER(pamhd::mhd::athena::get_flux_hlld);
			break;
		case pamhd::mhd::Solver::roe_athena:
			std::tie(flux, max_vel) = SOLVER(pamhd::mhd::athena::get_flux_roe);
			break;
		default:
			abort();
		}
		#undef SOLVER
		flux[mas] *= face_area*dt;
		flux[mom] *= face_area*dt;
		flux[nrj] *= face_area*dt;

		max_dt = std::min(max_dt, cell_size / max_vel);

		if (cell_i > 0) {
			// positive flux flows neg->pos, i.e. out of current cell
			Mas_f(cell) -= flux[mas];
			Mom_f(cell) -= flux[mom];
			Nrj_f(cell) -= flux[nrj];
		}
		if (cell_i < std::tuple_size<Grid>::value - 2) {
			Mas_f(neighbor) += flux[mas];
			Mom_f(neighbor) += flux[mom];
			Nrj_f(neighbor) += flux[nrj];
		}

		// upwind E = -VxB
		const auto
			v_c = (Mom(cell) / Mas(cell)).eval(),
			v_n = (Mom(neighbor) / Mas(neighbor)).eval();
		Ele(cell)[0] = 0;
		if (v_c[0] + v_n[0] > 0) {
			Ele(cell)[1] = v_c[0]*Mag(cell)[2] - v_c[2]*Mag(cell)[0];
			Ele(cell)[2] = v_c[1]*Mag(cell)[0] - v_c[0]*Mag(cell)[1];
		} else {
			Ele(cell)[1] = v_n[0]*Mag(neighbor)[2] - v_n[2]*Mag(cell)[0];
			Ele(cell)[2] = v_n[1]*Mag(cell)[0] - v_n[0]*Mag(neighbor)[1];
		}

		Mag_f(cell)[0] = 0;
		Mag_f(cell)[1] = Ele(cell)[2];
		Mag_f(cell)[2] = -Ele(cell)[1];
		if (cell_i > 0) {
			Mag_f(cell)[1] -= Ele(grid[cell_i-1])[2];
			Mag_f(cell)[2] += Ele(grid[cell_i-1])[1];
		}
		Mag_f(cell) *= dt;
	}

	// apply fluxes
	const double inverse_volume = 1.0 / (face_area * cell_size);

	for (auto& cell: grid) {
		pamhd::mhd::apply_fluxes(
			cell,
			inverse_volume,
			adiabatic_index,
			vacuum_permeability,
			Mas, Mom, Nrj, Mag,
			Mas_f, Mom_f, Nrj_f, Mag_f
		);

		Mas_f(cell)    =
		Mom_f(cell)[0] =
		Mom_f(cell)[1] =
		Mom_f(cell)[2] =
		Nrj_f(cell)    =
		Mag_f(cell)[0] =
		Mag_f(cell)[1] =
		Mag_f(cell)[2] = 0;
	}

	return max_dt;
}

/*!
Writes given advection simulation to a file plottable with gnuplot.

Returns the name of the file written which
is derived from given simulation time.
*/
std::string plot_mhd(
	Grid& grid,
	const double simulation_time,
	const double adiabatic_index,
	const double vacuum_permeability
) {
	std::ostringstream time_string, normalization_string;
	time_string
		<< std::setw(4)
		<< std::setfill('0')
		<< size_t(simulation_time * 1000);

	const std::string
		gnuplot_file_name(
			"mhd_s_"
			+ time_string.str()
			+ "_ms.dat"
		),
		plot_file_name(
			"mhd_s_"
			+ time_string.str()
			+ "_ms.png"
		),
		B_plot_file_name(
			"mhd_s_B_"
			+ time_string.str()
			+ "_ms.png"
		),
		v_plot_file_name(
			"mhd_s_v_"
			+ time_string.str()
			+ "_ms.png"
		);

	std::ofstream gnuplot_file(gnuplot_file_name);

	gnuplot_file
		<< "set term png enhanced\nset output '"
		<< plot_file_name
		<< "'\nset xlabel 'position'\n"
		   "set format y '%1.2e'\n"
		   "set format y2 '%1.2e'\n"
		   "set ylabel 'density/pressure'\n"
		   "set ytics nomirror\n"
		   "set y2tics nomirror\n"
		   "set key horizontal outside bottom\n"
		   "plot "
		     "'-' using 1:2 with line linewidth 2 title 'density', "
		     "'-' u 1:2 axis x1y2 w l lw 2 t 'pressure'\n"
		     ;

	// mass density
	for (size_t cell_i = 0; cell_i < std::tuple_size<Grid>::value; cell_i++) {
		auto& cell_data = grid[cell_i];
		const double x
			= get_cell_center(cell_i)
			/ (get_grid_end() - get_grid_start());
		gnuplot_file << x << " " << Mas(cell_data) << "\n";
	}
	gnuplot_file << "end\n";
	// pressure
	for (size_t cell_i = 0; cell_i < std::tuple_size<Grid>::value; cell_i++) {
		auto& cell_data = grid[cell_i];
		const double x = get_cell_center(cell_i) / (get_grid_end() - get_grid_start());
		gnuplot_file
			<< x << " "
			<< pamhd::mhd::get_pressure(
				Mas(cell_data),
				Mom(cell_data),
				Nrj(cell_data),
				Mag(cell_data),
				adiabatic_index,
				vacuum_permeability
			)
			<< "\n";
	}
	gnuplot_file << "end\n";

	gnuplot_file
		<< "set term png enhanced\nset output '"
		<< v_plot_file_name
		<< "'\nset xlabel 'position'\n"
		   "set ylabel 'velocity'\n"
		   "unset y2tics\n"
		   "set key horizontal outside bottom\n"
		   "plot "
		     "'-' u 1:2 w l lw 2 t 'v_x', "
		     "'-' u 1:2 w l lw 2 t 'v_y', "
		     "'-' u 1:2 w l lw 2 t 'v_z'\n"
		;
	// vx
	for (size_t cell_i = 0; cell_i < std::tuple_size<Grid>::value; cell_i++) {
		const double
			x = get_cell_center(cell_i) / (get_grid_end() - get_grid_start()),
			vx = Mom(grid[cell_i])[0] / Mas(grid[cell_i]);
		gnuplot_file << x << " " << vx << "\n";
	}
	gnuplot_file << "end\n";
	// vy
	for (size_t cell_i = 0; cell_i < std::tuple_size<Grid>::value; cell_i++) {
		const double
			x = get_cell_center(cell_i) / (get_grid_end() - get_grid_start()),
			vy = Mom(grid[cell_i])[1] / Mas(grid[cell_i]);
		gnuplot_file << x << " " << vy << "\n";
	}
	gnuplot_file << "end\n";
	// vz
	for (size_t cell_i = 0; cell_i < std::tuple_size<Grid>::value; cell_i++) {
		const double
			x = get_cell_center(cell_i) / (get_grid_end() - get_grid_start()),
			vz = Mom(grid[cell_i])[2] / Mas(grid[cell_i]);
		gnuplot_file << x << " " << vz << "\n";
	}
	gnuplot_file << "end\n";

	gnuplot_file
		<< "set term png enhanced\nset output '"
		<< B_plot_file_name
		<< "'\nset xlabel 'position'\n"
		   "set ylabel 'magnetic field'\n"
		   "set key horizontal outside bottom\n"
		   "plot "
		     "'-' u 1:2 w l lw 2 t 'B_x', "
		     "'-' u 1:2 w l lw 2 t 'B_y', "
		     "'-' u 1:2 w l lw 2 t 'B_z'\n"
		     ;
	// Bx
	for (size_t cell_i = 0; cell_i < std::tuple_size<Grid>::value; cell_i++) {
		const double x
			= get_cell_center(cell_i)
			/ (get_grid_end() - get_grid_start());
		gnuplot_file << x << " " << Mag(grid[cell_i])[0] + Bg_Mag(grid[cell_i])[0] << "\n";
	}
	gnuplot_file << "end\n";
	// By
	for (size_t cell_i = 0; cell_i < std::tuple_size<Grid>::value; cell_i++) {
		const double x
			= get_cell_center(cell_i)
			/ (get_grid_end() - get_grid_start());
		gnuplot_file << x << " " << Mag(grid[cell_i])[1] + Bg_Mag(grid[cell_i])[1] << "\n";
	}
	gnuplot_file << "end\n";
	// Bz
	for (size_t cell_i = 0; cell_i < std::tuple_size<Grid>::value; cell_i++) {
		const double x
			= get_cell_center(cell_i)
			/ (get_grid_end() - get_grid_start());
		gnuplot_file << x << " " << Mag(grid[cell_i])[2] + Bg_Mag(grid[cell_i])[2] << "\n";
	}
	gnuplot_file << "end\n";

	return gnuplot_file_name;
}

/*!
Writes given advection simulation to an ascii file.

Returns the name of the file written which
is derived from given solver name.
*/
void save_mhd(
	const std::string& solver,
	Grid& grid
) {
	const std::string file_name("mhd_staggered_" + solver + ".dat");

	std::ofstream outfile(file_name);
	outfile << std::setprecision(16) << std::scientific;

	for (auto& cell_data: grid) {
		outfile
			<< Mas(cell_data) << " "
			<< Mom(cell_data)[0] << " "
			<< Mom(cell_data)[1] << " "
			<< Mom(cell_data)[2] << " "
			<< Nrj(cell_data) << " "
			<< Mag(cell_data)[0] + Bg_Mag(cell_data)[0] << " "
			<< Mag(cell_data)[1] + Bg_Mag(cell_data)[1] << " "
			<< Mag(cell_data)[2] + Bg_Mag(cell_data)[2] << "\n";
	}
	outfile << std::endl;
}


template <class T> T get_relative_error(const T a, const T b)
{
	if (a == T(0) && b == T(0)) {
		return {0};
	}

	return {std::fabs(a - b) / std::max(std::fabs(a), std::fabs(b))};
}


void verify_mhd(
	Grid& grid,
	const std::string& file_name
) {
	std::ifstream infile(file_name);
	if (not infile.good()) {
		// try a subdirectory in case run from repository root
		infile.open("tests/mhd/" + file_name);

		if (not infile.good()) {
			std::cerr <<  __FILE__ << ":" << __LINE__
				<< " Couldn't open file " << file_name
				<< std::endl;
			abort();
		}
	}

	constexpr double
		// maximum allowed relative difference from reference
		max_error = 1e-6,
		// values smaller than this are assumed correct
		min_value = 1e-25;

	for (size_t cell_i = 0; cell_i < std::tuple_size<Grid>::value; cell_i++) {
		auto& cell = grid[cell_i];
		const auto rho = Mas(cell);
		const auto mom = Mom(cell);
		const auto nrj = Nrj(cell);
		auto mag = Mag(cell);
		mag += Bg_Mag(cell);

		auto ref_rho = rho;
		auto ref_mom = mom;
		auto ref_nrj = nrj;
		auto ref_mag = mag;
		infile
			>> ref_rho
			>> ref_mom[0]
			>> ref_mom[1]
			>> ref_mom[2]
			>> ref_nrj
			>> ref_mag[0]
			>> ref_mag[1]
			>> ref_mag[2];

		if (
			(rho > min_value or ref_rho > min_value)
			and get_relative_error(rho, ref_rho) > max_error
		) {
			std::cerr <<  __FILE__ << ":" << __LINE__
				<< " density " << rho << " != " << ref_rho
				<< " at cell " << cell_i
				<< std::endl;
			abort();
		}

		if (
			(mom[0] > min_value or ref_mom[0] > min_value)
			and get_relative_error(mom[0], ref_mom[0]) > max_error
		) {
			std::cerr <<  __FILE__ << ":" << __LINE__
				<< " x momentum " << mom[0] << " != " << ref_mom[0]
				<< " at cell " << cell_i
				<< std::endl;
			abort();
		}
		if (
			(mom[1] > min_value or ref_mom[1] > min_value)
			and get_relative_error(mom[1], ref_mom[1]) > max_error
		) {
			std::cerr <<  __FILE__ << ":" << __LINE__
				<< " y momentum " << mom[1] << " != " << ref_mom[1]
				<< " at cell " << cell_i
				<< std::endl;
			abort();
		}
		if (
			(mom[2] > min_value or ref_mom[2] > min_value)
			and get_relative_error(mom[2], ref_mom[2]) > max_error
		) {
			std::cerr <<  __FILE__ << ":" << __LINE__
				<< " z momentum " << mom[2] << " != " << ref_mom[2]
				<< " at cell " << cell_i
				<< std::endl;
			abort();
		}

		if (
			(nrj > min_value or ref_nrj > min_value)
			and get_relative_error(nrj, ref_nrj) > max_error
		) {
			std::cerr <<  __FILE__ << ":" << __LINE__
				<< " energy " << nrj << " != " << ref_nrj
				<< " at cell " << cell_i
				<< std::endl;
			abort();
		}

		if (
			(mag[0] > min_value or ref_mag[0] > min_value)
			and get_relative_error(mag[0], ref_mag[0]) > max_error
		) {
			std::cerr <<  __FILE__ << ":" << __LINE__
				<< " x magnetic field " << mag[0] << " != " << ref_mag[0]
				<< " at cell " << cell_i
				<< std::endl;
			abort();
		}
		if (
			(mag[1] > min_value or ref_mag[1] > min_value)
			and get_relative_error(mag[1], ref_mag[1]) > max_error
		) {
			std::cerr <<  __FILE__ << ":" << __LINE__
				<< " y magnetic field " << mag[1] << " != " << ref_mag[1]
				<< " at cell " << cell_i
				<< std::endl;
			abort();
		}
		if (
			(mag[2] > min_value or ref_mag[2] > min_value)
			and get_relative_error(mag[2], ref_mag[2]) > max_error
		) {
			std::cerr <<  __FILE__ << ":" << __LINE__
				<< " z magnetic field " << mag[2] << " != " << ref_mag[2]
				<< " at cell " << cell_i
				<< std::endl;
			abort();
		}
	}
}


int main(int argc, char* argv[])
{
	bool save = false, plot = false, no_verify = false, verbose = false;
	std::string solver_str("rusanov");
	boost::program_options::options_description options(
		"Usage: program_name [options], where options are:"
	);
	options.add_options()
		("help", "Print this help message")
		("solver",
			boost::program_options::value<std::string>(&solver_str)
				->default_value(solver_str),
			"Solver to use, available: rusanov")
		("save", "Save end result to ascii file")
		("plot", "Plot results using gnuplot")
		("no-verify", "Do not verify against reference results")
		("verbose", "Print run time information");

	// read options from command line
	boost::program_options::variables_map option_variables;
	boost::program_options::store(boost::program_options::parse_command_line(argc, argv, options), option_variables);
	boost::program_options::notify(option_variables);

	if (option_variables.count("help") > 0) {
		cout << options << endl;
		return EXIT_SUCCESS;
	}

	if (option_variables.count("save") > 0) {
		save = true;
	}

	if (option_variables.count("plot") > 0) {
		plot = true;
	}

	if (option_variables.count("no-verify") > 0) {
		no_verify = true;
	}

	if (option_variables.count("verbose") > 0) {
		verbose = true;
	}

	const pamhd::mhd::Solver solver
		= [&solver_str](){
			if (solver_str == "rusanov") {
				return pamhd::mhd::Solver::rusanov;
			} else if (solver_str == "hll_athena") {
				return pamhd::mhd::Solver::hll_athena;
			} else if (solver_str == "hlld_athena") {
				return pamhd::mhd::Solver::hlld_athena;
			} else if (solver_str == "roe_athena") {
				return pamhd::mhd::Solver::roe_athena;
			} else {
				std::cerr <<  __FILE__ << "(" << __LINE__ << ") Invalid solver: "
					<< solver_str << ", use --help to list available solvers"
					<< std::endl;
				abort();
			}
		}();

	Grid grid;

	if (verbose) {
		cout << "Initializing MHD" << endl;
	}
	initialize_mhd(
		solver_str,
		grid,
		adiabatic_index,
		vacuum_permeability,
		Mas, Mom, Nrj, Mag, Ele,
		Mas_f, Mom_f, Nrj_f, Mag_f
	);

	if (plot) {
		const std::string mhd_gnuplot_file_name
			= plot_mhd(grid, 0, adiabatic_index, vacuum_permeability);
		system(("gnuplot " + mhd_gnuplot_file_name).c_str());
		if (verbose) {
			cout << "Initial state plotted from file " << mhd_gnuplot_file_name << endl;
		}
	}

	const double
		simulation_duration = 1,
		mhd_plot_interval = simulation_duration / 5;
	double
		max_dt = 0,
		simulation_time = 0,
		mhd_next_plot = mhd_plot_interval;
	while (simulation_time < simulation_duration) {

		const double
			CFL = 0.5,
			// don't step over the final simulation time
			until_end = simulation_duration - simulation_time,
			time_step = std::min(CFL * max_dt, until_end);

		max_dt = std::numeric_limits<double>::max();

		if (verbose) {
			cout << "Solving MHD at time " << simulation_time
				<< " s with time step " << time_step << " s" << endl;
		}

		max_dt = std::min(
			max_dt,
			solve_mhd(
				solver,
				grid,
				time_step,
				adiabatic_index,
				vacuum_permeability,
				Mas, Mom, Nrj, Mag, Ele,
				Mas_f, Mom_f, Nrj_f, Mag_f
			)
		);

		simulation_time += time_step;

		if (plot && mhd_next_plot <= simulation_time) {
			mhd_next_plot += mhd_plot_interval;

			if (verbose) {
				cout << "Plotting simulation at time " << simulation_time << endl;
			}
			const std::string
				mhd_gnuplot_file_name
					= plot_mhd(
						grid,
						simulation_time,
						adiabatic_index,
						vacuum_permeability
					);
			system(("gnuplot " + mhd_gnuplot_file_name).c_str());
		}
	}

	if (save) {
		if (verbose) {
			cout << "Saving MHD at time " << simulation_time << endl;
		}
		save_mhd(solver_str, grid);
	}

	if (not no_verify) {
		const std::string reference_name("mhd_staggered_" + solver_str + ".ref");

		if (verbose) {
			cout << "Verifying result against file " << reference_name << endl;
		}

		verify_mhd(grid, reference_name);
	}

	if (verbose) {
		cout << "Simulation finished at time " << simulation_time << endl;
	}

	return EXIT_SUCCESS;
}
