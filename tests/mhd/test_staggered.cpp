/*
MHD test program of PAMHD.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2018, 2019, 2022 Finnish Meteorological Institute
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
#include "cmath"
#include "cstdlib"
#include "fstream"
#include "iostream"
#include "limits"
#include "streambuf"
#include "string"
#include "vector"

#include "boost/filesystem.hpp"
#include "boost/lexical_cast.hpp"
#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "Eigen/Core" // must be included before gensimcell.hpp
#include "mpi.h" // must be included before gensimcell.hpp
#include "gensimcell.hpp"
#include "rapidjson/document.h"
#include "rapidjson/error/en.h"

#include "background_magnetic_field.hpp"
#include "boundaries/geometries.hpp"
#include "boundaries/multivariable_boundaries.hpp"
#include "boundaries/multivariable_initial_conditions.hpp"
#include "grid_options.hpp"
#include "math/staggered.hpp"
#include "mhd/boundaries.hpp"
#include "mhd/common.hpp"
#include "mhd/hll_athena.hpp"
#include "mhd/hlld_athena.hpp"
#include "mhd/initialize.hpp"
#include "mhd/options.hpp"
#include "mhd/roe_athena.hpp"
#include "mhd/rusanov.hpp"
#include "mhd/save.hpp"
#include "mhd/solve_staggered.hpp"
#include "mhd/variables.hpp"
#include "simulation_options.hpp"
#include "variables.hpp"


using namespace std;

// data stored in every cell of simulation grid
using Cell = pamhd::mhd::Cell_Staggered;
using Grid = dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry>;

// returns reference to background magnetic field at +X face of given cell
const auto Bg_B_Pos_X = [](Cell& cell_data)->auto& {
		return cell_data[pamhd::Bg_Magnetic_Field_Pos_X()];
	};
// reference to +Y face background magnetic field
const auto Bg_B_Pos_Y = [](Cell& cell_data)->auto& {
		return cell_data[pamhd::Bg_Magnetic_Field_Pos_Y()];
	};
// ref to +Z face bg B
const auto Bg_B_Pos_Z = [](Cell& cell_data)->auto& {
		return cell_data[pamhd::Bg_Magnetic_Field_Pos_Z()];
	};

// returns reference to total mass density in given cell
const auto Mas = [](Cell& cell_data)->auto& {
		return cell_data[pamhd::mhd::MHD_State_Conservative()][pamhd::mhd::Mass_Density()];
	};
const auto Mom = [](Cell& cell_data)->auto& {
		return cell_data[pamhd::mhd::MHD_State_Conservative()][pamhd::mhd::Momentum_Density()];
	};
const auto Nrj = [](Cell& cell_data)->auto& {
		return cell_data[pamhd::mhd::MHD_State_Conservative()][pamhd::mhd::Total_Energy_Density()];
	};
const auto Mag = [](Cell& cell_data)->auto& {
		return cell_data[pamhd::mhd::MHD_State_Conservative()][pamhd::Magnetic_Field()];
	};
const auto Face_B = [](Cell& cell_data)->auto& {
		return cell_data[pamhd::Face_Magnetic_Field()];
	};
// divergence of magnetic field
const auto Mag_div = [](Cell& cell_data)->auto&{
		return cell_data[pamhd::Magnetic_Field_Divergence()];
	};
const auto Edge_E = [](Cell& cell_data)->auto& {
		return cell_data[pamhd::Edge_Electric_Field()];
	};

// solver info variable for boundary logic
const auto Sol_Info = [](Cell& cell_data)->auto& {
		return cell_data[pamhd::mhd::Solver_Info()];
	};
// flux of mass density through positive x face of cell
const auto Mas_fx = [](Cell& cell_data)->auto& {
		return cell_data[pamhd::mhd::MHD_Flux_Pos_X()][pamhd::mhd::Mass_Density()];
	};
const auto Mas_fy = [](Cell& cell_data)->auto& {
		return cell_data[pamhd::mhd::MHD_Flux_Pos_Y()][pamhd::mhd::Mass_Density()];
	};
const auto Mas_fz = [](Cell& cell_data)->auto& {
		return cell_data[pamhd::mhd::MHD_Flux_Pos_Z()][pamhd::mhd::Mass_Density()];
	};
const auto Mom_fx = [](Cell& cell_data)->auto& {
		return cell_data[pamhd::mhd::MHD_Flux_Pos_X()][pamhd::mhd::Momentum_Density()];
	};
const auto Mom_fy = [](Cell& cell_data)->auto& {
		return cell_data[pamhd::mhd::MHD_Flux_Pos_Y()][pamhd::mhd::Momentum_Density()];
	};
const auto Mom_fz = [](Cell& cell_data)->auto& {
		return cell_data[pamhd::mhd::MHD_Flux_Pos_Z()][pamhd::mhd::Momentum_Density()];
	};
const auto Nrj_fx = [](Cell& cell_data)->auto& {
		return cell_data[pamhd::mhd::MHD_Flux_Pos_X()][pamhd::mhd::Total_Energy_Density()];
	};
const auto Nrj_fy = [](Cell& cell_data)->auto& {
		return cell_data[pamhd::mhd::MHD_Flux_Pos_Y()][pamhd::mhd::Total_Energy_Density()];
	};
const auto Nrj_fz = [](Cell& cell_data)->auto& {
		return cell_data[pamhd::mhd::MHD_Flux_Pos_Z()][pamhd::mhd::Total_Energy_Density()];
	};
const auto Mag_fx = [](Cell& cell_data)->auto& {
		return cell_data[pamhd::mhd::MHD_Flux_Pos_X()][pamhd::Magnetic_Field()];
	};
const auto Mag_fy = [](Cell& cell_data)->auto& {
		return cell_data[pamhd::mhd::MHD_Flux_Pos_Y()][pamhd::Magnetic_Field()];
	};
const auto Mag_fz = [](Cell& cell_data)->auto& {
		return cell_data[pamhd::mhd::MHD_Flux_Pos_Z()][pamhd::Magnetic_Field()];
	};


int main(int argc, char* argv[])
{
	using std::asin;
	using std::atan2;
	using std::get;
	using std::min;
	using std::sqrt;

	/*
	Initialize MPI
	*/

	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		std::cerr << "Couldn't initialize MPI." << std::endl;
		abort();
	}

	MPI_Comm comm = MPI_COMM_WORLD;

	int rank = 0, comm_size = 0;
	if (MPI_Comm_rank(comm, &rank) != MPI_SUCCESS) {
		std::cerr << "Couldn't obtain MPI rank." << std::endl;
		abort();
	}
	if (MPI_Comm_size(comm, &comm_size) != MPI_SUCCESS) {
		std::cerr << "Couldn't obtain size of MPI communicator." << std::endl;
		abort();
	}

	// intialize Zoltan
	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
		std::cerr << "Zoltan_Initialize failed." << std::endl;
		abort();
	}

	/*
	Parse configuration file
	*/

	if (argc != 2) {
		if (argc < 2 and rank == 0) {
			std::cerr
				<< "Name of configuration file required."
				<< std::endl;
		}
		if (argc > 2 and rank == 0) {
			std::cerr
				<< "Too many arguments given to " << argv[0]
				<< ": " << argc - 1 << ", should be 1"
				<< std::endl;
		}
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	std::ifstream json_file(argv[1]);
	if (not json_file.good()) {
		if (rank == 0) {
			std::cerr << "Couldn't open configuration file " << argv[1] << std::endl;
		}
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	std::string json{
		std::istreambuf_iterator<char>(json_file),
		std::istreambuf_iterator<char>()
	};

	rapidjson::Document document;
	document.Parse(json.c_str());
	if (document.HasParseError()) {
		std::cerr << "Couldn't parse json data in file " << argv[1]
			<< " at character position " << document.GetErrorOffset()
			<< ": " << rapidjson::GetParseError_En(document.GetParseError())
			<< std::endl;
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	pamhd::Options options_sim{document};
	pamhd::grid::Options options_grid{document};
	pamhd::mhd::Options options_mhd{document};

	if (rank == 0 and options_sim.output_directory != "") {
		std::cout << "Saving results into directory " << options_sim.output_directory << std::endl;
		try {
			boost::filesystem::create_directories(options_sim.output_directory);
		} catch (const boost::filesystem::filesystem_error& e) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << ") "
				"Couldn't create output directory "
				<< options_sim.output_directory << ": "
				<< e.what()
				<< std::endl;
			abort();
		}
	}

	using geometry_id_t = unsigned int;

	pamhd::boundaries::Geometries<
		geometry_id_t,
		std::array<double, 3>,
		double,
		uint64_t
	> geometries;
	geometries.set(document);

	pamhd::boundaries::Multivariable_Initial_Conditions<
		geometry_id_t,
		pamhd::mhd::Number_Density,
		pamhd::mhd::Velocity,
		pamhd::mhd::Pressure,
		pamhd::Magnetic_Field
	> initial_conditions;
	initial_conditions.set(document);

	pamhd::boundaries::Multivariable_Boundaries<
		uint64_t,
		geometry_id_t,
		pamhd::mhd::Number_Density,
		pamhd::mhd::Velocity,
		pamhd::mhd::Pressure,
		pamhd::Magnetic_Field
	> boundaries;
	boundaries.set(document);

	pamhd::Background_Magnetic_Field<
		double,
		pamhd::Magnetic_Field::data_type
	> background_B;
	background_B.set(document);

	const auto mhd_solver
		= [&options_mhd, &background_B, &rank](){
			if (options_mhd.solver == "rusanov") {
				return pamhd::mhd::Solver::rusanov;
			} else if (options_mhd.solver == "hll-athena") {
				return pamhd::mhd::Solver::hll_athena;
			} else if (options_mhd.solver == "hlld-athena") {
				if (background_B.exists() and rank == 0) {
					std::cout << "NOTE: background magnetic field ignored by hlld-athena solver." << std::endl;
				}
				return pamhd::mhd::Solver::hlld_athena;
			} else if (options_mhd.solver == "roe-athena") {
				if (background_B.exists() and rank == 0) {
					std::cout << "NOTE: background magnetic field ignored by roe-athena solver." << std::endl;
				}
				return pamhd::mhd::Solver::roe_athena;
			} else {
				std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
					<< "Unsupported solver: " << options_mhd.solver
					<< std::endl;
				abort();
			}
		}();

	/*
	Initialize simulation grid
	*/
	const unsigned int neighborhood_size = 1;
	const auto& number_of_cells = options_grid.get_number_of_cells();
	const auto& periodic = options_grid.get_periodic();

	Grid grid;
	grid
		.set_initial_length(number_of_cells)
		.set_neighborhood_length(neighborhood_size)
		.set_periodic(periodic[0], periodic[1], periodic[2])
		.set_load_balancing_method(options_sim.lb_name)
		.set_maximum_refinement_level(0)
		.initialize(comm)
		.balance_load();

	// set grid geometry
	const std::array<double, 3>
		simulation_volume
			= options_grid.get_volume(),
		cell_volume{
			simulation_volume[0] / number_of_cells[0],
			simulation_volume[1] / number_of_cells[1],
			simulation_volume[2] / number_of_cells[2]
		};

	dccrg::Cartesian_Geometry::Parameters geom_params;
	geom_params.start = options_grid.get_start();
	geom_params.level_0_cell_length = cell_volume;

	try {
		grid.set_geometry(geom_params);
	} catch (...) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Couldn't set grid geometry."
			<< std::endl;
		abort();
	}

	// update owner process of cells for saving into file
	for (const auto& cell: grid.local_cells()) {
		(*cell.data)[pamhd::MPI_Rank()] = rank;
	}

	// assign cells into boundary geometries
	for (const auto& cell: grid.local_cells()) {
		const auto
			start = grid.geometry.get_min(cell.id),
			end = grid.geometry.get_max(cell.id);
		geometries.overlaps(start, end, cell.id);
	}

	/*
	Simulate
	*/

	const double time_end = options_sim.time_start + options_sim.time_length;
	double
		max_dt_mhd = 0,
		simulation_time = options_sim.time_start,
		next_mhd_save = options_mhd.save_n;

	// initialize MHD
	if (rank == 0) {
		cout << "Initializing MHD... " << endl;
	}

	pamhd::mhd::initialize_magnetic_field<pamhd::Magnetic_Field>(
		geometries,
		initial_conditions,
		background_B,
		grid,
		simulation_time,
		options_sim.vacuum_permeability,
		Face_B, Mag_fx,
		Bg_B_Pos_X, Bg_B_Pos_Y, Bg_B_Pos_Z
	);
	for (const auto& cell: grid.local_cells()) {
		Mag_fy(*cell.data) =
		Mag_fz(*cell.data) = {0, 0, 0};
	}

	// update background field between processes
	// update face B for calculating vol B
	Cell::set_transfer_all(
		true,
		pamhd::Face_Magnetic_Field(),
		pamhd::Bg_Magnetic_Field_Pos_X(),
		pamhd::Bg_Magnetic_Field_Pos_Y(),
		pamhd::Bg_Magnetic_Field_Pos_Z()
	);
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(
		false,
		pamhd::Face_Magnetic_Field(),
		pamhd::Bg_Magnetic_Field_Pos_X(),
		pamhd::Bg_Magnetic_Field_Pos_Y(),
		pamhd::Bg_Magnetic_Field_Pos_Z()
	);

	pamhd::mhd::average_magnetic_field<pamhd::mhd::Solver_Info>(
		grid.local_cells(),
		Mas, Mom, Nrj, Mag, Face_B,
		Sol_Info,
		options_sim.adiabatic_index,
		options_sim.vacuum_permeability,
		false
	);

	// update vol B for calculating fluid pressure
	Cell::set_transfer_all(true, pamhd::mhd::MHD_State_Conservative());
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false, pamhd::mhd::MHD_State_Conservative());

	pamhd::mhd::initialize_fluid(
		geometries,
		initial_conditions,
		grid,
		simulation_time,
		options_sim.adiabatic_index,
		options_sim.vacuum_permeability,
		options_sim.proton_mass,
		true,
		Mas, Mom, Nrj, Mag,
		Mas_fx, Mom_fx, Nrj_fx
	);
	for (const auto& cell: grid.local_cells()) {
		Mas_fy(*cell.data) =
		Mas_fz(*cell.data) =
		Nrj_fy(*cell.data) =
		Nrj_fz(*cell.data) = 0;
		Mom_fy(*cell.data) =
		Mom_fz(*cell.data) =
		Mag_fy(*cell.data) =
		Mag_fz(*cell.data) = {0, 0, 0};
	}

	pamhd::mhd::apply_magnetic_field_boundaries(
		grid,
		boundaries,
		geometries,
		simulation_time,
		Face_B
	);

	Cell::set_transfer_all(true, pamhd::Face_Magnetic_Field());
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false, pamhd::Face_Magnetic_Field());

	pamhd::mhd::average_magnetic_field<pamhd::mhd::Solver_Info>(
		grid.local_cells(),
		Mas, Mom, Nrj, Mag, Face_B,
		Sol_Info,
		options_sim.adiabatic_index,
		options_sim.vacuum_permeability,
		true
	);

	Cell::set_transfer_all(true, pamhd::mhd::MHD_State_Conservative());
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false, pamhd::mhd::MHD_State_Conservative());

	pamhd::mhd::apply_fluid_boundaries(
		grid,
		boundaries,
		geometries,
		simulation_time,
		Mas, Mom, Nrj, Mag,
		options_sim.proton_mass,
		options_sim.adiabatic_index,
		options_sim.vacuum_permeability
	);

	if (rank == 0) {
		cout << "Done initializing MHD" << endl;
	}

	/*
	Classify cells into normal, boundary and dont_solve
	*/

	Cell::set_transfer_all(true, pamhd::mhd::Solver_Info());
	pamhd::mhd::set_solver_info<pamhd::mhd::Solver_Info>(
		grid, boundaries, geometries, Sol_Info
	);
	Cell::set_transfer_all(false, pamhd::mhd::Solver_Info());

	size_t simulation_step = 0;
	while (simulation_time < time_end) {
		simulation_step++;

		/*
		Get maximum allowed time step
		*/
		double
			// don't step over the final simulation time
			until_end = time_end - simulation_time,
			local_time_step = min(options_mhd.time_step_factor * max_dt_mhd, until_end),
			time_step = -1;

		if (
			MPI_Allreduce(
				&local_time_step,
				&time_step,
				1,
				MPI_DOUBLE,
				MPI_MIN,
				comm
			) != MPI_SUCCESS
		) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< ": Couldn't set reduce time step."
				<< std::endl;
			abort();
		}

		max_dt_mhd = std::numeric_limits<double>::max();

		if (rank == 0) {
			cout << "Solving MHD at time " << simulation_time
				<< " s with time step " << time_step << " s";
		}

		/*
		Solve
		*/

		Cell::set_transfer_all(true, pamhd::mhd::MHD_State_Conservative());
		grid.start_remote_neighbor_copy_updates();

		double solve_max_dt = pamhd::mhd::solve_staggered<pamhd::mhd::Solver_Info>(
			mhd_solver,
			grid.inner_cells(),
			grid,
			time_step,
			options_sim.adiabatic_index,
			options_sim.vacuum_permeability,
			Mas, Mom, Nrj, Mag,
			Bg_B_Pos_X, Bg_B_Pos_Y, Bg_B_Pos_Z,
			Mas_fx, Mas_fy, Mas_fz,
			Mom_fx, Mom_fy, Mom_fz,
			Nrj_fx, Nrj_fy, Nrj_fz,
			Mag_fx, Mag_fy, Mag_fz,
			Sol_Info
		);
		max_dt_mhd = min(solve_max_dt, max_dt_mhd);

		grid.wait_remote_neighbor_copy_update_receives();

		solve_max_dt = pamhd::mhd::solve_staggered<pamhd::mhd::Solver_Info>(
			mhd_solver,
			grid.outer_cells(),
			grid,
			time_step,
			options_sim.adiabatic_index,
			options_sim.vacuum_permeability,
			Mas, Mom, Nrj, Mag,
			Bg_B_Pos_X, Bg_B_Pos_Y, Bg_B_Pos_Z,
			Mas_fx, Mas_fy, Mas_fz,
			Mom_fx, Mom_fy, Mom_fz,
			Nrj_fx, Nrj_fy, Nrj_fz,
			Mag_fx, Mag_fy, Mag_fz,
			Sol_Info
		);
		max_dt_mhd = min(solve_max_dt, max_dt_mhd);

		grid.wait_remote_neighbor_copy_update_sends();
		Cell::set_transfer_all(false, pamhd::mhd::MHD_State_Conservative());


		/*
		Apply solution
		*/

		Cell::set_transfer_all(
			true,
			pamhd::mhd::MHD_Flux_Pos_X(),
			pamhd::mhd::MHD_Flux_Pos_Y(),
			pamhd::mhd::MHD_Flux_Pos_Z()
		);
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(
			false,
			pamhd::mhd::MHD_Flux_Pos_X(),
			pamhd::mhd::MHD_Flux_Pos_Y(),
			pamhd::mhd::MHD_Flux_Pos_Z()
		);
		// TODO: split into inner and outer cells
		pamhd::mhd::apply_fluxes_staggered<pamhd::mhd::Solver_Info>(
			grid, time_step,
			Mas, Mom, Nrj, Mag, Edge_E,
			Mas_fx, Mas_fy, Mas_fz,
			Mom_fx, Mom_fy, Mom_fz,
			Nrj_fx, Nrj_fy, Nrj_fz,
			Mag_fx, Mag_fy, Mag_fz,
			Sol_Info
		);
		Cell::set_transfer_all(true, pamhd::Edge_Electric_Field());
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false, pamhd::Edge_Electric_Field());
		pamhd::mhd::solve_B<pamhd::mhd::Solver_Info>(
			grid.local_cells(),
			grid,
			time_step,
			Face_B,
			Edge_E,
			Sol_Info
		);

		// constant thermal pressure when updating vol B after solution
		Cell::set_transfer_all(true, pamhd::Face_Magnetic_Field());
		grid.start_remote_neighbor_copy_updates();

		pamhd::mhd::average_magnetic_field<pamhd::mhd::Solver_Info>(
			grid.inner_cells(),
			Mas, Mom, Nrj, Mag, Face_B,
			Sol_Info,
			options_sim.adiabatic_index,
			options_sim.vacuum_permeability,
			true
		);

		grid.wait_remote_neighbor_copy_update_receives();

		pamhd::mhd::average_magnetic_field<pamhd::mhd::Solver_Info>(
			grid.outer_cells(),
			Mas, Mom, Nrj, Mag, Face_B,
			Sol_Info,
			options_sim.adiabatic_index,
			options_sim.vacuum_permeability,
			true
		);

		grid.wait_remote_neighbor_copy_update_sends();
		Cell::set_transfer_all(false, pamhd::Face_Magnetic_Field());

		simulation_time += time_step;


		/*
		Update boundaries
		*/

		// TODO: split into inner and outer cells
		pamhd::mhd::apply_magnetic_field_boundaries(
			grid,
			boundaries,
			geometries,
			simulation_time,
			Face_B
		);


		Cell::set_transfer_all(true, pamhd::Face_Magnetic_Field());
		grid.start_remote_neighbor_copy_updates();

		pamhd::mhd::average_magnetic_field<pamhd::mhd::Solver_Info>(
			grid.inner_cells(),
			Mas, Mom, Nrj, Mag, Face_B,
			Sol_Info,
			options_sim.adiabatic_index,
			options_sim.vacuum_permeability,
			false
		);

		grid.wait_remote_neighbor_copy_update_receives();

		pamhd::mhd::average_magnetic_field<pamhd::mhd::Solver_Info>(
			grid.outer_cells(),
			Mas, Mom, Nrj, Mag, Face_B,
			Sol_Info,
			options_sim.adiabatic_index,
			options_sim.vacuum_permeability,
			false
		);

		grid.wait_remote_neighbor_copy_update_sends();
		Cell::set_transfer_all(false, pamhd::Face_Magnetic_Field());

		const auto total_div = pamhd::math::get_divergence_staggered(
			grid.local_cells(),
			grid,
			Face_B,
			Mag_div,
			Sol_Info
		);
		if (rank == 0) {
			cout << " total divergence " << total_div << endl;
		}
		Cell::set_transfer_all(true, pamhd::Magnetic_Field_Divergence());
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false, pamhd::Magnetic_Field_Divergence());

		Cell::set_transfer_all(true, pamhd::mhd::MHD_State_Conservative());
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false, pamhd::mhd::MHD_State_Conservative());

		// TODO: split into inner and outer cells
		pamhd::mhd::apply_fluid_boundaries(
			grid,
			boundaries,
			geometries,
			simulation_time,
			Mas, Mom, Nrj, Mag,
			options_sim.proton_mass,
			options_sim.adiabatic_index,
			options_sim.vacuum_permeability
		);

		/*
		Save simulation to disk
		*/

		if (
			(
				options_mhd.save_n >= 0
				and (
					simulation_time == options_sim.time_start
					or simulation_time >= time_end
				)
			) or (options_mhd.save_n > 0 and simulation_time >= next_mhd_save)
		) {
			if (next_mhd_save <= simulation_time) {
				next_mhd_save
					+= options_mhd.save_n
					* ceil((simulation_time - next_mhd_save) / options_mhd.save_n);
			}

			if (rank == 0) {
				cout << "Saving MHD at time " << simulation_time << endl;
			}

			constexpr uint64_t file_version = 2;
			if (
				not pamhd::mhd::save(
					boost::filesystem::canonical(
						boost::filesystem::path(options_sim.output_directory)
					).append("mhd_staggered_").generic_string(),
					grid,
					file_version,
					simulation_step,
					simulation_time,
					options_sim.adiabatic_index,
					options_sim.proton_mass,
					options_sim.vacuum_permeability,
					pamhd::mhd::MHD_State_Conservative(),
					pamhd::mhd::Solver_Info(),
					pamhd::MPI_Rank(),
					pamhd::Face_Magnetic_Field(),
					pamhd::Edge_Electric_Field(),
					pamhd::Bg_Magnetic_Field_Pos_X(),
					pamhd::Bg_Magnetic_Field_Pos_Y(),
					pamhd::Bg_Magnetic_Field_Pos_Z(),
					pamhd::Magnetic_Field_Divergence()
				)
			) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
					"Couldn't save mhd result."
					<< std::endl;
				abort();
			}
		}
	}

	if (rank == 0) {
		cout << "Simulation finished at time " << simulation_time << endl;
	}
	MPI_Finalize();

	return EXIT_SUCCESS;
}
