/*
Test program of PAMHD with one planet in solar wind.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2018, 2019, 2022, 2023, 2024 Finnish Meteorological Institute
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


Author(s): Ilja Honkonen
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
#include "solar_wind_box_options.hpp"
#include "boundaries/geometries.hpp"
#include "boundaries/multivariable_boundaries.hpp"
#include "boundaries/multivariable_initial_conditions.hpp"
#include "grid/amr.hpp"
#include "grid/options.hpp"
#include "grid/solar_wind_box.hpp"
#include "grid/variables.hpp"
#include "math/staggered.hpp"
#include "mhd/amr.hpp"
#include "mhd/boundaries.hpp"
#include "mhd/common.hpp"
#include "mhd/hll_athena.hpp"
#include "mhd/hlld_athena.hpp"
#include "mhd/initialize.hpp"
#include "mhd/initialize_staggered.hpp"
#include "mhd/options.hpp"
#include "mhd/roe_athena.hpp"
#include "mhd/rusanov.hpp"
#include "mhd/save.hpp"
#include "mhd/solar_wind_box.hpp"
#include "mhd/solve_staggered.hpp"
#include "mhd/variables.hpp"
#include "simulation_options.hpp"
#include "variables.hpp"


// data stored in every cell of simulation grid
using Cell = pamhd::mhd::Cell_Staggered;
using Grid = dccrg::Dccrg<
	Cell,
	dccrg::Cartesian_Geometry,
	std::tuple<pamhd::grid::Cell_Is_Local>,
	std::tuple<
		pamhd::grid::Face_Neighbor,
		pamhd::grid::Edge_Neighbor,
		pamhd::grid::Relative_Size,
		pamhd::grid::Neighbor_Is_Local>
>;

// returns reference to background magnetic field on cell faces
const auto Bg_B = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::Bg_Magnetic_Field()];
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
const auto Face_dB = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::Face_dB()];
};
// divergence of magnetic field
const auto Mag_div = [](Cell& cell_data)->auto&{
	return cell_data[pamhd::Magnetic_Field_Divergence()];
};

// solver info variable for boundary logic
const auto SInfo = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::Solver_Info()];
};

const auto Substep = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::Substepping_Period()];
};
const auto Substep_Min = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::Substep_Min()];
};
const auto Substep_Max = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::Substep_Max()];
};

const auto Timestep = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::Timestep()];
};

const auto Max_v = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::Max_Velocity()];
};

// flux of mass density through positive x face of cell
const auto Mas_pfx = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](0, 1)[pamhd::mhd::Mass_Density()];
};
const auto Mas_nfx = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](0, 0)[pamhd::mhd::Mass_Density()];
};
const auto Mas_pfy = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](1, 1)[pamhd::mhd::Mass_Density()];
};
const auto Mas_nfy = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](1, 0)[pamhd::mhd::Mass_Density()];
};
const auto Mas_pfz = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](2, 1)[pamhd::mhd::Mass_Density()];
};
const auto Mas_nfz = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](2, 0)[pamhd::mhd::Mass_Density()];
};
const auto Mom_pfx = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](0, 1)[pamhd::mhd::Momentum_Density()];
};
const auto Mom_nfx = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](0, 0)[pamhd::mhd::Momentum_Density()];
};
const auto Mom_pfy = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](1, 1)[pamhd::mhd::Momentum_Density()];
};
const auto Mom_nfy = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](1, 0)[pamhd::mhd::Momentum_Density()];
};
const auto Mom_pfz = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](2, 1)[pamhd::mhd::Momentum_Density()];
};
const auto Mom_nfz = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](2, 0)[pamhd::mhd::Momentum_Density()];
};
const auto Nrj_pfx = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](0, 1)[pamhd::mhd::Total_Energy_Density()];
};
const auto Nrj_nfx = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](0, 0)[pamhd::mhd::Total_Energy_Density()];
};
const auto Nrj_pfy = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](1, 1)[pamhd::mhd::Total_Energy_Density()];
};
const auto Nrj_nfy = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](1, 0)[pamhd::mhd::Total_Energy_Density()];
};
const auto Nrj_pfz = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](2, 1)[pamhd::mhd::Total_Energy_Density()];
};
const auto Nrj_nfz = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](2, 0)[pamhd::mhd::Total_Energy_Density()];
};
const auto Mag_pfx = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](0, 1)[pamhd::Magnetic_Field()];
};
const auto Mag_nfx = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](0, 0)[pamhd::Magnetic_Field()];
};
const auto Mag_pfy = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](1, 1)[pamhd::Magnetic_Field()];
};
const auto Mag_nfy = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](1, 0)[pamhd::Magnetic_Field()];
};
const auto Mag_pfz = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](2, 1)[pamhd::Magnetic_Field()];
};
const auto Mag_nfz = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::mhd::MHD_Flux()](2, 0)[pamhd::Magnetic_Field()];
};

// collections of above to shorten function arguments
const auto Mas_fs = std::make_tuple(
	Mas_nfx, Mas_pfx, Mas_nfy, Mas_pfy, Mas_nfz, Mas_pfz
);
const auto Mom_fs = std::make_tuple(
	Mom_nfx, Mom_pfx, Mom_nfy, Mom_pfy, Mom_nfz, Mom_pfz
);
const auto Nrj_fs = std::make_tuple(
	Nrj_nfx, Nrj_pfx, Nrj_nfy, Nrj_pfy, Nrj_nfz, Nrj_pfz
);
const auto Mag_fs = std::make_tuple(
	Mag_nfx, Mag_pfx, Mag_nfy, Mag_pfy, Mag_nfz, Mag_pfz
);

const auto Ref_min = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::grid::Target_Refinement_Level_Min()];
};
const auto Ref_max = [](Cell& cell_data)->auto& {
	return cell_data[pamhd::grid::Target_Refinement_Level_Max()];
};


template<
	class JSON,
	class Cells
> void apply_solar_wind_boundaries(
	const JSON& json,
	const Cells& solar_wind_cells,
	const int& solar_wind_dir,
	const double& sim_time,
	const double& adiabatic_index,
	const double& vacuum_permeability,
	const double& proton_mass
) try {
	const auto [
		mass, pressure, velocity, magnetic_field
	] = pamhd::mhd::get_solar_wind_parameters(
		json, sim_time, proton_mass
	);
	for (const auto& cell: solar_wind_cells) {
		Mas(*cell.data) = mass;
		Mom(*cell.data) = {
			mass*velocity[0],
			mass*velocity[1],
			mass*velocity[2]
		};
		// prevent div(B) at solar wind boundary
		if (solar_wind_dir != abs(1)) {
			Mag(*cell.data)[0]     =
			Face_B(*cell.data)(-1) =
			Face_B(*cell.data)(+1) = magnetic_field[0];
		}
		if (solar_wind_dir != abs(2)) {
			Mag(*cell.data)[1]     =
			Face_B(*cell.data)(-2) =
			Face_B(*cell.data)(+2) = magnetic_field[1];
		}
		if (solar_wind_dir != abs(3)) {
			Mag(*cell.data)[2]     =
			Face_B(*cell.data)(-3) =
			Face_B(*cell.data)(+3) = magnetic_field[2];
		}
		Nrj(*cell.data) = pamhd::mhd::get_total_energy_density(
			mass, velocity, pressure, magnetic_field,
			adiabatic_index, vacuum_permeability
		);
	}
} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


template<
	class Grid,
	class JSON,
	class SW_Cells,
	class Face_Cells,
	class Edge_Cells,
	class Vert_Cells,
	class Planet_Cells
> void apply_boundaries(
	Grid& grid,
	const double& sim_time,
	const JSON& json,
	const int& solar_wind_dir,
	const SW_Cells& solar_wind_cells,
	const Face_Cells& face_cells,
	const Edge_Cells& edge_cells,
	const Vert_Cells& vert_cells,
	const Planet_Cells& planet_cells,
	const double& adiabatic_index,
	const double& vacuum_permeability,
	const double& proton_mass
) try {
	// boundary and normal cell share face
	for (int dir: {-3,-2,-1,+1,+2,+3}) {
		for (const auto& cell: face_cells(dir)) {
			for (const auto& neighbor: cell.neighbors_of) {
				const auto& fn = neighbor.face_neighbor;
				if (fn != -dir) continue;
				Mas(*cell.data) = Mas(*neighbor.data);
				Mom(*cell.data) = Mom(*neighbor.data);
				Nrj(*cell.data) = Nrj(*neighbor.data);
				Mag(*cell.data) = Mag(*neighbor.data);
				for (int dir2: {-3,-2,-1,+1,+2,+3}) {
					if (dir2 == fn or dir2 == -fn) {
						Face_B(*cell.data)(dir2) = Face_B(*neighbor.data)(-fn);
					} else {
						Face_B(*cell.data)(dir2) = Face_B(*neighbor.data)(dir2);
					}
				}
			}
		}
	}

	// boundary and normal cell share edge
	for (int dim: {0, 1, 2})
	for (int dir1: {-1, +1})
	for (int dir2: {-1, +1}) {
		for (const auto& cell: edge_cells(dim, dir1, dir2)) {
			for (const auto& neighbor: cell.neighbors_of) {
				const auto& en = neighbor.edge_neighbor;
				if (en[0] != dim or en[1] != -dir1 or en[2] != -dir2) continue;
				Mas(*cell.data) = Mas(*neighbor.data);
				Mom(*cell.data) = Mom(*neighbor.data);
				Nrj(*cell.data) = Nrj(*neighbor.data);
				Mag(*cell.data) = Mag(*neighbor.data);
				if (dim == 0) {
					Face_B(*cell.data)(-1) = Face_B(*neighbor.data)(-1);
					Face_B(*cell.data)(+1) = Face_B(*neighbor.data)(+1);
					if (dir1 < 0) {
						Face_B(*cell.data)(-2) =
						Face_B(*cell.data)(+2) = Face_B(*neighbor.data)(-2);
					} else {
						Face_B(*cell.data)(-2) =
						Face_B(*cell.data)(+2) = Face_B(*neighbor.data)(+2);
					}
					if (dir2 < 0) {
						Face_B(*cell.data)(-3) =
						Face_B(*cell.data)(+3) = Face_B(*neighbor.data)(-3);
					} else {
						Face_B(*cell.data)(-3) =
						Face_B(*cell.data)(+3) = Face_B(*neighbor.data)(+3);
					}
				}
				if (dim == 1) {
					Face_B(*cell.data)(-2) = Face_B(*neighbor.data)(-2);
					Face_B(*cell.data)(+2) = Face_B(*neighbor.data)(+2);
					if (dir1 < 0) {
						Face_B(*cell.data)(-1) =
						Face_B(*cell.data)(+1) = Face_B(*neighbor.data)(-1);
					} else {
						Face_B(*cell.data)(-1) =
						Face_B(*cell.data)(+1) = Face_B(*neighbor.data)(+1);
					}
					if (dir2 < 0) {
						Face_B(*cell.data)(-3) =
						Face_B(*cell.data)(+3) = Face_B(*neighbor.data)(-3);
					} else {
						Face_B(*cell.data)(-3) =
						Face_B(*cell.data)(+3) = Face_B(*neighbor.data)(+3);
					}
				}
				if (dim == 2) {
					Face_B(*cell.data)(-3) = Face_B(*neighbor.data)(-3);
					Face_B(*cell.data)(+3) = Face_B(*neighbor.data)(+3);
					if (dir1 < 0) {
						Face_B(*cell.data)(-1) =
						Face_B(*cell.data)(+1) = Face_B(*neighbor.data)(-1);
					} else {
						Face_B(*cell.data)(-1) =
						Face_B(*cell.data)(+1) = Face_B(*neighbor.data)(+1);
					}
					if (dir2 < 0) {
						Face_B(*cell.data)(-2) =
						Face_B(*cell.data)(+2) = Face_B(*neighbor.data)(-2);
					} else {
						Face_B(*cell.data)(-2) =
						Face_B(*cell.data)(+2) = Face_B(*neighbor.data)(+2);
					}
				}
			}
		}
	}

	// boundary and normal cell share vertex
	for (const auto& cell: vert_cells) {
		for (const auto& neighbor: cell.neighbors_of) {
			const auto& fn = neighbor.face_neighbor;
			if (fn != 0) continue;
			const auto& en = neighbor.edge_neighbor;
			if (en[0] >= 0) continue;

			Mas(*cell.data) = Mas(*neighbor.data);
			Mom(*cell.data) = Mom(*neighbor.data);
			Nrj(*cell.data) = Nrj(*neighbor.data);
			Mag(*cell.data) = Mag(*neighbor.data);
			if (neighbor.x < 0) {
				Face_B(*cell.data)(-1) =
				Face_B(*cell.data)(+1) = Face_B(*neighbor.data)(+1);
			} else {
				Face_B(*cell.data)(-1) =
				Face_B(*cell.data)(+1) = Face_B(*neighbor.data)(-1);
			}
			if (neighbor.y < 0) {
				Face_B(*cell.data)(-2) =
				Face_B(*cell.data)(+2) = Face_B(*neighbor.data)(+2);
			} else {
				Face_B(*cell.data)(-2) =
				Face_B(*cell.data)(+2) = Face_B(*neighbor.data)(-2);
			}
			if (neighbor.z < 0) {
				Face_B(*cell.data)(-3) =
				Face_B(*cell.data)(+3) = Face_B(*neighbor.data)(+3);
			} else {
				Face_B(*cell.data)(-3) =
				Face_B(*cell.data)(+3) = Face_B(*neighbor.data)(-3);
			}
		}
	}

	// planetary boundary cells
	for (const auto& cell: planet_cells) {
		Mas(*cell.data) = proton_mass * 1e9;
		Mom(*cell.data) = {0, 0, 0};
		Mag(*cell.data) = {0, 0, 0};
		Nrj(*cell.data) = pamhd::mhd::get_total_energy_density(
			Mas(*cell.data),
			Mom(*cell.data),
			1e-11,
			Mag(*cell.data),
			adiabatic_index,
			vacuum_permeability
		);
		Face_B(*cell.data)(-1) =
		Face_B(*cell.data)(+1) =
		Face_B(*cell.data)(-2) =
		Face_B(*cell.data)(+2) =
		Face_B(*cell.data)(-3) =
		Face_B(*cell.data)(+3) = 0;

		pamhd::grid::Face_Type<bool> have_value{false, false, false, false, false, false};
		// corrections to Face_B from normal neighbors
		for (const auto& neighbor: cell.neighbors_of) {
			if (SInfo(*neighbor.data) < 1) continue;

			const auto& fn = neighbor.face_neighbor;
			const auto& en = neighbor.edge_neighbor;

			if (fn != 0) {
				Face_B(*cell.data)(fn) = Face_B(*neighbor.data)(-fn);
				have_value(fn) = true;
			} else if (en[0] >= 0) {
				// TODO
			} else {
				// TODO
			}
		}
		for (int dir: {-3,-2,-1,+1,+2,+3}) {
			if (not have_value(dir)) {
				if (have_value(-dir)) {
					Face_B(*cell.data)(dir) = Face_B(*cell.data)(-dir);
				} else {
				}
			}
		}
		Mag(*cell.data) = {
			0.5*(Face_B(*cell.data)(-1) + Face_B(*cell.data)(+1)),
			0.5*(Face_B(*cell.data)(-2) + Face_B(*cell.data)(+2)),
			0.5*(Face_B(*cell.data)(-3) + Face_B(*cell.data)(+3))
		};
		Nrj(*cell.data) = pamhd::mhd::get_total_energy_density(
			Mas(*cell.data),
			Mom(*cell.data),
			1e-11,
			Mag(*cell.data),
			adiabatic_index,
			vacuum_permeability
		);
	}

	apply_solar_wind_boundaries(
		json, solar_wind_cells,
		solar_wind_dir, sim_time,
		adiabatic_index,
		vacuum_permeability,
		proton_mass
	);
	grid.update_copies_of_remote_neighbors();

} catch (const std::exception& e) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + "): " + e.what());
} catch (...) {
	throw std::runtime_error(__FILE__ "(" + std::to_string(__LINE__) + ")");
}


int main(int argc, char* argv[]) {
	using std::abs;
	using std::ceil;
	using std::cerr;
	using std::cout;
	using std::endl;
	using std::flush;
	using std::max;
	using std::min;
	using std::runtime_error;
	using std::to_string;

	/*
	Initialize MPI
	*/

	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		cerr << "Couldn't initialize MPI." << endl;
		abort();
	}

	MPI_Comm comm = MPI_COMM_WORLD;

	int rank = 0, comm_size = 0;
	if (MPI_Comm_rank(comm, &rank) != MPI_SUCCESS) {
		cerr << "Couldn't obtain MPI rank." << endl;
		abort();
	}
	if (MPI_Comm_size(comm, &comm_size) != MPI_SUCCESS) {
		cerr << "Couldn't obtain size of MPI communicator." << endl;
		abort();
	}

	// initialize Zoltan
	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
		cerr << "Zoltan_Initialize failed." << endl;
		abort();
	}

	/*
	Parse configuration file
	*/

	if (argc != 2) {
		if (argc < 2 and rank == 0) {
			cerr
				<< "Name of configuration file required."
				<< endl;
		}
		if (argc > 2 and rank == 0) {
			cerr
				<< "Too many arguments given to " << argv[0]
				<< ": " << argc - 1 << ", should be 1"
				<< endl;
		}
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	std::ifstream json_file(argv[1]);
	if (not json_file.good()) {
		if (rank == 0) {
			cerr << "Couldn't open configuration file " << argv[1] << endl;
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
		if (rank == 0) {
			cerr << "Couldn't parse json data in file "
				<< argv[1] << " at character position "
				<< document.GetErrorOffset() << ": "
				<< rapidjson::GetParseError_En(document.GetParseError())
				<< endl;
		}
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	pamhd::Options options_sim{document};
	pamhd::grid::Options options_grid{document};
	pamhd::mhd::Options options_mhd{document};
	pamhd::Solar_Wind_Box_Options options_sw{document};
	const int solar_wind_dir = [&options_sw](){
		if (options_sw.sw_dir == "-x") return -1;
		else if (options_sw.sw_dir == "+x") return +1;
		else if (options_sw.sw_dir == "-y") return -2;
		else if (options_sw.sw_dir == "+y") return +2;
		else if (options_sw.sw_dir == "-z") return -3;
		else if (options_sw.sw_dir == "+z") return +3;
		else {
			std::cerr << "Invalid solar wind boundary direction: "
				<< options_sw.sw_dir << std::endl;
			MPI_Finalize();
			return EXIT_FAILURE;
		}
	}();

	if (rank == 0 and options_sim.output_directory != "") {
		cout << "Saving results into directory " << options_sim.output_directory << endl;
		try {
			boost::filesystem::create_directories(options_sim.output_directory);
		} catch (const boost::filesystem::filesystem_error& e) {
			cerr <<  __FILE__ << "(" << __LINE__ << ") "
				"Couldn't create output directory "
				<< options_sim.output_directory << ": "
				<< e.what()
				<< endl;
			abort();
		}
	}

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
					cout << "NOTE: background magnetic field ignored by hlld-athena solver." << endl;
				}
				return pamhd::mhd::Solver::hlld_athena;
			} else if (options_mhd.solver == "roe-athena") {
				if (background_B.exists() and rank == 0) {
					cout << "NOTE: background magnetic field ignored by roe-athena solver." << endl;
				}
				return pamhd::mhd::Solver::roe_athena;
			} else {
				cerr <<  __FILE__ << "(" << __LINE__ << "): "
					<< "Unsupported solver: " << options_mhd.solver
					<< endl;
				abort();
			}
		}();

	/*
	Initialize simulation grid
	*/
	const unsigned int neighborhood_size = 3;
	const auto& number_of_cells = options_grid.get_number_of_cells();
	const size_t min_cell0_count = 5 + 2*options_grid.get_max_ref_lvl();
	for (auto dim: {0, 1, 2}) {
		if (number_of_cells[dim] < min_cell0_count) {
			std::cout << "Number of initial cells in dimension " << dim
				<< " must be at least " << min_cell0_count
				<< " but " << number_of_cells[dim] << " given" << std::endl;
			abort();
		}
	}

	Grid grid; grid
		.set_initial_length(number_of_cells)
		.set_neighborhood_length(neighborhood_size)
		.set_periodic(false, false, false)
		.set_load_balancing_method(options_sim.lb_name)
		.set_maximum_refinement_level(options_grid.get_max_ref_lvl())
		.initialize(comm);
	const auto mrlvl = grid.get_maximum_refinement_level();

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
		cerr << __FILE__ << ":" << __LINE__
			<< ": Couldn't set grid geometry."
			<< endl;
		abort();
	}

	if (rank == 0) {
		cout << "Adapting and balancing grid at time "
			<< options_sim.time_start << "...  " << flush;
	}

	Cell::set_transfer_all(true,
		pamhd::grid::Target_Refinement_Level_Min(),
		pamhd::grid::Target_Refinement_Level_Max(),
		pamhd::mhd::Solver_Info()
	);
	prepare_grid(
		options_sw.inner_bdy_radius,
		grid, SInfo, Ref_max, Ref_min
	);
	pamhd::grid::set_minmax_refinement_level(
		grid.local_cells(), grid, options_grid,
		options_sim.time_start, Ref_min, Ref_max, true
	);
	pamhd::grid::adapt_grid(grid, Ref_min, Ref_max, SInfo);
	grid.balance_load();
	Cell::set_transfer_all(false,
		pamhd::grid::Target_Refinement_Level_Min(),
		pamhd::grid::Target_Refinement_Level_Max(),
		pamhd::mhd::Solver_Info()
	);
	if (rank == 0) {
		cout << "done" << endl;
	}

	auto solar_wind_cells
		= get_solar_wind_cells(solar_wind_dir, grid);
	for (const auto& cell: solar_wind_cells) {
		if (SInfo(*cell.data) != 0) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		if (Ref_max(*cell.data) != 0) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	}

	auto [
		face_cells, edge_cells, vert_cells
	] = get_outflow_cells(grid);
	for (int dir: {-3,-2,-1,+1,+2,+3}) {
		for (const auto& cell: face_cells(dir)) {
			if (SInfo(*cell.data) != 0) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			if (Ref_max(*cell.data) != 0) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		}
	}
	for (int dim: {0, 1, 2})
	for (int dir1: {-1, +1})
	for (int dir2: {-1, +1}) {
		for (const auto& cell: edge_cells(dim, dir1, dir2)) {
			if (SInfo(*cell.data) != 0) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
			if (Ref_max(*cell.data) != 0) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		}
	}
	for (const auto& cell: vert_cells) {
		if (SInfo(*cell.data) != 0) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
		if (Ref_max(*cell.data) != 0) throw runtime_error(__FILE__"(" + to_string(__LINE__) + ")");
	}

	auto planet_cells = pamhd::grid::get_planet_cells(
		options_sw.inner_bdy_radius, grid, SInfo
	);
	for (const auto& cell: planet_cells) {
		const auto [inside, outside]
			= pamhd::grid::at_inner_boundary(options_sw.inner_bdy_radius, cell.id, grid);
		if (not (inside and outside)) throw runtime_error(__FILE__"(" + to_string(__LINE__) + "): " + to_string(cell.id));
		if (SInfo(*cell.data) == 0) {
			bool have_solve = false;
			for (const auto& neighbor: cell.neighbors_of) {
				if (
					neighbor.face_neighbor == 0
					and neighbor.edge_neighbor[0] < 0
				) {
					continue;
				}
				if (SInfo(*neighbor.data) == 1) {
					have_solve = true;
					break;
				}
			}
			if (not have_solve) throw runtime_error(__FILE__"(" + to_string(__LINE__) + "): " + to_string(cell.id));
		}
		if (SInfo(*cell.data) == 1) {
			throw runtime_error(__FILE__"(" + to_string(__LINE__) + "): " + to_string(cell.id) + " " + to_string(SInfo(*cell.data)));
		}
	}

	for (const auto& cell: grid.local_cells()) {
		(*cell.data)[pamhd::MPI_Rank()] = rank;
		Substep(*cell.data) = 1;
		Max_v(*cell.data) = {-1, -1, -1, -1, -1, -1};
	}
	pamhd::mhd::set_minmax_substepping_period(
		options_sim.time_start, grid,
		options_mhd, Substep_Min, Substep_Max);
	Cell::set_transfer_all(true,
		pamhd::mhd::Max_Velocity(),
		pamhd::mhd::Substepping_Period(),
		pamhd::mhd::Substep_Min(),
		pamhd::mhd::Substep_Max(),
		pamhd::MPI_Rank()
	);
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false,
		pamhd::mhd::Max_Velocity(),
		pamhd::mhd::Substepping_Period(),
		pamhd::mhd::Substep_Min(),
		pamhd::mhd::Substep_Max(),
		pamhd::MPI_Rank()
	);

	/*
	Simulate
	*/

	const double time_end = options_sim.time_start + options_sim.time_length;
	double
		simulation_time = options_sim.time_start,
		next_mhd_save = options_mhd.save_n;

	// initialize MHD
	if (rank == 0) {
		cout << "Initializing MHD... " << endl;
	}

	Cell::set_transfer_all(true,
		pamhd::Face_Magnetic_Field(),
		pamhd::mhd::MHD_State_Conservative(),
		pamhd::Bg_Magnetic_Field(),
		pamhd::mhd::Solver_Info()
	);
	initialize_plasma(
		grid, simulation_time,
		document, background_B,
		Mas, Mom, Nrj, Mag, Face_B,
		Face_dB, Bg_B,
		options_sim.adiabatic_index,
		options_sim.vacuum_permeability,
		options_sim.proton_mass
	);
	grid.update_copies_of_remote_neighbors();
	pamhd::mhd::update_B_consistency(
		0, grid.local_cells(),
		Mas, Mom, Nrj, Mag, Face_B,
		SInfo, Substep,
		options_sim.adiabatic_index,
		options_sim.vacuum_permeability,
		true
	);
	grid.update_copies_of_remote_neighbors();

	apply_boundaries(
		grid, simulation_time, document,
		solar_wind_dir, solar_wind_cells,
		face_cells, edge_cells, vert_cells,
		planet_cells,
		options_sim.adiabatic_index,
		options_sim.vacuum_permeability,
		options_sim.proton_mass
	);
	grid.update_copies_of_remote_neighbors();

	pamhd::mhd::update_B_consistency(
		0, grid.local_cells(),
		Mas, Mom, Nrj, Mag, Face_B,
		SInfo, Substep,
		options_sim.adiabatic_index,
		options_sim.vacuum_permeability,
		true
	);
	grid.update_copies_of_remote_neighbors();
	Cell::set_transfer_all(false,
		pamhd::Face_Magnetic_Field(),
		pamhd::mhd::MHD_State_Conservative(),
		pamhd::Bg_Magnetic_Field(),
		pamhd::mhd::Solver_Info()
	);

	if (rank == 0) {
		cout << "Done initializing MHD" << endl;
	}

	// final init with timestep of 0
	pamhd::mhd::timestep(
		mhd_solver, grid, options_mhd, options_sim.time_start,
		0, options_mhd.time_step_factor,
		options_sim.adiabatic_index,
		options_sim.vacuum_permeability,
		Mas, Mom, Nrj, Mag, Face_B, Face_dB, Bg_B,
		Mas_fs, Mom_fs, Nrj_fs, Mag_fs, SInfo,
		Timestep, Substep, Substep_Min, Substep_Max, Max_v
	);
	size_t simulation_step = 0;
	constexpr uint64_t file_version = 3;
	if (options_mhd.save_n >= 0) {
		if (rank == 0) {
			cout << "Saving MHD at time " << simulation_time << endl;
		}
		if (
			not pamhd::mhd::save_staggered(
				boost::filesystem::canonical(
					boost::filesystem::path(options_sim.output_directory)
				).append("mhd_staggered_").generic_string(),
				grid,
				file_version,
				simulation_step,
				simulation_time,
				options_sim.adiabatic_index,
				options_sim.proton_mass,
				options_sim.vacuum_permeability
			)
		) {
			cerr <<  __FILE__ << "(" << __LINE__ << "): "
				"Couldn't save mhd result."
				<< endl;
			abort();
		}
	}

	while (simulation_time < time_end) {
		simulation_step++;

		// don't step over the final simulation time
		double until_end = time_end - simulation_time;
		const double dt = pamhd::mhd::timestep(
			mhd_solver, grid, options_mhd, simulation_time,
			until_end, options_mhd.time_step_factor,
			options_sim.adiabatic_index,
			options_sim.vacuum_permeability,
			Mas, Mom, Nrj, Mag, Face_B, Face_dB, Bg_B,
			Mas_fs, Mom_fs, Nrj_fs, Mag_fs, SInfo,
			Timestep, Substep, Substep_Min, Substep_Max, Max_v
		);
		if (rank == 0) {
			cout << "Solved MHD at time " << simulation_time
				<< " s with time step " << dt << " s" << flush;
		}
		simulation_time += dt;

		const auto avg_div = pamhd::math::get_divergence_staggered(
			grid.local_cells(), grid,
			Face_B, Mag_div, SInfo
		);
		if (rank == 0) {
			cout << " average divergence " << avg_div << endl;
		}
		Cell::set_transfer_all(true, pamhd::Magnetic_Field_Divergence());
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false, pamhd::Magnetic_Field_Divergence());

		Cell::set_transfer_all(true,
			pamhd::Face_Magnetic_Field(),
			pamhd::mhd::MHD_State_Conservative(),
			pamhd::Bg_Magnetic_Field(),
			pamhd::mhd::Solver_Info()
		);
		apply_boundaries(
			grid, simulation_time, document,
			solar_wind_dir, solar_wind_cells,
			face_cells, edge_cells, vert_cells,
			planet_cells,
			options_sim.adiabatic_index,
			options_sim.vacuum_permeability,
			options_sim.proton_mass
		);
		grid.update_copies_of_remote_neighbors();
		pamhd::mhd::update_B_consistency(
			0, grid.local_cells(),
			Mas, Mom, Nrj, Mag, Face_B,
			SInfo, Substep,
			options_sim.adiabatic_index,
			options_sim.vacuum_permeability,
			true
		);
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false,
			pamhd::Face_Magnetic_Field(),
			pamhd::mhd::MHD_State_Conservative(),
			pamhd::Bg_Magnetic_Field(),
			pamhd::mhd::Solver_Info()
		);

		if (
			(options_mhd.save_n >= 0 and simulation_time >= time_end)
			or (options_mhd.save_n > 0 and simulation_time >= next_mhd_save)
		) {
			if (rank == 0) {
				cout << "Saving MHD at time " << simulation_time << endl;
			}
			if (next_mhd_save <= simulation_time) {
				next_mhd_save
					+= options_mhd.save_n
					* ceil(max(options_mhd.save_n, simulation_time - next_mhd_save) / options_mhd.save_n);
			}
			if (
				not pamhd::mhd::save_staggered(
					boost::filesystem::canonical(
						boost::filesystem::path(options_sim.output_directory)
					).append("mhd_staggered_").generic_string(),
					grid,
					file_version,
					simulation_step,
					simulation_time,
					options_sim.adiabatic_index,
					options_sim.proton_mass,
					options_sim.vacuum_permeability
				)
			) {
				cerr <<  __FILE__ << "(" << __LINE__ << "): "
					"Couldn't save mhd result."
					<< endl;
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
