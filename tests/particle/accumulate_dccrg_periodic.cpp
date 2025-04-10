/*
Tests particle accumulator of PAMHD with periodic grid.

Copyright 2015, 2016, 2017 Ilja Honkonen
Copyright 2019, 2025 Finnish Meteorological Institute
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
*/

#include "algorithm"
#include "cmath"
#include "cstdlib"
#include "exception"
#include "iostream"
#include "random"
#include "utility"
#include "vector"

#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "Eigen/Core" // must be included before gensimcell.hpp
#include "mpi.h" // must be included before gensimcell.hpp
#include "gensimcell.hpp"

#include "particle/accumulate_dccrg.hpp"
#include "particle/accumulation_variables.hpp"
#include "particle/variables.hpp"
#include "variable_getter.hpp"


double mass_density()
{
	return 3.0;
}

struct Mass_Density {
	using data_type = double;
};

using Accumulated_To_Cell = pamhd::particle::Accumulated_To_Cell_T<Mass_Density>;
using Accumulated_To_Cells = pamhd::particle::Accumulated_To_Cells_T<Accumulated_To_Cell>;

using Cell = gensimcell::Cell<
	gensimcell::Optional_Transfer,
	pamhd::particle::Particles_Internal,
	Mass_Density,
	pamhd::Solver_Info,
	// only following two are transferred between processes
	pamhd::particle::Nr_Accumulated_To_Cells,
	Accumulated_To_Cells
>;
using Grid = dccrg::Dccrg<
	Cell,
	dccrg::Cartesian_Geometry,
	std::tuple<>,
	std::tuple<pamhd::particle::Is_Local>
>;


/*
Creates evenly spaced particles in given cells.
*/
void create_particles(
	const size_t values_per_cell,
	Grid& grid
) {
	const pamhd::particle::Velocity Vel{};
	const pamhd::particle::Position Pos{};

	std::mt19937_64 random_source;
	for (const auto& cell: grid.local_cells()) {
		random_source.seed(cell.id);

		const auto
			cell_min = grid.geometry.get_min(cell.id),
			cell_max = grid.geometry.get_max(cell.id),
			cell_length = grid.geometry.get_length(cell.id);
		const auto volume = cell_length[0] * cell_length[1] * cell_length[2];

		std::uniform_real_distribution<>
			position_generator_x(cell_min[0], cell_max[0]),
			position_generator_y(cell_min[1], cell_max[1]),
			position_generator_z(cell_min[2], cell_max[2]);

		for (size_t i = 0; i < values_per_cell; i++) {
			pamhd::particle::Particle_Internal new_particle;
			new_particle[Vel][0] =
			new_particle[Vel][1] =
			new_particle[Vel][2] =
			new_particle[pamhd::particle::Charge_Mass_Ratio()] = 0;

			new_particle[Pos][0] = position_generator_x(random_source);
			new_particle[Pos][1] = position_generator_y(random_source);
			new_particle[Pos][2] = position_generator_z(random_source);

			new_particle[pamhd::particle::Mass()]
				= mass_density() * volume / values_per_cell;

			(*cell.data)[pamhd::particle::Particles_Internal()].push_back(new_particle);
		}
	}
}


// returns a reference to data accumulated from particles in given cell
const auto bulk_value_getter = pamhd::Variable_Getter<Mass_Density>();

// returns a reference to data accumulated from particles in accumulation list
const auto list_bulk_value_getter
	= [](Accumulated_To_Cell& accu_item)->Mass_Density::data_type&{
		return accu_item[Mass_Density()];
	};

// returns a reference to target of accumulated values in accumulation list
const auto list_target_getter
	= [](Accumulated_To_Cell& accu_item)->pamhd::particle::Target::data_type&{
		return accu_item[pamhd::particle::Target()];
	};

/*
Returns a reference to given cell's list of data accumulated
from particles in given cell into other cells
*/
const auto accumulation_list_getter
	= [](Cell& cell_data)->Accumulated_To_Cells::data_type&{
		return cell_data[Accumulated_To_Cells()];
	};

// returns a reference to length of list above incoming from other processes
const auto accumulation_list_length_getter
	= [](Cell& cell_data)->pamhd::particle::Nr_Accumulated_To_Cells::data_type&{
		return cell_data[pamhd::particle::Nr_Accumulated_To_Cells()];
	};

const auto solver_info_getter = pamhd::Variable_Getter<pamhd::Solver_Info>();

const auto accumulate_from_remote_neighbors
	= [](Grid& grid){
		pamhd::particle::accumulate_from_remote_neighbors(
			grid,
			bulk_value_getter,
			list_bulk_value_getter,
			list_target_getter,
			accumulation_list_getter,
			solver_info_getter
		);
	};

const auto allocate_accumulation_lists
	= [](Grid& grid){
		pamhd::particle::allocate_accumulation_lists(
			grid,
			accumulation_list_getter,
			accumulation_list_length_getter
		);
	};

// returns infinite norm between analytic results and that in given cells
double get_norm(
	const Grid& grid,
	MPI_Comm& comm
) {
	double
		density_local = 0,
		density_global = 0;
	for (const auto& cell: grid.local_cells()) {
		density_local += (*cell.data)[Mass_Density()];
	}
	if (
		MPI_Allreduce(
			&density_local,
			&density_global,
			1,
			MPI_DOUBLE,
			MPI_SUM,
			comm
		) != MPI_SUCCESS
	) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Couldn't set reduce norm."
			<< std::endl;
		abort();
	}

	size_t
		cells_local = std::distance(grid.local_cells().begin(), grid.local_cells().end()),
		cells_global = 0;
	if (
		MPI_Allreduce(
			&cells_local,
			&cells_global,
			1,
			MPI_UINT64_T,
			MPI_SUM,
			comm
		) != MPI_SUCCESS
	) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< ": Couldn't set reduce norm."
			<< std::endl;
		abort();
	}

	density_global /= cells_global;

	return std::fabs(density_global - mass_density());
}


int main(int argc, char* argv[])
{
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

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
		std::cerr << "Zoltan_Initialize failed." << std::endl;
		abort();
	}


	constexpr size_t
		nr_of_values = 100,
		max_nr_of_cells = 512;

	const unsigned int neighborhood_size = 1;
	const std::array<bool, 3> periodic = {{true, true, true}};

	dccrg::Cartesian_Geometry::Parameters geom_params;
	geom_params.start = {{-3,  -5,  -7}};
	geom_params.level_0_cell_length = {{7,   5,   3}};

	for (size_t nr_cells = 1; nr_cells <= max_nr_of_cells; nr_cells *= 2) {
		const std::array<uint64_t, 3> nr_of_cells{{nr_cells, 1, 1}};
		Grid grid; grid
			.set_neighborhood_length(neighborhood_size)
			.set_maximum_refinement_level(0)
			.set_load_balancing_method("RANDOM")
			.set_periodic(periodic[0], periodic[1], periodic[2])
			.set_initial_length(nr_of_cells)
			.initialize(comm)
			.balance_load();

		grid.set_geometry(geom_params);

		grid.update_copies_of_remote_neighbors();

		create_particles(nr_of_values, grid);

		for (const auto& cell: grid.local_cells()) {
			bulk_value_getter.data(*cell.data) = 0;
			solver_info_getter.data(*cell.data) = 0;
		}
		pamhd::particle::accumulate(
			grid.local_cells(),
			grid,
			[](Cell& cell)->pamhd::particle::Particles_Internal::data_type&{
				return cell[pamhd::particle::Particles_Internal()];
			},
			[](pamhd::particle::Particle_Internal& particle)
				->pamhd::particle::Position::data_type&
			{
				return particle[pamhd::particle::Position()];
			},
			[](Cell&, pamhd::particle::Particle_Internal& particle)
				->pamhd::particle::Mass::data_type&
			{
				return particle[pamhd::particle::Mass()];
			},
			bulk_value_getter,
			list_bulk_value_getter,
			list_target_getter,
			accumulation_list_length_getter,
			accumulation_list_getter,
			solver_info_getter
		);

		Cell::set_transfer_all(true, pamhd::particle::Nr_Accumulated_To_Cells());
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false, pamhd::particle::Nr_Accumulated_To_Cells());

		allocate_accumulation_lists(grid);

		// transfer accumulated values between processes
		Cell::set_transfer_all(true, Accumulated_To_Cells());
		grid.update_copies_of_remote_neighbors();
		Cell::set_transfer_all(false, Accumulated_To_Cells());

		accumulate_from_remote_neighbors(grid);

		// transform accumulated mass to mass density
		for (const auto& cell: grid.local_cells()) {
			const auto cell_length = grid.geometry.get_length(cell.id);
			const auto volume = cell_length[0] * cell_length[1] * cell_length[2];
			(*cell.data)[Mass_Density()] /= volume;
		}

		const double norm = get_norm(grid, comm);

		if (norm > 1e-10) {
			if (rank == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< ": Norm is too large: " << norm
					<< " with " << nr_cells << " cell(s)"
					<< std::endl;
			}
			MPI_Finalize();
			return EXIT_FAILURE;
		}
	}

	MPI_Finalize();
	return EXIT_SUCCESS;
}
