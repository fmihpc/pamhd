/*
Particle variables and cell class of PAMHD.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2024, 2025 Finnish Meteorological Institute
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

#ifndef PAMHD_PARTICLE_VARIABLES_HPP
#define PAMHD_PARTICLE_VARIABLES_HPP


#include "vector"

#ifndef DONT_USE_MPI
#include "mpi.h" // must be included before gensimcell.hpp
#endif
#include "Eigen/Core" // must be included before gensimcell.hpp
#include "gensimcell.hpp"

#include "mhd/variables.hpp"
#include "particle/accumulation_variables.hpp"


namespace pamhd {
namespace particle {


/*
Variables used by particle solver
*/

struct Position {
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"position"}; }
	static const std::string get_option_name() { return {"position"}; }
	static const std::string get_option_help() { return {"Particle position"}; }
};

struct Velocity {
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"velocity"}; }
	static const std::string get_option_name() { return {"velocity"}; }
	static const std::string get_option_help() { return {"Particle velocity"}; }
};

//! Test particles have 0 mass.
struct Mass {
	using data_type = double;
	static const std::string get_name() { return {"mass"}; }
	static const std::string get_option_name() { return {"mass"}; }
	static const std::string get_option_help() { return {"Particle mass)"}; }
};

//! Mass of particle's species
struct Species_Mass {
	using data_type = double;
	static const std::string get_name() { return {"species mass"}; }
	static const std::string get_option_name() { return {"species-mass"}; }
	static const std::string get_option_help() { return {"Mass of particle's species)"}; }
};

//! Represents the ratio of particle mass and particle charge
struct Charge_Mass_Ratio {
	using data_type = double;
	static const std::string get_name() { return {"charge mass ratio"}; }
	static const std::string get_option_name() { return {"charge-mass-ratio"}; }
	static const std::string get_option_help() { return {"Particle charge to mass ratio"}; }
};


/*!
Template for a particle type that stores basic particle data
and given optional variables.
*/
template<
	class... Extra_Variables
> using Particle_T = gensimcell::Cell<
	gensimcell::Always_Transfer,
	Position,
	Velocity,
	Mass,
	Species_Mass,
	Charge_Mass_Ratio,
	Extra_Variables...
>;


//! Unique id of each particle
struct Particle_ID {
	using data_type = unsigned long long int;
};

//! A particle not moving between cells
using Particle_Internal = Particle_T<Particle_ID>;

struct Particles_Internal {
	using data_type = std::vector<Particle_Internal>;
};

//! Represents number of particles not moving between cells.
struct Nr_Particles_Internal {
	using data_type = unsigned long long int;
};


//! Represents destination cell of particles moving between processes
struct Destination_Cell {
	using data_type = unsigned long long int;
};

using Particle_External = Particle_T<Particle_ID, Destination_Cell>;

//! Represents particles moving between cells
struct Particles_External {
	static bool is_stale;
	using data_type = std::vector<Particle_External>;
};

//! Represents number of particles moving between cells.
struct Nr_Particles_External {
	static bool is_stale;
	using data_type = unsigned long long int;
};


//! of particles
struct Max_Spatial_Velocity {
	static bool is_stale;
	using data_type = double;
};

//! of particles
struct Max_Angular_Velocity {
	static bool is_stale;
	using data_type = double;
};

//! max linear acceleration of particles within cell
struct Max_Acceleration {
	using data_type = double;
};

struct Max_Gyrofrequency {
	using data_type = double;
};

struct Electric_Field {
	static bool is_stale;
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"electric field"}; }
	static const std::string get_option_name() { return {"electric-field"}; }
	static const std::string get_option_help() { return {"Cell-centered electric field"}; }
};

// number of macroparticles in simulation cell in initial & boundary conditions
struct Bdy_Nr_Particles_In_Cell {
	using data_type = double;
	static std::string get_name() { return std::string("Nr particles"); }
	static std::string get_option_name() { return std::string("nr-particles"); }
	static std::string get_option_help() { return std::string(""); }
};

// temperature of particles in initial & boundary conditions
struct Bdy_Temperature {
	using data_type = double;
	static std::string get_name() { return std::string("temperature"); }
	static std::string get_option_name() { return std::string("temperature"); }
	static std::string get_option_help() { return std::string(""); }
};

// number density of particles in initial & boundary conditions
struct Bdy_Number_Density {
	using data_type = double;
	static std::string get_name() { return std::string("number density"); }
	static std::string get_option_name() { return std::string("number-density"); }
	static std::string get_option_help() { return std::string(""); }
};

// velocity of particles in initial & boundary conditions
struct Bdy_Velocity {
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"velocity"}; }
	static const std::string get_option_name() { return {"velocity"}; }
	static const std::string get_option_help() { return {""}; }
};

// mass of particle species
struct Bdy_Species_Mass {
	using data_type = double;
	static const std::string get_name() { return {"species mass"}; }
	static const std::string get_option_name() { return {"species-mass"}; }
	static const std::string get_option_help() { return {""}; }
};

// charge to mass ratio of particles
struct Bdy_Charge_Mass_Ratio {
	using data_type = double;
	static const std::string get_name() { return {"charge mass ratio"}; }
	static const std::string get_option_name() { return {"charge-mass-ratio"}; }
	static const std::string get_option_help() { return {""}; }
};

struct Bulk_Mass {
	static bool is_stale;
	using data_type = double;
	static const std::string get_name() { return {"bulk mass"}; }
	static const std::string get_option_name() { return {"bulk-mass"}; }
	static const std::string get_option_help() { return {"Accumulated mass of particles"}; }
};

struct Bulk_Momentum {
	static bool is_stale;
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"bulk momentum"}; }
	static const std::string get_option_name() { return {"bulk-momentum"}; }
	static const std::string get_option_help() { return {"bulk mass * bulk velocity"}; }
};

struct Bulk_Velocity {
	//! second value used for tracking total weight of particles in cell
	static bool is_stale;
	using data_type = std::pair<Eigen::Vector3d, double>; // TODO: another cell type?
	static const std::string get_name() { return {"bulk velocity"}; }
	static const std::string get_option_name() { return {"bulk-velocity"}; }
	static const std::string get_option_help() { return {"Accumulated velocity of particles"}; }
};

struct Number_Of_Particles {
	static bool is_stale;
	using data_type = double;
	static const std::string get_name() { return {"number of species particles"}; }
	static const std::string get_option_name() { return {"number-of-species-particles"}; }
	static const std::string get_option_help() { return {"Accumulated number of particles (mass / species mass)"}; }
};

struct Bulk_Relative_Velocity2 {
	static bool is_stale;
	using data_type = double;
	static const std::string get_name() { return {"bulk relative velocity2"}; }
	static const std::string get_option_name() { return {"bulk-relative-velocity2"}; }
	static const std::string get_option_help() { return {"Accumulated square of particle velocity relative to bulk velocity"}; }
};

using Accumulated_To_Cell
	= Accumulated_To_Cell_T<
		Number_Of_Particles,
		Bulk_Mass,
		Bulk_Velocity,
		Bulk_Relative_Velocity2
	>;
using Accumulated_To_Cells = Accumulated_To_Cells_T<Accumulated_To_Cell>;


struct Current_Minus_Velocity {
	static bool is_stale;
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"J-V"}; }
	static const std::string get_option_name() { return {"J-V"}; }
	static const std::string get_option_help() { return {"Current minus Velocity for interpolating electric field to particle position"}; }
};

/*!
Cell type for test (massless) particle model of PAMHD.

Fixed size variables are first so that they
are saved at a known position in the file by dccrg.
gensimcell puts the variables in an MPI datatype
in the same order as they are given here, which
dccrg uses to save the file.
*/
using Cell_test_particle = gensimcell::Cell<
	gensimcell::Optional_Transfer,
	pamhd::Electric_Current_Density, // output compatible with regular model
	pamhd::Solver_Info,
	pamhd::MPI_Rank,
	pamhd::Magnetic_Field,
	pamhd::particle::Max_Spatial_Velocity,
	pamhd::particle::Max_Angular_Velocity,
	pamhd::particle::Electric_Field,
	pamhd::particle::Number_Of_Particles,
	pamhd::particle::Bdy_Number_Density,
	pamhd::particle::Bdy_Velocity,
	pamhd::particle::Bdy_Temperature,
	pamhd::particle::Bdy_Species_Mass,
	pamhd::particle::Bdy_Charge_Mass_Ratio,
	pamhd::particle::Bdy_Nr_Particles_In_Cell,
	pamhd::particle::Nr_Particles_Internal,
	pamhd::particle::Nr_Particles_External,
	pamhd::particle::Particles_Internal,
	pamhd::particle::Particles_External
>;

/*!
Cell type for particle model of PAMHD.

See Cell_test_particle for info.
*/
using Cell_hyb_particle = gensimcell::Cell<
	gensimcell::Optional_Transfer,
	pamhd::mhd::HD_State_Conservative,
	pamhd::Electric_Current_Density,
	pamhd::Solver_Info,
	pamhd::MPI_Rank,
	pamhd::Resistivity,
	pamhd::Magnetic_Field,
	pamhd::Bg_Magnetic_Field,
	pamhd::Magnetic_Field_Resistive,
	pamhd::Magnetic_Field_Temp,
	pamhd::Magnetic_Field_Divergence,
	pamhd::particle::Electric_Field,
	pamhd::particle::Number_Of_Particles,
	pamhd::particle::Bdy_Number_Density,
	pamhd::particle::Bdy_Velocity,
	pamhd::particle::Bdy_Temperature,
	pamhd::particle::Bdy_Species_Mass,
	pamhd::particle::Bdy_Charge_Mass_Ratio,
	pamhd::particle::Bdy_Nr_Particles_In_Cell,
	pamhd::particle::Bulk_Mass,
	pamhd::particle::Bulk_Momentum,
	pamhd::particle::Bulk_Velocity,
	pamhd::particle::Current_Minus_Velocity,
	pamhd::particle::Bulk_Relative_Velocity2,
	pamhd::particle::Nr_Particles_Internal,
	pamhd::particle::Nr_Particles_External,
	pamhd::particle::Nr_Accumulated_To_Cells,
	pamhd::particle::Particles_Internal,
	pamhd::particle::Particles_External,
	pamhd::particle::Accumulated_To_Cells,
	pamhd::mhd::HD_Flux_Conservative,
	pamhd::Magnetic_Field_Flux
>;

using Cell_hyb_particle_staggered = gensimcell::Cell<
	gensimcell::Optional_Transfer,
	pamhd::MPI_Rank,
	pamhd::Bg_Magnetic_Field,
	pamhd::Magnetic_Field_Divergence,
	pamhd::Face_Magnetic_Field,
	pamhd::Face_dB,
	pamhd::Electric_Current_Density,
	pamhd::grid::Target_Refinement_Level_Max,
	pamhd::grid::Target_Refinement_Level_Min,
	pamhd::Solver_Info,
	pamhd::mhd::Mass_Density,
	pamhd::mhd::Momentum_Density,
	pamhd::mhd::Total_Energy_Density,
	pamhd::Magnetic_Field,
	pamhd::mhd::Face_Boundary_Type,
	pamhd::mhd::Timestep,
	pamhd::mhd::Substepping_Period,
	pamhd::mhd::Substep_Min,
	pamhd::mhd::Substep_Max,
	pamhd::mhd::Max_Velocity,
	pamhd::mhd::MHD_Flux,
	pamhd::particle::Electric_Field,
	pamhd::particle::Number_Of_Particles,
	pamhd::particle::Max_Spatial_Velocity,
	pamhd::particle::Max_Angular_Velocity,
	pamhd::particle::Bdy_Number_Density,
	pamhd::particle::Bdy_Velocity,
	pamhd::particle::Bdy_Temperature,
	pamhd::particle::Bdy_Species_Mass,
	pamhd::particle::Bdy_Charge_Mass_Ratio,
	pamhd::particle::Bdy_Nr_Particles_In_Cell,
	pamhd::particle::Bulk_Mass,
	pamhd::particle::Bulk_Momentum,
	pamhd::particle::Bulk_Velocity,
	pamhd::particle::Current_Minus_Velocity,
	pamhd::particle::Bulk_Relative_Velocity2,
	pamhd::particle::Nr_Particles_Internal,
	pamhd::particle::Nr_Particles_External,
	pamhd::particle::Nr_Accumulated_To_Cells,
	pamhd::particle::Particles_Internal,
	pamhd::particle::Particles_External,
	pamhd::particle::Accumulated_To_Cells
>;

// Save this info in every neighbor iterator of every cell.
// Updated automatically by dccrg.
struct Is_Local {
	bool is_local = false;
	template<
		class Grid, class Cell_Item, class Neighbor_Item
	> void update(
		const Grid& grid, const Cell_Item&, const Neighbor_Item& neighbor, const int&, const Is_Local&
	) {
		is_local = grid.is_local(neighbor.id);
	}
};


}} // namespaces

#endif // ifndef PAMHD_PARTICLE_VARIABLES_HPP
