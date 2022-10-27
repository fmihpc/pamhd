/*
MHD variables and cell class of PAMHD.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

* Neither the name of copyright holders nor the names of their contributors
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

#ifndef PAMHD_MHD_VARIABLES_HPP
#define PAMHD_MHD_VARIABLES_HPP

#include "Eigen/Core" // must be included before gensimcell.hpp
#include "gensimcell.hpp"

#include "../variables.hpp"


namespace pamhd {
namespace mhd {


/*
Variables used by magnetohydrodynamic (MHD) solver
*/

struct Mass_Density {
	using data_type = double;
	static const std::string get_name() { return {"mass density"}; }
	static const std::string get_option_name() { return {"mass-density"}; }
	static const std::string get_option_help() { return {"Plasma mass density (kg / m^3)"}; }
};

struct Momentum_Density {
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"momentum density"}; }
	static const std::string get_option_name() { return {"momentum-density"}; }
	static const std::string get_option_help() { return {"Plasma momentum density (kg / m^3 * m / s)"}; }
};

struct Total_Energy_Density {
	using data_type = double;
	static const std::string get_name() { return {"total energy density"}; }
	static const std::string get_option_name() { return {"total-energy-density"}; }
	static const std::string get_option_help() { return {"Plasma total energy density (J / m^3)"}; }
};

struct Velocity {
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"velocity"}; }
	static const std::string get_option_name() { return {"velocity"}; }
	static const std::string get_option_help() { return {"Plasma velocity (m / s)"}; }
};

// velocity for initial/boundary conditions of two population version of test program
struct Velocity2 {
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"velocity2"}; }
	static const std::string get_option_name() { return {"velocity2"}; }
	static const std::string get_option_help() { return {"Plasma velocity, second population (m / s)"}; }
};

struct Pressure {
	using data_type = double;
	static const std::string get_name() { return {"pressure"}; }
	static const std::string get_option_name() { return {"pressure"}; }
	static const std::string get_option_help() { return {"Plasma thermal pressure (Pa)"}; }
};

// pressure for initial/boundary conditions of two population version of test program
struct Pressure2 {
	using data_type = double;
	static const std::string get_name() { return {"pressure2"}; }
	static const std::string get_option_name() { return {"pressure2"}; }
	static const std::string get_option_help() { return {"Plasma thermal pressure, second population (Pa)"}; }
};

struct Number_Density {
	using data_type = double;
	static const std::string get_name() { return {"number density"}; }
	static const std::string get_option_name() { return {"number-density"}; }
	static const std::string get_option_help() { return {"Plasma number density (protons / m^3)"}; }
};

// density for initial/boundary conditions of two population version of test program
struct Number_Density2 {
	using data_type = double;
	static const std::string get_name() { return {"number density2"}; }
	static const std::string get_option_name() { return {"number-density2"}; }
	static const std::string get_option_help() { return {"Plasma number density, second population (protons / m^3)"}; }
};

/*!
Information for MHD solver on how to handle a simulation cell.

Dont solve cells are ignored by MHD solver. In other cells
fluxes are calculated normally but the flux of a variable
is only applied if the cell isn't a value or copy boundary
of that variable.
*/
struct Solver_Info {
	using data_type = unsigned int;

	/*!
	Corresponding bit of this variable is set when a cell is
	of a boundary type for a variable.
	*/
	static const unsigned int
		dont_solve = 1 << 0,
		mass_density_bdy = 1 << 1,
		velocity_bdy = 1 << 2,
		pressure_bdy = 1 << 3,
		magnetic_field_bdy = 1 << 4,
		mass_density2_bdy = 1 << 5,
		velocity2_bdy = 1 << 6,
		pressure2_bdy = 1 << 7;
};

//! Conservative hydrodynamic variables
using HD_Conservative = gensimcell::Cell<
	gensimcell::Always_Transfer,
	Mass_Density,
	Momentum_Density,
	Total_Energy_Density
>;

//! Primitive HD variables
using HD_Primitive = gensimcell::Cell<
	gensimcell::Never_Transfer,
	Mass_Density,
	Velocity,
	Pressure
>;

//! Current state of simulation cell
struct HD_State_Conservative {
	using data_type = HD_Conservative;
	static const std::string get_name() { return {"conservative HD variables"}; }
	static const std::string get_option_name() { return {"HD-conservative"}; }
	static const std::string get_option_help() {
		return {"Conservative HD variables"};
	}
};

//! Change of simulation state in cell for next step
struct HD_Flux_Conservative {
	using data_type = HD_Conservative;
	static const std::string get_name() { return {"conservative HD flux variables"}; }
	static const std::string get_option_name() { return {"HD-flux-conservative"}; }
	static const std::string get_option_help() {
		return {"Conservative HD flux variables"};
	}
};

// cell type for MHD test program
using Cell = gensimcell::Cell<
	gensimcell::Optional_Transfer,
	pamhd::mhd::HD_State_Conservative,
	pamhd::Electric_Current_Density,
	pamhd::mhd::Solver_Info,
	pamhd::MPI_Rank,
	pamhd::Resistivity,
	pamhd::Magnetic_Field,
	// TODO: move background fields to inner cell?
	pamhd::Bg_Magnetic_Field_Pos_X,
	pamhd::Bg_Magnetic_Field_Pos_Y,
	pamhd::Bg_Magnetic_Field_Pos_Z,
	pamhd::Magnetic_Field_Resistive,
	pamhd::Magnetic_Field_Temp,
	pamhd::Magnetic_Field_Divergence,
	pamhd::Scalar_Potential_Gradient,
	pamhd::mhd::HD_Flux_Conservative,
	pamhd::Magnetic_Field_Flux,
	pamhd::Face_Magnetic_Field,
	pamhd::Edge_Electric_Field
>;


struct HD2_State_Conservative { using data_type = HD_Conservative; };
struct HD2_Flux_Conservative { using data_type = HD_Conservative; };

// cell type for multipopulation MHD test program
using Cell2 = gensimcell::Cell<
	gensimcell::Optional_Transfer,
	pamhd::mhd::HD_State_Conservative,
	pamhd::mhd::HD2_State_Conservative,
	pamhd::Electric_Current_Density,
	pamhd::mhd::Solver_Info,
	pamhd::MPI_Rank,
	pamhd::Resistivity,
	pamhd::Magnetic_Field,
	pamhd::Bg_Magnetic_Field_Pos_X,
	pamhd::Bg_Magnetic_Field_Pos_Y,
	pamhd::Bg_Magnetic_Field_Pos_Z,
	pamhd::Magnetic_Field_Resistive,
	pamhd::Magnetic_Field_Temp,
	pamhd::Magnetic_Field_Divergence,
	pamhd::Scalar_Potential_Gradient,
	pamhd::mhd::HD_Flux_Conservative,
	pamhd::mhd::HD2_Flux_Conservative,
	pamhd::Magnetic_Field_Flux
>;


namespace detail {

// internal data type used by MHD solvers
using MHD = gensimcell::Cell<
	gensimcell::Never_Transfer,
	pamhd::mhd::Mass_Density,
	pamhd::mhd::Momentum_Density,
	pamhd::mhd::Total_Energy_Density,
	pamhd::Magnetic_Field
>;

}}} // namespaces

#endif // ifndef PAMHD_MHD_VARIABLES_HPP
