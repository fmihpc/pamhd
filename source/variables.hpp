/*
Variables of PAMHD common to several test programs.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2023, 2024, 2025 Finnish Meteorological Institute
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


Author(s): Ilja Honkonen
*/

#ifndef PAMHD_VARIABLES_HPP
#define PAMHD_VARIABLES_HPP


#include "Eigen/Core"

#include "grid/amr.hpp"


namespace pamhd {


/*
Variables used by magnetohydrodynamic (MHD) solver
*/

struct Magnetic_Field {
	static bool is_stale;
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"magnetic field"}; }
	static const std::string get_option_name() { return {"magnetic-field"}; }
	static const std::string get_option_help() { return {"Plasma magnetic field (T)"}; }
};

struct Magnetic_Field_Flux {
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"magnetic field flux"}; }
	static const std::string get_option_name() { return {"magnetic-field-flux"}; }
	static const std::string get_option_help() { return {"Flux of magnetic field (T)"}; }
};

//! stores B before divergence removal so B can be restored after failed removal
struct Magnetic_Field_Temp {
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"temporary magnetic field"}; }
	static const std::string get_option_name() { return {"temporary-magnetic-field"}; }
	static const std::string get_option_help() { return {"Temporary value of magnetic field in plasma (T)"}; }
};

//! stores change in B due to resistivity
struct Magnetic_Field_Resistive {
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"resistive magnetic field"}; }
	static const std::string get_option_name() { return {"resistive-magnetic-field"}; }
	static const std::string get_option_help() { return {"Change in magnetic field due to resistivity"}; }
};

//! Background magnetic field vector at cell faces
struct Bg_Magnetic_Field {
	static bool is_stale;
	using data_type = pamhd::grid::Face_Type<Eigen::Vector3d>;
	static const std::string get_name() { return {"face background magnetic fields"}; }
	static const std::string get_option_name() { return {"bg-b"}; }
	static const std::string get_option_help() { return {"background magnetic field vector on cell faces"}; }
};

/*! Magnetic field component at cell faces

Magnetic field stored on cell's faces.
Component of magnetic field normal to the face is stored for each face.
*/
struct Face_Magnetic_Field {
	static bool is_stale;
	using data_type = pamhd::grid::Face_Type<double>;
	static const std::string get_name() { return {"magnetic field on faces"}; }
	static const std::string get_option_name() { return {"face-b"}; }
	static const std::string get_option_help() { return {"magnetic field on cell faces normal to cell faces"}; }
};

struct Face_dB {
	using data_type = pamhd::grid::Face_Type<double>;
};

//! Electric field along each cell edge
struct Edge_Electric_Field {
	using data_type = pamhd::grid::Edge_Type<double>;
	static const std::string get_name() { return {"electric field on edges"}; }
	static const std::string get_option_name() { return {"edge-e"}; }
	static const std::string get_option_help() { return {"electric field on cell edges parallel to cell edges"}; }
};

struct MPI_Rank {
	using data_type = int;
	static const std::string get_name() { return {"MPI rank"}; }
	static const std::string get_option_name() { return {"mpi-rank"}; }
	static const std::string get_option_help() { return {"Owner (MPI process) of cell"}; }
};

struct Magnetic_Field_Divergence {
	using data_type = double;
	static const std::string get_name() { return {"magnetic field divergence"}; }
	static const std::string get_option_name() { return {"magnetic-field-divergence"}; }
	static const std::string get_option_help() { return {"Divergence of plasma magnetic field (T/m)"}; }
};

//! J in J = ∇×B
struct Electric_Current_Density {
	static bool is_stale;
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"current density"}; }
	static const std::string get_option_name() { return {"current-density"}; }
	static const std::string get_option_help() { return {"Density of electric current"}; }
};

//! Electrical resistivity of plasma (n in E = -VxB + nJ + ...)
struct Resistivity {
	using data_type = double;
	static const std::string get_name() { return {"electrical resistivity"}; }
	static const std::string get_option_name() { return {"resistivity"}; }
	static const std::string get_option_help() { return {"Electrical resistivity"}; }
};


/*!
Information for solver(s) on how to handle a simulation cell.

-1 means dont_solve cell that's neither read nor written
0 means boundary read-only cell
1 means normal read-write cell
*/
struct Solver_Info {
	static bool is_stale;
	using data_type = int;
};


} // namespace

#endif // ifndef PAMHD_VARIABLES_HPP
