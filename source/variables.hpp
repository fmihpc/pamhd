/*
Variables of PAMHD common to several test programs.

Copyright 2014, 2015, 2016, 2017 Ilja Honkonen
Copyright 2023 Finnish Meteorological Institute
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

#ifndef PAMHD_VARIABLES_HPP
#define PAMHD_VARIABLES_HPP


#include "Eigen/Core"


namespace pamhd {


/*
Variables used by magnetohydrodynamic (MHD) solver
*/

struct Magnetic_Field {
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
	struct Bg_B_type {
		std::array<Eigen::Vector3d, 6> bg_b;

		/*
		dim_i = dimension of face, 0 = x, 1 = y, 2 = z
		side_i = side of face w.r.t. cell center, 0 = negative side
		*/
		const Eigen::Vector3d& operator()(
			const size_t dim_i,
			const size_t side_i
		) const {
			if (dim_i > 2) {
				throw std::domain_error("Dimension > 2");
			}
			if (side_i > 1) {
				throw std::domain_error("Side > 1");
			}
			return this->bg_b[dim_i*2 + side_i];
		}

		// https://stackoverflow.com/a/123995
		Eigen::Vector3d& operator()(
			const size_t dim_i,
			const size_t side_i
		) {
			return const_cast<Eigen::Vector3d&>(static_cast<const Bg_B_type&>(*this).operator()(dim_i, side_i));
		}

		#ifdef MPI_VERSION
		std::tuple<void*, int, MPI_Datatype> get_mpi_datatype() const {
			return std::make_tuple(
				(void*) this->bg_b.data(),
				this->bg_b.size(),
				MPI_DOUBLE
			);
		}
		#endif
	};
	using data_type = Bg_B_type;
	static const std::string get_name() { return {"face background magnetic fields"}; }
	static const std::string get_option_name() { return {"bg-b"}; }
	static const std::string get_option_help() { return {"background magnetic field vector on cell faces"}; }
};

/*! Magnetic field component at cell faces

Magnetic field stored on cell's faces on positive side from cell's center.
Component of magnetic field normal to the face is stored for each face,
e.g. Bx is located on the yz face, By on xz face and Bz on xy face.
*/
struct Face_Magnetic_Field {
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"magnetic field on positive faces"}; }
	static const std::string get_option_name() { return {"face-b"}; }
	static const std::string get_option_help() { return {"magnetic field normal component at cell faces on positive sides of cell"}; }
};

//! Magnetic field on negative cell faces
struct Face_Magnetic_Field_Neg {
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"magnetic field on negative faces"}; }
	static const std::string get_option_name() { return {"face-b-neg"}; }
	static const std::string get_option_help() { return {"magnetic field normal component at cell faces on negative sides of cell"}; }
};

//! Electric field along each cell edge
struct Edge_Electric_Field {
	struct Edge_E_type {
		std::array<double, 12> edge_e;

		/*
		par_dim_i = dimension to which electric field is parallel to,
		0 = x, 1 = y, 2 = z.
		first_perp_dim_i = negative or positive side of cell in
		lexically earlier dimension perpendicular to electric field,
		e.g. if par_dim_i = 0, first_perp_dim_i = 1 means positive side of
		cell in y dimension.
		second_perp_dim = neg or pos side in lexically later dimension,
		if par_dim_i = 2, second_perp_dim_i = 0 means negative side of
		cell in y dimension.
		par_dim_i | first_perp_dim_i | second..._i | E on edge of cell
		    0     |         0        |      0      | x directed: -y, -z
		    0     |         1        |      1      | x dir:      +y, +z
		...
		    2     |         1        |      1      | z dir:      +x, +y
		*/
		const double& operator()(
			const size_t par_dim_i,
			const size_t first_perp_dim_i,
			const size_t second_perp_dim_i
		) const {
			if (par_dim_i > 2) {
				throw std::domain_error("Parallel dimension > 2");
			}
			if (first_perp_dim_i > 1) {
				throw std::domain_error("First perpendicular dimension > 1");
			}
			if (second_perp_dim_i > 1) {
				throw std::domain_error("Second perpendicular dimension > 1");
			}
			return this->edge_e[par_dim_i*2*2 + first_perp_dim_i*2 + second_perp_dim_i];
		}

		// https://stackoverflow.com/a/123995
		double& operator()(
			const size_t par_dim_i,
			const size_t first_perp_dim_i,
			const size_t second_perp_dim_i
		) {
			return const_cast<double&>(static_cast<const Edge_E_type&>(*this).operator()(par_dim_i, first_perp_dim_i, second_perp_dim_i));
		}

		#ifdef MPI_VERSION
		std::tuple<void*, int, MPI_Datatype> get_mpi_datatype() const {
			return std::make_tuple(
				(void*) this->edge_e.data(),
				this->edge_e.size(),
				MPI_DOUBLE
			);
		}
		#endif
	};
	using data_type = Edge_E_type;
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

struct Scalar_Potential_Gradient {
	using data_type = Eigen::Vector3d;
	static const std::string get_name() { return {"scalar potential gradient"}; }
	static const std::string get_option_name() { return {"scalar-potential-gradient"}; }
	static const std::string get_option_help() { return {"Gradient of scalar potential from Poisson's equation"}; }
};

//! J in J = ∇×B
struct Electric_Current_Density {
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

} // namespace

#endif // ifndef PAMHD_VARIABLES_HPP
