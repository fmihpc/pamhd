/*
Roe MHD solver of PAMHD adapted from Athena (https://trac.princeton.edu/Athena/wiki)
Copyright 2003 James M. Stone
Copyright 2003 Thomas A. Gardiner
Copyright 2003 Peter J. Teuben
Copyright 2003 John F. Hawley
Copyright 2014, 2015, 2016, 2017 Ilja Honkonen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PAMHD_MHD_ROE_ATHENA_HPP
#define PAMHD_MHD_ROE_ATHENA_HPP


#include "algorithm"
#include "array"
#include "cmath"
#include "exception"
#include "limits"
#include "string"
#include "tuple"
#include "utility"

#include "gensimcell.hpp"

#include "mhd/common.hpp"
#include "mhd/variables.hpp"


namespace pamhd {
namespace mhd {
namespace athena {

/*!
\brief ADIABATIC MHD

- Input: d,v1,v2,v3,h,b1,b2,b3=Roe averaged density, velocities, enthalpy, B
         x,y = numerical factors (see eqn XX)
- Output: eigenvalues[7], right_eigenmatrix[7,7], left_eigenmatrix[7,7];
*/
void esys_roe_adb_mhd(const double d, const double v1, const double v2, const double v3,
  const double h, const double b1, const double b2, const double b3, 
  const double x, const double y,
  const double adiabatic_index,
  std::array<double, 7>& eigenvalues,
  std::array<std::array<double, 7>, 7>& right_eigenmatrix,
  std::array<std::array<double, 7>, 7>& left_eigenmatrix)
{
  constexpr double TINY_NUMBER = std::numeric_limits<double>::epsilon();
  const double
  	Gamma_1 = adiabatic_index - 1,
  	Gamma_2 = adiabatic_index - 2;
  double di,vsq,btsq,bt_starsq,vaxsq,hp,twid_asq,cfsq,cf,cssq,cs;
  double bt,bt_star,bet2,bet3,bet2_star,bet3_star,bet_starsq,vbet,alpha_f,alpha_s;
  double isqrtd,sqrtd,s,twid_a,qf,qs,af_prime,as_prime,afpbb,aspbb,vax;
  double norm,cff,css,af,as,afpb,aspb,q2_star,q3_star,vqstr;
  double ct2,tsum,tdif,cf2_cs2;
  double qa,qb,qc,qd;
  di = 1.0/d;
  vsq = v1*v1 + v2*v2 + v3*v3;
  btsq = b2*b2 + b3*b3;
  bt_starsq = (Gamma_1 - Gamma_2*y)*btsq;
  vaxsq = b1*b1*di;
  hp = h - (vaxsq + btsq*di);
  twid_asq = std::max((Gamma_1*(hp-0.5*vsq)-Gamma_2*x), TINY_NUMBER);

/* Compute fast- and slow-magnetosonic speeds (eq. B18) */

  ct2 = bt_starsq*di;
  tsum = vaxsq + ct2 + twid_asq;
  tdif = vaxsq + ct2 - twid_asq;
  cf2_cs2 = std::sqrt(tdif*tdif + 4.0*twid_asq*ct2);

  cfsq = 0.5*(tsum + cf2_cs2);
  cf = std::sqrt(cfsq);

  cssq = twid_asq*vaxsq/cfsq;
  cs = std::sqrt(cssq);

/* Compute beta(s) (eqs. A17, B20, B28) */

  bt = std::sqrt(btsq);
  bt_star = std::sqrt(bt_starsq);
  if (bt == 0.0) {
    bet2 = 1.0;
    bet3 = 0.0;
  } else {
    bet2 = b2/bt;
    bet3 = b3/bt;
  }
  bet2_star = bet2/std::sqrt(Gamma_1 - Gamma_2*y);
  bet3_star = bet3/std::sqrt(Gamma_1 - Gamma_2*y);
  bet_starsq = bet2_star*bet2_star + bet3_star*bet3_star;
  vbet = v2*bet2_star + v3*bet3_star;

/* Compute alpha(s) (eq. A16) */

  if ((cfsq-cssq) == 0.0) {
    alpha_f = 1.0;
    alpha_s = 0.0;
  } else if ( (twid_asq - cssq) <= 0.0) {
    alpha_f = 0.0;
    alpha_s = 1.0;
  } else if ( (cfsq - twid_asq) <= 0.0) {
    alpha_f = 1.0;
    alpha_s = 0.0;
  } else {
    alpha_f = std::sqrt((twid_asq - cssq)/(cfsq - cssq));
    alpha_s = std::sqrt((cfsq - twid_asq)/(cfsq - cssq));
  }

/* Compute Q(s) and A(s) (eq. A14-15), etc. */

  sqrtd = std::sqrt(d);
  isqrtd = 1.0/sqrtd;
  s = (b1 > 0) ? 1 : -1;
  twid_a = std::sqrt(twid_asq);
  qf = cf*alpha_f*s;
  qs = cs*alpha_s*s;
  af_prime = twid_a*alpha_f*isqrtd;
  as_prime = twid_a*alpha_s*isqrtd;
  afpbb = af_prime*bt_star*bet_starsq;
  aspbb = as_prime*bt_star*bet_starsq;

/* Compute eigenvalues (eq. B17) */

  vax = std::sqrt(vaxsq);
  eigenvalues[0] = v1 - cf;
  eigenvalues[1] = v1 - vax;
  eigenvalues[2] = v1 - cs;
  eigenvalues[3] = v1;
  eigenvalues[4] = v1 + cs;
  eigenvalues[5] = v1 + vax;
  eigenvalues[6] = v1 + cf;

/* Right-eigenvectors, stored as COLUMNS (eq. B21) */
/* Note statements are grouped in ROWS for optimization, even though rem[*][n]
 * is the nth right eigenvector */

  right_eigenmatrix[0][0] = alpha_f;
  right_eigenmatrix[0][1] = 0.0;
  right_eigenmatrix[0][2] = alpha_s;
  right_eigenmatrix[0][3] = 1.0;
  right_eigenmatrix[0][4] = alpha_s;
  right_eigenmatrix[0][5] = 0.0;
  right_eigenmatrix[0][6] = alpha_f;

  right_eigenmatrix[1][0] = alpha_f*eigenvalues[0];
  right_eigenmatrix[1][1] = 0.0;
  right_eigenmatrix[1][2] = alpha_s*eigenvalues[2];
  right_eigenmatrix[1][3] = v1;
  right_eigenmatrix[1][4] = alpha_s*eigenvalues[4];
  right_eigenmatrix[1][5] = 0.0;
  right_eigenmatrix[1][6] = alpha_f*eigenvalues[6];

  qa = alpha_f*v2;
  qb = alpha_s*v2;
  qc = qs*bet2_star;
  qd = qf*bet2_star;
  right_eigenmatrix[2][0] = qa + qc;
  right_eigenmatrix[2][1] = -bet3;
  right_eigenmatrix[2][2] = qb - qd;
  right_eigenmatrix[2][3] = v2;
  right_eigenmatrix[2][4] = qb + qd;
  right_eigenmatrix[2][5] = bet3;
  right_eigenmatrix[2][6] = qa - qc;

  qa = alpha_f*v3;
  qb = alpha_s*v3;
  qc = qs*bet3_star;
  qd = qf*bet3_star;
  right_eigenmatrix[3][0] = qa + qc;
  right_eigenmatrix[3][1] = bet2;
  right_eigenmatrix[3][2] = qb - qd;
  right_eigenmatrix[3][3] = v3;
  right_eigenmatrix[3][4] = qb + qd;
  right_eigenmatrix[3][5] = -bet2;
  right_eigenmatrix[3][6] = qa - qc;

  right_eigenmatrix[4][0] = alpha_f*(hp - v1*cf) + qs*vbet + aspbb;
  right_eigenmatrix[4][1] = -(v2*bet3 - v3*bet2);
  right_eigenmatrix[4][2] = alpha_s*(hp - v1*cs) - qf*vbet - afpbb;
  right_eigenmatrix[4][3] = 0.5*vsq + Gamma_2*x/Gamma_1;
  right_eigenmatrix[4][4] = alpha_s*(hp + v1*cs) + qf*vbet - afpbb;
  right_eigenmatrix[4][5] = -right_eigenmatrix[4][1];
  right_eigenmatrix[4][6] = alpha_f*(hp + v1*cf) - qs*vbet + aspbb;

  right_eigenmatrix[5][0] = as_prime*bet2_star;
  right_eigenmatrix[5][1] = -bet3*s*isqrtd;
  right_eigenmatrix[5][2] = -af_prime*bet2_star;
  right_eigenmatrix[5][3] = 0.0;
  right_eigenmatrix[5][4] = right_eigenmatrix[5][2];
  right_eigenmatrix[5][5] = right_eigenmatrix[5][1];
  right_eigenmatrix[5][6] = right_eigenmatrix[5][0];

  right_eigenmatrix[6][0] = as_prime*bet3_star;
  right_eigenmatrix[6][1] = bet2*s*isqrtd;
  right_eigenmatrix[6][2] = -af_prime*bet3_star;
  right_eigenmatrix[6][3] = 0.0;
  right_eigenmatrix[6][4] = right_eigenmatrix[6][2];
  right_eigenmatrix[6][5] = right_eigenmatrix[6][1];
  right_eigenmatrix[6][6] = right_eigenmatrix[6][0];

/* Left-eigenvectors, stored as ROWS (eq. B29) */

/* Normalize by 1/2a^{2}: quantities denoted by \hat{f} */
  norm = 0.5/twid_asq;
  cff = norm*alpha_f*cf;
  css = norm*alpha_s*cs;
  qf *= norm;
  qs *= norm;
  af = norm*af_prime*d;
  as = norm*as_prime*d;
  afpb = norm*af_prime*bt_star;
  aspb = norm*as_prime*bt_star;

/* Normalize by (gamma-1)/2a^{2}: quantities denoted by \bar{f} */
  norm *= Gamma_1;
  alpha_f *= norm;
  alpha_s *= norm;
  q2_star = bet2_star/bet_starsq;
  q3_star = bet3_star/bet_starsq;
  vqstr = (v2*q2_star + v3*q3_star);
  norm *= 2.0;

  left_eigenmatrix[0][0] = alpha_f*(vsq-hp) + cff*(cf+v1) - qs*vqstr - aspb;
  left_eigenmatrix[0][1] = -alpha_f*v1 - cff;
  left_eigenmatrix[0][2] = -alpha_f*v2 + qs*q2_star;
  left_eigenmatrix[0][3] = -alpha_f*v3 + qs*q3_star;
  left_eigenmatrix[0][4] = alpha_f;
  left_eigenmatrix[0][5] = as*q2_star - alpha_f*b2;
  left_eigenmatrix[0][6] = as*q3_star - alpha_f*b3;

  left_eigenmatrix[1][0] = 0.5*(v2*bet3 - v3*bet2);
  left_eigenmatrix[1][1] = 0.0;
  left_eigenmatrix[1][2] = -0.5*bet3;
  left_eigenmatrix[1][3] = 0.5*bet2;
  left_eigenmatrix[1][4] = 0.0;
  left_eigenmatrix[1][5] = -0.5*sqrtd*bet3*s;
  left_eigenmatrix[1][6] = 0.5*sqrtd*bet2*s;

  left_eigenmatrix[2][0] = alpha_s*(vsq-hp) + css*(cs+v1) + qf*vqstr + afpb;
  left_eigenmatrix[2][1] = -alpha_s*v1 - css;
  left_eigenmatrix[2][2] = -alpha_s*v2 - qf*q2_star;
  left_eigenmatrix[2][3] = -alpha_s*v3 - qf*q3_star;
  left_eigenmatrix[2][4] = alpha_s;
  left_eigenmatrix[2][5] = -af*q2_star - alpha_s*b2;
  left_eigenmatrix[2][6] = -af*q3_star - alpha_s*b3;

  left_eigenmatrix[3][0] = 1.0 - norm*(0.5*vsq - Gamma_2*x/Gamma_1); 
  left_eigenmatrix[3][1] = norm*v1;
  left_eigenmatrix[3][2] = norm*v2;
  left_eigenmatrix[3][3] = norm*v3;
  left_eigenmatrix[3][4] = -norm;
  left_eigenmatrix[3][5] = norm*b2;
  left_eigenmatrix[3][6] = norm*b3;

  left_eigenmatrix[4][0] = alpha_s*(vsq-hp) + css*(cs-v1) - qf*vqstr + afpb;
  left_eigenmatrix[4][1] = -alpha_s*v1 + css;
  left_eigenmatrix[4][2] = -alpha_s*v2 + qf*q2_star;
  left_eigenmatrix[4][3] = -alpha_s*v3 + qf*q3_star;
  left_eigenmatrix[4][4] = alpha_s;
  left_eigenmatrix[4][5] = left_eigenmatrix[2][5];
  left_eigenmatrix[4][6] = left_eigenmatrix[2][6];

  left_eigenmatrix[5][0] = -left_eigenmatrix[1][0];
  left_eigenmatrix[5][1] = 0.0;
  left_eigenmatrix[5][2] = -left_eigenmatrix[1][2];
  left_eigenmatrix[5][3] = -left_eigenmatrix[1][3];
  left_eigenmatrix[5][4] = 0.0;
  left_eigenmatrix[5][5] = left_eigenmatrix[1][5];
  left_eigenmatrix[5][6] = left_eigenmatrix[1][6];

  left_eigenmatrix[6][0] = alpha_f*(vsq-hp) + cff*(cf-v1) + qs*vqstr - aspb;
  left_eigenmatrix[6][1] = -alpha_f*v1 + cff;
  left_eigenmatrix[6][2] = -alpha_f*v2 - qs*q2_star;
  left_eigenmatrix[6][3] = -alpha_f*v3 - qs*q3_star;
  left_eigenmatrix[6][4] = alpha_f;
  left_eigenmatrix[6][5] = left_eigenmatrix[0][5];
  left_eigenmatrix[6][6] = left_eigenmatrix[0][6];
}




/*!
 *  \brief Conserved variables in 1D (does not contain Bx).
 *  IMPORTANT!! The order of the elements in Cons1DS CANNOT be changed.
 */
struct Cons1DS {
  double d;
  double Mx;
  double My;
  double Mz;
  double E;
  double By;
  double Bz;
};

/*!
 *  \brief Primitive variables in 1D (does not contain Bx).
 *  IMPORTANT!! The order of the elements in Prim1DS CANNOT be changed.
 */
struct Prim1DS{
  double d;
  double Vx;
  double Vy;
  double Vz;			/*!< velocity in Z-direction */
  double P;			/*!< pressure */
  double By;			/*!< cell centered magnetic fields in Y-dir */
  double Bz;			/*!< cell centered magnetic fields in Z-dir */
};



/*----------------------------------------------------------------------------*/
/*! \fn void fluxes(const Cons1DS Ul, const Cons1DS Ur,
 *            const Prim1DS Wl, const Prim1DS Wr,
 *            const double Bxi, ...)
 *  \brief Computes 1D fluxes
 *   Input Arguments:
 *   - Bxi = B normal to cell interface
 *   - Ul,Ur = L/R-states of CONSERVED variables at cell interface
 *   Returns flux contribution from negative and positive state respectively
 */
Cons1DS athena_roe_fluxes(
	const Cons1DS& Ul,
	const Cons1DS& Ur,
	const Prim1DS& Wl,
	const Prim1DS& Wr,
	const double Bxi,
	const double adiabatic_index
) {
  constexpr int NWAVE = 7;

  constexpr double etah = 0.0;

  double sqrtdl,sqrtdr,isdlpdr,droe,v1roe,v2roe,v3roe,pbl=0.0,pbr=0.0;
  double hroe;
  double b2roe,b3roe,x,y;
  std::array<double, NWAVE> coeff, ev, dU, a, u_inter;
  std::array<std::array<double, NWAVE>, NWAVE> rem,lem;
  double p_inter=0.0;

  for (int n=0; n<NWAVE; n++) {
    for (int m=0; m<NWAVE; m++) {
      rem[n][m] = 0.0;
      lem[n][m] = 0.0;
    }
  }

/*--- Step 2. ------------------------------------------------------------------
 * Compute Roe-averaged data from left- and right-states
 */

  sqrtdl = std::sqrt(Wl.d);
  sqrtdr = std::sqrt(Wr.d);
  isdlpdr = 1.0/(sqrtdl + sqrtdr);

  droe  = sqrtdl*sqrtdr;
  v1roe = (sqrtdl*Wl.Vx + sqrtdr*Wr.Vx)*isdlpdr;
  v2roe = (sqrtdl*Wl.Vy + sqrtdr*Wr.Vy)*isdlpdr;
  v3roe = (sqrtdl*Wl.Vz + sqrtdr*Wr.Vz)*isdlpdr;

/* The Roe average of the magnetic field is defined differently  */

  b2roe = (sqrtdr*Wl.By + sqrtdl*Wr.By)*isdlpdr;
  b3roe = (sqrtdr*Wl.Bz + sqrtdl*Wr.Bz)*isdlpdr;
  x = 0.5*(std::pow(Wl.By - Wr.By, 2) + std::pow(Wl.Bz - Wr.Bz, 2))/std::pow(sqrtdl + sqrtdr, 2);
  y = 0.5*(Wl.d + Wr.d)/droe;
  pbl = 0.5*(std::pow(Bxi, 2) + std::pow(Wl.By, 2) + std::pow(Wl.Bz, 2));
  pbr = 0.5*(std::pow(Bxi, 2) + std::pow(Wr.By, 2) + std::pow(Wr.Bz, 2));

/*
 * Following Roe(1981), the enthalpy H=(E+P)/d is averaged for adiabatic flows,
 * rather than E or P directly.  sqrtdl*hl = sqrtdl*(el+pl)/dl = (el+pl)/sqrtdl
 */

  hroe  = ((Ul.E + Wl.P + pbl)/sqrtdl + (Ur.E + Wr.P + pbr)/sqrtdr)*isdlpdr;

/*--- Step 3. ------------------------------------------------------------------
 * Compute eigenvalues and eigenmatrices using Roe-averaged values
 */

  esys_roe_adb_mhd(
    droe,
    v1roe,
    v2roe,
    v3roe,
    hroe,
    Bxi,
    b2roe,
    b3roe,
    x,
    y,
    adiabatic_index,
    ev,
    rem,
    lem
  );

/*--- Step 4. ------------------------------------------------------------------
 * Compute L/R fluxes 
 */

  Cons1DS Fl,Fr;

  Fl.d  = Ul.Mx;
  Fr.d  = Ur.Mx;

  Fl.Mx = Ul.Mx*Wl.Vx;
  Fr.Mx = Ur.Mx*Wr.Vx;

  Fl.My = Ul.Mx*Wl.Vy;
  Fr.My = Ur.Mx*Wr.Vy;

  Fl.Mz = Ul.Mx*Wl.Vz;
  Fr.Mz = Ur.Mx*Wr.Vz;

  Fl.Mx += Wl.P;
  Fr.Mx += Wr.P;

  Fl.E  = (Ul.E + Wl.P)*Wl.Vx;
  Fr.E  = (Ur.E + Wr.P)*Wr.Vx;

  Fl.Mx -= 0.5*(Bxi*Bxi - std::pow(Wl.By, 2) - std::pow(Wl.Bz, 2));
  Fr.Mx -= 0.5*(Bxi*Bxi - std::pow(Wr.By, 2) - std::pow(Wr.Bz, 2));

  Fl.My -= Bxi*Wl.By;
  Fr.My -= Bxi*Wr.By;

  Fl.Mz -= Bxi*Wl.Bz;
  Fr.Mz -= Bxi*Wr.Bz;

  Fl.E += pbl*Wl.Vx - Bxi*(Bxi*Wl.Vx + Wl.By*Wl.Vy + Wl.Bz*Wl.Vz);
  Fr.E += pbr*Wr.Vx - Bxi*(Bxi*Wr.Vx + Wr.By*Wr.Vy + Wr.Bz*Wr.Vz);

  Fl.By = Wl.By*Wl.Vx - Bxi*Wl.Vy;
  Fr.By = Wr.By*Wr.Vx - Bxi*Wr.Vy;

  Fl.Bz = Wl.Bz*Wl.Vx - Bxi*Wl.Vz;
  Fr.Bz = Wr.Bz*Wr.Vx - Bxi*Wr.Vz;

/*--- Step 5. ------------------------------------------------------------------
 * Return upwind flux if flow is supersonic
 */

  if(ev[0] >= 0.0) {
    return Fl;
  }

  if(ev[NWAVE-1] <= 0.0) {
    return Fr;
  }

/*--- Step 6. ------------------------------------------------------------------
 * Compute projection of dU onto L eigenvectors ("vector A")
 */

  dU[0] = Ur.d - Ul.d;
  dU[1] = Ur.Mx - Ul.Mx;
  dU[2] = Ur.My - Ul.My;
  dU[3] = Ur.Mz - Ul.Mz;
  dU[4] = Ur.E - Ul.E;
  dU[5] = Ur.By - Ul.By;
  dU[6] = Ur.Bz - Ul.Bz;

  for (int n=0; n<NWAVE; n++) {
    a[n] = 0.0;
    for (int m=0; m<NWAVE; m++) {
      a[n] += lem[n][m]*dU[m];
    }
  }

/*--- Step 7. ------------------------------------------------------------------
 * Check that the density and pressure in the intermediate states are positive.
 */

  u_inter[0] = Ul.d;
  u_inter[1] = Ul.Mx;
  u_inter[2] = Ul.My;
  u_inter[3] = Ul.Mz;
  u_inter[4] = Ul.E;
  u_inter[5] = Ul.By;
  u_inter[6] = Ul.Bz;

  for (int n=0; n<NWAVE-1; n++) {
    for (int m=0; m<NWAVE; m++) u_inter[m] += a[n]*rem[m][n];
    if(ev[n+1] > ev[n]) {
      if (u_inter[0] <= 0.0) {
        throw std::domain_error(
          std::string("Non-physical density in intermediate state of ")
          + __func__
        );
      }
      p_inter
        = u_inter[4]
        - 0.5 * (
          std::pow(u_inter[1], 2)
          + std::pow(u_inter[2], 2)
          + std::pow(u_inter[3], 2)
        ) / u_inter[0];
      p_inter -= 0.5 * (
        std::pow(u_inter[NWAVE-2], 2)
        + std::pow(u_inter[NWAVE-1], 2)
        + std::pow(Bxi, 2)
      );
      if (p_inter < 0.0) {
        throw std::domain_error(
          std::string("Non-physical pressure in intermediate state of ")
          + __func__
        );
      }
    }
  }

/*--- Step 8. ------------------------------------------------------------------
 * Compute Roe flux */

  for (int m=0; m<NWAVE; m++) {
    coeff[m] = 0.5 * a[m] * std::max(std::fabs(ev[m]), etah);
  }

  Cons1DS flux;
  flux.d = 0.5 * (Fl.d + Fr.d);
  flux.Mx = 0.5 * (Fl.Mx + Fr.Mx);
  flux.My = 0.5 * (Fl.My + Fr.My);
  flux.Mz = 0.5 * (Fl.Mz + Fr.Mz);
  flux.E = 0.5 * (Fl.E + Fr.E);
  flux.By = 0.5 * (Fl.By + Fr.By);
  flux.Bz = 0.5 * (Fl.Bz + Fr.Bz);

  for (size_t m = 0; m < NWAVE; m++) {
    flux.d -= coeff[m] * rem[0][m];
    flux.Mx -= coeff[m] * rem[1][m];
    flux.My -= coeff[m] * rem[2][m];
    flux.Mz -= coeff[m] * rem[3][m];
    flux.E -= coeff[m] * rem[4][m];
    flux.By -= coeff[m] * rem[5][m];
    flux.Bz -= coeff[m] * rem[6][m];
  }

  return flux;
}


/*!
See get_flux_hll() in hll_athena.hpp

Ignores background field.
*/
template <
	class Mass_Density_T,
	class Momentum_Density_T,
	class Total_Energy_Density_T,
	class Magnetic_Field_T,
	class MHD,
	class Vector,
	class Scalar
> std::tuple<MHD, Scalar> get_flux_roe(
	MHD& state_neg,
	MHD& state_pos,
	const Vector& /*bg_face_magnetic_field*/,
	const Scalar& adiabatic_index,
	const Scalar& vacuum_permeability
) {
	using std::to_string;

	const Vector bg_face_magnetic_field{0, 0, 0};

	// shorthand notation for simulation variables
	const Mass_Density_T Mas{};
	const Momentum_Density_T Mom{};
	const Total_Energy_Density_T Nrj{};
	const Magnetic_Field_T Mag{};
	const Velocity Vel{};
	const Pressure Pre{};

	// getter functions for required variables
	const auto Mas_g
		= [](MHD& data)->typename Mass_Density_T::data_type&{
			return data[Mass_Density_T()];
		};
	const auto Mom_g
		= [](MHD& data)->typename Momentum_Density_T::data_type&{
			return data[Momentum_Density_T()];
		};
	const auto Nrj_g
		= [](MHD& data)->typename Total_Energy_Density_T::data_type&{
			return data[Total_Energy_Density_T()];
		};
	const auto Mas_g_p
		= [](HD_Primitive& data)->typename Mass_Density_T::data_type&{
			return data[Mass_Density_T()];
		};
	const auto Vel_g
		= [](HD_Primitive& data)->typename Velocity::data_type&{
			return data[Velocity()];
		};
	const auto Pre_g
		= [](HD_Primitive& data)->typename Pressure::data_type&{
			return data[Pressure()];
		};

	if (state_neg[Mas] <= 0) {
		throw std::domain_error(
			std::string("Non-positive mass density on negative side given to ")
			+ __func__
			+ std::string(": ")
			+ to_string(state_neg[Mas])
		);
	}
	if (state_pos[Mas] <= 0) {
		throw std::domain_error(
			std::string("Non-positive mass density on positive side given to ")
			+ __func__
			+ std::string(": ")
			+ to_string(state_pos[Mas])
		);
	}

	const auto pressure_neg
		= get_pressure(
			state_neg[Mas],
			state_neg[Mom],
			state_neg[Nrj],
			state_neg[Mag],
			adiabatic_index,
			vacuum_permeability
		);
	if (pressure_neg <= 0) {
		throw std::domain_error(
			std::string("Non-positive pressure on negative side given to ")
			+ __func__
			+ std::string(": ")
			+ to_string(pressure_neg)
			+ " with adiabatic index "
			+ to_string(adiabatic_index)
			+ " and vacuum permeability "
			+ to_string(vacuum_permeability)
		);
	}

	const auto pressure_pos
		= get_pressure(
			state_pos[Mas],
			state_pos[Mom],
			state_pos[Nrj],
			state_pos[Mag],
			adiabatic_index,
			vacuum_permeability
		);
	if (pressure_pos <= 0) {
		throw std::domain_error(
			std::string("Non-positive pressure on positive side given to ")
			+ __func__
			+ std::string(": ")
			+ to_string(pressure_pos)
			+ " with adiabatic index "
			+ to_string(adiabatic_index)
			+ " and vacuum permeability "
			+ to_string(vacuum_permeability)
		);
	}

	const double Bxi
		= (state_neg[Mag][0] + state_pos[Mag][0])
		/ (2.0 * std::sqrt(vacuum_permeability));

	Cons1DS
		Ul{
			state_neg[Mas],
			state_neg[Mom][0],
			state_neg[Mom][1],
			state_neg[Mom][2],
			state_neg[Nrj],
			state_neg[Mag][1] / std::sqrt(vacuum_permeability),
			state_neg[Mag][2] / std::sqrt(vacuum_permeability)
		},
		Ur{
			state_pos[Mas],
			state_pos[Mom][0],
			state_pos[Mom][1],
			state_pos[Mom][2],
			state_pos[Nrj],
			state_pos[Mag][1] / std::sqrt(vacuum_permeability),
			state_pos[Mag][2] / std::sqrt(vacuum_permeability)
		};

	const auto
		prim_neg_temp = get_primitive<HD_Primitive>(
			state_neg,
			state_neg[Mag],
			adiabatic_index,
			1.0,
			Mas_g_p, Vel_g, Pre_g,
			Mas_g, Mom_g, Nrj_g
		),
		prim_pos_temp = get_primitive<HD_Primitive>(
			state_pos,
			state_pos[Mag],
			adiabatic_index,
			1.0,
			Mas_g_p, Vel_g, Pre_g,
			Mas_g, Mom_g, Nrj_g
		);

	const Prim1DS
		Wl{
			prim_neg_temp[Mas],
			prim_neg_temp[Vel][0],
			prim_neg_temp[Vel][1],
			prim_neg_temp[Vel][2],
			prim_neg_temp[Pre],
			state_neg[Mag][1] / std::sqrt(vacuum_permeability),
			state_neg[Mag][2] / std::sqrt(vacuum_permeability)
		},
		Wr{
			prim_pos_temp[Mas],
			prim_pos_temp[Vel][0],
			prim_pos_temp[Vel][1],
			prim_pos_temp[Vel][2],
			prim_pos_temp[Pre],
			state_pos[Mag][1] / std::sqrt(vacuum_permeability),
			state_pos[Mag][2] / std::sqrt(vacuum_permeability)
		};

	const auto flux_temp
		= athena_roe_fluxes(
			Ul,
			Ur,
			Wl,
			Wr,
			Bxi,
			adiabatic_index
		);

	MHD flux;
	flux[Mas] = flux_temp.d;
	flux[Mom][0] = flux_temp.Mx;
	flux[Mom][1] = flux_temp.My;
	flux[Mom][2] = flux_temp.Mz;
	flux[Nrj] = flux_temp.E;
	flux[Mag][0] = 0;
	flux[Mag][1] = flux_temp.By * std::sqrt(vacuum_permeability);
	flux[Mag][2] = flux_temp.Bz * std::sqrt(vacuum_permeability);

	// get maximum signal speed
	const auto
		fast_magnetosonic_neg
			= get_fast_magnetosonic_speed(
				state_neg[Mas],
				state_neg[Mom],
				state_neg[Nrj],
				state_neg[Mag],
				bg_face_magnetic_field,
				adiabatic_index,
				vacuum_permeability
			),

		fast_magnetosonic_pos
			= get_fast_magnetosonic_speed(
				state_pos[Mas],
				state_pos[Mom],
				state_pos[Nrj],
				state_pos[Mag],
				bg_face_magnetic_field,
				adiabatic_index,
				vacuum_permeability
			),

		max_signal = std::max(fast_magnetosonic_neg, fast_magnetosonic_pos);

	const auto
		flow_v_neg
			= state_neg[Mom]
			/ state_neg[Mas],

		flow_v_pos
			= state_pos[Mom]
			/ state_pos[Mas];

	const auto
		max_signal_neg
			= (flow_v_neg[0] <= flow_v_pos[0])
			? flow_v_neg[0] - max_signal
			: flow_v_pos[0] - max_signal,

		max_signal_pos
			= (flow_v_neg[0] <= flow_v_pos[0])
			? flow_v_pos[0] + max_signal
			: flow_v_neg[0] + max_signal;

	double ret_signal_speed = 0;
	if (max_signal_neg >= 0.0) {
		ret_signal_speed = std::fabs(max_signal_pos);
	} else if (max_signal_pos <= 0.0) {
		ret_signal_speed = std::fabs(max_signal_neg);
	} else {
		ret_signal_speed
			= std::max(
				std::fabs(max_signal_neg),
				std::fabs(max_signal_pos)
			);
	}

	return std::make_tuple(flux, ret_signal_speed);
}

}}} // namespaces

#endif // ifndef PAMHD_MHD_ROE_ATHENA_HPP
