/*****************************************************************************
 *
 * libdiffpy         by DANSE Diffraction group
 *                   Simon J. L. Billinge
 *                   (c) 2009 The Trustees of Columbia University
 *                   in the City of New York.  All rights reserved.
 *
 * File coded by:    Pavol Juhas
 *
 * See AUTHORS.rst for a list of people who contributed.
 * See LICENSE_DANSE.rst for license information.
 *
 ******************************************************************************
 *
 * class Lattice -- vector and matrix conversions between general and
 *     Cartesian coordinate systems
 *
 *****************************************************************************/

#include <array>
#include <cassert>
#include <stdexcept>

#include <diffpy/srreal/Lattice.hpp>

using namespace std;
using namespace diffpy::srreal;

/////////////////////////////////////////////////////////////////////////////
// class Lattice
/////////////////////////////////////////////////////////////////////////////

// Constructors --------------------------------------------------------------

Lattice::Lattice() {
  mbaserot = R3::identity();
  this->assignLatPar(1.0, 1.0, 1.0, 90.0, 90.0, 90.0);
}

Lattice::Lattice(double a0, double b0, double c0, double alpha0, double beta0,
                 double gamma0) {
  mbaserot = R3::identity();
  this->assignLatPar(a0, b0, c0, alpha0, beta0, gamma0);
}

// Public Methods ------------------------------------------------------------

// The setters leave this lattice unchanged when they throw.

void Lattice::setLatPar(double a0, double b0, double c0, double alpha0,
                        double beta0, double gamma0) {
  Lattice L(*this);
  L.assignLatPar(a0, b0, c0, alpha0, beta0, gamma0);
  *this = L;
}

void Lattice::setLatBase(const R3::Vector& va0, const R3::Vector& vb0,
                         const R3::Vector& vc0) {
  Lattice L(*this);
  L.assignLatBase(va0, vb0, vc0);
  *this = L;
}

void Lattice::assignLatPar(double a0, double b0, double c0, double alpha0,
                           double beta0, double gamma0) {
  using diffpy::mathutils::acosd;
  using diffpy::mathutils::cosd;
  using diffpy::mathutils::sind;
  ma = a0;
  mb = b0;
  mc = c0;
  malpha = alpha0;
  mbeta = beta0;
  mgamma = gamma0;
  mcosa = cosd(malpha);
  msina = sind(malpha);
  mcosb = cosd(mbeta);
  msinb = sind(mbeta);
  mcosg = cosd(mgamma);
  msing = sind(mgamma);
  // Vnorm is a volume of unit cell with a=b=c=1
  double Vnorm = this->volumeNormal();
  // reciprocal lattice
  mar = msina / (ma * Vnorm);
  mbr = msinb / (mb * Vnorm);
  mcr = msing / (mc * Vnorm);
  mcosar = (mcosb * mcosg - mcosa) / (msinb * msing);
  mcosbr = (mcosa * mcosg - mcosb) / (msina * msing);
  mcosgr = (mcosa * mcosb - mcosg) / (msina * msinb);
  msinar = sqrt(1.0 - mcosar * mcosar);
  msinbr = sqrt(1.0 - mcosbr * mcosbr);
  msingr = sqrt(1.0 - mcosgr * mcosgr);
  malphar = acosd(mcosar);
  mbetar = acosd(mcosbr);
  mgammar = acosd(mcosgr);
  // metric tensor
  this->updateMetrics();
  // standard cartesian coordinates of lattice vectors
  this->updateStandardBase();
  // calculate unit cell rotation matrix baserot, base = stdbase * baserot
  mbase = R3::prod(mstdbase, mbaserot);
  mva = R3::row(mbase, 0);
  mvb = R3::row(mbase, 1);
  mvc = R3::row(mbase, 2);
  mrecbase = R3::inverse(mbase);
  mvar = R3::column(mrecbase, 0);
  mvbr = R3::column(mrecbase, 1);
  mvcr = R3::column(mrecbase, 2);
  R3::setrow(mnormbase, 0, R3::row(mbase, 0) * mar);
  R3::setrow(mnormbase, 1, R3::row(mbase, 1) * mbr);
  R3::setrow(mnormbase, 2, R3::row(mbase, 2) * mcr);
  R3::setcolumn(mrecnormbase, 0, R3::column(mrecbase, 0) / mar);
  R3::setcolumn(mrecnormbase, 1, R3::column(mrecbase, 1) / mbr);
  R3::setcolumn(mrecnormbase, 2, R3::column(mrecbase, 2) / mcr);
}

void Lattice::assignLatBase(const R3::Vector& va0, const R3::Vector& vb0,
                            const R3::Vector& vc0) {
  using diffpy::mathutils::acosd;
  R3::setrow(mbase, 0, va0);
  R3::setrow(mbase, 1, vb0);
  R3::setrow(mbase, 2, vc0);
  double detbase = R3::determinant(mbase);
  if (fabs(detbase) < 1.0e-8) {
    throw invalid_argument("base vectors are degenerate");
  } else if (detbase < 0.0) {
    throw invalid_argument("base is not right-handed");
  }
  mva = va0;
  mvb = vb0;
  mvc = vc0;
  ma = R3::norm(mva);
  mb = R3::norm(mvb);
  mc = R3::norm(mvc);
  mcosa = R3::dot(mvb, mvc) / (mb * mc);
  mcosb = R3::dot(mva, mvc) / (ma * mc);
  mcosg = R3::dot(mva, mvb) / (ma * mb);
  msina = sqrt(1.0 - mcosa * mcosa);
  msinb = sqrt(1.0 - mcosb * mcosb);
  msing = sqrt(1.0 - mcosg * mcosg);
  malpha = acosd(mcosa);
  mbeta = acosd(mcosb);
  mgamma = acosd(mcosg);
  // Vnorm is a volume of unit cell with a=b=c=1
  double Vnorm = this->volumeNormal();
  // reciprocal lattice
  mar = msina / (ma * Vnorm);
  mbr = msinb / (mb * Vnorm);
  mcr = msing / (mc * Vnorm);
  mcosar = (mcosb * mcosg - mcosa) / (msinb * msing);
  mcosbr = (mcosa * mcosg - mcosb) / (msina * msing);
  mcosgr = (mcosa * mcosb - mcosg) / (msina * msinb);
  msinar = sqrt(1.0 - mcosar * mcosar);
  msinbr = sqrt(1.0 - mcosbr * mcosbr);
  msingr = sqrt(1.0 - mcosgr * mcosgr);
  malphar = acosd(mcosar);
  mbetar = acosd(mcosbr);
  mgammar = acosd(mcosgr);
  // standard orientation of lattice vectors
  this->updateStandardBase();
  // calculate unit cell rotation matrix,  base = stdbase*baserot
  mbaserot = R3::prod(R3::inverse(mstdbase), mbase);
  mrecbase = R3::inverse(mbase);
  mvar = R3::column(mrecbase, 0);
  mvbr = R3::column(mrecbase, 1);
  mvcr = R3::column(mrecbase, 2);
  // bases normalized to unit reciprocal vectors
  R3::setrow(mnormbase, 0, R3::row(mbase, 0) * mar);
  R3::setrow(mnormbase, 1, R3::row(mbase, 1) * mbr);
  R3::setrow(mnormbase, 2, R3::row(mbase, 2) * mcr);
  R3::setcolumn(mrecnormbase, 0, R3::column(mrecbase, 0) / mar);
  R3::setcolumn(mrecnormbase, 1, R3::column(mrecbase, 1) / mbr);
  R3::setcolumn(mrecnormbase, 2, R3::column(mrecbase, 2) / mcr);
  this->updateMetrics();
}

// volumeNormal is a volume of unit cell with a=b=c=1
double Lattice::volumeNormal() const {
  double rv = sqrt(1.0 + 2.0 * mcosa * mcosb * mcosg - mcosa * mcosa -
                   mcosb * mcosb - mcosg * mcosg);
  return rv;
}

double Lattice::volume() const {
  double rv = this->volumeNormal() * this->a() * this->b() * this->c();
  return rv;
}

R3::Vector Lattice::cartesian(const R3::Vector& lv) const {
  return R3::mxvecproduct(lv, mbase);
}

R3::Vector Lattice::fractional(const R3::Vector& cv) const {
  return R3::mxvecproduct(cv, mrecbase);
}

R3::Vector Lattice::ucvCartesian(const R3::Vector& cv) const {
  return cartesian(ucvFractional(fractional(cv)));
}

R3::Vector Lattice::ucvFractional(const R3::Vector& lv) const {
  using mathutils::eps_eq;
  R3::Vector res = lv - floor(lv);
  if (eps_eq(res[0], 1.0)) res[0] = 0.0;
  if (eps_eq(res[1], 1.0)) res[1] = 0.0;
  if (eps_eq(res[2], 1.0)) res[2] = 0.0;
  return res;
}

R3::Matrix Lattice::cartesianMatrix(const R3::Matrix& Ml) const {
  return R3::prod(R3::trans(mnormbase), R3::prod(Ml, mnormbase));
}

R3::Matrix Lattice::fractionalMatrix(const R3::Matrix& Mc) const {
  return R3::prod(R3::trans(mrecnormbase), R3::prod(Mc, mrecnormbase));
}

const R3::Vector& Lattice::ucMaxDiagonal() const {
  static constexpr array<R3::Vector, 4> ucdiagonals = {{
    {+1, +1, +1},
    {-1, +1, +1},
    {+1, -1, +1},
    {+1, +1, -1},
  }};
  double maxnorm = -1;
  const R3::Vector* maxucd = nullptr;
  for (const R3::Vector& ucd : ucdiagonals) {
    double normucd = this->norm(ucd);
    if (normucd > maxnorm) {
      maxnorm = normucd;
      maxucd = &ucd;
    }
  }
  assert(maxucd);
  return *maxucd;
}

double Lattice::ucMaxDiagonalLength() const {
  const R3::Vector& ucmd = this->ucMaxDiagonal();
  double res = this->norm(ucmd);
  return res;
}

bool Lattice::operator==(const Lattice& other) const {
  bool rv = (this->base() == other.base());
  return rv;
}

bool Lattice::operator!=(const Lattice& other) const {
  return !(*this == other);
}

// Private Methods -----------------------------------------------------------

void Lattice::updateMetrics() {
  mmetrics(0, 0) = ma * ma;
  mmetrics(1, 1) = mb * mb;
  mmetrics(2, 2) = mc * mc;
  mmetrics(0, 1) = mmetrics(1, 0) = ma * mb * mcosg;
  mmetrics(0, 2) = mmetrics(2, 0) = ma * mc * mcosb;
  mmetrics(1, 2) = mmetrics(2, 1) = mb * mc * mcosa;
}

void Lattice::updateStandardBase() {
  // a
  mstdbase(0, 0) = 1.0 / mar;
  mstdbase(1, 0) = 0.0;
  mstdbase(2, 0) = 0.0;
  // b
  mstdbase(0, 1) = -mcosgr / msingr / mar;
  mstdbase(1, 1) = mb * msina;
  mstdbase(2, 1) = 0.0;
  // c
  mstdbase(0, 2) = mcosb * ma;
  mstdbase(1, 2) = mb * mcosa;
  mstdbase(2, 2) = mc;
}

// End of file
