/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* R3linalg -- vector and matrix types and linar algebra
*     operations in R3 space
*
*****************************************************************************/

#ifndef R3LINALG_HPP_INCLUDED
#define R3LINALG_HPP_INCLUDED

#include <algorithm>
#include <boost/serialization/base_object.hpp>
#include <blitz/array.h>
#include <diffpy/mathutils.hpp>

namespace diffpy {
namespace srreal {
namespace R3 {

// Declarations --------------------------------------------------------------

// Constants

const int Ndim = 3;
using ::diffpy::mathutils::SQRT_DOUBLE_EPS;

// Types

typedef blitz::TinyMatrix<double,Ndim,Ndim> Matrix;
typedef blitz::TinyVector<double,Ndim> Vector;

// Functions

const Matrix& identity();
const Matrix& zeros();
double determinant(const Matrix& A);
Matrix inverse(const Matrix& A);
Matrix transpose(const Matrix& A);
const Matrix& product(const Matrix&, const Matrix&);

template <class V> double norm(const V&);
template <class V> double distance(const V& u, const V& v);
template <class V> double dot(const V& u, const V& v);
template <class V> Vector cross(const V& u, const V& v);
template <class V> const Vector& mxvecproduct(const Matrix&, const V&);
template <class V> const Vector& mxvecproduct(const V&, const Matrix&);

// Template functions --------------------------------------------------------

template <class V>
double norm(const V& u)
{
    return sqrt(R3::dot(u, u));
}


template <class V>
double distance(const V& u, const V& v)
{
    static R3::Vector duv;
    duv[0] = u[0] - v[0];
    duv[1] = u[1] - v[1];
    duv[2] = u[2] - v[2];
    return R3::norm(duv);
}


template <class V>
double dot(const V& u, const V& v)
{
    return (u[0]*v[0] + u[1]*v[1] + u[2]*v[2]);
}


template <class V>
Vector cross(const V& u, const V& v)
{
    Vector res;
    res[0] = u[1]*v[2] - u[2]*v[1];
    res[1] = u[2]*v[0] - u[0]*v[2];
    res[2] = u[0]*v[1] - u[1]*v[0];
    return res;
}


template <class V>
const Vector& mxvecproduct(const Matrix& M, const V& u)
{
    static Vector res;
    res[0] = M(0,0)*u[0] + M(0,1)*u[1] + M(0,2)*u[2];
    res[1] = M(1,0)*u[0] + M(1,1)*u[1] + M(1,2)*u[2];
    res[2] = M(2,0)*u[0] + M(2,1)*u[1] + M(2,2)*u[2];
    return res;
}


template <class V>
const Vector& mxvecproduct(const V& u, const Matrix& M)
{
    static Vector res;
    res[0] = u[0]*M(0,0) + u[1]*M(1,0) + u[2]*M(2,0);
    res[1] = u[0]*M(0,1) + u[1]*M(1,1) + u[2]*M(2,1);
    res[2] = u[0]*M(0,2) + u[1]*M(1,2) + u[2]*M(2,2);
    return res;
}

}   // namespace R3
}   // namespace srreal

namespace mathutils {

// Template specializations --------------------------------------------------

template<>
bool EpsilonLess::operator()<srreal::R3::Matrix, srreal::R3::Matrix>(
        const srreal::R3::Matrix& A, const srreal::R3::Matrix& B) const;

template<>
bool EpsilonEqual::operator()<srreal::R3::Vector, srreal::R3::Vector>(
        const srreal::R3::Vector& u, const srreal::R3::Vector& v) const;

template<>
bool EpsilonEqual::operator()<srreal::R3::Matrix, srreal::R3::Matrix>(
        const srreal::R3::Matrix& A, const srreal::R3::Matrix& B) const;

}   // namespace mathutils
}   // namespace diffpy

namespace blitz {

// Equality ------------------------------------------------------------------

inline
bool operator==(const diffpy::srreal::R3::Vector& u,
        const diffpy::srreal::R3::Vector& v)
{
    diffpy::srreal::R3::Vector::const_iterator pu = u.begin();
    diffpy::srreal::R3::Vector::const_iterator pv = v.begin();
    bool rv = (pu == pv) || (
            *pu++ == *pv++ &&
            *pu++ == *pv++ &&
            *pu++ == *pv++
            );
    return rv;
}


inline
bool operator==(const diffpy::srreal::R3::Matrix& A,
        const diffpy::srreal::R3::Matrix& B)
{
    using diffpy::srreal::R3::Ndim;
    const size_t sz = Ndim * Ndim;
    bool rv = (&A == &B) || std::equal(A.data(), A.data() + sz, B.data());
    return rv;
}


// Hashing -------------------------------------------------------------------

size_t hash_value(const diffpy::srreal::R3::Vector& v);
size_t hash_value(const diffpy::srreal::R3::Matrix& A);

}   // namespace blitz

// Serialization -------------------------------------------------------------

namespace boost {
namespace serialization {

template<class Archive>
void serialize(Archive& ar,
        diffpy::srreal::R3::Vector& v, const unsigned int version)
{
    ar & v[0] & v[1] & v[2];
}


template<class Archive>
void serialize(Archive& ar,
        diffpy::srreal::R3::Matrix& A, const unsigned int version)
{
    using namespace diffpy::srreal;
    R3::Matrix::T_numtype* p = A.data();
    R3::Matrix::T_numtype* plast = p + R3::Ndim * R3::Ndim;
    for (; p != plast; ++p)     ar & (*p);
}


} // namespace serialization
} // namespace boost


#endif  // R3LINALG_HPP_INCLUDED
