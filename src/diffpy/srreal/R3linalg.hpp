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
* $Id$
*
*****************************************************************************/

#ifndef R3LINALG_HPP_INCLUDED
#define R3LINALG_HPP_INCLUDED

#include <algorithm>
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>
#include <diffpy/mathutils.hpp>

namespace diffpy {
namespace srreal {
namespace R3 {

// Declarations --------------------------------------------------------------

// Constants

const int Ndim = 3;
using blitz::product;
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

template <class V> double norm(const V&);
template <class V> double distance(const V& u, const V& v);
template <class V> double dot(const V& u, const V& v);
template <class V> Vector cross(const V& u, const V& v);
template <class V> const Vector& mxvecproduct(const Matrix&, const V&);
template <class V> const Vector& mxvecproduct(const V&, const Matrix&);

class EpsCompare;

template <class T>
    bool EpsEqual(const T& A, const T& B, double eps=SQRT_DOUBLE_EPS);

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


class EpsCompare
{
    public:

        // constructor
        EpsCompare(double eps) : mepscmp(eps)  { }

        bool operator()(const Vector& u, const Vector& v) const
        {
            bool rv = std::lexicographical_compare(
                    &(u[0]), &(u[0]) + Ndim,
                    &(v[0]), &(v[0]) + Ndim, mepscmp);
            return rv;
        }

        bool operator()(const Matrix& A, const Matrix& B) const
        {
            bool rv = std::lexicographical_compare(
                    A.data(), A.data() + Ndim * Ndim,
                    B.data(), B.data() + Ndim * Ndim, mepscmp);
            return rv;
        }

    private:

        diffpy::mathutils::EpsilonCompare mepscmp;
};


template <class T>
bool EpsEqual(const T& A, const T& B, double eps)
{
    EpsCompare epscmp(eps);
    return !epscmp(A, B) && !epscmp(B, A);
}


}   // namespace R3
}   // namespace srreal
}   // namespace diffpy

#endif  // R3LINALG_HPP_INCLUDED
