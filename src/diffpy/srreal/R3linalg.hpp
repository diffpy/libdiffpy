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
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <diffpy/mathutils.hpp>

namespace diffpy {
namespace srreal {
namespace R3 {

// Declarations --------------------------------------------------------------

namespace ublas = boost::numeric::ublas;
using ublas::prod;
using ublas::row;
using ublas::column;
using ublas::trans;

// Constants

const int Ndim = 3;
using ::diffpy::mathutils::SQRT_DOUBLE_EPS;
const ublas::zero_vector<double> zerovector(Ndim);

// Classes

class Vector : public ublas::vector<double, ublas::bounded_array<double,3> >
{
        typedef ublas::vector<double, ublas::bounded_array<double,3> >
            BaseVector;

    public:

	// constructors
	Vector() : BaseVector(3)  { }

	Vector(const double& x, const double& y, const double& z)
            : BaseVector(3)
        {
            Vector::iterator xi = this->begin();
            *(xi++) = x;
            *(xi++) = y;
            *(xi++) = z;
        }

	template <class T>
            Vector(const ublas::vector_expression<T>& r) : BaseVector(r)
        { }

	template <class T>
            void operator=(const ublas::vector_expression<T>& r)
        {
            BaseVector::operator=(r);
        }

	template <class T>
            void operator=(const BaseVector& r)
        {
            BaseVector::operator=(r);
        }

    private:

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            using boost::serialization::base_object;
            ar & base_object<BaseVector>(*this);
        }
};


class Matrix : public ublas::matrix<double,
    ublas::row_major, ublas::bounded_array<double,9> >
{
    typedef ublas::matrix<double, ublas::row_major,
        ublas::bounded_array<double,9> >  BaseMatrix;

    public:

	// constructors
	Matrix() : BaseMatrix(3, 3)  { }

	Matrix(const double& x0, const double& x1, const double& x2,
               const double& x3, const double& x4, const double& x5,
               const double& x6, const double& x7, const double& x8) :
            BaseMatrix(3, 3)
        {
            Matrix::array_type::iterator xi = this->data().begin();
            *(xi++) = x0; *(xi++) = x1; *(xi++) = x2;
            *(xi++) = x3; *(xi++) = x4; *(xi++) = x5;
            *(xi++) = x6; *(xi++) = x7; *(xi++) = x8;
        }

	template <class T>
            Matrix(const ublas::matrix_expression<T>& r) : BaseMatrix(r)
        { }

	template <class T>
            void operator=(const ublas::matrix_expression<T>& r)
        {
            BaseMatrix::operator=(r);
        }

	template <class T>
            void operator=(const BaseMatrix& r)
        {
            BaseMatrix::operator=(r);
        }

    private:

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            using boost::serialization::base_object;
            ar & base_object<BaseMatrix>(*this);
        }
};

// Functions

const Matrix& identity();
const Matrix& zeromatrix();
double determinant(const Matrix& A);
const Matrix& inverse(const Matrix& A);

const Vector& floor(const Vector&);
template <class V> double norm(const V&);
template <class V> double distance(const V& u, const V& v);
template <class V> double dot(const V& u, const V& v);
template <class V> Vector cross(const V& u, const V& v);
template <class V> const Vector& mxvecproduct(const Matrix&, const V&);
template <class V> const Vector& mxvecproduct(const V&, const Matrix&);

// Equality ------------------------------------------------------------------

inline
bool operator==(const Vector& u, const Vector& v)
{
    bool rv = std::equal(u.begin(), u.end(), v.begin());
    return rv;
}


inline
bool operator==(const Matrix& A, const Matrix& B)
{
    bool rv = (&A == &B) ||
        std::equal(A.data().begin(), A.data().end(), B.data().begin());
    return rv;
}

// Hashing -------------------------------------------------------------------

size_t hash_value(const Vector& v);
size_t hash_value(const Matrix& A);

// Inlined functions ---------------------------------------------------------

inline
const Vector& floor(const Vector& v)
{
    static Vector res;
    Vector::const_iterator xi = v.begin();
    Vector::iterator xo = res.begin();
    for (; xi != v.end(); ++xi, ++xo)  *xo = std::floor(*xi);
    return res;
}

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

#endif  // R3LINALG_HPP_INCLUDED
