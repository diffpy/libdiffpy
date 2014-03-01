/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* R3linalg -- vector and matrix types and linar algebra operations in R3 space
*
*****************************************************************************/

#include <boost/functional/hash.hpp>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <diffpy/srreal/R3linalg.hpp>

namespace diffpy {
namespace srreal {
namespace R3 {

const Matrix& identity()
{
    static Matrix mx = ublas::identity_matrix<double>(Ndim);
    return mx;
}


const Matrix& zeromatrix()
{
    static Matrix mx = ublas::zero_matrix<double>(Ndim, Ndim);
    return mx;
}


double determinant(const Matrix& A)
{
    gsl_matrix* gA = gsl_matrix_alloc(Ndim, Ndim);
    for (int i = 0; i != Ndim; ++i)
    {
        for (int j = 0; j != Ndim; ++j)
        {
            gsl_matrix_set(gA, i, j, A(i,j));
        }
    }
    gsl_permutation* gP = gsl_permutation_alloc(Ndim);
    int signum;
    gsl_linalg_LU_decomp(gA, gP, &signum);
    double det = gsl_linalg_LU_det(gA, signum);
    gsl_permutation_free(gP);
    gsl_matrix_free(gA);
    return det;
}


const Matrix& inverse(const Matrix& A)
{
    static Matrix B;
    gsl_matrix* gA = gsl_matrix_alloc(Ndim, Ndim);
    for (int i = 0; i != Ndim; ++i)
    {
        for (int j = 0; j != Ndim; ++j)
        {
            gsl_matrix_set(gA, i, j, A(i,j));
        }
    }
    gsl_permutation* gP = gsl_permutation_alloc(Ndim);
    int signum;
    gsl_linalg_LU_decomp(gA, gP, &signum);
    double* bdata = &(B.data()[0]);
    gsl_matrix_view gB = gsl_matrix_view_array(bdata, Ndim, Ndim);
    gsl_linalg_LU_invert(gA, gP, &gB.matrix);
    gsl_permutation_free(gP);
    gsl_matrix_free(gA);
    return B;
}


size_t hash_value(const Vector& v)
{
    return boost::hash_range(v.begin(), v.end());
}


size_t hash_value(const Matrix& A)
{
    return boost::hash_range(A.data().begin(), A.data().end());
}

}   // namespace R3
}   // namespace srreal

namespace mathutils {

// EpsilonLess specialization ------------------------------------------------

template<>
bool EpsilonLess::operator()<srreal::R3::Matrix, srreal::R3::Matrix>(
        const srreal::R3::Matrix& A, const srreal::R3::Matrix& B) const
{
    bool rv = std::lexicographical_compare(
            A.data().begin(), A.data().end(),
            B.data().begin(), B.data().end(), *this);
    return rv;
}

// EpsilonEqual specializations ----------------------------------------------

template<>
bool EpsilonEqual::operator()<srreal::R3::Vector, srreal::R3::Vector>(
        const srreal::R3::Vector& u, const srreal::R3::Vector& v) const
{
    bool rv = std::equal(u.begin(), u.end(), v.begin(), *this);
    return rv;
}


template<>
bool EpsilonEqual::operator()<srreal::R3::Matrix, srreal::R3::Matrix>(
        const srreal::R3::Matrix& A, const srreal::R3::Matrix& B) const
{
    bool rv = std::equal(A.data().begin(), A.data().end(),
            B.data().begin(), *this);
    return rv;
}

}   // namespace mathutils
}   // namespace diffpy

// End of file
