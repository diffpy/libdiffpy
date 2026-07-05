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
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* R3linalg -- vector and matrix types and linar algebra operations in R3 space
*
*****************************************************************************/

#include <algorithm>
#include <cmath>
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


void eigen_solve_3x3(const Matrix& A, Vector& w, Matrix& V)
{
    V = identity();
    Matrix m = A;

    const int max_iter = 50;
    const double eps = 1e-10;

    for (int iter = 0; iter < max_iter; ++iter)
    {
        double max_off_diag = 0.0;
        int p = 0;
        int q = 1;

        for (int i = 0; i < Ndim; ++i)
        {
            for (int j = i + 1; j < Ndim; ++j)
            {
                if (std::abs(m(i, j)) > max_off_diag)
                {
                    max_off_diag = std::abs(m(i, j));
                    p = i;
                    q = j;
                }
            }
        }

        if (max_off_diag < eps)  break;

        double phi = 0.5 * std::atan2(
                2.0 * m(p, q), m(q, q) - m(p, p));
        double c = std::cos(phi);
        double s = std::sin(phi);

        double m_pp = m(p, p);
        double m_qq = m(q, q);
        double m_pq = m(p, q);

        m(p, p) = c * c * m_pp - 2.0 * s * c * m_pq + s * s * m_qq;
        m(q, q) = s * s * m_pp + 2.0 * s * c * m_pq + c * c * m_qq;
        m(p, q) = 0.0;
        m(q, p) = 0.0;

        for (int i = 0; i < Ndim; ++i)
        {
            if (i == p || i == q)  continue;
            double m_ip = m(i, p);
            double m_iq = m(i, q);
            m(i, p) = c * m_ip - s * m_iq;
            m(p, i) = m(i, p);
            m(i, q) = s * m_ip + c * m_iq;
            m(q, i) = m(i, q);
        }

        for (int i = 0; i < Ndim; ++i)
        {
            double v_ip = V(i, p);
            double v_iq = V(i, q);
            V(i, p) = c * v_ip - s * v_iq;
            V(i, q) = s * v_ip + c * v_iq;
        }
    }

    w[0] = m(0, 0);
    w[1] = m(1, 1);
    w[2] = m(2, 2);

    for (int i = 0; i < Ndim - 1; ++i)
    {
        for (int j = 0; j < Ndim - 1 - i; ++j)
        {
            if (w[j] <= w[j + 1])  continue;
            std::swap(w[j], w[j + 1]);
            for (int k = 0; k < Ndim; ++k)
            {
                std::swap(V(k, j), V(k, j + 1));
            }
        }
    }
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
