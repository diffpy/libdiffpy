/***********************************************************************
* Short Title: linear algebra functions on R3
*
* Comments: defininitions of linear algebra functions for
*     blitz::TinyVector  and  blitz::TinyMatrix
*
* <license text>
***********************************************************************/

#include <boost/functional/hash.hpp>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <diffpy/srreal/R3linalg.hpp>

namespace diffpy {
namespace srreal {

const R3::Matrix& R3::identity()
{
    static R3::Matrix mx;
    static bool mx_ready = false;
    if (!mx_ready)
    {
        mx = 1.0, 0.0, 0.0,
             0.0, 1.0, 0.0,
             0.0, 0.0, 1.0;
        mx_ready = true;
    }
    return mx;
}


const R3::Matrix& R3::zeros()
{
    static R3::Matrix mx;
    static bool mx_ready = false;
    if (!mx_ready)  mx = 0.0;
    return mx;
}


double R3::determinant(const R3::Matrix& A)
{
    gsl_matrix* gA = gsl_matrix_alloc(R3::Ndim, R3::Ndim);
    for (int i = 0; i != R3::Ndim; ++i)
    {
        for (int j = 0; j != R3::Ndim; ++j)
        {
            gsl_matrix_set(gA, i, j, A(i,j));
        }
    }
    gsl_permutation* gP = gsl_permutation_alloc(R3::Ndim);
    int signum;
    gsl_linalg_LU_decomp(gA, gP, &signum);
    double det = gsl_linalg_LU_det(gA, signum);
    gsl_permutation_free(gP);
    gsl_matrix_free(gA);
    return det;
}


R3::Matrix R3::inverse(const R3::Matrix& A)
{
    gsl_matrix* gA = gsl_matrix_alloc(R3::Ndim, R3::Ndim);
    for (int i = 0; i != R3::Ndim; ++i)
    {
        for (int j = 0; j != R3::Ndim; ++j)
        {
            gsl_matrix_set(gA, i, j, A(i,j));
        }
    }
    gsl_permutation* gP = gsl_permutation_alloc(R3::Ndim);
    int signum;
    gsl_linalg_LU_decomp(gA, gP, &signum);
    R3::Matrix B;
    gsl_matrix_view gB = gsl_matrix_view_array(B.data(), R3::Ndim, R3::Ndim);
    gsl_linalg_LU_invert(gA, gP, &gB.matrix);
    gsl_permutation_free(gP);
    gsl_matrix_free(gA);
    return B;
}


R3::Matrix R3::transpose(const R3::Matrix& A)
{
    R3::Matrix res;
    res = A(0,0), A(1,0), A(2,0),
          A(0,1), A(1,1), A(2,1),
          A(0,2), A(1,2), A(2,2);
    return res;
}


}   // namespace srreal

namespace mathutils {

// EpsilonLess specialization ------------------------------------------------

template<>
bool EpsilonLess::operator()<srreal::R3::Matrix, srreal::R3::Matrix>(
        const srreal::R3::Matrix& A, const srreal::R3::Matrix& B) const
{
    using srreal::R3::Ndim;
    bool rv = std::lexicographical_compare(
            A.data(), A.data() + Ndim * Ndim,
            B.data(), B.data() + Ndim * Ndim, *this);
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
    using srreal::R3::Ndim;
    bool rv = std::equal(A.data(), A.data() + Ndim * Ndim, B.data(), *this);
    return rv;
}

}   // namespace mathutils
}   // namespace diffpy

namespace blitz {

using namespace diffpy::srreal;

size_t hash_value(const R3::Vector& v)
{
    return boost::hash_range(v.begin(), v.end());
}


size_t hash_value(const R3::Matrix& A)
{
    const size_t sz = R3::Ndim * R3::Ndim;
    return boost::hash_range(A.data(), A.data() + sz);
}

}   // namespace blitz

// End of file
