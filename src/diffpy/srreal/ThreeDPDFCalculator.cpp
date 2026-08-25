#include <diffpy/srreal/ThreeDPDFCalculator.hpp>
#include <diffpy/srreal/BaseBondGenerator.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <diffpy/serialization.hpp>
#include <diffpy/srreal/R3linalg.hpp>
#include <diffpy/srreal/PQEvaluator.hpp>
#include <gsl/gsl_fft_complex.h>
#include <fstream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

namespace diffpy {
namespace srreal {

// Static constants
constexpr double ThreeDPDFCalculator::DEFAULT_RMAX_3D;
constexpr double ThreeDPDFCalculator::DEFAULT_GRID_STEP;

// Constructor ---------------------------------------------------------------

ThreeDPDFCalculator::ThreeDPDFCalculator() :
    mdr(DEFAULT_GRID_STEP),
    mnbins(0),
    maccumblocksize(32),
    mapplyrho0background3d(true),
    musecqwindow3d(true),
    mcalculationmode3d(0),
    mhistogramweightmode3d(0),
    menablennDelta3d(false),
    mnnDelta3d(0.0),
    mnnDelta3dUpperBound(0.0),
    mnnDeltaPositiveEta3d(0.9),
    mdistanceDelta1_3d(0.0),
    mdistanceDelta2_3d(0.0),
    mdeltaPairA3d("*"),
    mdeltaPairB3d("*"),
    mdeltaShellIndex3d(1),
    mdeltaShellTolerance3d(0.05),
    mdeltaKeyTolerance3d(1.0e-6),
    museadpscaleSensitivity3d(false),
    madpScale3d(1.0),
    mrho0backgroundscale3d(1.0)
{
    this->registerDoubleAttribute("enable_nn_delta3d",
            this,
            &ThreeDPDFCalculator::getEnableNNDelta3DAttr,
            &ThreeDPDFCalculator::setEnableNNDelta3DAttr);
    this->registerDoubleAttribute("nn_delta3d",
            this,
            &ThreeDPDFCalculator::getNNDelta3D,
            &ThreeDPDFCalculator::setNNDelta3D);
    this->registerDoubleAttribute("nn_delta3d_upper_bound",
            this,
            &ThreeDPDFCalculator::getNNDelta3DUpperBound);
    this->registerDoubleAttribute("nn_delta_positive_eta3d",
            this,
            &ThreeDPDFCalculator::getNNDeltaPositiveEta3D,
            &ThreeDPDFCalculator::setNNDeltaPositiveEta3D);
    this->registerDoubleAttribute("delta1_3d",
            this,
            &ThreeDPDFCalculator::getDistanceDelta1_3D,
            &ThreeDPDFCalculator::setDistanceDelta1_3D);
    this->registerDoubleAttribute("delta2_3d",
            this,
            &ThreeDPDFCalculator::getDistanceDelta2_3D,
            &ThreeDPDFCalculator::setDistanceDelta2_3D);
    this->registerDoubleAttribute("delta_shell_index3d",
            this,
            &ThreeDPDFCalculator::getDeltaShellIndex3DAttr,
            &ThreeDPDFCalculator::setDeltaShellIndex3DAttr);
    this->registerDoubleAttribute("delta_shell_tolerance3d",
            this,
            &ThreeDPDFCalculator::getDeltaShellTolerance3D,
            &ThreeDPDFCalculator::setDeltaShellTolerance3D);
    this->registerDoubleAttribute("delta_key_tolerance3d",
            this,
            &ThreeDPDFCalculator::getDeltaKeyTolerance3D,
            &ThreeDPDFCalculator::setDeltaKeyTolerance3D);
    this->registerDoubleAttribute("use_adp_scale_sensitivity3d",
            this,
            &ThreeDPDFCalculator::getUseADPScaleSensitivity3DAttr,
            &ThreeDPDFCalculator::setUseADPScaleSensitivity3DAttr);
    this->registerDoubleAttribute("adp_scale3d",
            this,
            &ThreeDPDFCalculator::getADPScale3D,
            &ThreeDPDFCalculator::setADPScale3D);
    this->registerDoubleAttribute("rho0_background_scale3d",
            this,
            &ThreeDPDFCalculator::getRho0BackgroundScale3D,
            &ThreeDPDFCalculator::setRho0BackgroundScale3D);
    this->registerDoubleAttribute("calculation_mode3d",
            this,
            &ThreeDPDFCalculator::getCalculationMode3DAttr,
            &ThreeDPDFCalculator::setCalculationMode3DAttr);
    this->registerDoubleAttribute("histogram_weight_mode3d",
            this,
            &ThreeDPDFCalculator::getHistogramWeightMode3DAttr,
            &ThreeDPDFCalculator::setHistogramWeightMode3DAttr);

    this->setRmax(DEFAULT_RMAX_3D);
    this->setRstep(mdr);
    this->setQmax(12.0);
}

// Public Methods ------------------------------------------------------------

void ThreeDPDFCalculator::setGridStep(double dr)
{
    if (dr <= 0) throw std::invalid_argument("Grid step must be positive.");
    if (dr != mdr)
    {
        mdr = dr;
        this->setRstep(dr);
        this->resetValue(); // triggers re-allocation
    }
}

double ThreeDPDFCalculator::getGridStep() const
{
    return mdr;
}

void ThreeDPDFCalculator::setAccumBlockSize(int bs)
{
    if (bs <= 0) throw std::invalid_argument("Accumulation block size must be positive.");
    maccumblocksize = bs;
}

int ThreeDPDFCalculator::getAccumBlockSize() const
{
    return maccumblocksize;
}

void ThreeDPDFCalculator::setApplyRho0Background3D(bool v)
{
    mapplyrho0background3d = v;
}

bool ThreeDPDFCalculator::getApplyRho0Background3D() const
{
    return mapplyrho0background3d;
}

void ThreeDPDFCalculator::setUseCQWindow3D(bool v)
{
    musecqwindow3d = v;
}

bool ThreeDPDFCalculator::getUseCQWindow3D() const
{
    return musecqwindow3d;
}

void ThreeDPDFCalculator::setCalculationMode3D(int mode)
{
    if (mode != 0 && mode != 1)
        throw std::invalid_argument("3D calculation mode must be 0 (ADP PDF) or 1 (vector histogram).");
    if (mcalculationmode3d != mode)  mticker.click();
    mcalculationmode3d = mode;
}

int ThreeDPDFCalculator::getCalculationMode3D() const
{
    return mcalculationmode3d;
}

double ThreeDPDFCalculator::getCalculationMode3DAttr() const
{
    return static_cast<double>(mcalculationmode3d);
}

void ThreeDPDFCalculator::setCalculationMode3DAttr(double v)
{
    const int mode = static_cast<int>(std::lround(v));
    if (std::fabs(v - static_cast<double>(mode)) > 1.0e-8)
        throw std::invalid_argument("3D calculation mode must be an integer value.");
    this->setCalculationMode3D(mode);
}

void ThreeDPDFCalculator::setHistogramWeightMode3D(int mode)
{
    if (mode != 0 && mode != 1)
        throw std::invalid_argument("3D histogram weight mode must be 0 (scattering) or 1 (count).");
    if (mhistogramweightmode3d != mode)  mticker.click();
    mhistogramweightmode3d = mode;
}

int ThreeDPDFCalculator::getHistogramWeightMode3D() const
{
    return mhistogramweightmode3d;
}

double ThreeDPDFCalculator::getHistogramWeightMode3DAttr() const
{
    return static_cast<double>(mhistogramweightmode3d);
}

void ThreeDPDFCalculator::setHistogramWeightMode3DAttr(double v)
{
    const int mode = static_cast<int>(std::lround(v));
    if (std::fabs(v - static_cast<double>(mode)) > 1.0e-8)
        throw std::invalid_argument("3D histogram weight mode must be an integer value.");
    this->setHistogramWeightMode3D(mode);
}

void ThreeDPDFCalculator::setEnableNNDelta3D(bool v)
{
    if (menablennDelta3d != v)  mticker.click();
    menablennDelta3d = v;
}

bool ThreeDPDFCalculator::getEnableNNDelta3D() const
{
    return menablennDelta3d;
}

double ThreeDPDFCalculator::getEnableNNDelta3DAttr() const
{
    return menablennDelta3d ? 1.0 : 0.0;
}

void ThreeDPDFCalculator::setEnableNNDelta3DAttr(double v)
{
    this->setEnableNNDelta3D(v != 0.0);
}

void ThreeDPDFCalculator::setNNDelta3D(double v)
{
    if (v < 0.0) throw std::invalid_argument("3D NN delta must be non-negative.");
    const double bound = this->currentNNDeltaUpperBound();
    if (!mdeltaShellRecords.empty())
    {
        if (bound <= 0.0 && v > 0.0)
            throw std::invalid_argument("3D NN delta cannot be positive without a positive-definite shell bound.");
        if (bound > 0.0 && v >= bound)
            throw std::invalid_argument("3D NN delta exceeds the positive-definite shell bound.");
    }
    if (mnnDelta3d != v)  mticker.click();
    mnnDelta3d = v;
}

const double& ThreeDPDFCalculator::getNNDelta3D() const
{
    return mnnDelta3d;
}

double ThreeDPDFCalculator::getNNDelta3DUpperBound() const
{
    return this->currentNNDeltaUpperBound();
}

const double& ThreeDPDFCalculator::getNNDeltaPositiveEta3D() const
{
    return mnnDeltaPositiveEta3d;
}

void ThreeDPDFCalculator::setNNDeltaPositiveEta3D(double v)
{
    if (v <= 0.0 || v >= 1.0)
        throw std::invalid_argument("Positive-definite eta must be between 0 and 1.");
    const double minproj = (mnnDeltaPositiveEta3d > 0.0) ?
        (mnnDelta3dUpperBound / mnnDeltaPositiveEta3d) : 0.0;
    if (mnnDeltaPositiveEta3d != v)  mticker.click();
    mnnDeltaPositiveEta3d = v;
    if (minproj > 0.0)  mnnDelta3dUpperBound = mnnDeltaPositiveEta3d * minproj;
    this->validateNNDelta3D();
    this->validateDistanceDecayDelta3D();
}

void ThreeDPDFCalculator::setDistanceDelta1_3D(double v)
{
    if (v < 0.0) throw std::invalid_argument("3D distance-delta delta1 must be non-negative.");
    if (mdistanceDelta1_3d != v)  mticker.click();
    mdistanceDelta1_3d = v;
    this->validateDistanceDecayDelta3D();
}

const double& ThreeDPDFCalculator::getDistanceDelta1_3D() const
{
    return mdistanceDelta1_3d;
}

void ThreeDPDFCalculator::setDistanceDelta2_3D(double v)
{
    if (v < 0.0) throw std::invalid_argument("3D distance-delta delta2 must be non-negative.");
    if (mdistanceDelta2_3d != v)  mticker.click();
    mdistanceDelta2_3d = v;
    this->validateDistanceDecayDelta3D();
}

const double& ThreeDPDFCalculator::getDistanceDelta2_3D() const
{
    return mdistanceDelta2_3d;
}

void ThreeDPDFCalculator::setDeltaPairTypes3D(const std::string& a, const std::string& b)
{
    if (a.empty() || b.empty())
        throw std::invalid_argument("Delta pair atom types must be non-empty.");
    if (mdeltaPairA3d != a || mdeltaPairB3d != b)  mticker.click();
    mdeltaPairA3d = a;
    mdeltaPairB3d = b;
}

const std::string& ThreeDPDFCalculator::getDeltaPairA3D() const
{
    return mdeltaPairA3d;
}

const std::string& ThreeDPDFCalculator::getDeltaPairB3D() const
{
    return mdeltaPairB3d;
}

void ThreeDPDFCalculator::setDeltaShellIndex3D(int v)
{
    if (v <= 0) throw std::invalid_argument("Delta shell index must be positive.");
    if (mdeltaShellIndex3d != v)  mticker.click();
    mdeltaShellIndex3d = v;
}

int ThreeDPDFCalculator::getDeltaShellIndex3D() const
{
    return mdeltaShellIndex3d;
}

double ThreeDPDFCalculator::getDeltaShellIndex3DAttr() const
{
    return static_cast<double>(mdeltaShellIndex3d);
}

void ThreeDPDFCalculator::setDeltaShellIndex3DAttr(double v)
{
    this->setDeltaShellIndex3D(static_cast<int>(std::floor(v + 0.5)));
}

const double& ThreeDPDFCalculator::getDeltaShellTolerance3D() const
{
    return mdeltaShellTolerance3d;
}

void ThreeDPDFCalculator::setDeltaShellTolerance3D(double v)
{
    if (v <= 0.0) throw std::invalid_argument("Delta shell tolerance must be positive.");
    if (mdeltaShellTolerance3d != v)  mticker.click();
    mdeltaShellTolerance3d = v;
}

const double& ThreeDPDFCalculator::getDeltaKeyTolerance3D() const
{
    return mdeltaKeyTolerance3d;
}

void ThreeDPDFCalculator::setDeltaKeyTolerance3D(double v)
{
    if (v <= 0.0) throw std::invalid_argument("Delta key tolerance must be positive.");
    if (mdeltaKeyTolerance3d != v)  mticker.click();
    mdeltaKeyTolerance3d = v;
}

int ThreeDPDFCalculator::getDeltaEligiblePairCount3D() const
{
    return static_cast<int>(mdeltaEligibleBondKeys.size());
}

void ThreeDPDFCalculator::setUseADPScaleSensitivity3D(bool v)
{
    if (museadpscaleSensitivity3d != v)  mticker.click();
    museadpscaleSensitivity3d = v;
    this->validateNNDelta3D();
}

bool ThreeDPDFCalculator::getUseADPScaleSensitivity3D() const
{
    return museadpscaleSensitivity3d;
}

double ThreeDPDFCalculator::getUseADPScaleSensitivity3DAttr() const
{
    return museadpscaleSensitivity3d ? 1.0 : 0.0;
}

void ThreeDPDFCalculator::setUseADPScaleSensitivity3DAttr(double v)
{
    this->setUseADPScaleSensitivity3D(v != 0.0);
}

void ThreeDPDFCalculator::setADPScale3D(double v)
{
    if (v <= 0.0) throw std::invalid_argument("3D ADP scale must be positive.");
    if (madpScale3d != v)  mticker.click();
    madpScale3d = v;
    this->validateNNDelta3D();
}

const double& ThreeDPDFCalculator::getADPScale3D() const
{
    return madpScale3d;
}

void ThreeDPDFCalculator::setRho0BackgroundScale3D(double v)
{
    if (v < 0.0) throw std::invalid_argument("3D rho0 background scale must be non-negative.");
    if (mrho0backgroundscale3d != v)  mticker.click();
    mrho0backgroundscale3d = v;
}

const double& ThreeDPDFCalculator::getRho0BackgroundScale3D() const
{
    return mrho0backgroundscale3d;
}

QuantityType ThreeDPDFCalculator::getThreeDPDF() const
{
    QuantityType result;
    std::vector<double> grid = mgrid3d;
    if (!this->isVectorHistogramMode3D())
    {
        this->applyPostProcessing3D(grid);
    }

    size_t nonzero_count = 0;
    for (double val : grid) {
        if (val != 0.0) ++nonzero_count;
    }
    result.reserve(nonzero_count * 4);

    for (size_t i = 0; i < grid.size(); ++i)
    {
        if (grid[i] != 0.0)
        {
            double x, y, z;
            indexToCoord(i, x, y, z);

            result.push_back(x);
            result.push_back(y);
            result.push_back(z);
            result.push_back(grid[i]);
        }
    }
    return result;
}

QuantityType ThreeDPDFCalculator::getRadialHistogram3D() const
{
    QuantityType result;
    result.reserve(mradialhistogram3d.size() * 2);
    for (size_t i = 0; i < mradialhistogram3d.size(); ++i)
    {
        result.push_back(static_cast<double>(i) * mdr);
        result.push_back(mradialhistogram3d[i]);
    }
    return result;
}

void ThreeDPDFCalculator::exportGrid3DBinary(const std::string& path, bool usefloat32, bool applypost) const
{
    std::vector<double> grid = mgrid3d;
    if (applypost)
    {
        this->applyPostProcessing3D(grid);
    }

    std::ofstream ofs(path.c_str(), std::ios::binary | std::ios::trunc);
    if (!ofs) throw std::runtime_error("Failed to open output file: " + path);

    const size_t nxy = static_cast<size_t>(mnbins) * mnbins;
    if (usefloat32)
    {
        std::vector<float> row(nxy);
        for (int iz = 0; iz < mnbins; ++iz)
        {
            const size_t base = static_cast<size_t>(iz) * nxy;
            for (size_t i = 0; i < nxy; ++i) row[i] = static_cast<float>(grid[base + i]);
            ofs.write(reinterpret_cast<const char*>(row.data()), static_cast<std::streamsize>(nxy * sizeof(float)));
        }
    }
    else
    {
        for (int iz = 0; iz < mnbins; ++iz)
        {
            const size_t base = static_cast<size_t>(iz) * nxy;
            ofs.write(reinterpret_cast<const char*>(&grid[base]), static_cast<std::streamsize>(nxy * sizeof(double)));
        }
    }
    if (!ofs) throw std::runtime_error("Failed while writing output file: " + path);
}

void ThreeDPDFCalculator::applyPostProcessing3D(std::vector<double>& grid) const
{
    if (musecqwindow3d) applyQWindow3D(grid);

    const double rdf_scale = this->getRDFScale();
    if (rdf_scale != 1.0)
    {
        for (double& val : grid) val *= rdf_scale;
    }

    if (mapplyrho0background3d)
    {
        const double rho0_bg = computeRho0Background();
        if (rho0_bg != 0.0)
        {
            for (double& val : grid) val -= rho0_bg;
        }
    }

    double qdamp = 0.0;
    try
    {
        qdamp = this->getEnvelopeByType("qresolution")->getDoubleAttr("qdamp");
    }
    catch (...)
    {
        qdamp = 0.0;
    }
    if (qdamp > 0.0)
    {
        for (size_t i = 0; i < grid.size(); ++i)
        {
            if (grid[i] == 0.0) continue;
            double x, y, z;
            indexToCoord(i, x, y, z);
            const double r = sqrt(x * x + y * y + z * z);
            grid[i] *= exp(-0.5 * (r * qdamp) * (r * qdamp));
        }
    }
}

double ThreeDPDFCalculator::computeRho0Background() const
{
    const StructureAdapterPtr& structure = this->getStructure();
    if (!structure)  return 0.0;
    const double partialpdfscale = this->getPartialPDFScale();
    return mrho0backgroundscale3d * partialpdfscale * structure->numberDensity();
}

void ThreeDPDFCalculator::applyQWindow3D(std::vector<double>& grid) const
{
    if (grid.empty())  return;

    const double qmin = this->getQmin();
    const double qmax = this->getQmax();
    if (qmin <= 0.0 && qmax <= 0.0)  return;

    const int n = mnbins;
    int npad = 1;
    while (npad < n)  npad <<= 1;

    const size_t n3 = static_cast<size_t>(npad) * npad * npad;
    std::vector<double> data(2 * n3, 0.0);

    const int center = n / 2;
    for (int iz = 0; iz < n; ++iz)
    {
        const int sz = (iz - center + npad) % npad;
        for (int iy = 0; iy < n; ++iy)
        {
            const int sy = (iy - center + npad) % npad;
            const size_t base_src = static_cast<size_t>(iz * n + iy) * n;
            const size_t base_dst = static_cast<size_t>(sz * npad + sy) * npad;
            for (int ix = 0; ix < n; ++ix)
            {
                const int sx = (ix - center + npad) % npad;
                const size_t src = base_src + ix;
                const size_t dst = base_dst + sx;
                data[2 * dst] = grid[src];
            }
        }
    }

    for (int iz = 0; iz < npad; ++iz)
    {
        for (int iy = 0; iy < npad; ++iy)
        {
            double* row = &data[2 * ((iz * npad + iy) * npad)];
            gsl_fft_complex_radix2_forward(row, 1, npad);
        }
    }

    for (int iz = 0; iz < npad; ++iz)
    {
        for (int ix = 0; ix < npad; ++ix)
        {
            double* col = &data[2 * (iz * npad * npad + ix)];
            gsl_fft_complex_radix2_forward(col, npad, npad);
        }
    }

    for (int iy = 0; iy < npad; ++iy)
    {
        for (int ix = 0; ix < npad; ++ix)
        {
            double* line = &data[2 * (iy * npad + ix)];
            gsl_fft_complex_radix2_forward(line, npad * npad, npad);
        }
    }

    const double qstep = (npad > 0) ? (2.0 * M_PI / (npad * mdr)) : 0.0;
    for (int iz = 0; iz < npad; ++iz)
    {
        const int kz = (iz <= npad / 2) ? iz : iz - npad;
        const double qz = kz * qstep;
        for (int iy = 0; iy < npad; ++iy)
        {
            const int ky = (iy <= npad / 2) ? iy : iy - npad;
            const double qy = ky * qstep;
            for (int ix = 0; ix < npad; ++ix)
            {
                const int kx = (ix <= npad / 2) ? ix : ix - npad;
                const double qx = kx * qstep;
                const double q = sqrt(qx * qx + qy * qy + qz * qz);
                if ((qmax > 0.0 && q > qmax) || (qmin > 0.0 && q < qmin))
                {
                    const size_t idx = (static_cast<size_t>(iz) * npad + iy) * npad + ix;
                    data[2 * idx] = 0.0;
                    data[2 * idx + 1] = 0.0;
                }
            }
        }
    }

    for (int iy = 0; iy < npad; ++iy)
    {
        for (int ix = 0; ix < npad; ++ix)
        {
            double* line = &data[2 * (iy * npad + ix)];
            gsl_fft_complex_radix2_inverse(line, npad * npad, npad);
        }
    }

    for (int iz = 0; iz < npad; ++iz)
    {
        for (int ix = 0; ix < npad; ++ix)
        {
            double* col = &data[2 * (iz * npad * npad + ix)];
            gsl_fft_complex_radix2_inverse(col, npad, npad);
        }
    }

    for (int iz = 0; iz < npad; ++iz)
    {
        for (int iy = 0; iy < npad; ++iy)
        {
            double* row = &data[2 * ((iz * npad + iy) * npad)];
            gsl_fft_complex_radix2_inverse(row, 1, npad);
        }
    }

    for (int iz = 0; iz < n; ++iz)
    {
        const int sz = (iz - center + npad) % npad;
        for (int iy = 0; iy < n; ++iy)
        {
            const int sy = (iy - center + npad) % npad;
            const size_t base_src = static_cast<size_t>(sz * npad + sy) * npad;
            const size_t base_dst = static_cast<size_t>(iz * n + iy) * n;
            for (int ix = 0; ix < n; ++ix)
            {
                const int sx = (ix - center + npad) % npad;
                const size_t src = base_src + sx;
                const size_t dst = base_dst + ix;
                grid[dst] = data[2 * src];
            }
        }
    }
}

// Protected Methods ---------------------------------------------------------

void ThreeDPDFCalculator::resetValue()
{
    // Ensure odd number of bins so origin is centered.
    mnbins = static_cast<int>(2 * ceil(this->getRmax() / mdr)) + 1;
    size_t total = static_cast<size_t>(mnbins) * mnbins * mnbins;
    mgrid3d.assign(total, 0.0);
    const size_t nr = static_cast<size_t>(ceil(this->getRmax() / mdr)) + 1;
    mradialhistogram3d.assign(nr, 0.0);

    PDFCalculator::resetValue();
    if (!this->isVectorHistogramMode3D())
    {
        this->buildDeltaEligibleShellTable();
        this->validateNNDelta3D();
        this->validateDistanceDecayDelta3D();
    }
    if (mevaluator) mevaluator->setFlag(USEFULLSUM, true);
}

void ThreeDPDFCalculator::addPairContribution(const BaseBondGenerator& bnds, int summationscale)
{
    if (bnds.distance() == 0.0) return;

    int i0 = bnds.site0();
    int i1 = bnds.site1();
    int cntsites = this->countSites();

    if (i0 >= cntsites || i1 >= cntsites)
        return;

    const R3::Vector& rvec = bnds.r01();
    const double pairscale = bnds.multiplicity() * static_cast<double>(summationscale);
    double sfprod = this->sfSite(i0) * this->sfSite(i1) * pairscale;

    if (this->isVectorHistogramMode3D())
    {
        const double weight = (mhistogramweightmode3d == 1) ? pairscale : sfprod;
        addVectorHistogramToGrid(rvec, bnds.distance(), weight);
        return;
    }

    const R3::Matrix& U_i = bnds.Ucartesian0();
    const R3::Matrix& U_j = bnds.Ucartesian1();
    R3::Matrix Sigma = this->effectivePairCovariance(bnds, U_i, U_j);

    addAnisotropicGaussianToGrid(rvec, Sigma, sfprod);
}

void ThreeDPDFCalculator::addVectorHistogramToGrid(const R3::Vector& r_ij, double distance, double weight)
{
    if (weight == 0.0)  return;

    const size_t idx = this->coordToIndex(r_ij[0], r_ij[1], r_ij[2]);
    if (idx < mgrid3d.size())
    {
#ifdef _OPENMP
#pragma omp atomic
#endif
        mgrid3d[idx] += weight;
    }

    if (distance >= 0.0 && !mradialhistogram3d.empty())
    {
        const int ir = static_cast<int>(std::lround(distance / mdr));
        if (ir >= 0 && ir < static_cast<int>(mradialhistogram3d.size()))
        {
#ifdef _OPENMP
#pragma omp atomic
#endif
            mradialhistogram3d[static_cast<size_t>(ir)] += weight;
        }
    }
}

void ThreeDPDFCalculator::addAnisotropicGaussianToGrid(const R3::Vector& r_ij, const R3::Matrix& Sigma, double sfprod)
{
    // Eigenvalue decomposition
    R3::Vector eigenvalues;
    R3::Matrix eigenvectors;

    R3::eigen_solve_3x3(Sigma, eigenvalues, eigenvectors);

    // Check for positive definiteness (eigenvalues are sorted ascending)
    if (eigenvalues[0] <= 1e-8) {
        return;
    }

    // 3. Invert the covariance matrix
    R3::Matrix Sigma_inv = R3::inverse(Sigma);
    double det_Sigma = R3::determinant(Sigma);

    // 4. Normalization factor
    const double two_pi = 2.0 * M_PI;
    double norm_factor = sfprod / (pow(two_pi, 1.5) * sqrt(det_Sigma));

    // 5. Determine sampling bounding box
    double max_sigma = 0.0;
    for(int k=0; k<3; ++k) {
        max_sigma = std::max(max_sigma, sqrt(eigenvalues[k]));
    }

    // 4.0 sigma cutoff
    double cutoff_radius = 4.0 * max_sigma;

    // Safety clamps
    if (cutoff_radius < mdr) cutoff_radius = mdr;
    if (cutoff_radius > 10.0) cutoff_radius = 10.0;

    // 6. Iterate over the local grid indices
    double halfspan = (mnbins / 2) * mdr;

    auto get_index_range = [&](double center_val) -> std::pair<int, int> {
        double min_val = center_val - cutoff_radius;
        double max_val = center_val + cutoff_radius;

        // Node-centered grid: points at k * md
        int start = static_cast<int>(ceil((min_val + halfspan) / mdr));
        int end   = static_cast<int>(floor((max_val + halfspan) / mdr));

        start = std::max(0, start);
        end   = std::min(mnbins - 1, end);

        return {start, end};
    };

    std::pair<int, int> xr = get_index_range(r_ij[0]);
    std::pair<int, int> yr = get_index_range(r_ij[1]);
    std::pair<int, int> zr = get_index_range(r_ij[2]);

    double cutoff_sq = 16.0;

    const int bs = std::max(1, maccumblocksize);
    const int nbz = (zr.second - zr.first + bs) / bs;
    const int nby = (yr.second - yr.first + bs) / bs;
    const int nbx = (xr.second - xr.first + bs) / bs;

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::vector<double> local;

#ifdef _OPENMP
#pragma omp for collapse(3) schedule(dynamic, 1)
#endif
        for (int tbz = 0; tbz < nbz; ++tbz)
        {
            for (int tby = 0; tby < nby; ++tby)
            {
                for (int tbx = 0; tbx < nbx; ++tbx)
                {
                    const int bz = zr.first + tbz * bs;
                    const int by = yr.first + tby * bs;
                    const int bx = xr.first + tbx * bs;

                    const int izhi = std::min(zr.second, bz + bs - 1);
                    const int iyhi = std::min(yr.second, by + bs - 1);
                    const int ixhi = std::min(xr.second, bx + bs - 1);

                    const int tz = izhi - bz + 1;
                    const int ty = iyhi - by + 1;
                    const int tx = ixhi - bx + 1;
                    const size_t nloc = static_cast<size_t>(tz) * ty * tx;
                    local.assign(nloc, 0.0);

                    bool any = false;
                    for (int iz = bz; iz <= izhi; ++iz)
                    {
                        const double z_grid = iz * mdr - halfspan;
                        const double dz = z_grid - r_ij[2];
                        for (int iy = by; iy <= iyhi; ++iy)
                        {
                            const double y_grid = iy * mdr - halfspan;
                            const double dy = y_grid - r_ij[1];
                            for (int ix = bx; ix <= ixhi; ++ix)
                            {
                                const double x_grid = ix * mdr - halfspan;
                                const double dx = x_grid - r_ij[0];
                                R3::Vector delta(dx, dy, dz);
                                R3::Vector tmp = R3::mxvecproduct(Sigma_inv, delta);
                                const double mahalanobis_sq = R3::dot(delta, tmp);
                                if (mahalanobis_sq > cutoff_sq) continue;

                                const double val = norm_factor * exp(-0.5 * mahalanobis_sq);
                                const int lz = iz - bz;
                                const int ly = iy - by;
                                const int lx = ix - bx;
                                const size_t lidx = (static_cast<size_t>(lz) * ty + ly) * tx + lx;
                                local[lidx] += val;
                                any = true;
                            }
                        }
                    }

                    if (!any) continue;

#ifdef _OPENMP
#pragma omp critical(threedpdf_tile_reduce)
#endif
                    {
                        for (int lz = 0; lz < tz; ++lz)
                        {
                            const int iz = bz + lz;
                            for (int ly = 0; ly < ty; ++ly)
                            {
                                const int iy = by + ly;
                                for (int lx = 0; lx < tx; ++lx)
                                {
                                    const size_t lidx = (static_cast<size_t>(lz) * ty + ly) * tx + lx;
                                    const double v = local[lidx];
                                    if (v == 0.0) continue;
                                    const int ix = bx + lx;
                                    const size_t gidx = (static_cast<size_t>(iz) * mnbins + iy) * mnbins + ix;
                                    mgrid3d[gidx] += v;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// Private Helpers -----------------------------------------------------------

bool ThreeDPDFCalculator::DeltaBondKey::operator==(const DeltaBondKey& other) const
{
    return site0 == other.site0 && site1 == other.site1 &&
        rx == other.rx && ry == other.ry && rz == other.rz;
}

size_t ThreeDPDFCalculator::DeltaBondKeyHash::operator()(const DeltaBondKey& key) const
{
    size_t seed = 0;
    seed ^= std::hash<int>()(key.site0) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= std::hash<int>()(key.site1) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= std::hash<long>()(key.rx) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= std::hash<long>()(key.ry) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= std::hash<long>()(key.rz) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
}

bool ThreeDPDFCalculator::matchesDeltaPairTypes(
        const std::string& atom0, const std::string& atom1) const
{
    const bool allpairs = (mdeltaPairA3d == "*" && mdeltaPairB3d == "*");
    if (allpairs)  return true;
    const bool forward =
        (mdeltaPairA3d == "*" || atom0 == mdeltaPairA3d) &&
        (mdeltaPairB3d == "*" || atom1 == mdeltaPairB3d);
    const bool reverse =
        (mdeltaPairA3d == "*" || atom1 == mdeltaPairA3d) &&
        (mdeltaPairB3d == "*" || atom0 == mdeltaPairB3d);
    if (forward || reverse)  return true;
    return (atom0 == mdeltaPairA3d && atom1 == mdeltaPairB3d) ||
        (atom0 == mdeltaPairB3d && atom1 == mdeltaPairA3d);
}

ThreeDPDFCalculator::DeltaBondKey
ThreeDPDFCalculator::makeDeltaBondKey(const BaseBondGenerator& bnds) const
{
    const R3::Vector& r = bnds.r01();
    DeltaBondKey key;
    key.site0 = bnds.site0();
    key.site1 = bnds.site1();
    key.rx = static_cast<long>(std::lround(r[0] / mdeltaKeyTolerance3d));
    key.ry = static_cast<long>(std::lround(r[1] / mdeltaKeyTolerance3d));
    key.rz = static_cast<long>(std::lround(r[2] / mdeltaKeyTolerance3d));
    return key;
}

bool ThreeDPDFCalculator::isDeltaEligiblePair(const BaseBondGenerator& bnds) const
{
    if (!menablennDelta3d || mnnDelta3d == 0.0)  return false;
    return mdeltaEligibleBondKeys.find(this->makeDeltaBondKey(bnds)) !=
        mdeltaEligibleBondKeys.end();
}

double ThreeDPDFCalculator::projectedSigmaAlongBond(
        const R3::Matrix& Sigma,
        const R3::Vector& r,
        double distance) const
{
    if (distance <= 0.0)  return 0.0;
    R3::Vector e(r[0] / distance, r[1] / distance, r[2] / distance);
    R3::Vector tmp = R3::mxvecproduct(Sigma, e);
    return R3::dot(e, tmp);
}

double ThreeDPDFCalculator::positiveDeltaBoundAlongBond(
        const R3::Matrix& Sigma,
        const R3::Vector& r,
        double distance) const
{
    if (distance <= 0.0)  return 0.0;
    R3::Vector e(r[0] / distance, r[1] / distance, r[2] / distance);

    R3::Vector eigenvalues;
    R3::Matrix eigenvectors;
    R3::eigen_solve_3x3(Sigma, eigenvalues, eigenvectors);
    if (eigenvalues[0] <= 1.0e-10)
        return 0.0;

    try
    {
        R3::Matrix Sigma_inv = R3::inverse(Sigma);
        R3::Vector tmp = R3::mxvecproduct(Sigma_inv, e);
        const double denom = R3::dot(e, tmp);
        return (std::isfinite(denom) && denom > 0.0) ? (1.0 / denom) : 0.0;
    }
    catch (...)
    {
        return 0.0;
    }
}

bool ThreeDPDFCalculator::isDistanceDecayDeltaActive() const
{
    return menablennDelta3d &&
        (mdistanceDelta1_3d > 0.0 || mdistanceDelta2_3d > 0.0);
}

double ThreeDPDFCalculator::distanceDecayDeltaFraction(double distance) const
{
    if (distance <= 0.0)  return 0.0;
    return mdistanceDelta1_3d / distance +
        mdistanceDelta2_3d / (distance * distance);
}

double ThreeDPDFCalculator::currentNNDeltaUpperBound() const
{
    if (mnnDelta3dUpperBound <= 0.0)  return 0.0;
    return museadpscaleSensitivity3d ?
        (mnnDelta3dUpperBound * madpScale3d) : mnnDelta3dUpperBound;
}

bool ThreeDPDFCalculator::isVectorHistogramMode3D() const
{
    return mcalculationmode3d == 1;
}

void ThreeDPDFCalculator::validateNNDelta3D() const
{
    if (this->isVectorHistogramMode3D())  return;
    if (mnnDelta3d < 0.0)
        throw std::invalid_argument("3D NN delta must be non-negative.");
    if (mnnDelta3d == 0.0 || mdeltaShellRecords.empty())  return;
    const double bound = this->currentNNDeltaUpperBound();
    if (bound <= 0.0)
        throw std::invalid_argument("3D NN delta cannot be positive without a positive-definite shell bound.");
    if (mnnDelta3d >= bound)
        throw std::invalid_argument("3D NN delta exceeds the positive-definite shell bound.");
}

void ThreeDPDFCalculator::validateDistanceDecayDelta3D() const
{
    if (this->isVectorHistogramMode3D())  return;
    if (mdistanceDelta1_3d < 0.0 || mdistanceDelta2_3d < 0.0)
        throw std::invalid_argument("3D distance-decay delta parameters must be non-negative.");
    if (!this->isDistanceDecayDeltaActive() || mdeltaShellRecords.empty())  return;

    for (const DeltaShellRecord& record : mdeltaShellRecords)
    {
        if (record.distance <= 0.0)  continue;
        const double f = this->distanceDecayDeltaFraction(record.distance);
        if (f <= 0.0)  continue;
        if (record.projectedSigma <= 0.0 || record.positiveDeltaBound <= 0.0)
            continue;
        const double delta = record.projectedSigma * f;
        const double bound = mnnDeltaPositiveEta3d * record.positiveDeltaBound;
        if (bound <= 0.0 || delta >= bound)
            throw std::invalid_argument("3D distance-decay delta violates the pair positive-definite covariance bound.");
    }
}

void ThreeDPDFCalculator::buildDeltaEligibleShellTable()
{
    mdeltaEligibleBondKeys.clear();
    mdeltaShellRecords.clear();
    mnnDelta3dUpperBound = 0.0;

    const StructureAdapterPtr& structure = this->getStructure();
    if (!structure)  return;

    BaseBondGeneratorPtr bnds = structure->createBondGenerator();
    this->configureBondGenerator(*bnds);
    const int cntsites = structure->countSites();

    for (int i0 = 0; i0 < cntsites; ++i0)
    {
        bnds->selectAnchorSite(i0);
        bnds->selectSiteRange(0, cntsites);
        for (bnds->rewind(); !bnds->finished(); bnds->next())
        {
            if (bnds->distance() == 0.0)  continue;
            const std::string& atom0 = structure->siteAtomType(bnds->site0());
            const std::string& atom1 = structure->siteAtomType(bnds->site1());
            if (!this->matchesDeltaPairTypes(atom0, atom1))  continue;

            const R3::Vector& r = bnds->r01();
            R3::Matrix Sigma = bnds->Ucartesian0() + bnds->Ucartesian1();
            if (museadpscaleSensitivity3d)
            {
                for (int row = 0; row < 3; ++row)
                    for (int col = 0; col < 3; ++col)
                        Sigma(row, col) *= madpScale3d;
            }
            const double projected = this->projectedSigmaAlongBond(
                    Sigma, r, bnds->distance());
            const double positive_bound = this->positiveDeltaBoundAlongBond(
                    Sigma, r, bnds->distance());

            DeltaShellRecord record;
            record.site0 = bnds->site0();
            record.site1 = bnds->site1();
            record.r01x = r[0];
            record.r01y = r[1];
            record.r01z = r[2];
            record.distance = bnds->distance();
            record.projectedSigma = projected;
            record.positiveDeltaBound = positive_bound;
            mdeltaShellRecords.push_back(record);
        }
    }

    std::vector<size_t> order(mdeltaShellRecords.size());
    for (size_t i = 0; i < order.size(); ++i)  order[i] = i;
    std::sort(order.begin(), order.end(),
            [&](size_t a, size_t b) {
                return mdeltaShellRecords[a].distance < mdeltaShellRecords[b].distance;
            });

    int shell = 0;
    double last_distance = 0.0;
    bool first = true;
    double min_positive_bound = std::numeric_limits<double>::infinity();
    for (size_t idx : order)
    {
        const double d = mdeltaShellRecords[idx].distance;
        if (first || d - last_distance > mdeltaShellTolerance3d)
        {
            ++shell;
            first = false;
        }
        last_distance = d;
        if (shell == mdeltaShellIndex3d)
        {
            if (mdeltaShellRecords[idx].positiveDeltaBound < min_positive_bound)
                min_positive_bound = mdeltaShellRecords[idx].positiveDeltaBound;
            DeltaBondKey key;
            key.site0 = mdeltaShellRecords[idx].site0;
            key.site1 = mdeltaShellRecords[idx].site1;
            key.rx = static_cast<long>(std::lround(
                        mdeltaShellRecords[idx].r01x / mdeltaKeyTolerance3d));
            key.ry = static_cast<long>(std::lround(
                        mdeltaShellRecords[idx].r01y / mdeltaKeyTolerance3d));
            key.rz = static_cast<long>(std::lround(
                        mdeltaShellRecords[idx].r01z / mdeltaKeyTolerance3d));
            mdeltaEligibleBondKeys.insert(key);
        }
    }

    if (min_positive_bound < std::numeric_limits<double>::infinity() &&
        min_positive_bound > 0.0)
    {
        mnnDelta3dUpperBound = mnnDeltaPositiveEta3d * min_positive_bound;
    }
}

R3::Matrix ThreeDPDFCalculator::effectivePairCovariance(
        const BaseBondGenerator& bnds,
        const R3::Matrix& U_i,
        const R3::Matrix& U_j) const
{
    R3::Matrix Sigma = U_i + U_j;
    if (museadpscaleSensitivity3d)
    {
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col)
                Sigma(row, col) *= madpScale3d;
    }

    const double d = bnds.distance();
    if (d <= 0.0)  return Sigma;

    R3::Vector e(bnds.r01()[0] / d, bnds.r01()[1] / d, bnds.r01()[2] / d);
    double delta_to_subtract = 0.0;
    if (this->isDistanceDecayDeltaActive())
    {
        const StructureAdapterPtr& structure = this->getStructure();
        if (!structure)  return Sigma;
        const std::string& atom0 = structure->siteAtomType(bnds.site0());
        const std::string& atom1 = structure->siteAtomType(bnds.site1());
        if (!this->matchesDeltaPairTypes(atom0, atom1))  return Sigma;

        const double fraction = this->distanceDecayDeltaFraction(d);
        if (fraction <= 0.0)  return Sigma;
        const double projected = this->projectedSigmaAlongBond(Sigma, bnds.r01(), d);
        delta_to_subtract = projected * fraction;
    }
    else
    {
        if (!this->isDeltaEligiblePair(bnds))  return Sigma;
        delta_to_subtract = mnnDelta3d;
    }

    if (delta_to_subtract <= 0.0)  return Sigma;

    const double positive_bound = this->positiveDeltaBoundAlongBond(Sigma, bnds.r01(), d);
    if (positive_bound <= 0.0)  return Sigma;

    const double local_bound = mnnDeltaPositiveEta3d * positive_bound;
    if (local_bound <= 0.0 || delta_to_subtract >= local_bound)
        throw std::invalid_argument("3D delta violates the pair positive-definite covariance bound.");

    for (int row = 0; row < 3; ++row)
    {
        for (int col = 0; col < 3; ++col)
        {
            Sigma(row, col) -= delta_to_subtract * e[row] * e[col];
        }
    }
    return Sigma;
}

size_t ThreeDPDFCalculator::coordToIndex(double x, double y, double z) const
{
    double halfspan = (mnbins / 2) * mdr;
    if (fabs(x) > halfspan || fabs(y) > halfspan || fabs(z) > halfspan)
        return mgrid3d.size(); // out of bounds

    int ix = static_cast<int>(std::lround((x + halfspan) / mdr));
    int iy = static_cast<int>(std::lround((y + halfspan) / mdr));
    int iz = static_cast<int>(std::lround((z + halfspan) / mdr));

    ix = std::max(0, std::min(ix, mnbins - 1));
    iy = std::max(0, std::min(iy, mnbins - 1));
    iz = std::max(0, std::min(iz, mnbins - 1));

    return static_cast<size_t>((iz * mnbins + iy) * mnbins + ix);
}

void ThreeDPDFCalculator::indexToCoord(size_t idx, double& x, double& y, double& z) const
{
    size_t iz = idx / (static_cast<size_t>(mnbins) * mnbins);
    size_t rem = idx % (static_cast<size_t>(mnbins) * mnbins);
    size_t iy = rem / mnbins;
    size_t ix = rem % mnbins;

    double halfspan = (mnbins / 2) * mdr;
    x = ix * mdr - halfspan;
    y = iy * mdr - halfspan;
    z = iz * mdr - halfspan;
}

}   // namespace srreal
}   // namespace diffpy

#include <diffpy/serialization.ipp>
DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::ThreeDPDFCalculator)
