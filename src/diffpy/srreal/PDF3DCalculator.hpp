#ifndef PDF3DCALCULATOR_HPP_INCLUDED
#define PDF3DCALCULATOR_HPP_INCLUDED

#include <diffpy/srreal/PDFCalculator.hpp>
#include <vector>
#include <string>
#include <unordered_set>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <diffpy/srreal/R3linalg.hpp>

namespace diffpy {
namespace srreal {

class PDF3DCalculator : public PDFCalculator
{
public:
    // constructor
    PDF3DCalculator();

    // Public interface to retrieve 3D PDF data
    // Returns a flat vector of [x0, y0, z0, G0, x1, y1, z1, G1, ...]
    QuantityType get3DPDF() const;

    // Export dense 3D grid (nz, ny, nx) directly to binary file.
    // usefloat32=true writes float32, otherwise float64.
    // applypost=true applies q-window/rdf_scale/rho0/qdamp before export.
    void exportGrid3DBinary(const std::string& path, bool usefloat32 = true, bool applypost = true) const;
    QuantityType getRadialHistogram3D() const;

    // Grid configuration
    void setGridStep(double dr);
    double getGridStep() const;

    // Blocked accumulation configuration
    void setAccumBlockSize(int bs);
    int getAccumBlockSize() const;

    void setApplyRho0Background3D(bool);
    bool getApplyRho0Background3D() const;

    void setUseCQWindow3D(bool v);
    bool getUseCQWindow3D() const;

    // Nearest-neighbor correlated-motion correction for anisotropic 3D PDF models.
    // Boolean switches are also exposed as double attributes with 0/1
    // values because libdiffpy's generic attribute system is double-only.
    void setEnableNNDelta3D(bool);
    bool getEnableNNDelta3D() const;
    double getEnableNNDelta3DAttr() const;
    void setEnableNNDelta3DAttr(double);

    void setNNDelta3D(double);
    const double& getNNDelta3D() const;
    double getNNDelta3DUpperBound() const;
    const double& getNNDeltaPositiveEta3D() const;
    void setNNDeltaPositiveEta3D(double);
    void setDistanceDelta1_3D(double);
    const double& getDistanceDelta1_3D() const;
    void setDistanceDelta2_3D(double);
    const double& getDistanceDelta2_3D() const;

    void setDeltaPairTypes3D(const std::string&, const std::string&);
    const std::string& getDeltaPairA3D() const;
    const std::string& getDeltaPairB3D() const;
    void setDeltaShellIndex3D(int);
    int getDeltaShellIndex3D() const;
    double getDeltaShellIndex3DAttr() const;
    void setDeltaShellIndex3DAttr(double);
    const double& getDeltaShellTolerance3D() const;
    void setDeltaShellTolerance3D(double);
    const double& getDeltaKeyTolerance3D() const;
    void setDeltaKeyTolerance3D(double);
    int getDeltaEligiblePairCount3D() const;

    void setUseADPScaleSensitivity3D(bool);
    bool getUseADPScaleSensitivity3D() const;
    double getUseADPScaleSensitivity3DAttr() const;
    void setUseADPScaleSensitivity3DAttr(double);
    void setADPScale3D(double);
    const double& getADPScale3D() const;

    void setRho0BackgroundScale3D(double);
    const double& getRho0BackgroundScale3D() const;

    void setCalculationMode3D(int);
    int getCalculationMode3D() const;
    double getCalculationMode3DAttr() const;
    void setCalculationMode3DAttr(double);

    void setHistogramWeightMode3D(int);
    int getHistogramWeightMode3D() const;
    double getHistogramWeightMode3DAttr() const;
    void setHistogramWeightMode3DAttr(double);

protected:
    // Override PairQuantity virtual methods
    virtual void resetValue() override;
    virtual void addPairContribution(const BaseBondGenerator& bnds, int) override;

private:
    struct DeltaBondKey
    {
        int site0;
        int site1;
        long rx;
        long ry;
        long rz;
        bool operator==(const DeltaBondKey&) const;
    };

    struct DeltaBondKeyHash
    {
        size_t operator()(const DeltaBondKey&) const;
    };

    struct DeltaShellRecord
    {
        int site0;
        int site1;
        double r01x;
        double r01y;
        double r01z;
        double distance;
        double projectedSigma;
        double positiveDeltaBound;
    };

    // Helper to add anisotropic Gaussian to the grid
    void addAnisotropicGaussianToGrid(const R3::Vector& r_ij, const R3::Matrix& Sigma, double sfprod);
    void addVectorHistogramToGrid(const R3::Vector& r_ij, double distance, double weight);

    // Helper: map 3D position to linear index
    size_t coordToIndex(double x, double y, double z) const;
    void indexToCoord(size_t idx, double& x, double& y, double& z) const;

    void applyPostProcessing3D(std::vector<double>& grid) const;
    void applyQWindow3D(std::vector<double>& grid) const;
    double computeRho0Background() const;
    void buildDeltaEligibleShellTable();
    bool isDeltaEligiblePair(const BaseBondGenerator&) const;
    bool matchesDeltaPairTypes(const std::string&, const std::string&) const;
    DeltaBondKey makeDeltaBondKey(const BaseBondGenerator&) const;
    R3::Matrix effectivePairCovariance(
            const BaseBondGenerator&,
            const R3::Matrix&,
            const R3::Matrix&) const;
    double projectedSigmaAlongBond(
            const R3::Matrix&,
            const R3::Vector&,
            double distance) const;
    double positiveDeltaBoundAlongBond(
            const R3::Matrix&,
            const R3::Vector&,
            double distance) const;
    bool isDistanceDecayDeltaActive() const;
    double distanceDecayDeltaFraction(double distance) const;
    void validateDistanceDecayDelta3D() const;
    void validateNNDelta3D() const;
    double currentNNDeltaUpperBound() const;
    bool isVectorHistogramMode3D() const;

    // Data members
    double mdr;                      // grid spacing (assumes cubic grid centered at origin)
    int mnbins;                      // number of bins per dimension (odd number, center at 0)
    std::vector<double> mgrid3d;     // flattened 3D histogram (size = mnbins^3)
    std::vector<double> mradialhistogram3d;
    int maccumblocksize;              // block size for tiled accumulation loops

    bool mapplyrho0background3d;
    bool musecqwindow3d;
    int mcalculationmode3d;
    int mhistogramweightmode3d;

    bool menablennDelta3d;
    double mnnDelta3d;
    double mnnDelta3dUpperBound;
    double mnnDeltaPositiveEta3d;
    double mdistanceDelta1_3d;
    double mdistanceDelta2_3d;
    std::string mdeltaPairA3d;
    std::string mdeltaPairB3d;
    int mdeltaShellIndex3d;
    double mdeltaShellTolerance3d;
    double mdeltaKeyTolerance3d;
    bool museadpscaleSensitivity3d;
    double madpScale3d;
    double mrho0backgroundscale3d;
    std::unordered_set<DeltaBondKey, DeltaBondKeyHash> mdeltaEligibleBondKeys;
    std::vector<DeltaShellRecord> mdeltaShellRecords;

    // Internal constants
    static constexpr double DEFAULT_RMAX_3D = 10.0;  // Angstrom
    static constexpr double DEFAULT_GRID_STEP = 0.1; // Angstrom

    // serialization
    friend class boost::serialization::access;
    template<class Archive>
        void serialize(Archive& ar, const unsigned int version)
    {
        using boost::serialization::base_object;
        ar & base_object<PDFCalculator>(*this);
        ar & mdr;
        ar & mnbins;
        ar & mgrid3d;
        ar & mradialhistogram3d;
        ar & maccumblocksize;
        ar & mapplyrho0background3d;
        ar & musecqwindow3d;
        ar & mcalculationmode3d;
        ar & mhistogramweightmode3d;
        ar & menablennDelta3d;
        ar & mnnDelta3d;
        ar & mnnDelta3dUpperBound;
        ar & mnnDeltaPositiveEta3d;
        ar & mdistanceDelta1_3d;
        ar & mdistanceDelta2_3d;
        ar & mdeltaPairA3d;
        ar & mdeltaPairB3d;
        ar & mdeltaShellIndex3d;
        ar & mdeltaShellTolerance3d;
        ar & mdeltaKeyTolerance3d;
        ar & museadpscaleSensitivity3d;
        ar & madpScale3d;
        ar & mrho0backgroundscale3d;
    }
};

}   // namespace srreal
}   // namespace diffpy

#endif  // PDF3DCALCULATOR_HPP_INCLUDED
