/*****************************************************************************
 *
 * libdiffpy         by DANSE Diffraction group
 *                   Simon J. L. Billinge
 *                   (c) 2010 The Trustees of Columbia University
 *                   in the City of New York.  All rights reserved.
 *
 * File coded by:    Pavol Juhas
 *
 * See AUTHORS.rst for a list of people who contributed.
 * See LICENSE_DANSE.rst for license information.
 *
 ******************************************************************************
 *
 * class BaseDebyeSum -- base class for concrete Debye sum calculators
 *
 *****************************************************************************/

#ifndef BASEDEBYESUM_HPP_INCLUDED
#define BASEDEBYESUM_HPP_INCLUDED

#include <cstdint>
#include <map>
#include <string>
#include <utility>

#include <diffpy/srreal/PairQuantity.hpp>
#include <diffpy/srreal/PeakWidthModel.hpp>
#include <diffpy/srreal/PDFUtils.hpp>

#include <diffpy/Export.hpp>

namespace diffpy {
namespace srreal {

class DLL_EXPORT BaseDebyeSum : public PairQuantity,
                                public PeakWidthModelOwner {
 public:
  // constructor
  BaseDebyeSum();

  // PairQuantity overloads
  virtual eventticker::EventTicker& ticker() const;

  // results
  /// F values on a full Q-grid starting at 0
  QuantityType getF() const;

  // Q-range methods
  /// Full Q-grid starting at 0
  QuantityType getQgrid() const;
  // Q-range configuration
  void setQmin(double);
  const double& getQmin() const;
  void setQmax(double);
  const double& getQmax() const;
  void setQstep(double);
  const double& getQstep() const;

  // Summation cutoff due to Q-dependent pair scaling or large distance
  /// set relative cutoff value for Debye sum contribution
  void setDebyePrecision(double);
  /// return relative cutoff value for Debye sum contribution
  const double& getDebyePrecision() const;

 protected:
  // PairQuantity overloads
  virtual void resetValue();
  virtual void addPairContribution(const BaseBondGenerator&, int);
  virtual bool addAllPairContributions();
  /// addAllPairContributions for derived classes that keep the
  /// addPairContribution of this class.
  bool addAllPairsOfDebyeSum();
  // support for PQEvaluatorOptimized
  virtual void stashPartialValue();
  virtual void restorePartialValue();

  // own methods
  virtual double sfSiteAtQ(int, const double& Q) const;

 private:
  // methods
  double sfAverageAtkQ(int kq) const;
  /// cache structure factors data for a quick access during summation
  void cacheStructureData();

  // data
  // configuration
  double mqmin;
  double mqmax;
  double mqstep;
  double mdebyeprecision;
  struct SiteTypes {
    std::vector<int> typeofsite;
    // A scattering type names an elemental curve and a site occupancy.
    std::vector<int> curve;
    std::vector<double> occupancy;
  };
  struct {
    SiteTypes types;
    QuantityType factors;  // consecutive curves on the Q grid
    int nq = 0;
    QuantityType sfaverageatkq;
    double totaloccupancy;
  } mstructure_cache;
  QuantityType mdbsumstash;
  std::map<std::string, int> melements;
  eventticker::EventTicker msfticker;

  using SiteKey = std::pair<std::string, std::uint64_t>;
  // Representative curves identify changes to virtual, site-specific factors.
  std::map<SiteKey, int> mcustomTypes;
  bool mcustomAmbiguous = false;
  void cacheCustomFactors(SiteTypes&, int nq);
  void cacheElementFactors(SiteTypes&, int nq);
  std::vector<QuantityType> archiveFactors() const;
  void loadFactors(const std::vector<QuantityType>& curves);

  // serialization
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    using boost::serialization::base_object;
    ar& base_object<PairQuantity>(*this);
    ar& base_object<PeakWidthModelOwner>(*this);
    ar & mqmin;
    ar & mqmax;
    ar & mqstep;
    ar & mdebyeprecision;
    ar & mstructure_cache.types.typeofsite;
    // Preserve the archive representation used before elemental curves.
    std::vector<QuantityType> curves;
    if constexpr (Archive::is_saving::value) curves = this->archiveFactors();
    ar & curves;
    if constexpr (Archive::is_loading::value) this->loadFactors(curves);
    ar & mstructure_cache.sfaverageatkq;
    ar & mstructure_cache.totaloccupancy;
  }

};  // class BaseDebyeSum

}  // namespace srreal
}  // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_CLASS_EXPORT_KEY(diffpy::srreal::BaseDebyeSum)

#endif  // BASEDEBYESUM_HPP_INCLUDED
