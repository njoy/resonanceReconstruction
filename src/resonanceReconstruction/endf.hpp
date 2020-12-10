#ifndef NJOY_R2_ENDF
#define NJOY_R2_ENDF

// system includes

// other includes
#include "ENDFtk/section/2/151.hpp"

namespace endf {

  using ScatteringRadius = njoy::ENDFtk::section::Type< 2, 151 >::ScatteringRadius;
  using ResonanceRange = njoy::ENDFtk::section::Type< 2, 151 >::ResonanceRange;
  using BreitWignerLValue = njoy::ENDFtk::section::Type< 2, 151 >::BreitWignerLValue;
  using ReichMooreLValue = njoy::ENDFtk::section::Type< 2, 151 >::ReichMooreLValue;
  using SingleLevelBreitWigner = njoy::ENDFtk::section::Type< 2, 151 >::SingleLevelBreitWigner;
  using MultiLevelBreitWigner = njoy::ENDFtk::section::Type< 2, 151 >::MultiLevelBreitWigner;
  using ReichMoore = njoy::ENDFtk::section::Type< 2, 151 >::ReichMoore;
  using RMatrixLimited = njoy::ENDFtk::section::Type< 2, 151 >::RMatrixLimited;
  using UnresolvedEnergyIndependent = njoy::ENDFtk::section::Type< 2, 151 >::UnresolvedEnergyIndependent;
  using UnresolvedEnergyDependentFissionWidths = njoy::ENDFtk::section::Type< 2, 151 >::UnresolvedEnergyDependentFissionWidths;
  using UnresolvedEnergyDependent = njoy::ENDFtk::section::Type< 2, 151 >::UnresolvedEnergyDependent;

} // endf namespace

#endif
