// physical quantities used in resonance reconstruction
#include "resonanceReconstruction/Quantity.hpp"

#ifndef NJOY_RESONANCE_RECONSTRUCTION
#define NJOY_RESONANCE_RECONSTRUCTION

#include <complex>
#include <iterator>
#include <memory>

#include "interpolation.hpp"
#include "elementary.hpp"

#include "range/v3/distance.hpp"
#include "range/v3/iterator_range.hpp"
#include "range/v3/to_container.hpp"
#include "range/v3/action/sort.hpp"
#include "range/v3/action/unique.hpp"
#include "range/v3/algorithm/count_if.hpp"
#include "range/v3/algorithm/for_each.hpp"
#include "range/v3/algorithm/find_if.hpp"
#include "range/v3/numeric/accumulate.hpp"
#include "range/v3/utility/iterator.hpp"
#include "range/v3/view/all.hpp"
#include "range/v3/view/cartesian_product.hpp"
#include "range/v3/view/chunk.hpp"
#include "range/v3/view/concat.hpp"
#include "range/v3/view/filter.hpp"
#include "range/v3/view/group_by.hpp"
#include "range/v3/view/indices.hpp"
#include "range/v3/view/iota.hpp"
#include "range/v3/view/join.hpp"
#include "range/v3/view/repeat_n.hpp"
#include "range/v3/view/single.hpp"
#include "range/v3/view/transform.hpp"
#include "range/v3/view/zip.hpp"
#include "range/v3/view/zip_with.hpp"

// ENDF components
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
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wextra"
#include <Eigen/Dense>
#pragma GCC diagnostic pop

namespace njoy {
namespace resonanceReconstruction {

using Matrix3x3 = Eigen::Matrix3cd;

template< int i >
using Integer = std::integral_constant< int, i >;

#include "resonanceReconstruction/src/radius.hpp"
#include "resonanceReconstruction/src/channelRadius.hpp"
#include "resonanceReconstruction/src/root.hpp"
#include "resonanceReconstruction/src/neutronWaveNumber.hpp"
#include "resonanceReconstruction/src/penetrationShift.hpp"
#include "resonanceReconstruction/src/phaseShift.hpp"

#include "resonanceReconstruction/pack.hpp"
#include "resonanceReconstruction/EnergyRange.hpp"
#include "resonanceReconstruction/ZeroWidth.hpp"

#include "resonanceReconstruction/breitWigner.hpp"
#include "resonanceReconstruction/reichMoore.hpp"
}
}

#endif

#include "resonanceReconstruction/rmatrix.hpp"
