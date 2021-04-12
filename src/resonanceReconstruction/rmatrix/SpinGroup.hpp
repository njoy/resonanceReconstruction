#ifndef NJOY_R2_RMATRIX_SPINGROUP
#define NJOY_R2_RMATRIX_SPINGROUP

// system includes
#include <algorithm>
#include <complex>
#include <vector>

// other includes
#include "Log.hpp"
#include "range/v3/range/conversion.hpp"
#include "range/v3/iterator/operations.hpp"
#include "range/v3/action/sort.hpp"
#include "range/v3/action/unique.hpp"
#include "range/v3/algorithm/count_if.hpp"
#include "range/v3/algorithm/for_each.hpp"
#include "range/v3/algorithm/find_if.hpp"
#include "range/v3/view/subrange.hpp"
#include "range/v3/numeric/accumulate.hpp"
#include "range/v3/view/all.hpp"
#include "range/v3/view/concat.hpp"
#include "range/v3/view/filter.hpp"
#include "range/v3/view/repeat_n.hpp"
#include "range/v3/view/single.hpp"
#include "range/v3/view/transform.hpp"
#include "range/v3/view/zip_with.hpp"
#include "range/v3/view/zip.hpp"
#include "resonanceReconstruction/Quantity.hpp"
#include "resonanceReconstruction/rmatrix/Map.hpp"
#include "resonanceReconstruction/rmatrix/ReactionChannelID.hpp"
#include "resonanceReconstruction/rmatrix/Formalism.hpp"
#include "resonanceReconstruction/rmatrix/BoundaryOption.hpp"
#include "resonanceReconstruction/rmatrix/ParticleChannel.hpp"
#include "resonanceReconstruction/rmatrix/ParticleChannelData.hpp"
#include "resonanceReconstruction/rmatrix/ResonanceTable.hpp"
#include "resonanceReconstruction/rmatrix/RLMatrixCalculator.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {

/**
 *  @class
 *  @brief A spin group corresponding to a Jpi quantum number set
 */
template < typename Formalism, typename BoundaryOption >
class SpinGroup {

  /* fields */
  RLMatrixCalculator< Formalism, BoundaryOption > rlmatrix_;

  std::vector< ReactionID > reactions_;
  std::vector< unsigned int > incident_;
  std::vector< ParticleChannel > channels_;
  ResonanceTable parameters_;

  /* auxiliary functions */
  #include "resonanceReconstruction/rmatrix/SpinGroup/src/makeChannels.hpp"
  #include "resonanceReconstruction/rmatrix/SpinGroup/src/makeResonanceTable.hpp"

  #include "resonanceReconstruction/rmatrix/SpinGroup/src/makeReactionIdentifiers.hpp"
  #include "resonanceReconstruction/rmatrix/SpinGroup/src/determineIncidentChannels.hpp"

  #include "resonanceReconstruction/rmatrix/SpinGroup/src/penetrabilities.hpp"
  #include "resonanceReconstruction/rmatrix/SpinGroup/src/phaseShifts.hpp"
  #include "resonanceReconstruction/rmatrix/SpinGroup/src/coulombShifts.hpp"
  #include "resonanceReconstruction/rmatrix/SpinGroup/src/sqrtPenetrabilities.hpp"
  #include "resonanceReconstruction/rmatrix/SpinGroup/src/omegas.hpp"

  #include "resonanceReconstruction/rmatrix/SpinGroup/src/verifyChannels.hpp"
  #include "resonanceReconstruction/rmatrix/SpinGroup/src/verifyIncidentChannels.hpp"
  #include "resonanceReconstruction/rmatrix/SpinGroup/src/verifyResonanceChannels.hpp"

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/SpinGroup/src/ctor.hpp"

  /**
   *  @brief Return the channels in the spin group
   */
  auto channels() const {

    return ranges::cpp20::views::all( this->channels_ );
  }

  /**
   *  @brief Return the current incident channels in the spin group
   */
  auto incidentChannels() const {

    return ranges::cpp20::views::all( this->incident_ )
             | ranges::views::transform( [&] ( const unsigned int i )
                                             { return this->channels_[i]; } ); }

  /**
   *  @brief Return the current incident particle pair
   */
  const ParticlePair incidentPair() const {

    auto incidentParticlePair = [] ( const auto& channel )
                                   { return channel.incidentParticlePair(); };

    return std::visit( incidentParticlePair, this->channels_.front() );
  }

  /**
   *  @brief Return the channel identifiers associated to each channel
   *
   *  These channel identifiers are unique: they reference a given particle pair
   *  (which may or may not be unique) and the quantum numbers of the channel.
   *  The order in which these are given equals the order of the channels in the
   *  spin group.
   */
  auto channelIDs() const {

    auto channelID = [] ( const auto& channel )
                        { return channel.channelID(); };

    return this->channels()
             | ranges::views::transform(
                   [=] ( const auto& channel )
                       { return std::visit( channelID, channel ); } );
  }

  /**
   *  @brief Return the reaction identifiers associated to each channel
   *
   *  These reaction identifiers are not unique: a given reaction identifier
   *  may appear multiple times because multiple channels can contribute to the
   *  same reactions (e.g. multiple fission channels). The order in which these
   *  are given equals the order of the channels in the spin group.
   */
  auto reactionIDs() const {

    return ranges::cpp20::views::all( this->reactions_ );
  }

  /**
   *  @brief Return the resonance table
   */
  const ResonanceTable& resonanceTable() const { return this->parameters_; }

  //#include "resonanceReconstruction/rmatrix/SpinGroup/src/switchIncidentPair.hpp"
  #include "resonanceReconstruction/rmatrix/SpinGroup/src/evaluate.hpp"
  #include "resonanceReconstruction/rmatrix/SpinGroup/src/evaluateTMatrix.hpp"
  #include "resonanceReconstruction/rmatrix/SpinGroup/src/grid.hpp"
};

} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
