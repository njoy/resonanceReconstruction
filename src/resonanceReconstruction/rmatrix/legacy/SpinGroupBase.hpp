#ifndef NJOY_R2_RMATRIX_LEGACY_SPINGROUPBASE
#define NJOY_R2_RMATRIX_LEGACY_SPINGROUPBASE

// system includes

// other includes
#include "resonanceReconstruction/rmatrix/ChannelType.hpp"
#include "resonanceReconstruction/rmatrix/Channel.hpp"
#include "resonanceReconstruction/rmatrix/ReactionID.hpp"
#include "range/v3/algorithm/count_if.hpp"
#include "range/v3/view/all.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {
namespace legacy {

/**
 *  @class
 *  @brief Base interface for a legacy l,J spin group
 */
template < typename ResonanceTableType > class SpinGroupBase {

  /* fields */
  Channel< Neutron > incident_;
  ResonanceTableType table_;
  std::vector< ReactionID > reactions_;

  /* auxiliary functions */
  #include "resonanceReconstruction/rmatrix/legacy/SpinGroupBase/src/makeReactionIdentifiers.hpp"

protected:

  const ReactionID& elasticID() const { return this->reactions_[0]; }
  const ReactionID& captureID() const { return this->reactions_[1]; }
  const ReactionID& fissionID() const { return this->reactions_[2]; }

public:

  /* constructor */
  #include "resonanceReconstruction/rmatrix/legacy/SpinGroupBase/src/ctor.hpp"

  /* methods */

  /**
   *  @brief Return the incident channel data
   *
   *  REMARK: not all information in the channel will be correct (channel spin,
   *          target parity, target charge, etc. since this data is not
   *          required for the legacy resolved/unresolved resonances)
   */
  const Channel< Neutron >& incidentChannel() const { return this->incident_; }

  /**
   *  @brief Return the orbital angular momentum l for this l,J pair
   */
  const OrbitalAngularMomentum& orbitalAngularMomentum() const {

    return this->incidentChannel().quantumNumbers().orbitalAngularMomentum();
  }

  /**
   *  @brief Return the total angular momentum J for this l,J pair
   */
  const TotalAngularMomentum& totalAngularMomentum() const {

    return this->incidentChannel().quantumNumbers().totalAngularMomentum();
  }

  /**
   *  @brief Return the resonance table
   */
  const ResonanceTableType& resonanceTable() const { return this->table_; }

  /**
   *  @brief Return the reactions defined here
   */
  auto reactionIDs() const {

    return ranges::cpp20::views::all( this->reactions_ ); 
  }

  /**
   *  @brief Return whether or not the spin group has fission or not
   */
  bool hasFission() const { return this->reactions_.size() == 3; }
};

} // legacy namespace
} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
