#ifndef NJOY_R2_RMATRIX_LEGACY_RESONANCETABLEBASE
#define NJOY_R2_RMATRIX_LEGACY_RESONANCETABLEBASE

// system includes
#include <vector>

// other includes
#include "range/v3/view/transform.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {
namespace legacy {

/**
 *  @class
 *  @brief Base interface for a table of resonances for a specific l,J value
 */
template < typename ResonanceType > class ResonanceTableBase {

  /* fields */
  std::vector< ResonanceType > resonances_;

public:

  /* constructor */

  /**
   *  @brief Constructor
   *
   *  @param[in] resonances   the resonances (ne values)
   */
  ResonanceTableBase( std::vector< ResonanceType >&& resonances ) :
      resonances_( resonances ) {}

  /* methods */

  /**
   *  @brief Return the number of resonances
   */
  unsigned int numberResonances() const { return this->resonances_.size(); }

  /**
   *  @brief Return the resonances
   */
  const std::vector< ResonanceType >& resonances() const {

    return this->resonances_;
  }

  /**
   *  @brief Return the resonance energies
   */
  auto energies() const {

    return this->resonances()
             | ranges::cpp20::views::transform(
                   [] ( const auto& resonance ) -> decltype(auto)
                      { return resonance.energy(); } );
  }
};

} // legacy namespace
} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
