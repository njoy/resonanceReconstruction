#ifndef NJOY_R2_RMATRIX_LEGACY_UNRESOLVED_RESONANCETABLE
#define NJOY_R2_RMATRIX_LEGACY_UNRESOLVED_RESONANCETABLE

// system includes
#include <vector>

// other includes
#include "Log.hpp"
#include "range/v3/range/conversion.hpp"
#include "resonanceReconstruction/Quantity.hpp"
#include "resonanceReconstruction/rmatrix/Table.hpp"
#include "resonanceReconstruction/rmatrix/legacy/ResonanceTableBase.hpp"
#include "resonanceReconstruction/rmatrix/legacy/unresolved/Degrees.hpp"
#include "resonanceReconstruction/rmatrix/legacy/unresolved/Resonance.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {
namespace legacy {
namespace unresolved {

/**
 *  @class
 *  @brief Unresolved resonance parameters for a specific l,J value
 *
 *  Since it is required to interpolate on the resonance parameters, this
 *  class also stores interpolation tables to provide unresolved resonances
 *  at any energy in the unresolved resonance region.
 *
 *  @todo since C++17 does not allow us to decltype() something using lambdas,
 *        a ranges::views::transformed table cannot be created and stored as
 *        a member of the class. It could be done in the call function but then
 *        the interpolation table would be created every time.
 *
 *  @todo we currently only assume linear interpolation on the parameters
 */
class ResonanceTable : protected ResonanceTableBase< Resonance > {

  /* aliases */
  using LevelSpacingTable = LinLinTable< Energy, LevelSpacing >;
  using ReducedWidthTable = LinLinTable< Energy, ReducedWidth >;
  using WidthTable = LinLinTable< Energy, Width >;

  /* fields */
  Degrees degrees_;
  LevelSpacingTable level_spacing_table_;
  ReducedWidthTable elastic_table_;
  WidthTable capture_table_;
  WidthTable fission_table_;
  WidthTable competition_table_;

  /* auxiliary functions */
  #include "resonanceReconstruction/rmatrix/legacy/unresolved/ResonanceTable/src/verifyTable.hpp"
  #include "resonanceReconstruction/rmatrix/legacy/unresolved/ResonanceTable/src/make.hpp"

public:

  /* methods */

  using ResonanceTableBase::numberResonances;
  using ResonanceTableBase::resonances;
  using ResonanceTableBase::energies;

  /**
   *  @brief Return the degrees of freedom for each channel
   */
  const Degrees& degreesOfFreedom() const { return this->degrees_; }

  #include "resonanceReconstruction/rmatrix/legacy/unresolved/ResonanceTable/src/call.hpp"

  /* constructor */
  #include "resonanceReconstruction/rmatrix/legacy/unresolved/ResonanceTable/src/ctor.hpp"
};

} // unresolved namespace
} // legacy namespace
} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
