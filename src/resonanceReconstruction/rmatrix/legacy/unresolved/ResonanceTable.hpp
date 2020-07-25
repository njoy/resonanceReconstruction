/**
 *  @class
 *  @brief Unresolved resonance parameters for a specific J,pi or l,J value
 *
 *  Since it is required to interpolate on the resonance parameters, this
 *  class also stores interpolation tables to provide unresolved resonances
 *  at any energy in the unresolved resonance region.
 *
 *  @todo since C++17 does not allow us to decltype() something using lambdas,
 *        a ranges::view::transformed table cannot be created and stored as
 *        a member of the class. It could be done in the call function but then
 *        the interpolation table would be created every time.
 *
 *  @todo we currently only assume linear interpolation on the parameters
 */
class ResonanceTable {

  /* aliases */
  template < typename XType, typename YType >
  using Table = interpolation::Table<
                  interpolation::table::Type<
                    interpolation::LinearLinear,
                    interpolation::table::search::Binary,
			              interpolation::table::discontinuity::TakeLeft,
                		std::vector< XType >, std::vector< YType > >,
                  interpolation::table::right::interval::Throws,
                  interpolation::table::left::interval::Throws >;
  using LevelSpacingTable = Table< Energy, LevelSpacing >;
  using ReducedWidthTable = Table< Energy, ReducedWidth >;
  using WidthTable = Table< Energy, Width >;

  /* fields */
  std::vector< Resonance > widths_;
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

  /**
   *  @brief Return the number of resonances
   */
  unsigned int numberResonances() const { return this->widths_.size(); }

  /**
   *  @brief Return the resonances
   */
  auto resonances() const { return ranges::view::all( this->widths_ ); }

  /**
   *  @brief Return the resonance energies
   */
  auto energies() const {

    return this->resonances()
             | ranges::view::transform( [] ( const auto& resonance )
                                           { return resonance.energy(); } );
  }

  /**
   *  @brief Return the degrees of freedom for each channel
   */
  const Degrees& degreesOfFreedom() const { return this->degrees_; }

  #include "resonanceReconstruction/rmatrix/legacy/unresolved/ResonanceTable/src/call.hpp"

  /* constructor */
  #include "resonanceReconstruction/rmatrix/legacy/unresolved/ResonanceTable/src/ctor.hpp"
};
