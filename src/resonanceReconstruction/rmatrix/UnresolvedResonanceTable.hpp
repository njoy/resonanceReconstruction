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
class UnresolvedResonanceTable :
  protected BaseResonanceTable< UnresolvedResonance > {

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
  using LevelSpacingTable = Table< Energy, Energy >;
  using WidthTable = Table< Energy, ReducedWidth >;

  /* fields */
  std::vector< unsigned int > degrees_;
  LevelSpacingTable level_spacing_table_;
  std::vector< WidthTable > width_tables_;

  /* auxiliary functions */
  #include "resonanceReconstruction/rmatrix/UnresolvedResonanceTable/src/verifyTable.hpp"
  #include "resonanceReconstruction/rmatrix/UnresolvedResonanceTable/src/make.hpp"

public:

  /* methods */
  using BaseResonanceTable::numberChannels;
  using BaseResonanceTable::numberResonances;
  using BaseResonanceTable::channels;
  using BaseResonanceTable::resonances;
  using BaseResonanceTable::energies;

  /**
   *  @brief Return the degrees of freedom for each channel
   */
  auto degreesOfFreedom() const { return ranges::view::all( this->degrees_ ); }

  #include "resonanceReconstruction/rmatrix/UnresolvedResonanceTable/src/call.hpp"

  /* constructor */
  #include "resonanceReconstruction/rmatrix/UnresolvedResonanceTable/src/ctor.hpp"
};
