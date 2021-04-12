#ifndef NJOY_R2_RMATRIX_RECONSTRUCTOR
#define NJOY_R2_RMATRIX_RECONSTRUCTOR

// system includes
#include <variant>
#include <vector>

// other includes
#include "range/v3/view/filter.hpp"
#include "resonanceReconstruction/Quantity.hpp"
#include "resonanceReconstruction/rmatrix/Map.hpp"
#include "resonanceReconstruction/rmatrix/ReactionID.hpp"
#include "resonanceReconstruction/rmatrix/Formalism.hpp"
#include "resonanceReconstruction/rmatrix/BoundaryOption.hpp"
#include "resonanceReconstruction/rmatrix/CompoundSystem.hpp"
#include "resonanceReconstruction/rmatrix/legacy/resolved/CompoundSystem.hpp"
#include "resonanceReconstruction/rmatrix/legacy/unresolved/CompoundSystem.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {

using CompoundSystemVariant =
    std::variant< legacy::resolved::CompoundSystem< SingleLevelBreitWigner >,
                  legacy::resolved::CompoundSystem< MultiLevelBreitWigner >,
                  CompoundSystem< ReichMoore, ShiftFactor >,
                  CompoundSystem< ReichMoore, Constant >,
                  //CompoundSystem< GeneralRMatrix, ShiftFactor >,
                  //CompoundSystem< GeneralRMatrix, Constant >,
                  legacy::unresolved::CompoundSystem >;

/**
 *  @class
 *  @brief Class used to reconstruct cross sections from ENDF resonances
 */
class Reconstructor {

  /* fields */
  Energy lower_;
  Energy upper_;
  CompoundSystemVariant system_;

public:

  /* constructor */
  Reconstructor( const Energy& lower, const Energy& upper,
                 const CompoundSystemVariant& reconstructor ) :
      lower_( lower ), upper_( upper ),
      system_( reconstructor ) {}

  /**
   *  @brief Return the compound system
   */
  CompoundSystemVariant& compoundSystem() { return this->system_; }

  /**
   *  @brief Return whether or not the reconstructor is for resolved resonances
   */
  bool isResolved() {

    return this->system_.index() !=
           std::variant_size_v< CompoundSystemVariant > - 1;
  }

  /**
   *  @brief Return whether or not the reconstructor is for unresolved resonances
   */
  bool isUnresolved() {

    return not this->isResolved();
  }

  /**
   *  @brief Return the lower energy boundary
   */
  const Energy& lowerEnergy() const { return this->lower_; }

  /**
   *  @brief Return the lower energy boundary
   */
  const Energy& upperEnergy() const { return this->upper_; }

  /**
   *  @brief Return the reaction identifiers for this reconstructor
   */
  auto reactionIDs() {

    return std::visit( [] ( const auto& system )
                          { return system.reactionIDs(); },
                       this->compoundSystem() );
  }

  /**
   *  @brief Return the minimal energy grid derived from the resonance
   *         parameters
   */
  std::vector< Energy > grid() const {

    auto filter = [&] ( const auto& energy ) {

      return ( this->lowerEnergy() <= energy ) and
             ( energy <= this->upperEnergy() );
    };

    auto energies = std::visit( [&] ( auto& system ) { return system.grid(); },
                                this->system_ );
    return energies | ranges::cpp20::views::filter( filter );
  }

  /**
   *  @brief Reconstruct the cross sections at the given energy
   *
   *  @param[in] energy   the energy for which cross sections are to be
   *                      calculated
   */
  Map< ReactionID, CrossSection > operator()( const Energy& energy ) {

    Map< ReactionID, CrossSection > result;
    if ( ( energy >= this->lowerEnergy() ) and
         ( energy <= this->upperEnergy() ) ) {

      std::visit( [&] ( auto& system )
                      { system.evaluate( energy, result ); },
                  this->system_ );
    }
    return result;
  }

  /**
   *  @brief Return the interpolation scheme
   */
  std::optional< int > interpolation() const {

    return std::visit(
             njoy::utility::overload{
                 [&] ( const legacy::unresolved::CompoundSystem& system )
                     -> std::optional< int >
                     { return std::make_optional( system.interpolation() ); },
                 [&] ( const auto& )
                     -> std::optional< int >
                     { return std::nullopt; } },
             this->system_ );
  }
};

} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
