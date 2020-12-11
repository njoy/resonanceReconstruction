#ifndef NJOY_R2_RMATRIX_FROMENDF
#define NJOY_R2_RMATRIX_FROMENDF

// system includes

// other includes
#include "resonanceReconstruction/endf.hpp"
#include "resonanceReconstruction/Quantity.hpp"
#include "resonanceReconstruction/rmatrix/Formalism.hpp"
#include "resonanceReconstruction/rmatrix/BoundaryOption.hpp"
#include "resonanceReconstruction/rmatrix/CompoundSystem.hpp"
#include "resonanceReconstruction/rmatrix/legacy/resolved/CompoundSystem.hpp"
#include "resonanceReconstruction/rmatrix/legacy/unresolved/CompoundSystem.hpp"
#include "resonanceReconstruction/rmatrix/possibleChannelSpinValues.hpp"
#include "resonanceReconstruction/rmatrix/possibleChannelTotalAngularMomentumValues.hpp"
#include "resonanceReconstruction/rmatrix/Reconstructor.hpp"

namespace njoy {
namespace resonanceReconstruction {
namespace rmatrix {

#include "resonanceReconstruction/rmatrix/src/makeQuantumNumbers.hpp"
#include "resonanceReconstruction/rmatrix/src/makeParticlePairs.hpp"
#include "resonanceReconstruction/rmatrix/src/makeParticleChannels.hpp"
#include "resonanceReconstruction/rmatrix/src/makeParticleChannelData.hpp"
#include "resonanceReconstruction/rmatrix/src/makeChannelRadiusTable.hpp"
#include "resonanceReconstruction/rmatrix/src/makeChannelRadii.hpp"
#include "resonanceReconstruction/rmatrix/src/makeCompoundSystem.hpp"
#include "resonanceReconstruction/rmatrix/src/makeReichMooreChannelData.hpp"
#include "resonanceReconstruction/rmatrix/src/makeReichMooreCompoundSystem.hpp"
#include "resonanceReconstruction/rmatrix/src/makeLegacyBreitWignerSpinGroups.hpp"
#include "resonanceReconstruction/rmatrix/src/makeLegacyBreitWignerCompoundSystem.hpp"
#include "resonanceReconstruction/rmatrix/src/makeLegacyUnresolvedResonanceTable.hpp"
#include "resonanceReconstruction/rmatrix/src/makeLegacyUnresolvedSpinGroups.hpp"
#include "resonanceReconstruction/rmatrix/src/makeLegacyUnresolvedCompoundSystem.hpp"

Reconstructor
fromENDF( const endf::ResonanceRange& endfResonanceRange,
          const AtomicMass& neutronMass,
          const ElectricalCharge& elementaryCharge,
          const ParticleID& incident,
          const ParticleID& target ) {

  const auto lower = endfResonanceRange.lowerEnergy();
  const auto upper = endfResonanceRange.upperEnergy();
  const auto naps = endfResonanceRange.scatteringRadiusCalculationOption();
  std::optional< ChannelRadiusTable > nro =
    makeChannelRadiusTable( endfResonanceRange.scatteringRadius() );

  switch ( endfResonanceRange.type() ) {

    // resolved resonances
    case 1 : {

      switch ( endfResonanceRange.representation() ) {

        // SLBW
        case 1 : {

          auto endfSLBW =
            std::get< endf::SingleLevelBreitWigner >( endfResonanceRange.parameters() );

          return Reconstructor(
                     lower * electronVolt,
                     upper * electronVolt,
                     makeLegacyBreitWignerCompoundSystem(
                         endfSLBW, neutronMass,
                         incident, target, nro, naps,
                         SingleLevelBreitWigner() ) );
        }
        // MLBW
        case 2 : {

          auto endfMLBW =
            std::get< endf::MultiLevelBreitWigner >( endfResonanceRange.parameters() );

          return Reconstructor(
                     lower * electronVolt,
                     upper * electronVolt,
                     makeLegacyBreitWignerCompoundSystem(
                       endfMLBW, neutronMass,
                       incident, target, nro, naps,
                       MultiLevelBreitWigner() ) );
        }
        // ReichMoore
        case 3 : {

          auto endfReichMoore =
            std::get< endf::ReichMoore >( endfResonanceRange.parameters() );

          return Reconstructor(
                     lower * electronVolt,
                     upper * electronVolt,
                     makeReichMooreCompoundSystem(
                         endfReichMoore, neutronMass,
                         incident, target, nro, naps ) );
        }
        // R-matrix limited
        case 7 : {

          auto endfRMatrix =
            std::get< endf::RMatrixLimited >( endfResonanceRange.parameters() );
          bool shiftFactorBoundary =
              endfRMatrix.particlePairs().shiftFactorFlag().front() == 0 ? true : false;

          switch ( endfRMatrix.formalism() ) {

            // Reich-Moore
            case 3 : {

              if ( shiftFactorBoundary ) {

                return Reconstructor(
                           lower * electronVolt,
                           upper * electronVolt,
                           makeCompoundSystem( endfRMatrix,
                                               neutronMass, elementaryCharge,
                                               incident, target,
                                               ReichMoore(), ShiftFactor() ) );
              }
              else {

                return Reconstructor(
                           lower * electronVolt,
                           upper * electronVolt,
                           makeCompoundSystem( endfRMatrix,
                                               neutronMass, elementaryCharge,
                                               incident, target,
                                               ReichMoore(), Constant() ) );
              }
            }
            // general R-matrix
            case 4 : {

              /*if ( shiftFactorBoundary ) {

                return makeCompoundSystem( endfRMatrix,
                                           neutronMass, elementaryCharge,
                                           GeneralRMatrix(), ShiftFactor() );
              }
              else {

                return makeCompoundSystem( endfRMatrix,
                                           neutronMass, elementaryCharge,
                                           GeneralRMatrix(), Constant() );
              }*/
              throw std::runtime_error( "fromENDF is not implemented for "
                                        "general R-matrix in the R-matrix limited "
                                        "resolved resonances" );
            }
            // SLBW and MLBW
            default : {

              throw std::runtime_error( "fromENDF is not implemented for SLBW "
                                        "and MLBW in the R-matrix limited "
                                        "resolved resonances" );
            }
          }
        }
        // AA
        default : {

          throw std::runtime_error( "fromENDF is not implemented for the "
                                    "AA resolved resonances" );
        }
      }
    }
    // unresolved resonances
    case 2 : {

      switch ( endfResonanceRange.parameters().index() ) {

        case 5: {

          auto endfEnergyIndependent =
            std::get< endf::UnresolvedEnergyIndependent >( endfResonanceRange.parameters() );
          return Reconstructor(
                     lower * electronVolt,
                     upper * electronVolt,
                     makeLegacyUnresolvedCompoundSystem(
                         endfEnergyIndependent,
                         neutronMass, elementaryCharge,
                         incident, target, nro, naps, lower, upper ) );
        }
        case 6: {

          auto endfEnergyDependentFission =
            std::get< endf::UnresolvedEnergyDependentFissionWidths >( endfResonanceRange.parameters() );
          return Reconstructor(
                     lower * electronVolt,
                     upper * electronVolt,
                     makeLegacyUnresolvedCompoundSystem(
                         endfEnergyDependentFission,
                         neutronMass, elementaryCharge,
                         incident, target, nro, naps, lower, upper ) );
        }
        case 7: {

          auto endfEnergyDependent =
            std::get< endf::UnresolvedEnergyDependent >( endfResonanceRange.parameters() );
          return Reconstructor(
                     lower * electronVolt,
                     upper * electronVolt,
                     makeLegacyUnresolvedCompoundSystem(
                         endfEnergyDependent,
                         neutronMass, elementaryCharge,
                         incident, target, nro, naps, lower, upper ) );
        }
        default : {

          throw std::runtime_error( "You somehow reached unreachable code" );
        }
      }
    }
    // special case and unresolved resonances
    default : {

      throw std::runtime_error( "fromENDF is not implemented for the special "
                                "case" );
    }
  }
}

} // rmatrix namespace
} // resonanceReconstruction namespace
} // njoy namespace

#endif
