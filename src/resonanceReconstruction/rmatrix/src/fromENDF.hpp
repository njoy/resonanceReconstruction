Reconstructor
fromENDF( const ENDF::ResonanceRange& endfResonanceRange,
          const AtomicMass& neutronMass,
          const ElectricalCharge& elementaryCharge,
          const ParticleID& incident,
          const ParticleID& target ) {

  auto lower = endfResonanceRange.lowerEnergy();
  auto upper = endfResonanceRange.upperEnergy();
  auto nro = endfResonanceRange.energyDependentScatteringRadius();
  auto naps = endfResonanceRange.scatteringRadiusCalculationOption();

  switch ( endfResonanceRange.type() ) {

    // resolved resonances
    case 1 : {

      switch ( endfResonanceRange.representation() ) {

        // ReichMoore
        case 3 : {

          throw std::runtime_error( "fromENDF is not implemented for ENDF "
                                    "Reich-Moore resolved resonances (LRF = 3)" );
        }
        // R-matrix limited
        case 7 : {

          auto endfRMatrix =
            std::get< ENDF::resolved::RMatrixLimited >( endfResonanceRange.parameters() );
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
                                               ReichMoore(), ShiftFactor() ) );
              }
              else {

                return Reconstructor(
                           lower * electronVolt,
                           upper * electronVolt,
                           makeCompoundSystem( endfRMatrix,
                                               neutronMass, elementaryCharge,
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
        // SLBW, MLBW, AA
        default : {

          throw std::runtime_error( "fromENDF is not implemented for the SLBW, "
                                    "MLBW or AA resolved resonances" );
        }
      }
    }
    // unresolved resonances
    case 2 : {

      if ( nro ) {

        throw std::runtime_error( "Energy dependent scattering radii have not "
                                  "been implemented" );
      }

      switch ( endfResonanceRange.parameters().index() ) {

        case 5: {

          auto endfEnergyIndependent =
            std::get< ENDF::unresolved::EnergyIndependent >( endfResonanceRange.parameters() );
          return Reconstructor(
                     lower * electronVolt,
                     upper * electronVolt,
                     makeLegacyUnresolvedCompoundSystem(
                         endfEnergyIndependent,
                         neutronMass, elementaryCharge,
                         incident, target, naps, lower, upper ) );
        }
        case 6: {

          auto endfEnergyDependentFission =
            std::get< ENDF::unresolved::EnergyDependentFissionWidths >( endfResonanceRange.parameters() );
          return Reconstructor(
                     lower * electronVolt,
                     upper * electronVolt,
                     makeLegacyUnresolvedCompoundSystem(
                         endfEnergyDependentFission,
                         neutronMass, elementaryCharge,
                         incident, target, naps, lower, upper ) );
        }
        case 7: {

          auto endfEnergyDependent =
            std::get< ENDF::unresolved::EnergyDependent >( endfResonanceRange.parameters() );
          return Reconstructor(
                     lower * electronVolt,
                     upper * electronVolt,
                     makeLegacyUnresolvedCompoundSystem(
                         endfEnergyDependent,
                         neutronMass, elementaryCharge,
                         incident, target, naps, lower, upper ) );
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
