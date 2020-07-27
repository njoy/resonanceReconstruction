inline Reconstructor
fromENDF( const ENDF::ResonanceRange& endfResonanceRange,
          const AtomicMass& neutronMass,
          const ElectricalCharge& elementaryCharge ) {

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
                           endfResonanceRange.lowerEnergy() * electronVolt,
                           endfResonanceRange.upperEnergy() * electronVolt,
                           makeCompoundSystem( endfRMatrix,
                                               neutronMass, elementaryCharge,
                                               ReichMoore(), ShiftFactor() ) );
              }
              else {

                return Reconstructor(
                           endfResonanceRange.lowerEnergy() * electronVolt,
                           endfResonanceRange.upperEnergy() * electronVolt,
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
    // special case and unresolved resonances
    default : {

      throw std::runtime_error( "fromENDF is not implemented for the special "
                                "case or unresolved resonances" );
    }
  }
}
