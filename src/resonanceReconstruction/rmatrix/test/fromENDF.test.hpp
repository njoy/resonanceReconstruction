std::string Fe54();
std::string Ca40();
std::string Cl35();

SCENARIO( "fromENDF" ) {

  GIVEN( "valid ENDF data for Fe54" ) {

    std::string string = Fe54();
    auto begin = string.begin();
    auto end = string.end();
    long lineNumber = 1;

    njoy::ENDFtk::HeadRecord head( begin, end, lineNumber );
    njoy::ENDFtk::section::Type< 2, 151 > endf( head, begin, end, lineNumber, 2625 );
    ResonanceRange endfResonanceRange = endf.isotopes().front().resonanceRanges().front();

    auto resonances = fromENDF( endfResonanceRange, neutronMass, elementaryCharge, ParticleID( "n" ), ParticleID( "Fe54" ) );

    THEN( "the appropriate CompoundSystem is returned" ) {

      CHECK( true == resonances.isResolved() );
      CHECK( false == resonances.isUnresolved() );
      CHECK( 1e-5 == Approx( resonances.lowerEnergy().value ) );
      CHECK( 1.036e+6 == Approx( resonances.upperEnergy().value ) );

      auto compoundsystem = std::get< CompoundSystem< ReichMoore, ShiftFactor > >( resonances.compoundSystem() );

      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // content verification
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      // spin groups
      auto spingroups = compoundsystem.spinGroups();
      CHECK( 5 == spingroups.size() );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup0 = spingroups[0];

      // channels
      auto channels0 = spingroup0.channels();

      CHECK( 1 == channels0.size() ); // 1 normal channel + 1 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel00 = std::get< Channel< Neutron > >( channels0[0] );

      // incident particle pair
      const auto incident00 = channel00.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident00.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident00.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident00.particle().spin() ) );
      CHECK( +1 == incident00.particle().parity() );
      CHECK( 53.47624 * 1.008664 == Approx( incident00.residual().mass().value ) );
      CHECK( 26.0 * 1.602e-19 == Approx( incident00.residual().charge().value ) );
      CHECK( 0.0 == Approx( incident00.residual().spin() ) );
      CHECK( +1 == incident00.residual().parity() );
      CHECK( "n,Fe54" == incident00.pairID().symbol() );

      // particle pairs
      const auto pair00 = channel00.particlePair();
      CHECK( 1.008664 == Approx( pair00.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair00.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair00.particle().spin() ) );
      CHECK( +1 == pair00.particle().parity() );
      CHECK( 53.47624 * 1.008664 == Approx( pair00.residual().mass().value ) );
      CHECK( 26.0 * 1.602e-19 == Approx( pair00.residual().charge().value ) );
      CHECK( 0.0 == Approx( pair00.residual().spin() ) );
      CHECK( +1 == pair00.residual().parity() );
      CHECK( "n,Fe54" == pair00.pairID().symbol() );

      // quantum numbers
      const auto numbers00 = channel00.quantumNumbers();
      CHECK( 1 == numbers00.orbitalAngularMomentum() );
      CHECK( 0.5 == numbers00.spin() );
      CHECK( 0.5 == numbers00.totalAngularMomentum() );
      CHECK( -1 == numbers00.parity() );
      CHECK( "{1,1/2,1/2-}" == numbers00.toString() );

      // radii
      const auto radii00 = channel00.radii();
      CHECK( .54373 == Approx( radii00.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .54373 == Approx( radii00.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .54373 == Approx( radii00.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel00.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel00.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup1 = spingroups[1];

      // channels
      auto channels1 = spingroup1.channels();

      CHECK( 1 == channels1.size() ); // 1 normal channel + 1 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel10 = std::get< Channel< Neutron > >( channels1[0] );

      // incident particle pair
      const auto incident10 = channel10.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident10.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident10.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident10.particle().spin() ) );
      CHECK( +1 == incident10.particle().parity() );
      CHECK( 53.47624 * 1.008664 == Approx( incident10.residual().mass().value ) );
      CHECK( 26.0 * 1.602e-19 == Approx( incident10.residual().charge().value ) );
      CHECK( 0.0 == Approx( incident10.residual().spin() ) );
      CHECK( +1 == incident10.residual().parity() );
      CHECK( "n,Fe54" == incident10.pairID().symbol() );

      // particle pair
      const auto pair10 = channel10.particlePair();
      CHECK( 1.008664 == Approx( pair10.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair10.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair10.particle().spin() ) );
      CHECK( +1 == pair10.particle().parity() );
      CHECK( 53.47624 * 1.008664 == Approx( pair10.residual().mass().value ) );
      CHECK( 26.0 * 1.602e-19 == Approx( pair10.residual().charge().value ) );
      CHECK( 0.0 == Approx( pair10.residual().spin() ) );
      CHECK( +1 == pair10.residual().parity() );
      CHECK( "n,Fe54" == pair10.pairID().symbol() );

      // quantum numbers
      const auto numbers10 = channel10.quantumNumbers();
      CHECK( 0 == numbers10.orbitalAngularMomentum() );
      CHECK( 0.5 == numbers10.spin() );
      CHECK( 0.5 == numbers10.totalAngularMomentum() );
      CHECK( +1 == numbers10.parity() );
      CHECK( "{0,1/2,1/2+}" == numbers10.toString() );

      // radii
      const auto radii10 = channel10.radii();
      CHECK( .54373 == Approx( radii10.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .54373 == Approx( radii10.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .54373 == Approx( radii10.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel10.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel10.Q().value ) );

      // resonance table
      auto table1 = spingroup1.resonanceTable();

      CHECK( 1 == table1.numberChannels() ); // 1 normal channel + 1 eliminated
      CHECK( 148 == table1.numberResonances() );

      auto energies1 = table1.energies();
      CHECK( -1.223300e+6 == Approx( energies1.front().value ) );
      CHECK( 1.505510e+6 == Approx( energies1.back().value ) );

      auto resonances1 = table1.resonances();
      CHECK( -1.223300e+6 == Approx( resonances1.front().energy().value ) );
      CHECK( 1.505510e+6 == Approx( resonances1.back().energy().value ) );
      CHECK( 1 == resonances1.front().widths().size() );
      CHECK( 1 == resonances1.back().widths().size() );
      CHECK( std::sqrt( 961108.6 / 2. / channel10.penetrability( -1.223300e+6 * electronVolt ) )
             == Approx( resonances1.front().widths().front().value ) );
      CHECK( std::sqrt( 1601.084 / 2. / channel10.penetrability( 1.505510e+6 * electronVolt ) )
             == Approx( resonances1.back().widths().front().value ) );
      CHECK( std::sqrt( 1. / 2. ) == Approx( resonances1.front().eliminatedWidth().value ) );
      CHECK( std::sqrt( .64906 / 2. ) == Approx( resonances1.back().eliminatedWidth().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 2
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup2 = spingroups[2];

      // channels
      auto channels2 = spingroup2.channels();

      CHECK( 1 == channels2.size() ); // 1 normal channel + 1 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 2, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel20 = std::get< Channel< Neutron > >( channels2[0] );

      // incident particle pair
      const auto incident20 = channel20.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident20.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident20.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident20.particle().spin() ) );
      CHECK( +1 == incident20.particle().parity() );
      CHECK( 53.47624 * 1.008664 == Approx( incident20.residual().mass().value ) );
      CHECK( 26.0 * 1.602e-19 == Approx( incident20.residual().charge().value ) );
      CHECK( 0.0 == Approx( incident20.residual().spin() ) );
      CHECK( +1 == incident20.residual().parity() );
      CHECK( "n,Fe54" == incident20.pairID().symbol() );

      // particle pairs
      const auto pair20 = channel20.particlePair();
      CHECK( 1.008664 == Approx( pair20.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair20.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair20.particle().spin() ) );
      CHECK( +1 == pair20.particle().parity() );
      CHECK( 53.47624 * 1.008664 == Approx( pair20.residual().mass().value ) );
      CHECK( 26.0 * 1.602e-19 == Approx( pair20.residual().charge().value ) );
      CHECK( 0.0 == Approx( pair20.residual().spin() ) );
      CHECK( +1 == pair20.residual().parity() );
      CHECK( "n,Fe54" == pair20.pairID().symbol() );

      // quantum numbers
      const auto numbers20 = channel20.quantumNumbers();
      CHECK( 1 == numbers20.orbitalAngularMomentum() );
      CHECK( 0.5 == numbers20.spin() );
      CHECK( 1.5 == numbers20.totalAngularMomentum() );
      CHECK( -1 == numbers20.parity() );
      CHECK( "{1,1/2,3/2-}" == numbers20.toString() );

      // radii
      const auto radii20 = channel20.radii();
      CHECK( .54373 == Approx( radii20.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .54373 == Approx( radii20.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .54373 == Approx( radii20.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel20.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel20.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 3
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup3 = spingroups[3];

      // channels
      auto channels3 = spingroup3.channels();

      CHECK( 1 == channels3.size() ); // 1 normal channel + 1 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 3, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel30 = std::get< Channel< Neutron > >( channels3[0] );

      // incident particle pair
      const auto incident30 = channel30.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident30.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident30.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident30.particle().spin() ) );
      CHECK( +1 == incident30.particle().parity() );
      CHECK( 53.47624 * 1.008664 == Approx( incident30.residual().mass().value ) );
      CHECK( 26.0 * 1.602e-19 == Approx( incident30.residual().charge().value ) );
      CHECK( 0.0 == Approx( incident30.residual().spin() ) );
      CHECK( +1 == incident30.residual().parity() );
      CHECK( "n,Fe54" == incident30.pairID().symbol() );

      // particle pairs
      const auto pair30 = channel30.particlePair();
      CHECK( 1.008664 == Approx( pair30.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair30.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair30.particle().spin() ) );
      CHECK( +1 == pair30.particle().parity() );
      CHECK( 53.47624 * 1.008664 == Approx( pair30.residual().mass().value ) );
      CHECK( 26.0 * 1.602e-19 == Approx( pair30.residual().charge().value ) );
      CHECK( 0.0 == Approx( pair30.residual().spin() ) );
      CHECK( +1 == pair30.residual().parity() );
      CHECK( "n,Fe54" == pair30.pairID().symbol() );

      // quantum numbers
      const auto numbers30 = channel30.quantumNumbers();
      CHECK( 2 == numbers30.orbitalAngularMomentum() );
      CHECK( 0.5 == numbers30.spin() );
      CHECK( 1.5 == numbers30.totalAngularMomentum() );
      CHECK( +1 == numbers30.parity() );
      CHECK( "{2,1/2,3/2+}" == numbers30.toString() );

      // radii
      const auto radii30 = channel30.radii();
      CHECK( .54373 == Approx( radii30.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .54373 == Approx( radii30.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .54373 == Approx( radii30.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel30.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel30.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 4
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup4 = spingroups[4];

      // channels
      auto channels4 = spingroup4.channels();

      CHECK( 1 == channels4.size() ); // 1 normal channel + 1 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 4, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel40 = std::get< Channel< Neutron > >( channels4[0] );

      // incident particle pair
      const auto incident40 = channel40.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident40.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident40.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident40.particle().spin() ) );
      CHECK( +1 == incident40.particle().parity() );
      CHECK( 53.47624 * 1.008664 == Approx( incident40.residual().mass().value ) );
      CHECK( 26.0 * 1.602e-19 == Approx( incident40.residual().charge().value ) );
      CHECK( 0.0 == Approx( incident40.residual().spin() ) );
      CHECK( +1 == incident40.residual().parity() );
      CHECK( "n,Fe54" == incident40.pairID().symbol() );

      // particle pairs
      const auto pair40 = channel40.particlePair();
      CHECK( 1.008664 == Approx( pair40.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair40.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair40.particle().spin() ) );
      CHECK( +1 == pair40.particle().parity() );
      CHECK( 53.47624 * 1.008664 == Approx( pair40.residual().mass().value ) );
      CHECK( 26.0 * 1.602e-19 == Approx( pair40.residual().charge().value ) );
      CHECK( 0.0 == Approx( pair40.residual().spin() ) );
      CHECK( +1 == pair40.residual().parity() );
      CHECK( "n,Fe54" == pair40.pairID().symbol() );

      // quantum numbers
      const auto numbers40 = channel40.quantumNumbers();
      CHECK( 2 == numbers40.orbitalAngularMomentum() );
      CHECK( 0.5 == numbers40.spin() );
      CHECK( 2.5 == numbers40.totalAngularMomentum() );
      CHECK( +1 == numbers40.parity() );
      CHECK( "{2,1/2,5/2+}" == numbers40.toString() );

      // radii
      const auto radii40 = channel40.radii();
      CHECK( .54373 == Approx( radii40.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .54373 == Approx( radii40.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .54373 == Approx( radii40.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel40.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel40.Q().value ) );
    } // THEN

    THEN( "cross sections can be reconstructed" ) {

      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // resonance reconstruction verification
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      // values taken from NJOY2016 PENDF tape for ENDF/B-VIII.0 Fe54

      ReactionID elas( "n,Fe54->n,Fe54" );
      ReactionID capt( "n,Fe54->capture" );
      std::map< ReactionID, CrossSection > xs;

      xs = resonances( 1e-5 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 2.162315 == Approx( xs[ elas ].value ) );
      CHECK( 1.133306e+2 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 2.162315 == Approx( xs[ elas ].value ) );
      CHECK( 3.583829e+1 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 2.162314 == Approx( xs[ elas ].value ) );
      CHECK( 1.133307e+1 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-2 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 2.162307 == Approx( xs[ elas ].value ) );
      CHECK( 3.583834 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-1 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 2.162233 == Approx( xs[ elas ].value ) );
      CHECK( 1.133321 == Approx( xs[ capt ].value ) );

      xs = resonances( 1. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 2.161498 == Approx( xs[ elas ].value ) );
      CHECK( 3.584290e-1 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e+1 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 2.154145 == Approx( xs[ elas ].value ) );
      CHECK( 1.134767e-1 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e+2 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 2.080799 == Approx( xs[ elas ].value ) );
      CHECK( 3.631075e-2 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 1.366975 == Approx( xs[ elas ].value ) );
      CHECK( 1.321792e-2 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 46.15276 == Approx( xs[ elas ].value ) );
      CHECK( 2.555440e-2 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e+5 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 3.750047 == Approx( xs[ elas ].value ) );
      CHECK( 5.879488e-3 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e+6 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 4.481752 == Approx( xs[ elas ].value ) );
      CHECK( 1.369045e-3 == Approx( xs[ capt ].value ) );
    } // THEN
  } // GIVEN

  GIVEN( "valid ENDF data for Ca40" ) {

    std::string string = Ca40();
    auto begin = string.begin();
    auto end = string.end();
    long lineNumber = 1;

    njoy::ENDFtk::HeadRecord head( begin, end, lineNumber );
    njoy::ENDFtk::section::Type< 2, 151 > endf( head, begin, end, lineNumber, 2025 );
    ResonanceRange endfResonanceRange = endf.isotopes().front().resonanceRanges().front();

    auto resonances = fromENDF( endfResonanceRange, neutronMass, elementaryCharge, ParticleID( "n" ), ParticleID( "Ca40" ) );

    THEN( "the appropriate CompoundSystem is returned" ) {

      CHECK( true == resonances.isResolved() );
      CHECK( false == resonances.isUnresolved() );
      CHECK( 1e-5 == Approx( resonances.lowerEnergy().value ) );
      CHECK( 1.5e+6 == Approx( resonances.upperEnergy().value ) );

      auto compoundsystem = std::get< CompoundSystem< ReichMoore, ShiftFactor > >( resonances.compoundSystem() );

      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // content verification
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      // spin groups
      auto spingroups = compoundsystem.spinGroups();
      CHECK( 5 == spingroups.size() );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup0 = spingroups[0];

      // channels
      auto channels0 = spingroup0.channels();

      CHECK( 2 == channels0.size() ); // 2 normal channel + 1 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel00 = std::get< Channel< Neutron > >( channels0[0] );

      // incident particle pair
      const auto incident00 = channel00.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident00.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident00.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident00.particle().spin() ) );
      CHECK( +1 == incident00.particle().parity() );
      CHECK( 39.61929 * 1.008664 == Approx( incident00.residual().mass().value ) );
      CHECK( 20.0 * 1.602e-19 == Approx( incident00.residual().charge().value ) );
      CHECK( 0.0 == Approx( incident00.residual().spin() ) );
      CHECK( +1 == incident00.residual().parity() );
      CHECK( "n,Ca40" == incident00.pairID().symbol() );

      // particle pair
      const auto pair00 = channel00.particlePair();
      CHECK( 1.008664 == Approx( pair00.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair00.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair00.particle().spin() ) );
      CHECK( +1 == pair00.particle().parity() );
      CHECK( 39.61929 * 1.008664 == Approx( pair00.residual().mass().value ) );
      CHECK( 20.0 * 1.602e-19 == Approx( pair00.residual().charge().value ) );
      CHECK( 0.0 == Approx( pair00.residual().spin() ) );
      CHECK( +1 == pair00.residual().parity() );
      CHECK( "n,Ca40" == pair00.pairID().symbol() );

      // quantum numbers
      const auto numbers00 = channel00.quantumNumbers();
      CHECK( 1 == numbers00.orbitalAngularMomentum() );
      CHECK( 0.5 == numbers00.spin() );
      CHECK( 0.5 == numbers00.totalAngularMomentum() );
      CHECK( -1 == numbers00.parity() );
      CHECK( "{1,1/2,1/2-}" == numbers00.toString() );

      // radii
      const auto radii00 = channel00.radii();
      CHECK( 0.4993153 == Approx( radii00.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4993153 == Approx( radii00.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.5574746 == Approx( radii00.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel00.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel00.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0, channel 1: n,a
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel01 = std::get< Channel< ChargedParticle > >( channels0[1] );

      // incident particle pair
      const auto incident01 = channel01.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident01.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident01.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident01.particle().spin() ) );
      CHECK( +1 == incident01.particle().parity() );
      CHECK( 39.61929 * 1.008664 == Approx( incident01.residual().mass().value ) );
      CHECK( 20.0 * 1.602e-19 == Approx( incident01.residual().charge().value ) );
      CHECK( 0.0 == Approx( incident01.residual().spin() ) );
      CHECK( +1 == incident01.residual().parity() );
      CHECK( "n,Ca40" == incident01.pairID().symbol() );

      // particle pair
      const auto pair01 = channel01.particlePair();
      CHECK( 3.967131 * 1.008664 == Approx( pair01.particle().mass().value ) );
      CHECK( 2.0 * 1.602e-19 == Approx( pair01.particle().charge().value ) );
      CHECK( 0.0 == Approx( pair01.particle().spin() ) );
      CHECK( +1 == pair01.particle().parity() );
      CHECK( 36.64921 * 1.008664 == Approx( pair01.residual().mass().value ) );
      CHECK( 18.0 * 1.602e-19 == Approx( pair01.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair01.residual().spin() ) );
      CHECK( +1 == pair01.residual().parity() );
      CHECK( "he4,Ar37" == pair01.pairID().symbol() );

      // quantum numbers
      const auto numbers01 = channel01.quantumNumbers();
      CHECK( 1 == numbers01.orbitalAngularMomentum() );
      CHECK( 1.5 == numbers01.spin() );
      CHECK( 0.5 == numbers01.totalAngularMomentum() );
      CHECK( -1 == numbers01.parity() );
      CHECK( "{1,3/2,1/2-}" == numbers01.toString() );

      // radii
      const auto radii01 = channel01.radii();
      CHECK( 0.4885639 == Approx( radii01.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4885639 == Approx( radii01.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4885639 == Approx( radii01.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel01.boundaryCondition() );

      // Q value
      CHECK( 1.747660e+6 == Approx( channel01.Q().value ) );

      // resonance table
      auto table0 = spingroup0.resonanceTable();

      CHECK( 2 == table0.numberChannels() ); // 2 normal channel + 1 eliminated
      CHECK( 51 == table0.numberResonances() );

      auto energies0 = table0.energies();
      CHECK( 2.880782e+3 == Approx( energies0.front().value ) );
      CHECK( 1.496383e+6 == Approx( energies0.back().value ) );

      auto resonances0 = table0.resonances();
      CHECK( 2.880782e+3 == Approx( resonances0.front().energy().value ) );
      CHECK( 1.496383e+6 == Approx( resonances0.back().energy().value ) );
      CHECK( 2 == resonances0.front().widths().size() );
      CHECK( 2 == resonances0.back().widths().size() );
      CHECK( std::sqrt( 1.810972e-3 / 2. / channel00.penetrability( 2.880782e+3 * electronVolt ) )
             == Approx( resonances0.front().widths()[0].value ) );
      CHECK( std::sqrt( 1.000029e-3 / 2. / channel01.penetrability( 2.880782e+3 * electronVolt ) )
             == Approx( resonances0.front().widths()[1].value ) );
      CHECK( std::sqrt( 1.448484e+2 / 2. / channel00.penetrability( 1.496383e+6 * electronVolt ) )
             == Approx( resonances0.back().widths()[0].value ) );
      CHECK( std::sqrt( 1.012630e-3 / 2. / channel01.penetrability( 1.496383e+6 * electronVolt ) )
             == Approx( resonances0.back().widths()[1].value ) );
      CHECK( std::sqrt( 3.499989e-1 / 2. ) == Approx( resonances0.front().eliminatedWidth().value ) );
      CHECK( std::sqrt( 1.050878 / 2. ) == Approx( resonances0.back().eliminatedWidth().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup1 = spingroups[1];

      // channels
      auto channels1 = spingroup1.channels();

      CHECK( 2 == channels1.size() ); // 2 normal channel + 1 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel10 = std::get< Channel< Neutron > >( channels1[0] );

      // incident particle pair
      const auto incident10 = channel10.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident10.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident10.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident10.particle().spin() ) );
      CHECK( +1 == incident10.particle().parity() );
      CHECK( 39.61929 * 1.008664 == Approx( incident10.residual().mass().value ) );
      CHECK( 20.0 * 1.602e-19 == Approx( incident10.residual().charge().value ) );
      CHECK( 0.0 == Approx( incident10.residual().spin() ) );
      CHECK( +1 == incident10.residual().parity() );
      CHECK( "n,Ca40" == incident10.pairID().symbol() );

      // particle pair
      const auto pair10 = channel10.particlePair();
      CHECK( 1.008664 == Approx( pair10.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair10.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair10.particle().spin() ) );
      CHECK( +1 == pair10.particle().parity() );
      CHECK( 39.61929 * 1.008664 == Approx( pair10.residual().mass().value ) );
      CHECK( 20.0 * 1.602e-19 == Approx( pair10.residual().charge().value ) );
      CHECK( 0.0 == Approx( pair10.residual().spin() ) );
      CHECK( +1 == pair10.residual().parity() );
      CHECK( "n,Ca40" == pair10.pairID().symbol() );

      // quantum numbers
      const auto numbers10 = channel10.quantumNumbers();
      CHECK( 0 == numbers10.orbitalAngularMomentum() );
      CHECK( 0.5 == numbers10.spin() );
      CHECK( 0.5 == numbers10.totalAngularMomentum() );
      CHECK( +1 == numbers10.parity() );
      CHECK( "{0,1/2,1/2+}" == numbers10.toString() );

      // radii
      const auto radii10 = channel10.radii();
      CHECK( 0.4993153 == Approx( radii10.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4993153 == Approx( radii10.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.3828792 == Approx( radii10.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel10.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel10.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1, channel 1: n,a
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel11 = std::get< Channel< ChargedParticle > >( channels1[1] );

      // incident particle pair
      const auto incident11 = channel11.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident11.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident11.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident11.particle().spin() ) );
      CHECK( +1 == incident11.particle().parity() );
      CHECK( 39.61929 * 1.008664 == Approx( incident11.residual().mass().value ) );
      CHECK( 20.0 * 1.602e-19 == Approx( incident11.residual().charge().value ) );
      CHECK( 0.0 == Approx( incident11.residual().spin() ) );
      CHECK( +1 == incident11.residual().parity() );
      CHECK( "n,Ca40" == incident11.pairID().symbol() );

      // particle pair
      const auto pair11 = channel11.particlePair();
      CHECK( 3.967131 * 1.008664 == Approx( pair11.particle().mass().value ) );
      CHECK( 2.0 * 1.602e-19 == Approx( pair11.particle().charge().value ) );
      CHECK( 0.0 == Approx( pair11.particle().spin() ) );
      CHECK( +1 == pair11.particle().parity() );
      CHECK( 36.64921 * 1.008664 == Approx( pair11.residual().mass().value ) );
      CHECK( 18.0 * 1.602e-19 == Approx( pair11.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair11.residual().spin() ) );
      CHECK( +1 == pair11.residual().parity() );
      CHECK( "he4,Ar37" == pair11.pairID().symbol() );

      // quantum numbers
      const auto numbers11 = channel11.quantumNumbers();
      CHECK( 2 == numbers11.orbitalAngularMomentum() );
      CHECK( 1.5 == numbers11.spin() );
      CHECK( 0.5 == numbers11.totalAngularMomentum() );
      CHECK( +1 == numbers11.parity() );
      CHECK( "{2,3/2,1/2+}" == numbers11.toString() );

      // radii
      const auto radii11 = channel11.radii();
      CHECK( 0.4885639 == Approx( radii11.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4885639 == Approx( radii11.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4885639 == Approx( radii11.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel11.boundaryCondition() );

      // Q value
      CHECK( 1.747660e+6 == Approx( channel11.Q().value ) );

      // resonance table
      auto table1 = spingroup1.resonanceTable();

      CHECK( 2 == table1.numberChannels() ); // 2 normal channel + 1 eliminated
      CHECK( 39 == table1.numberResonances() );

      auto energies1 = table1.energies();
      CHECK( -4.586687e+5 == Approx( energies1.front().value ) );
      CHECK( 1.913996e+6 == Approx( energies1.back().value ) );

      auto resonances1 = table1.resonances();
      CHECK( -4.586687e+5 == Approx( resonances1.front().energy().value ) );
      CHECK( 1.913996e+6 == Approx( resonances1.back().energy().value ) );
      CHECK( 2 == resonances1.front().widths().size() );
      CHECK( 2 == resonances1.back().widths().size() );
      CHECK( std::sqrt( 9.809761e+2 / 2. / channel10.penetrability( -4.586687e+5 * electronVolt ) )
             == Approx( resonances1.front().widths()[0].value ) );
      CHECK( std::sqrt( 2.199782e-3 / 2. / channel11.penetrability( -4.586687e+5 * electronVolt ) )
             == Approx( resonances1.front().widths()[1].value ) );
      CHECK( std::sqrt( 1.878975e+5 / 2. / channel10.penetrability( 1.913996e+6 * electronVolt ) )
             == Approx( resonances1.back().widths()[0].value ) );
      CHECK( std::sqrt( 2.784131e-4 / 2. / channel11.penetrability( 1.913996e+6 * electronVolt ) )
             == Approx( resonances1.back().widths()[1].value ) );
      CHECK( std::sqrt( 1.000091 / 2. ) == Approx( resonances1.front().eliminatedWidth().value ) );
      CHECK( std::sqrt( 1. / 2. ) == Approx( resonances1.back().eliminatedWidth().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 2
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup2 = spingroups[2];

      // channels
      auto channels2 = spingroup2.channels();

      CHECK( 3 == channels2.size() ); // 3 normal channel + 1 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 2, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel20 = std::get< Channel< Neutron > >( channels2[0] );

      // incident particle pair
      const auto incident20 = channel20.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident20.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident20.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident20.particle().spin() ) );
      CHECK( +1 == incident20.particle().parity() );
      CHECK( 39.61929 * 1.008664 == Approx( incident20.residual().mass().value ) );
      CHECK( 20.0 * 1.602e-19 == Approx( incident20.residual().charge().value ) );
      CHECK( 0.0 == Approx( incident20.residual().spin() ) );
      CHECK( +1 == incident20.residual().parity() );
      CHECK( "n,Ca40" == incident20.pairID().symbol() );

      // particle pair
      const auto pair20 = channel20.particlePair();
      CHECK( 1.008664 == Approx( pair20.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair20.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair20.particle().spin() ) );
      CHECK( +1 == pair20.particle().parity() );
      CHECK( 39.61929 * 1.008664 == Approx( pair20.residual().mass().value ) );
      CHECK( 20.0 * 1.602e-19 == Approx( pair20.residual().charge().value ) );
      CHECK( 0.0 == Approx( pair20.residual().spin() ) );
      CHECK( +1 == pair20.residual().parity() );
      CHECK( "n,Ca40" == pair20.pairID().symbol() );

      // quantum numbers
      const auto numbers20 = channel20.quantumNumbers();
      CHECK( 1 == numbers20.orbitalAngularMomentum() );
      CHECK( 0.5 == numbers20.spin() );
      CHECK( 1.5 == numbers20.totalAngularMomentum() );
      CHECK( -1 == numbers20.parity() );
      CHECK( "{1,1/2,3/2-}" == numbers20.toString() );

      // radii
      const auto radii20 = channel20.radii();
      CHECK( 0.4993153 == Approx( radii20.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4993153 == Approx( radii20.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.5574746 == Approx( radii20.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel20.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel20.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 2, channel 1: n,p
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel21 = std::get< Channel< ChargedParticle > >( channels2[1] );

      // incident particle pair
      const auto incident21 = channel21.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident21.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident21.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident21.particle().spin() ) );
      CHECK( +1 == incident21.particle().parity() );
      CHECK( 39.61929 * 1.008664 == Approx( incident21.residual().mass().value ) );
      CHECK( 20.0 * 1.602e-19 == Approx( incident21.residual().charge().value ) );
      CHECK( 0.0 == Approx( incident21.residual().spin() ) );
      CHECK( +1 == incident21.residual().parity() );
      CHECK( "n,Ca40" == incident21.pairID().symbol() );

      // particle pair
      const auto pair21 = channel21.particlePair();
      CHECK( 9.986235e-1 * 1.008664 == Approx( pair21.particle().mass().value ) );
      CHECK( 1.0 * 1.602e-19 == Approx( pair21.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair21.particle().spin() ) );
      CHECK( +1 == pair21.particle().parity() );
      CHECK( 39.62069 * 1.008664 == Approx( pair21.residual().mass().value ) );
      CHECK( 19.0 * 1.602e-19 == Approx( pair21.residual().charge().value ) );
      CHECK( 4.0 == Approx( pair21.residual().spin() ) );
      CHECK( -1 == pair21.residual().parity() );
      CHECK( "p,K40" == pair21.pairID().symbol() );

      // quantum numbers
      const auto numbers21 = channel21.quantumNumbers();
      CHECK( 2 == numbers21.orbitalAngularMomentum() );
      CHECK( 3.5 == numbers21.spin() );
      CHECK( 1.5 == numbers21.totalAngularMomentum() );
      CHECK( -1 == numbers21.parity() );
      CHECK( "{2,7/2,3/2-}" == numbers21.toString() );

      // radii
      const auto radii21 = channel21.radii();
      CHECK( 0.4993202 == Approx( radii21.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4993202 == Approx( radii21.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4993202 == Approx( radii21.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel21.boundaryCondition() );

      // Q value
      CHECK( -5.285469e+5 == Approx( channel21.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 2, channel 2: n,a
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel22 = std::get< Channel< ChargedParticle > >( channels2[2] );

      // incident particle pair
      const auto incident22 = channel22.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident22.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident22.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident22.particle().spin() ) );
      CHECK( +1 == incident22.particle().parity() );
      CHECK( 39.61929 * 1.008664 == Approx( incident22.residual().mass().value ) );
      CHECK( 20.0 * 1.602e-19 == Approx( incident22.residual().charge().value ) );
      CHECK( 0.0 == Approx( incident22.residual().spin() ) );
      CHECK( +1 == incident22.residual().parity() );
      CHECK( "n,Ca40" == incident22.pairID().symbol() );

      // particle pair
      const auto pair22 = channel22.particlePair();
      CHECK( 3.967131 * 1.008664 == Approx( pair22.particle().mass().value ) );
      CHECK( 2.0 * 1.602e-19 == Approx( pair22.particle().charge().value ) );
      CHECK( 0.0 == Approx( pair22.particle().spin() ) );
      CHECK( +1 == pair22.particle().parity() );
      CHECK( 36.64921 * 1.008664 == Approx( pair22.residual().mass().value ) );
      CHECK( 18.0 * 1.602e-19 == Approx( pair22.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair22.residual().spin() ) );
      CHECK( +1 == pair22.residual().parity() );
      CHECK( "he4,Ar37" == pair22.pairID().symbol() );

      // quantum numbers
      const auto numbers22 = channel22.quantumNumbers();
      CHECK( 1 == numbers22.orbitalAngularMomentum() );
      CHECK( 1.5 == numbers22.spin() );
      CHECK( 1.5 == numbers22.totalAngularMomentum() );
      CHECK( -1 == numbers22.parity() );
      CHECK( "{1,3/2,3/2-}" == numbers22.toString() );

      // radii
      const auto radii22 = channel22.radii();
      CHECK( 0.4885639 == Approx( radii22.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4885639 == Approx( radii22.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4885639 == Approx( radii22.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel22.boundaryCondition() );

      // Q value
      CHECK( 1.747660e+6 == Approx( channel22.Q().value ) );

      // resonance table
      auto table2 = spingroup2.resonanceTable();

      CHECK( 3 == table2.numberChannels() ); // 3 normal channel + 1 eliminated
      CHECK( 52 == table2.numberResonances() );

      auto energies2 = table2.energies();
      CHECK( 1.084167e+4 == Approx( energies2.front().value ) );
      CHECK( 1.477528e+6 == Approx( energies2.back().value ) );

      auto resonances2 = table2.resonances();
      CHECK( 1.084167e+4 == Approx( resonances2.front().energy().value ) );
      CHECK( 1.477528e+6 == Approx( resonances2.back().energy().value ) );
      CHECK( 3 == resonances2.front().widths().size() );
      CHECK( 3 == resonances2.back().widths().size() );
      CHECK( std::sqrt( 5.834810e-1 / 2. / channel20.penetrability( 1.084167e+4 * electronVolt ) )
             == Approx( resonances2.front().widths()[0].value ) );
      CHECK( std::sqrt( 0.0 / 2. / channel21.penetrability( 1.084167e+4 * electronVolt ) )
             == Approx( resonances2.front().widths()[1].value ) );
      CHECK( std::sqrt( 9.999433e-4 / 2. / channel22.penetrability( 1.084167e+4 * electronVolt ) )
             == Approx( resonances2.front().widths()[2].value ) );
      CHECK( std::sqrt( 1.603824e+2 / 2. / channel20.penetrability( 1.477528e+6 * electronVolt ) )
             == Approx( resonances2.back().widths()[0].value ) );
      CHECK( std::sqrt( 1.500000e-4 / 2. / channel21.penetrability( 1.477528e+6 * electronVolt ) )
             == Approx( resonances2.back().widths()[1].value ) );
      CHECK( std::sqrt( 1.008688e-3 / 2. / channel22.penetrability( 1.477528e+6 * electronVolt ) )
             == Approx( resonances2.back().widths()[2].value ) );
      CHECK( std::sqrt( 0.6385209 / 2. ) == Approx( resonances2.front().eliminatedWidth().value ) );
      CHECK( std::sqrt( 0.9389508 / 2. ) == Approx( resonances2.back().eliminatedWidth().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 3
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup3 = spingroups[3];

      // channels
      auto channels3 = spingroup3.channels();

      CHECK( 3 == channels3.size() ); // 2 normal channel + 1 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 3, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel30 = std::get< Channel< Neutron > >( channels3[0] );

      // incident particle pair
      const auto incident30 = channel30.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident30.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident30.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident30.particle().spin() ) );
      CHECK( +1 == incident30.particle().parity() );
      CHECK( 39.61929 * 1.008664 == Approx( incident30.residual().mass().value ) );
      CHECK( 20.0 * 1.602e-19 == Approx( incident30.residual().charge().value ) );
      CHECK( 0.0 == Approx( incident30.residual().spin() ) );
      CHECK( +1 == incident30.residual().parity() );
      CHECK( "n,Ca40" == incident30.pairID().symbol() );

      // particle pair
      const auto pair30 = channel30.particlePair();
      CHECK( 1.008664 == Approx( pair30.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair30.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair30.particle().spin() ) );
      CHECK( +1 == pair30.particle().parity() );
      CHECK( 39.61929 * 1.008664 == Approx( pair30.residual().mass().value ) );
      CHECK( 20.0 * 1.602e-19 == Approx( pair30.residual().charge().value ) );
      CHECK( 0.0 == Approx( pair30.residual().spin() ) );
      CHECK( +1 == pair30.residual().parity() );
      CHECK( "n,Ca40" == pair30.pairID().symbol() );

      // quantum numbers
      const auto numbers30 = channel30.quantumNumbers();
      CHECK( 2 == numbers30.orbitalAngularMomentum() );
      CHECK( 0.5 == numbers30.spin() );
      CHECK( 1.5 == numbers30.totalAngularMomentum() );
      CHECK( +1 == numbers30.parity() );
      CHECK( "{2,1/2,3/2+}" == numbers30.toString() );

      // radii
      const auto radii30 = channel30.radii();
      CHECK( 0.4993153 == Approx( radii30.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4993153 == Approx( radii30.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.3833842 == Approx( radii30.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel30.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel30.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 3, channel 1: n,a
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel31 = std::get< Channel< ChargedParticle > >( channels3[1] );

      // incident particle pair
      const auto incident31 = channel31.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident31.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident31.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident31.particle().spin() ) );
      CHECK( +1 == incident31.particle().parity() );
      CHECK( 39.61929 * 1.008664 == Approx( incident31.residual().mass().value ) );
      CHECK( 20.0 * 1.602e-19 == Approx( incident31.residual().charge().value ) );
      CHECK( 0.0 == Approx( incident31.residual().spin() ) );
      CHECK( +1 == incident31.residual().parity() );
      CHECK( "n,Ca40" == incident31.pairID().symbol() );

      // particle pair
      const auto pair31 = channel31.particlePair();
      CHECK( 3.967131 * 1.008664 == Approx( pair31.particle().mass().value ) );
      CHECK( 2.0 * 1.602e-19 == Approx( pair31.particle().charge().value ) );
      CHECK( 0.0 == Approx( pair31.particle().spin() ) );
      CHECK( +1 == pair31.particle().parity() );
      CHECK( 36.64921 * 1.008664 == Approx( pair31.residual().mass().value ) );
      CHECK( 18.0 * 1.602e-19 == Approx( pair31.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair31.residual().spin() ) );
      CHECK( +1 == pair31.residual().parity() );
      CHECK( "he4,Ar37" == pair31.pairID().symbol() );

      // quantum numbers
      const auto numbers31 = channel31.quantumNumbers();
      CHECK( 0 == numbers31.orbitalAngularMomentum() );
      CHECK( 1.5 == numbers31.spin() );
      CHECK( 1.5 == numbers31.totalAngularMomentum() );
      CHECK( +1 == numbers31.parity() );
      CHECK( "{0,3/2,3/2+}" == numbers31.toString() );

      // radii
      const auto radii31 = channel31.radii();
      CHECK( 0.4885639 == Approx( radii31.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4885639 == Approx( radii31.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4885639 == Approx( radii31.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel31.boundaryCondition() );

      // Q value
      CHECK( 1.747660e+6 == Approx( channel31.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 3, channel 2: n,a
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel32 = std::get< Channel< ChargedParticle > >( channels3[2] );

      // incident particle pair
      const auto incident32 = channel32.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident32.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident32.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident32.particle().spin() ) );
      CHECK( +1 == incident32.particle().parity() );
      CHECK( 39.61929 * 1.008664 == Approx( incident32.residual().mass().value ) );
      CHECK( 20.0 * 1.602e-19 == Approx( incident32.residual().charge().value ) );
      CHECK( 0.0 == Approx( incident32.residual().spin() ) );
      CHECK( +1 == incident32.residual().parity() );
      CHECK( "n,Ca40" == incident32.pairID().symbol() );

      // particle pair
      const auto pair32 = channel32.particlePair();
      CHECK( 3.967131 * 1.008664 == Approx( pair32.particle().mass().value ) );
      CHECK( 2.0 * 1.602e-19 == Approx( pair32.particle().charge().value ) );
      CHECK( 0.0 == Approx( pair32.particle().spin() ) );
      CHECK( +1 == pair32.particle().parity() );
      CHECK( 36.64921 * 1.008664 == Approx( pair32.residual().mass().value ) );
      CHECK( 18.0 * 1.602e-19 == Approx( pair32.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair32.residual().spin() ) );
      CHECK( +1 == pair32.residual().parity() );
      CHECK( "he4,Ar37" == pair32.pairID().symbol() );

      // quantum numbers
      const auto numbers32 = channel32.quantumNumbers();
      CHECK( 2 == numbers32.orbitalAngularMomentum() );
      CHECK( 1.5 == numbers32.spin() );
      CHECK( 1.5 == numbers32.totalAngularMomentum() );
      CHECK( +1 == numbers32.parity() );
      CHECK( "{2,3/2,3/2+}" == numbers32.toString() );

      // radii
      const auto radii32 = channel32.radii();
      CHECK( 0.4885639 == Approx( radii32.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4885639 == Approx( radii32.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4885639 == Approx( radii32.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel32.boundaryCondition() );

      // Q value
      CHECK( 1.747660e+6 == Approx( channel32.Q().value ) );

      // resonance table
      auto table3 = spingroup3.resonanceTable();

      CHECK( 3 == table3.numberChannels() ); // 3 normal channel + 1 eliminated
      CHECK( 54 == table3.numberResonances() );

      auto energies3 = table3.energies();
      CHECK( 3.026620e+3 == Approx( energies3.front().value ) );
      CHECK( 1.515000e+6 == Approx( energies3.back().value ) );

      auto resonances3 = table3.resonances();
      CHECK( 3.026620e+3 == Approx( resonances3.front().energy().value ) );
      CHECK( 1.515000e+6 == Approx( resonances3.back().energy().value ) );
      CHECK( 3 == resonances3.front().widths().size() );
      CHECK( 3 == resonances3.back().widths().size() );
      CHECK( std::sqrt( 1.501024e-4 / 2. / channel30.penetrability( 3.026620e+3 * electronVolt ) )
             == Approx( resonances3.front().widths()[0].value ) );
      CHECK( std::sqrt( 1.000013e-3 / 2. / channel31.penetrability( 3.026620e+3 * electronVolt ) )
             == Approx( resonances3.front().widths()[1].value ) );
      CHECK( std::sqrt( 1.000013e-3 / 2. / channel32.penetrability( 3.026620e+3 * electronVolt ) )
             == Approx( resonances3.front().widths()[2].value ) );
      CHECK( std::sqrt( 25. / 2. / channel30.penetrability( 1.515000e+6 * electronVolt ) )
             == Approx( resonances3.back().widths()[0].value ) );
      CHECK( std::sqrt( 1.e-3 / 2. / channel31.penetrability( 1.515000e+6 * electronVolt ) )
             == Approx( resonances3.back().widths()[1].value ) );
      CHECK( std::sqrt( 1.e-3 / 2. / channel32.penetrability( 1.515000e+6 * electronVolt ) )
             == Approx( resonances3.back().widths()[2].value ) );
      CHECK( std::sqrt( 5.099998e-1 / 2. ) == Approx( resonances3.front().eliminatedWidth().value ) );
      CHECK( std::sqrt( 1. / 2. ) == Approx( resonances3.back().eliminatedWidth().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 4
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup4 = spingroups[4];

      // channels
      auto channels4 = spingroup4.channels();

      CHECK( 3 == channels4.size() ); // 3 normal channel + 1 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 4, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel40 = std::get< Channel< Neutron > >( channels4[0] );

      // incident particle pair
      const auto incident40 = channel40.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident40.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident40.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident40.particle().spin() ) );
      CHECK( +1 == incident40.particle().parity() );
      CHECK( 39.61929 * 1.008664 == Approx( incident40.residual().mass().value ) );
      CHECK( 20.0 * 1.602e-19 == Approx( incident40.residual().charge().value ) );
      CHECK( 0.0 == Approx( incident40.residual().spin() ) );
      CHECK( +1 == incident40.residual().parity() );
      CHECK( "n,Ca40" == incident40.pairID().symbol() );

      // particle pair
      const auto pair40 = channel40.particlePair();
      CHECK( 1.008664 == Approx( pair40.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair40.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair40.particle().spin() ) );
      CHECK( +1 == pair40.particle().parity() );
      CHECK( 39.61929 * 1.008664 == Approx( pair40.residual().mass().value ) );
      CHECK( 20.0 * 1.602e-19 == Approx( pair40.residual().charge().value ) );
      CHECK( 0.0 == Approx( pair40.residual().spin() ) );
      CHECK( +1 == pair40.residual().parity() );
      CHECK( "n,Ca40" == pair40.pairID().symbol() );

      // quantum numbers
      const auto numbers40 = channel40.quantumNumbers();
      CHECK( 2 == numbers40.orbitalAngularMomentum() );
      CHECK( 0.5 == numbers40.spin() );
      CHECK( 2.5 == numbers40.totalAngularMomentum() );
      CHECK( +1 == numbers40.parity() );
      CHECK( "{2,1/2,5/2+}" == numbers40.toString() );

      // radii
      const auto radii40 = channel40.radii();
      CHECK( 0.4993153 == Approx( radii40.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4993153 == Approx( radii40.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.3833842 == Approx( radii40.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel40.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel40.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 4, channel 1: n,p
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel41 = std::get< Channel< ChargedParticle > >( channels4[1] );

      // incident particle pair
      const auto incident41 = channel41.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident41.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident41.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident41.particle().spin() ) );
      CHECK( +1 == incident41.particle().parity() );
      CHECK( 39.61929 * 1.008664 == Approx( incident41.residual().mass().value ) );
      CHECK( 20.0 * 1.602e-19 == Approx( incident41.residual().charge().value ) );
      CHECK( 0.0 == Approx( incident41.residual().spin() ) );
      CHECK( +1 == incident41.residual().parity() );
      CHECK( "n,Ca40" == incident41.pairID().symbol() );

      // particle pair
      const auto pair41 = channel41.particlePair();
      CHECK( 9.986235e-1 * 1.008664 == Approx( pair41.particle().mass().value ) );
      CHECK( 1.0 * 1.602e-19 == Approx( pair41.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair41.particle().spin() ) );
      CHECK( +1 == pair41.particle().parity() );
      CHECK( 39.62069 * 1.008664 == Approx( pair41.residual().mass().value ) );
      CHECK( 19.0 * 1.602e-19 == Approx( pair41.residual().charge().value ) );
      CHECK( 4.0 == Approx( pair41.residual().spin() ) );
      CHECK( -1 == pair41.residual().parity() );
      CHECK( "p,K40" == pair41.pairID().symbol() );

      // quantum numbers
      const auto numbers41 = channel41.quantumNumbers();
      CHECK( 1 == numbers41.orbitalAngularMomentum() );
      CHECK( 3.5 == numbers41.spin() );
      CHECK( 2.5 == numbers41.totalAngularMomentum() );
      CHECK( +1 == numbers41.parity() );
      CHECK( "{1,7/2,5/2+}" == numbers41.toString() );

      // radii
      const auto radii41 = channel41.radii();
      CHECK( 0.4993202 == Approx( radii41.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4993202 == Approx( radii41.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4993202 == Approx( radii41.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel41.boundaryCondition() );

      // Q value
      CHECK( -5.285469e+5 == Approx( channel41.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 4, channel 2: n,a
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel42 = std::get< Channel< ChargedParticle > >( channels4[2] );

      // incident particle pair
      const auto incident42 = channel42.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident42.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident42.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident42.particle().spin() ) );
      CHECK( +1 == incident42.particle().parity() );
      CHECK( 39.61929 * 1.008664 == Approx( incident42.residual().mass().value ) );
      CHECK( 20.0 * 1.602e-19 == Approx( incident42.residual().charge().value ) );
      CHECK( 0.0 == Approx( incident42.residual().spin() ) );
      CHECK( +1 == incident42.residual().parity() );
      CHECK( "n,Ca40" == incident42.pairID().symbol() );

      // particle pair
      const auto pair42 = channel42.particlePair();
      CHECK( 3.967131 * 1.008664 == Approx( pair42.particle().mass().value ) );
      CHECK( 2.0 * 1.602e-19 == Approx( pair42.particle().charge().value ) );
      CHECK( 0.0 == Approx( pair42.particle().spin() ) );
      CHECK( +1 == pair42.particle().parity() );
      CHECK( 36.64921 * 1.008664 == Approx( pair42.residual().mass().value ) );
      CHECK( 18.0 * 1.602e-19 == Approx( pair42.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair42.residual().spin() ) );
      CHECK( +1 == pair42.residual().parity() );
      CHECK( "he4,Ar37" == pair42.pairID().symbol() );

      // quantum numbers
      const auto numbers42 = channel42.quantumNumbers();
      CHECK( 2 == numbers42.orbitalAngularMomentum() );
      CHECK( 1.5 == numbers42.spin() );
      CHECK( 2.5 == numbers42.totalAngularMomentum() );
      CHECK( +1 == numbers42.parity() );
      CHECK( "{2,3/2,5/2+}" == numbers42.toString() );

      // radii
      const auto radii42 = channel42.radii();
      CHECK( 0.4885639 == Approx( radii42.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4885639 == Approx( radii42.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4885639 == Approx( radii42.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel42.boundaryCondition() );

      // Q value
      CHECK( 1.747660e+6 == Approx( channel42.Q().value ) );

      // resonance table
      auto table4 = spingroup4.resonanceTable();

      CHECK( 3 == table4.numberChannels() ); // 3 normal channel + 1 eliminated
      CHECK( 50 == table4.numberResonances() );

      auto energies4 = table4.energies();
      CHECK( 3.453854e+4 == Approx( energies4.front().value ) );
      CHECK( 1.530000e+6 == Approx( energies4.back().value ) );

      auto resonances4 = table4.resonances();
      CHECK( 3.453854e+4 == Approx( resonances4.front().energy().value ) );
      CHECK( 1.530000e+6 == Approx( resonances4.back().energy().value ) );
      CHECK( 3 == resonances4.front().widths().size() );
      CHECK( 3 == resonances4.back().widths().size() );
      CHECK( std::sqrt( 3.772489e-1 / 2. / channel40.penetrability( 3.453854e+4 * electronVolt ) )
             == Approx( resonances4.front().widths()[0].value ) );
      CHECK( std::sqrt( 0.0 / 2. / channel41.penetrability( 3.453854e+4 * electronVolt ) )
             == Approx( resonances4.front().widths()[1].value ) );
      CHECK( std::sqrt( 1.000300e-3 / 2. / channel42.penetrability( 3.453854e+4 * electronVolt ) )
             == Approx( resonances4.front().widths()[2].value ) );
      CHECK( std::sqrt( 45. / 2. / channel40.penetrability( 1.530000e+6 * electronVolt ) )
             == Approx( resonances4.back().widths()[0].value ) );
      CHECK( std::sqrt( 1.500000e-4 / 2. / channel41.penetrability( 1.530000e+6 * electronVolt ) )
             == Approx( resonances4.back().widths()[1].value ) );
      CHECK( std::sqrt( 1.e-3 / 2. / channel42.penetrability( 1.530000e+6 * electronVolt ) )
             == Approx( resonances4.back().widths()[2].value ) );
      CHECK( std::sqrt( 1.306229e-1 / 2. ) == Approx( resonances4.front().eliminatedWidth().value ) );
      CHECK( std::sqrt( 1. / 2. ) == Approx( resonances4.back().eliminatedWidth().value ) );
    } // THEN

    THEN( "cross sections can be reconstructed" ) {

      ReactionID elas( "n,Ca40->n,Ca40" );
      ReactionID prot( "n,Ca40->h1,K40" );
      ReactionID alph( "n,Ca40->he4,Ar37" );
      ReactionID capt( "n,Ca40->capture" );
      std::map< ReactionID, CrossSection > xs;

      // this value previously produced NaN cross section values
      xs = resonances( 542087 * electronVolt );
      CHECK( 4 == xs.size() );
      CHECK( 0.5977843018 == Approx( xs[ elas ].value ) );
      CHECK( 0. == Approx( xs[ prot ].value ) );
      CHECK( 9.8241e-6 == Approx( xs[ alph ].value ) );
      CHECK( 3.54988e-5 == Approx( xs[ capt ].value ) );
    } // THEN
  } // GIVEN

  GIVEN( "valid ENDF data for Cl35" ) {

    std::string string = Cl35();
    auto begin = string.begin();
    auto end = string.end();
    long lineNumber = 1;

    njoy::ENDFtk::HeadRecord head( begin, end, lineNumber );
    njoy::ENDFtk::section::Type< 2, 151 > endf( head, begin, end, lineNumber, 1725 );
    ResonanceRange endfResonanceRange = endf.isotopes().front().resonanceRanges().front();

    auto resonances = fromENDF( endfResonanceRange, neutronMass, elementaryCharge, ParticleID( "n" ), ParticleID( "Cl35" ) );

    THEN( "the appropriate CompoundSystem is returned" ) {

      CHECK( true == resonances.isResolved() );
      CHECK( false == resonances.isUnresolved() );
      CHECK( 1e-5 == Approx( resonances.lowerEnergy().value ) );
      CHECK( 1.2e+6 == Approx( resonances.upperEnergy().value ) );

      auto compoundsystem = std::get< CompoundSystem< ReichMoore, ShiftFactor > >( resonances.compoundSystem() );

      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // content verification
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      // spin groups
      auto spingroups = compoundsystem.spinGroups();
      CHECK( 6 == spingroups.size() );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup0 = spingroups[0];

      // channels
      auto channels0 = spingroup0.channels();

      CHECK( 1 == channels0.size() ); // 2 normal channel + 1 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel00 = std::get< Channel< Neutron > >( channels0[0] );

      // incident particle pair
      const auto incident00 = channel00.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident00.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident00.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident00.particle().spin() ) );
      CHECK( +1 == incident00.particle().parity() );
      CHECK( 34.66845 * 1.008664 == Approx( incident00.residual().mass().value ) );
      CHECK( 17.0 * 1.602e-19 == Approx( incident00.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident00.residual().spin() ) );
      CHECK( +1 == incident00.residual().parity() );
      CHECK( "n,Cl35" == incident00.pairID().symbol() );

      // particle pair
      const auto pair00 = channel00.particlePair();
      CHECK( 1.008664 == Approx( pair00.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair00.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair00.particle().spin() ) );
      CHECK( +1 == pair00.particle().parity() );
      CHECK( 34.66845 * 1.008664 == Approx( pair00.residual().mass().value ) );
      CHECK( 17.0 * 1.602e-19 == Approx( pair00.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair00.residual().spin() ) );
      CHECK( +1 == pair00.residual().parity() );
      CHECK( "n,Cl35" == pair00.pairID().symbol() );

      // quantum numbers
      const auto numbers00 = channel00.quantumNumbers();
      CHECK( 1 == numbers00.orbitalAngularMomentum() );
      CHECK( 1. == numbers00.spin() );
      CHECK( 0. == numbers00.totalAngularMomentum() );
      CHECK( -1 == numbers00.parity() );
      CHECK( "{1,1,0-}" == numbers00.toString() );

      // radii
      const auto radii00 = channel00.radii();
      CHECK( 0.4822220 == Approx( radii00.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4822220 == Approx( radii00.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4888750 == Approx( radii00.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel00.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel00.Q().value ) );

      // resonance table
      auto table0 = spingroup0.resonanceTable();

      CHECK( 1 == table0.numberChannels() ); // 2 normal channel + 1 eliminated
      CHECK( 9 == table0.numberResonances() );

      auto energies0 = table0.energies();
      CHECK( 2.239640e+4 == Approx( energies0.front().value ) );
      CHECK( 5.478545e+5 == Approx( energies0.back().value ) );

      auto resonances0 = table0.resonances();
      CHECK( 2.239640e+4 == Approx( resonances0.front().energy().value ) );
      CHECK( 5.478545e+5 == Approx( resonances0.back().energy().value ) );
      CHECK( 1 == resonances0.front().widths().size() );
      CHECK( 1 == resonances0.back().widths().size() );
      CHECK( std::sqrt( .9663670 / 2. / channel00.penetrability( 2.239640e+4 * electronVolt ) )
             == Approx( resonances0.front().widths()[0].value ) );
      CHECK( std::sqrt( 7.640130e+2 / 2. / channel00.penetrability( 5.478545e+5 * electronVolt ) )
             == Approx( resonances0.back().widths()[0].value ) );
      CHECK( std::sqrt( 1.724800 / 2. ) == Approx( resonances0.front().eliminatedWidth().value ) );
      CHECK( std::sqrt( 0.86 / 2. ) == Approx( resonances0.back().eliminatedWidth().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup1 = spingroups[1];

      // channels
      auto channels1 = spingroup1.channels();

      CHECK( 4 == channels1.size() ); // 4 normal channel + 1 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel10 = std::get< Channel< Neutron > >( channels1[0] );

      // incident particle pair
      const auto incident10 = channel10.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident10.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident10.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident10.particle().spin() ) );
      CHECK( +1 == incident10.particle().parity() );
      CHECK( 34.66845 * 1.008664 == Approx( incident10.residual().mass().value ) );
      CHECK( 17.0 * 1.602e-19 == Approx( incident10.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident10.residual().spin() ) );
      CHECK( +1 == incident10.residual().parity() );
      CHECK( "n,Cl35" == incident10.pairID().symbol() );

      // particle pair
      const auto pair10 = channel10.particlePair();
      CHECK( 1.008664 == Approx( pair10.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair10.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair10.particle().spin() ) );
      CHECK( +1 == pair10.particle().parity() );
      CHECK( 34.66845 * 1.008664 == Approx( pair10.residual().mass().value ) );
      CHECK( 17.0 * 1.602e-19 == Approx( pair10.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair10.residual().spin() ) );
      CHECK( +1 == pair10.residual().parity() );
      CHECK( "n,Cl35" == pair10.pairID().symbol() );

      // quantum numbers
      const auto numbers10 = channel10.quantumNumbers();
      CHECK( 1 == numbers10.orbitalAngularMomentum() );
      CHECK( 1. == numbers10.spin() );
      CHECK( 1. == numbers10.totalAngularMomentum() );
      CHECK( -1 == numbers10.parity() );
      CHECK( "{1,1,1-}" == numbers10.toString() );

      // radii
      const auto radii10 = channel10.radii();
      CHECK( 0.4822220 == Approx( radii10.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4822220 == Approx( radii10.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4888750 == Approx( radii10.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel10.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel10.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1, channel 1: n,p
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel11 = std::get< Channel< ChargedParticle > >( channels1[1] );

      // incident particle pair
      const auto incident11 = channel11.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident11.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident11.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident11.particle().spin() ) );
      CHECK( +1 == incident11.particle().parity() );
      CHECK( 34.66845 * 1.008664 == Approx( incident11.residual().mass().value ) );
      CHECK( 17.0 * 1.602e-19 == Approx( incident11.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident11.residual().spin() ) );
      CHECK( +1 == incident11.residual().parity() );
      CHECK( "n,Cl35" == incident11.pairID().symbol() );

      // particle pair
      const auto pair11 = channel11.particlePair();
      CHECK( 0.9986235 * 1.008664 == Approx( pair11.particle().mass().value ) );
      CHECK( 1.0 * 1.602e-19 == Approx( pair11.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair11.particle().spin() ) );
      CHECK( +1 == pair11.particle().parity() );
      CHECK( 34.66863 * 1.008664 == Approx( pair11.residual().mass().value ) );
      CHECK( 16.0 * 1.602e-19 == Approx( pair11.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair11.residual().spin() ) );
      CHECK( +1 == pair11.residual().parity() );
      CHECK( "p,S35" == pair11.pairID().symbol() );

      // quantum numbers
      const auto numbers11 = channel11.quantumNumbers();
      CHECK( 1 == numbers11.orbitalAngularMomentum() );
      CHECK( 1. == numbers11.spin() );
      CHECK( 1. == numbers11.totalAngularMomentum() );
      CHECK( -1 == numbers11.parity() );
      CHECK( "{1,1,1-}" == numbers11.toString() );

      // radii
      const auto radii11 = channel11.radii();
      CHECK( 0.4822220 == Approx( radii11.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4822220 == Approx( radii11.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4888750 == Approx( radii11.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel11.boundaryCondition() );

      // Q value
      CHECK( 6.152200e+5 == Approx( channel11.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1, channel 2: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel12 = std::get< Channel< Neutron > >( channels1[2] );

      // incident particle pair
      const auto incident12 = channel12.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident12.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident12.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident12.particle().spin() ) );
      CHECK( +1 == incident12.particle().parity() );
      CHECK( 34.66845 * 1.008664 == Approx( incident12.residual().mass().value ) );
      CHECK( 17.0 * 1.602e-19 == Approx( incident12.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident12.residual().spin() ) );
      CHECK( +1 == incident12.residual().parity() );
      CHECK( "n,Cl35" == incident12.pairID().symbol() );

      // particle pair
      const auto pair12 = channel12.particlePair();
      CHECK( 1.008664 == Approx( pair12.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair12.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair12.particle().spin() ) );
      CHECK( +1 == pair12.particle().parity() );
      CHECK( 34.66845 * 1.008664 == Approx( pair12.residual().mass().value ) );
      CHECK( 17.0 * 1.602e-19 == Approx( pair12.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair12.residual().spin() ) );
      CHECK( +1 == pair12.residual().parity() );
      CHECK( "n,Cl35" == pair12.pairID().symbol() );

      // quantum numbers
      const auto numbers12 = channel12.quantumNumbers();
      CHECK( 1 == numbers12.orbitalAngularMomentum() );
      CHECK( 2. == numbers12.spin() );
      CHECK( 1. == numbers12.totalAngularMomentum() );
      CHECK( -1 == numbers12.parity() );
      CHECK( "{1,2,1-}" == numbers12.toString() );

      // radii
      const auto radii12 = channel12.radii();
      CHECK( 0.4822220 == Approx( radii12.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4822220 == Approx( radii12.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4888750 == Approx( radii12.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel12.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel12.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1, channel 3: n,p
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel13 = std::get< Channel< ChargedParticle > >( channels1[3] );

      // incident particle pair
      const auto incident13 = channel13.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident13.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident13.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident13.particle().spin() ) );
      CHECK( +1 == incident13.particle().parity() );
      CHECK( 34.66845 * 1.008664 == Approx( incident13.residual().mass().value ) );
      CHECK( 17.0 * 1.602e-19 == Approx( incident13.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident13.residual().spin() ) );
      CHECK( +1 == incident13.residual().parity() );
      CHECK( "n,Cl35" == incident13.pairID().symbol() );

      // particle pair
      const auto pair13 = channel13.particlePair();
      CHECK( 0.9986235 * 1.008664 == Approx( pair13.particle().mass().value ) );
      CHECK( 1.0 * 1.602e-19 == Approx( pair13.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair13.particle().spin() ) );
      CHECK( +1 == pair13.particle().parity() );
      CHECK( 34.66863 * 1.008664 == Approx( pair13.residual().mass().value ) );
      CHECK( 16.0 * 1.602e-19 == Approx( pair13.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair13.residual().spin() ) );
      CHECK( +1 == pair13.residual().parity() );
      CHECK( "p,S35" == pair13.pairID().symbol() );

      // quantum numbers
      const auto numbers13 = channel13.quantumNumbers();
      CHECK( 1 == numbers13.orbitalAngularMomentum() );
      CHECK( 2. == numbers13.spin() );
      CHECK( 1. == numbers13.totalAngularMomentum() );
      CHECK( -1 == numbers13.parity() );
      CHECK( "{1,2,1-}" == numbers13.toString() );

      // radii
      const auto radii13 = channel13.radii();
      CHECK( 0.4822220 == Approx( radii13.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4822220 == Approx( radii13.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4888750 == Approx( radii13.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel13.boundaryCondition() );

      // Q value
      CHECK( 6.152200e+5 == Approx( channel13.Q().value ) );

      // resonance table
      auto table1 = spingroup1.resonanceTable();

      CHECK( 4 == table1.numberChannels() ); // 2 normal channel + 1 eliminated
      CHECK( 38 + 18 == table1.numberResonances() ); // 38 for l=1, 18 for l=2

      auto energies1 = table1.energies();
      CHECK( 4250.762 == Approx( energies1.front().value ) );
      CHECK( 1.435502e+6 == Approx( energies1.back().value ) );

      auto resonances1 = table1.resonances();
      CHECK( 4250.762 == Approx( resonances1.front().energy().value ) );
      CHECK( 1.435502e+6 == Approx( resonances1.back().energy().value ) );
      CHECK( 4 == resonances1.front().widths().size() );
      CHECK( 4 == resonances1.back().widths().size() );
      CHECK( std::sqrt( .628 / 2. / channel10.penetrability( 4250.762 * electronVolt ) )
             == Approx( resonances1.front().widths()[0].value ) );
      CHECK( std::sqrt( .23 / 2. / channel11.penetrability( 4250.762 * electronVolt ) )
             == Approx( resonances1.front().widths()[1].value ) );
      CHECK( 0. == Approx( resonances1.front().widths()[2].value ) );
      CHECK( 0. == Approx( resonances1.front().widths()[3].value ) );
      CHECK( std::sqrt( 5.365630e+3 / 2. / channel10.penetrability( 1.435502e+6 * electronVolt ) )
             == Approx( resonances1.back().widths()[0].value ) );
      CHECK( 0. == Approx( resonances1.back().widths()[1].value ) );
      CHECK( 0. == Approx( resonances1.back().widths()[2].value ) );
      CHECK( 0. == Approx( resonances1.back().widths()[3].value ) );
      CHECK( std::sqrt( 0.472 / 2. ) == Approx( resonances1.front().eliminatedWidth().value ) );
      CHECK( std::sqrt( 0.86 / 2. ) == Approx( resonances1.back().eliminatedWidth().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 2
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup2 = spingroups[2];

      // channels
      auto channels2 = spingroup2.channels();

      CHECK( 2 == channels2.size() ); // 2 normal channel + 1 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 2, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel20 = std::get< Channel< Neutron > >( channels2[0] );

      // incident particle pair
      const auto incident20 = channel20.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident20.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident20.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident20.particle().spin() ) );
      CHECK( +1 == incident20.particle().parity() );
      CHECK( 34.66845 * 1.008664 == Approx( incident20.residual().mass().value ) );
      CHECK( 17.0 * 1.602e-19 == Approx( incident20.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident20.residual().spin() ) );
      CHECK( +1 == incident20.residual().parity() );
      CHECK( "n,Cl35" == incident20.pairID().symbol() );

      // particle pair
      const auto pair20 = channel20.particlePair();
      CHECK( 1.008664 == Approx( pair20.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair20.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair20.particle().spin() ) );
      CHECK( +1 == pair20.particle().parity() );
      CHECK( 34.66845 * 1.008664 == Approx( pair20.residual().mass().value ) );
      CHECK( 17.0 * 1.602e-19 == Approx( pair20.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair20.residual().spin() ) );
      CHECK( +1 == pair20.residual().parity() );
      CHECK( "n,Cl35" == pair20.pairID().symbol() );

      // quantum numbers
      const auto numbers20 = channel20.quantumNumbers();
      CHECK( 0 == numbers20.orbitalAngularMomentum() );
      CHECK( 1. == numbers20.spin() );
      CHECK( 1. == numbers20.totalAngularMomentum() );
      CHECK( +1 == numbers20.parity() );
      CHECK( "{0,1,1+}" == numbers20.toString() );

      // radii
      const auto radii20 = channel20.radii();
      CHECK( 0.4822220 == Approx( radii20.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4822220 == Approx( radii20.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.3667980 == Approx( radii20.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel20.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel20.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 2, channel 1: n,p
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel21 = std::get< Channel< ChargedParticle > >( channels2[1] );

      // incident particle pair
      const auto incident21 = channel21.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident21.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident21.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident21.particle().spin() ) );
      CHECK( +1 == incident21.particle().parity() );
      CHECK( 34.66845 * 1.008664 == Approx( incident21.residual().mass().value ) );
      CHECK( 17.0 * 1.602e-19 == Approx( incident21.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident21.residual().spin() ) );
      CHECK( +1 == incident21.residual().parity() );
      CHECK( "n,Cl35" == incident21.pairID().symbol() );

      // particle pair
      const auto pair21 = channel21.particlePair();
      CHECK( 0.9986235 * 1.008664 == Approx( pair21.particle().mass().value ) );
      CHECK( 1.0 * 1.602e-19 == Approx( pair21.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair21.particle().spin() ) );
      CHECK( +1 == pair21.particle().parity() );
      CHECK( 34.66863 * 1.008664 == Approx( pair21.residual().mass().value ) );
      CHECK( 16.0 * 1.602e-19 == Approx( pair21.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair21.residual().spin() ) );
      CHECK( +1 == pair21.residual().parity() );
      CHECK( "p,S35" == pair21.pairID().symbol() );

      // quantum numbers
      const auto numbers21 = channel21.quantumNumbers();
      CHECK( 0 == numbers21.orbitalAngularMomentum() );
      CHECK( 1. == numbers21.spin() );
      CHECK( 1. == numbers21.totalAngularMomentum() );
      CHECK( +1 == numbers21.parity() );
      CHECK( "{0,1,1+}" == numbers21.toString() );

      // radii
      const auto radii21 = channel21.radii();
      CHECK( 0.4822220 == Approx( radii21.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4822220 == Approx( radii21.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.3667980 == Approx( radii21.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel21.boundaryCondition() );

      // Q value
      CHECK( 6.152200e+5 == Approx( channel21.Q().value ) );

      // resonance table
      auto table2 = spingroup2.resonanceTable();

      CHECK( 2 == table2.numberChannels() ); // 2 normal channel + 1 eliminated
      CHECK( 23 == table2.numberResonances() );

      auto energies2 = table2.energies();
      CHECK( 5.493200e+4 == Approx( energies2.front().value ) );
      CHECK( 1.205687e+6 == Approx( energies2.back().value ) );

      auto resonances2 = table2.resonances();
      CHECK( 5.493200e+4 == Approx( resonances2.front().energy().value ) );
      CHECK( 1.205687e+6 == Approx( resonances2.back().energy().value ) );
      CHECK( 2 == resonances2.front().widths().size() );
      CHECK( 2 == resonances2.back().widths().size() );
      CHECK( std::sqrt( 46.44240 / 2. / channel20.penetrability( 5.493200e+4 * electronVolt ) )
             == Approx( resonances2.front().widths()[0].value ) );
      CHECK( 0. == Approx( resonances2.front().widths()[1].value ) );
      CHECK( std::sqrt( 2.179040e+2 / 2. / channel20.penetrability( 6.823616e+4 * electronVolt ) )
             == Approx( resonances2[1].widths()[0].value ) );
      CHECK( std::sqrt( 1e-5 / 2. / channel21.penetrability( 6.823616e+4 * electronVolt ) )
             == Approx( resonances2[1].widths()[1].value ) );
      CHECK( std::sqrt( 6.425840e+2 / 2. / channel20.penetrability( 1.205687e+6 * electronVolt ) )
             == Approx( resonances2.back().widths()[0].value ) );
      CHECK( 0. == Approx( resonances2.back().widths()[1].value ) );
      CHECK( std::sqrt( 0.36726 / 2. ) == Approx( resonances2.front().eliminatedWidth().value ) );
      CHECK( std::sqrt( 0.606 / 2. ) == Approx( resonances2.back().eliminatedWidth().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 3
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup3 = spingroups[3];

      // channels
      auto channels3 = spingroup3.channels();

      CHECK( 4 == channels1.size() ); // 4 normal channel + 1 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 3, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel30 = std::get< Channel< Neutron > >( channels3[0] );

      // incident particle pair
      const auto incident30 = channel30.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident30.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident30.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident30.particle().spin() ) );
      CHECK( +1 == incident30.particle().parity() );
      CHECK( 34.66845 * 1.008664 == Approx( incident30.residual().mass().value ) );
      CHECK( 17.0 * 1.602e-19 == Approx( incident30.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident30.residual().spin() ) );
      CHECK( +1 == incident30.residual().parity() );
      CHECK( "n,Cl35" == incident30.pairID().symbol() );

      // particle pair
      const auto pair30 = channel30.particlePair();
      CHECK( 1.008664 == Approx( pair30.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair30.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair30.particle().spin() ) );
      CHECK( +1 == pair30.particle().parity() );
      CHECK( 34.66845 * 1.008664 == Approx( pair30.residual().mass().value ) );
      CHECK( 17.0 * 1.602e-19 == Approx( pair30.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair30.residual().spin() ) );
      CHECK( +1 == pair30.residual().parity() );
      CHECK( "n,Cl35" == pair30.pairID().symbol() );

      // quantum numbers
      const auto numbers30 = channel30.quantumNumbers();
      CHECK( 1 == numbers30.orbitalAngularMomentum() );
      CHECK( 1. == numbers30.spin() );
      CHECK( 2. == numbers30.totalAngularMomentum() );
      CHECK( -1 == numbers30.parity() );
      CHECK( "{1,1,2-}" == numbers30.toString() );

      // radii
      const auto radii30 = channel30.radii();
      CHECK( 0.4822220 == Approx( radii30.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4822220 == Approx( radii30.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4888750 == Approx( radii30.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel30.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel30.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 3, channel 1: n,p
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel31 = std::get< Channel< ChargedParticle > >( channels3[1] );

      // incident particle pair
      const auto incident31 = channel31.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident31.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident31.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident31.particle().spin() ) );
      CHECK( +1 == incident31.particle().parity() );
      CHECK( 34.66845 * 1.008664 == Approx( incident31.residual().mass().value ) );
      CHECK( 17.0 * 1.602e-19 == Approx( incident31.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident31.residual().spin() ) );
      CHECK( +1 == incident31.residual().parity() );
      CHECK( "n,Cl35" == incident31.pairID().symbol() );

      // particle pair
      const auto pair31 = channel31.particlePair();
      CHECK( 0.9986235 * 1.008664 == Approx( pair31.particle().mass().value ) );
      CHECK( 1.0 * 1.602e-19 == Approx( pair31.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair31.particle().spin() ) );
      CHECK( +1 == pair31.particle().parity() );
      CHECK( 34.66863 * 1.008664 == Approx( pair31.residual().mass().value ) );
      CHECK( 16.0 * 1.602e-19 == Approx( pair31.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair31.residual().spin() ) );
      CHECK( +1 == pair31.residual().parity() );
      CHECK( "p,S35" == pair31.pairID().symbol() );

      // quantum numbers
      const auto numbers31 = channel31.quantumNumbers();
      CHECK( 1 == numbers31.orbitalAngularMomentum() );
      CHECK( 1. == numbers31.spin() );
      CHECK( 2. == numbers31.totalAngularMomentum() );
      CHECK( -1 == numbers31.parity() );
      CHECK( "{1,1,2-}" == numbers31.toString() );

      // radii
      const auto radii31 = channel31.radii();
      CHECK( 0.4822220 == Approx( radii31.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4822220 == Approx( radii31.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4888750 == Approx( radii31.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel31.boundaryCondition() );

      // Q value
      CHECK( 6.152200e+5 == Approx( channel31.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 3, channel 2: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel32 = std::get< Channel< Neutron > >( channels3[2] );

      // incident particle pair
      const auto incident32 = channel32.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident32.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident32.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident32.particle().spin() ) );
      CHECK( +1 == incident32.particle().parity() );
      CHECK( 34.66845 * 1.008664 == Approx( incident32.residual().mass().value ) );
      CHECK( 17.0 * 1.602e-19 == Approx( incident32.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident32.residual().spin() ) );
      CHECK( +1 == incident32.residual().parity() );
      CHECK( "n,Cl35" == incident32.pairID().symbol() );

      // particle pair
      const auto pair32 = channel32.particlePair();
      CHECK( 1.008664 == Approx( pair32.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair32.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair32.particle().spin() ) );
      CHECK( +1 == pair32.particle().parity() );
      CHECK( 34.66845 * 1.008664 == Approx( pair32.residual().mass().value ) );
      CHECK( 17.0 * 1.602e-19 == Approx( pair32.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair32.residual().spin() ) );
      CHECK( +1 == pair32.residual().parity() );
      CHECK( "n,Cl35" == pair32.pairID().symbol() );

      // quantum numbers
      const auto numbers32 = channel32.quantumNumbers();
      CHECK( 1 == numbers32.orbitalAngularMomentum() );
      CHECK( 2. == numbers32.spin() );
      CHECK( 2. == numbers32.totalAngularMomentum() );
      CHECK( -1 == numbers32.parity() );
      CHECK( "{1,2,2-}" == numbers32.toString() );

      // radii
      const auto radii32 = channel32.radii();
      CHECK( 0.4822220 == Approx( radii32.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4822220 == Approx( radii32.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4888750 == Approx( radii32.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel32.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel32.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 3, channel 3: n,p
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel33 = std::get< Channel< ChargedParticle > >( channels3[3] );

      // incident particle pair
      const auto incident33 = channel33.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident33.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident33.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident33.particle().spin() ) );
      CHECK( +1 == incident33.particle().parity() );
      CHECK( 34.66845 * 1.008664 == Approx( incident33.residual().mass().value ) );
      CHECK( 17.0 * 1.602e-19 == Approx( incident33.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident33.residual().spin() ) );
      CHECK( +1 == incident33.residual().parity() );
      CHECK( "n,Cl35" == incident33.pairID().symbol() );

      // particle pair
      const auto pair33 = channel33.particlePair();
      CHECK( 0.9986235 * 1.008664 == Approx( pair33.particle().mass().value ) );
      CHECK( 1.0 * 1.602e-19 == Approx( pair33.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair33.particle().spin() ) );
      CHECK( +1 == pair33.particle().parity() );
      CHECK( 34.66863 * 1.008664 == Approx( pair33.residual().mass().value ) );
      CHECK( 16.0 * 1.602e-19 == Approx( pair33.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair33.residual().spin() ) );
      CHECK( +1 == pair33.residual().parity() );
      CHECK( "p,S35" == pair33.pairID().symbol() );

      // quantum numbers
      const auto numbers33 = channel33.quantumNumbers();
      CHECK( 1 == numbers33.orbitalAngularMomentum() );
      CHECK( 2. == numbers33.spin() );
      CHECK( 2. == numbers33.totalAngularMomentum() );
      CHECK( -1 == numbers33.parity() );
      CHECK( "{1,2,2-}" == numbers33.toString() );

      // radii
      const auto radii33 = channel33.radii();
      CHECK( 0.4822220 == Approx( radii33.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4822220 == Approx( radii33.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4888750 == Approx( radii33.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel33.boundaryCondition() );

      // Q value
      CHECK( 6.152200e+5 == Approx( channel33.Q().value ) );

      // resonance table
      auto table3 = spingroup3.resonanceTable();

      CHECK( 4 == table3.numberChannels() ); // 4 normal channel + 1 eliminated
      CHECK( 67 + 28 == table3.numberResonances() );

      auto energies3 = table3.energies();
      CHECK( -3.369334e+5 == Approx( energies3.front().value ) );
      CHECK( 1.441365e+6 == Approx( energies3.back().value ) );

      auto resonances3 = table3.resonances();
      CHECK( -3.369334e+5 == Approx( resonances3.front().energy().value ) );
      CHECK( 1.441365e+6 == Approx( resonances3.back().energy().value ) );
      CHECK( 4 == resonances3.front().widths().size() );
      CHECK( 4 == resonances3.back().widths().size() );
      CHECK( std::sqrt( 3.820180e+4 / 2. / channel30.penetrability( -3.369334e+5 * electronVolt ) )
             == Approx( resonances3.front().widths()[0].value ) );
      CHECK( 0. == Approx( resonances3.front().widths()[1].value ) );
      CHECK( 0. == Approx( resonances3.front().widths()[2].value ) );
      CHECK( 0. == Approx( resonances3.front().widths()[3].value ) );
      CHECK( 0. == Approx( resonances3.back().widths()[0].value ) );
      CHECK( 0. == Approx( resonances3.back().widths()[1].value ) );
      CHECK( std::sqrt( 1.608740e+3 / 2. / channel32.penetrability( 1.441365e+6 * electronVolt ) )
             == Approx( resonances3.back().widths()[2].value ) );
      CHECK( 0. == Approx( resonances3.back().widths()[3].value ) );
      CHECK( std::sqrt( 0.53401 / 2. ) == Approx( resonances3.front().eliminatedWidth().value ) );
      CHECK( std::sqrt( 0.86 / 2. ) == Approx( resonances3.back().eliminatedWidth().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 4
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup4 = spingroups[4];

      // channels
      auto channels4 = spingroup4.channels();

      CHECK( 2 == channels4.size() ); // 2 normal channel + 1 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 4, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel40 = std::get< Channel< Neutron > >( channels4[0] );

      // incident particle pair
      const auto incident40 = channel40.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident40.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident40.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident40.particle().spin() ) );
      CHECK( +1 == incident40.particle().parity() );
      CHECK( 34.66845 * 1.008664 == Approx( incident40.residual().mass().value ) );
      CHECK( 17.0 * 1.602e-19 == Approx( incident40.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident40.residual().spin() ) );
      CHECK( +1 == incident40.residual().parity() );
      CHECK( "n,Cl35" == incident40.pairID().symbol() );

      // particle pair
      const auto pair40 = channel40.particlePair();
      CHECK( 1.008664 == Approx( pair40.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair40.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair40.particle().spin() ) );
      CHECK( +1 == pair40.particle().parity() );
      CHECK( 34.66845 * 1.008664 == Approx( pair40.residual().mass().value ) );
      CHECK( 17.0 * 1.602e-19 == Approx( pair40.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair40.residual().spin() ) );
      CHECK( +1 == pair40.residual().parity() );
      CHECK( "n,Cl35" == pair40.pairID().symbol() );

      // quantum numbers
      const auto numbers40 = channel40.quantumNumbers();
      CHECK( 0 == numbers40.orbitalAngularMomentum() );
      CHECK( 2. == numbers40.spin() );
      CHECK( 2. == numbers40.totalAngularMomentum() );
      CHECK( +1 == numbers40.parity() );
      CHECK( "{0,2,2+}" == numbers40.toString() );

      // radii
      const auto radii40 = channel40.radii();
      CHECK( 0.4822220 == Approx( radii40.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4822220 == Approx( radii40.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.3667980 == Approx( radii40.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel40.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel40.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 4, channel 1: n,p
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel41 = std::get< Channel< ChargedParticle > >( channels4[1] );

      // incident particle pair
      const auto incident41 = channel41.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident41.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident41.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident41.particle().spin() ) );
      CHECK( +1 == incident41.particle().parity() );
      CHECK( 34.66845 * 1.008664 == Approx( incident41.residual().mass().value ) );
      CHECK( 17.0 * 1.602e-19 == Approx( incident41.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident41.residual().spin() ) );
      CHECK( +1 == incident41.residual().parity() );
      CHECK( "n,Cl35" == incident41.pairID().symbol() );

      // particle pair
      const auto pair41 = channel41.particlePair();
      CHECK( 0.9986235 * 1.008664 == Approx( pair41.particle().mass().value ) );
      CHECK( 1.0 * 1.602e-19 == Approx( pair41.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair41.particle().spin() ) );
      CHECK( +1 == pair41.particle().parity() );
      CHECK( 34.66863 * 1.008664 == Approx( pair41.residual().mass().value ) );
      CHECK( 16.0 * 1.602e-19 == Approx( pair41.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair41.residual().spin() ) );
      CHECK( +1 == pair41.residual().parity() );
      CHECK( "p,S35" == pair41.pairID().symbol() );

      // quantum numbers
      const auto numbers41 = channel41.quantumNumbers();
      CHECK( 0 == numbers41.orbitalAngularMomentum() );
      CHECK( 2. == numbers41.spin() );
      CHECK( 2. == numbers41.totalAngularMomentum() );
      CHECK( +1 == numbers41.parity() );
      CHECK( "{0,2,2+}" == numbers41.toString() );

      // radii
      const auto radii41 = channel41.radii();
      CHECK( 0.4822220 == Approx( radii41.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4822220 == Approx( radii41.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.3667980 == Approx( radii41.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel41.boundaryCondition() );

      // Q value
      CHECK( 6.152200e+5 == Approx( channel41.Q().value ) );

      // resonance table
      auto table4 = spingroup4.resonanceTable();

      CHECK( 2 == table4.numberChannels() ); // 2 normal channel + 1 eliminated
      CHECK( 32 == table4.numberResonances() );

      auto energies4 = table4.energies();
      CHECK( -1.806500e+2 == Approx( energies4.front().value ) );
      CHECK( 7.563145e+6 == Approx( energies4.back().value ) );

      auto resonances4 = table4.resonances();
      CHECK( -1.806500e+2 == Approx( resonances4.front().energy().value ) );
      CHECK( 7.563145e+6 == Approx( resonances4.back().energy().value ) );
      CHECK( 2 == resonances4.front().widths().size() );
      CHECK( 2 == resonances4.back().widths().size() );
      CHECK( std::sqrt( 13.27700 / 2. / channel40.penetrability( -1.806500e+2 * electronVolt ) )
             == Approx( resonances4.front().widths()[0].value ) );
      CHECK( std::sqrt( 5.992300e-3 / 2. / channel41.penetrability( -1.806500e+2 * electronVolt ) )
             == Approx( resonances4.front().widths()[1].value ) );
      CHECK( std::sqrt( 6.219050e+5 / 2. / channel40.penetrability( 7.563145e+6 * electronVolt ) )
             == Approx( resonances4.back().widths()[0].value ) );
      CHECK( std::sqrt( 1000. / 2. / channel41.penetrability( 7.563145e+6 * electronVolt ) )
             == Approx( resonances4.back().widths()[1].value ) );
      CHECK( std::sqrt( 0.53015 / 2. ) == Approx( resonances4.front().eliminatedWidth().value ) );
      CHECK( std::sqrt( 0.38398 / 2. ) == Approx( resonances4.back().eliminatedWidth().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 5
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup5 = spingroups[5];

      // channels
      auto channels5 = spingroup5.channels();

      CHECK( 2 == channels5.size() ); // 2 normal channel + 1 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 5, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel50 = std::get< Channel< Neutron > >( channels5[0] );

      // incident particle pair
      const auto incident50 = channel50.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident50.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident50.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident50.particle().spin() ) );
      CHECK( +1 == incident50.particle().parity() );
      CHECK( 34.66845 * 1.008664 == Approx( incident50.residual().mass().value ) );
      CHECK( 17.0 * 1.602e-19 == Approx( incident50.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident50.residual().spin() ) );
      CHECK( +1 == incident50.residual().parity() );
      CHECK( "n,Cl35" == incident50.pairID().symbol() );

      // particle pair
      const auto pair50 = channel50.particlePair();
      CHECK( 1.008664 == Approx( pair50.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair50.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair50.particle().spin() ) );
      CHECK( +1 == pair50.particle().parity() );
      CHECK( 34.66845 * 1.008664 == Approx( pair50.residual().mass().value ) );
      CHECK( 17.0 * 1.602e-19 == Approx( pair50.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair50.residual().spin() ) );
      CHECK( +1 == pair50.residual().parity() );
      CHECK( "n,Cl35" == pair50.pairID().symbol() );

      // quantum numbers
      const auto numbers50 = channel50.quantumNumbers();
      CHECK( 1 == numbers50.orbitalAngularMomentum() );
      CHECK( 2. == numbers50.spin() );
      CHECK( 3. == numbers50.totalAngularMomentum() );
      CHECK( -1 == numbers50.parity() );
      CHECK( "{1,2,3-}" == numbers50.toString() );

      // radii
      const auto radii50 = channel50.radii();
      CHECK( 0.4822220 == Approx( radii50.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4822220 == Approx( radii50.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4888750 == Approx( radii50.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel50.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel50.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 5, channel 1: n,p
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel51 = std::get< Channel< ChargedParticle > >( channels5[1] );

      // incident particle pair
      const auto incident51 = channel51.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident51.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident51.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident51.particle().spin() ) );
      CHECK( +1 == incident51.particle().parity() );
      CHECK( 34.66845 * 1.008664 == Approx( incident51.residual().mass().value ) );
      CHECK( 17.0 * 1.602e-19 == Approx( incident51.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident51.residual().spin() ) );
      CHECK( +1 == incident51.residual().parity() );
      CHECK( "n,Cl35" == incident51.pairID().symbol() );

      // particle pair
      const auto pair51 = channel51.particlePair();
      CHECK( 0.9986235 * 1.008664 == Approx( pair51.particle().mass().value ) );
      CHECK( 1.0 * 1.602e-19 == Approx( pair51.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair51.particle().spin() ) );
      CHECK( +1 == pair51.particle().parity() );
      CHECK( 34.66863 * 1.008664 == Approx( pair51.residual().mass().value ) );
      CHECK( 16.0 * 1.602e-19 == Approx( pair51.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair51.residual().spin() ) );
      CHECK( +1 == pair51.residual().parity() );
      CHECK( "p,S35" == pair51.pairID().symbol() );

      // quantum numbers
      const auto numbers51 = channel51.quantumNumbers();
      CHECK( 1 == numbers51.orbitalAngularMomentum() );
      CHECK( 2. == numbers51.spin() );
      CHECK( 3. == numbers51.totalAngularMomentum() );
      CHECK( -1 == numbers51.parity() );
      CHECK( "{1,2,3-}" == numbers51.toString() );

      // radii
      const auto radii51 = channel51.radii();
      CHECK( 0.4822220 == Approx( radii51.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4822220 == Approx( radii51.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 0.4888750 == Approx( radii51.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel51.boundaryCondition() );

      // Q value
      CHECK( 6.152200e+5 == Approx( channel51.Q().value ) );

      // resonance table
      auto table5 = spingroup5.resonanceTable();

      CHECK( 2 == table5.numberChannels() ); // 2 normal channel + 1 eliminated
      CHECK( 57 == table5.numberResonances() );

      auto energies5 = table5.energies();
      CHECK( 1.635612e+4 == Approx( energies5.front().value ) );
      CHECK( 1.485128e+6 == Approx( energies5.back().value ) );

      auto resonances5 = table5.resonances();
      CHECK( 1.635612e+4 == Approx( resonances5.front().energy().value ) );
      CHECK( 1.485128e+6 == Approx( resonances5.back().energy().value ) );
      CHECK( 2 == resonances5.front().widths().size() );
      CHECK( 2 == resonances5.back().widths().size() );
      CHECK( std::sqrt( 5.981800 / 2. / channel50.penetrability( 1.635612e+4 * electronVolt ) )
             == Approx( resonances5.front().widths()[0].value ) );
      CHECK( std::sqrt( 0.164019 / 2. / channel51.penetrability( 1.635612e+4 * electronVolt ) )
             == Approx( resonances5.front().widths()[1].value ) );
      CHECK( std::sqrt( 1.054090e+4 / 2. / channel50.penetrability( 1.485128e+6 * electronVolt ) )
             == Approx( resonances5.back().widths()[0].value ) );
      CHECK( 0. == Approx( resonances5.back().widths()[1].value ) );
      CHECK( std::sqrt( 0.3865 / 2. ) == Approx( resonances5.front().eliminatedWidth().value ) );
      CHECK( std::sqrt( 0.86 / 2. ) == Approx( resonances5.back().eliminatedWidth().value ) );
    } // THEN

    THEN( "cross sections can be reconstructed" ) {

      ReactionID elas( "n,Cl35->n,Cl35" );
      ReactionID prot( "n,Cl35->h1,S35" );
      ReactionID capt( "n,Cl35->capture" );
      std::map< ReactionID, CrossSection > xs;

      xs = resonances( 1e-5 * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 20.6887411443 == Approx( xs[ elas ].value ) );
      CHECK( 0. == Approx( xs[ prot ].value ) );
      CHECK( 2185.6323459543 == Approx( xs[ capt ].value ) );
    } // THEN
  } // GIVEN
} // SCENARIO

std::string Fe54() {

  // Fe54 ENDF/B-VIII.0 LRF=7 resonance evaluation

  return
    " 2.605400+4 5.347624+1          0          0          1          02625 2151     \n"
    " 2.605400+4 1.000000+0          0          0          1          02625 2151     \n"
    " 1.000000-5 1.036000+6          1          7          0          12625 2151     \n"
    " 0.000000+0 0.000000+0          0          3          5          02625 2151     \n"
    " 0.000000+0 0.000000+0          2          0         24          42625 2151     \n"
    " 0.000000+0 5.446635+1 0.000000+0 2.600000+1 1.000000+0 0.000000+02625 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 1.020000+2 0.000000+0 0.000000+02625 2151     \n"
    " 1.000000+0 5.347624+1 0.000000+0 2.600000+1 5.000000-1 0.000000+02625 2151     \n"
    " 0.000000+0 1.000000+0 0.000000+0 2.000000+0 0.000000+0 1.000000+02625 2151     \n"
    " 5.000000-1 0.000000+0          0          0         12          22625 2151     \n"
    " 1.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.000000+0 0.000000+0 5.000000-1 0.000000+0 5.437300-1 5.437300-12625 2151     \n"
    " 0.000000+0 0.000000+0          0        148        888        1482625 2151     \n"
    "-1.223300+6 1.000000+0 9.611086+5 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    "-2.500000+4 1.798017+0 1.032662+4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.788000+3 1.455000+0 1.187354+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.287200+4 2.000000+0 2.000345+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.190500+4 2.000000+0 1.781791+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.881300+4 2.000000+0 5.551680+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.301700+5 3.000000+0 3.109384+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.475400+5 2.600000+0 3.083417+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.743000+5 2.000000+0 4.114302+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.887000+5 2.000000+0 3.684975+4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.067800+5 4.000000-1 1.329600+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.236200+5 6.000000-1 8.984901+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.441500+5 2.000000+0 1.990200+4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.573900+5 4.000000-1 4.010800+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.922000+5 1.100000+0 9.317101+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.072400+5 4.000000-1 6.319500+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.201400+5 4.000000-1 2.220800+4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.316300+5 4.000000-1 3.519400+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.637100+5 4.000000-1 3.811500+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.712300+5 4.000000-1 4.752100+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.100200+5 4.000000-1 1.010700+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.112100+5 4.000000-1 3.435300+4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.254800+5 4.000000-1 1.463000+4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.302300+5 4.000000-1 3.839600+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.556500+5 4.000000-1 6.794200+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.832300+5 4.000000-1 2.227600+4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.852200+5 4.000000-1 4.562300+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.041979+5 4.000000-1 3.716947+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.106047+5 4.000000-1 2.039219+4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.268393+5 4.000000-1 9.906387+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.413763+5 4.000000-1 7.350561+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.587263+5 4.000000-1 6.587169+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.606915+5 4.000000-1 3.444498+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.928310+5 4.000000-1 2.908803+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.089289+5 4.000000-1 2.814805+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.420913+5 4.000000-1 1.435638+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.612270+5 4.000000-1 7.269890+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.726386+5 4.000000-1 2.202403+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.893352+5 4.000000-1 7.337199+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.056951+5 4.000000-1 1.145771-2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.197259+5 4.000000-1 5.241748+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.365430+5 4.000000-1 1.096020+4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.484317+5 4.000000-1 1.583502+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.625934+5 4.000000-1 2.142637+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.934283+5 4.000000-1 3.224649+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.155311+5 4.000000-1 2.624618+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.191293+5 4.000000-1 3.436016+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.338111+5 4.000000-1 1.988010+4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.440180+5 4.000000-1 7.598592+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.610880+5 4.000000-1 3.754241+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.712498+5 4.000000-1 4.100152+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.757820+5 4.000000-1 1.118407+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.770124+5 4.000000-1 1.126387-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.945130+5 4.000000-1 1.872254+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.005395+5 4.000000-1 9.809761-2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.045392+5 4.000000-1 4.966268+4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.208361+5 4.000000-1 1.102448+4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.270522+5 4.000000-1 5.401460+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.521109+5 4.000000-1 1.430473+4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.601766+5 4.000000-1 3.057673+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.746035+5 4.000000-1 1.702876+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.836143+5 4.000000-1 3.736347-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.925091+5 4.000000-1 2.168022+4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.989940+5 4.000000-1 3.905539+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.006835+6 7.000000-1 3.493911+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.007882+6 7.000000-1 1.857533+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.007962+6 7.000000-1 1.261542+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.012074+6 7.000000-1 4.408087+4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.015771+6 7.000000-1 1.211998+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.016541+6 7.000000-1 7.275909+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.016807+6 7.000000-1 8.587253+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.023773+6 4.000000-1 1.382577+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.029100+6 4.000000-1 1.760486+4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.030373+6 7.000000-1 1.878123+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.036413+6 4.000000-1 1.705277+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.036500+6 4.000000-1 1.196904+5 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.044895+6 7.000000-1 2.387808+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.048513+6 7.000000-1 2.865925+4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.050702+6 7.000000-1 4.321734+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.053288+6 7.000000-1 1.668717+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.072431+6 7.000000-1 1.037008+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.073040+6 7.000000-1 4.756084+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.074304+6 7.000000-1 2.581614+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.076827+6 7.000000-1 2.561767+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.092165+6 7.000000-1 5.178618+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.094301+6 7.000000-1 4.002575+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.100703+6 7.000000-1 8.000866+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.101455+6 7.000000-1 1.629439+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.102090+6 7.000000-1 7.044857+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.105549+6 7.000000-1 2.826743+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.106805+6 7.000000-1 2.914327+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.113871+6 7.000000-1 6.341946+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.115013+6 7.000000-1 2.983193+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.116032+6 7.000000-1 2.011558+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.119116+6 7.000000-1 1.832122+4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.125409+6 7.000000-1 1.602570+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.130170+6 7.000000-1 1.063606+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.137036+6 7.000000-1 1.177170+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.139917+6 7.000000-1 7.116193+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.148549+6 7.000000-1 6.479679+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.149079+6 7.000000-1 1.791819+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.153116+6 7.000000-1 1.291133+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.160226+6 7.000000-1 2.301283+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.166220+6 7.000000-1 1.626226+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.176926+6 7.000000-1 2.721416+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.177563+6 7.000000-1 2.833032+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.180217+6 7.000000-1 2.960697+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.184935+6 7.000000-1 1.349241+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.193010+6 7.000000-1 1.346723+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.211500+6 7.000000-1 8.129232+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.211900+6 7.000000-1 2.232713+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.216800+6 7.000000-1 7.532561+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.217600+6 7.000000-1 1.499205+4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.219900+6 7.000000-1 1.368553+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.222600+6 7.000000-1 2.580568+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.228000+6 7.000000-1 1.045278+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.231200+6 7.000000-1 8.469479+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.235500+6 7.000000-1 1.296387+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.238200+6 7.000000-1 3.595342+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.261800+6 7.000000-1 5.152908+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.265000+6 7.000000-1 2.320658+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.270700+6 7.000000-1 5.740736+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.272800+6 7.000000-1 1.121145+4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.277700+6 7.000000-1 1.963529+4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.286400+6 7.000000-1 5.449450+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.288900+6 7.000000-1 4.610701+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.289300+6 7.000000-1 3.541559+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.294000+6 7.000000-1 1.513853+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.301464+6 7.000000-1 1.829046+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.303552+6 7.000000-1 6.868546+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.352640+6 7.000000-1 2.170071+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.382017+6 7.000000-1 8.166709-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.401250+6 7.000000-1 2.397060+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.407300+6 7.000000-1 1.597179+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.411100+6 7.000000-1 4.537008+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.412100+6 7.000000-1 4.967713+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.414000+6 7.000000-1 6.834148+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.434900+6 7.000000-1 1.972391+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.436800+6 7.000000-1 1.253945+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.459200+6 7.000000-1 4.672500+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.461300+6 7.000000-1 2.205592+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.472900+6 7.000000-1 4.017182+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.473500+6 7.000000-1 4.226838+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.475000+6 7.000000-1 5.630068+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.478600+6 7.000000-1 1.844172+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.501105+6 6.490600-1 1.598740+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.503308+6 6.490600-1 1.599913+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.505510+6 6.490600-1 1.601084+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    "-5.000000-1 0.000000+0          0          0         12          22625 2151     \n"
    " 1.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.000000+0 1.000000+0 5.000000-1 0.000000+0 5.437300-1 5.437300-12625 2151     \n"
    " 0.000000+0 0.000000+0          0        176       1056        1762625 2151     \n"
    " 5.152000+4 3.600000-1 1.600200+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.359000+4 1.500000+0 1.700000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.545900+4 5.600000-1 3.200000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.432000+4 4.000000-1 2.601400+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.348000+4 5.900000-1 8.600000-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.135500+4 4.000000-1 6.601200+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.041600+5 4.000000-1 3.906100+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.158500+5 4.000000-1 2.152200+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.198300+5 4.000000-1 2.398400+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.263300+5 4.000000-1 6.520300+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.357600+5 3.500000-1 8.148600+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.504200+5 4.000000-1 7.445300+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.527000+5 4.000000-1 1.925500+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.571200+5 4.000000-1 1.077800+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.592800+5 4.000000-1 1.176400+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.652900+5 4.000000-1 6.101800+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.740300+5 4.000000-1 1.974800+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.779100+5 4.000000-1 7.004800+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.075000+5 4.000000-1 6.698100+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.076100+5 4.000000-1 1.182800+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.138500+5 4.000000-1 1.518300+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.253300+5 4.000000-1 4.758500+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.312300+5 5.900000-1 2.216667+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.377500+5 4.000000-1 7.063000+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.420400+5 4.000000-1 1.341400+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.623500+5 4.000000-1 4.649100+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.636000+5 4.000000-1 8.453800+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.714900+5 4.000000-1 7.560900+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.809900+5 3.200000-1 1.976800+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.893500+5 4.000000-1 6.920200+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.928600+5 4.000000-1 2.074600+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.035800+5 4.000000-1 1.919700+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.038200+5 4.000000-1 1.009800+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.091600+5 4.000000-1 2.961000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.117800+5 4.000000-1 4.457000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.260900+5 4.000000-1 9.192000+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.327500+5 4.000000-1 3.685100+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.367500+5 4.000000-1 5.497300+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.438500+5 4.000000-1 6.212900+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.573800+5 4.000000-1 2.848100+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.610300+5 5.900000-1 5.730000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.627300+5 5.900000-1 6.000000+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.673800+5 4.000000-1 2.910200+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.685200+5 4.000000-1 1.159600+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.701800+5 4.000000-1 1.031900+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.751300+5 4.000000-1 2.678900+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.847700+5 5.900000-1 7.800000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.896800+5 5.900000-1 1.240000+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.978500+5 5.900000-1 5.080000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.063100+5 4.000000-1 4.208800+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.163000+5 4.000000-1 7.167600+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.276600+5 4.000000-1 2.120400+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.277000+5 4.000000-1 3.344400+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.435300+5 4.000000-1 1.633100+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.464600+5 4.000000-1 8.464100+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.615900+5 4.000000-1 6.324700+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.686700+5 4.000000-1 5.382800+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.740000+5 4.000000-1 6.103500+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.906000+5 4.000000-1 1.008000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.930000+5 4.000000-1 1.442400+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.988200+5 4.000000-1 2.717300+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.109800+5 4.000000-1 1.087134+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.118797+5 4.000000-1 1.079896+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.242641+5 4.000000-1 4.204802+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.244715+5 4.000000-1 1.942133+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.260984+5 4.000000-1 1.099708+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.392393+5 4.000000-1 1.999082+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.490489+5 4.000000-1 6.418378+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.510297+5 4.000000-1 1.267160+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.612238+5 4.000000-1 9.983347+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.893758+5 4.000000-1 3.267754+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.064872+5 4.000000-1 3.363521+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.304033+5 4.000000-1 1.303078+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.361519+5 4.000000-1 8.801258+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.411656+5 4.000000-1 3.879805+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.569479+5 4.000000-1 4.685032+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.654921+5 4.000000-1 3.842323+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.848858+5 4.000000-1 1.630389-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.855506+5 4.000000-1 2.036212+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.079449+5 4.000000-1 3.262791+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.243508+5 4.000000-1 4.622010+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.390339+5 4.000000-1 5.374703+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.519040+5 4.000000-1 2.715538+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.913333+5 4.000000-1 9.720917+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.021403+5 4.000000-1 2.317511+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.052483+5 4.000000-1 2.727453+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.209221+5 4.000000-1 4.991260+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.230299+5 4.000000-1 1.540980-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.281600+5 4.000000-1 4.679765+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.388523+5 4.000000-1 6.492580+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.544406+5 4.000000-1 2.247298+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.742371+5 4.000000-1 1.183769+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.791545+5 4.000000-1 1.729120+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.851337+5 4.000000-1 2.502694+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.038761+5 4.000000-1 1.278049+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.148105+5 4.000000-1 5.179936+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.261085+5 4.000000-1 1.082181+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.458108+5 4.000000-1 1.290737+4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.501513+5 4.000000-1 2.264682+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.561678+5 4.000000-1 2.404672+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.571613+5 4.000000-1 1.372520+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.706518+5 4.000000-1 3.140170+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.827969+5 4.000000-1 3.958198+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.918115+5 4.000000-1 4.349787-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.924344+5 4.000000-1 8.988154+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.943551+5 4.000000-1 5.400306+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.007936+6 7.000000-1 4.757654+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.010554+6 7.000000-1 1.788109+4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.013662+6 7.000000-1 4.526361+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.014236+6 7.000000-1 1.639979+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.018384+6 7.000000-1 1.445604+4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.020035+6 7.000000-1 2.270929+4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.022489+6 7.000000-1 1.585767+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.024088+6 7.000000-1 2.741289+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.028337+6 7.000000-1 2.688911+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.028730+6 7.000000-1 1.169575+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.035029+6 7.000000-1 5.049021+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.041483+6 7.000000-1 6.932586+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.049194+6 7.000000-1 5.853234-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.051941+6 7.000000-1 2.732861-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.054551+6 7.000000-1 2.067607+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.058661+6 7.000000-1 7.427173+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.067529+6 7.000000-1 2.610806+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.079005+6 7.000000-1 4.984630+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.080066+6 7.000000-1 4.202750+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.087100+6 7.000000-1 2.263758+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.093357+6 7.000000-1 1.181581+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.095347+6 7.000000-1 3.341273+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.114263+6 7.000000-1 3.951291-3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.129215+6 7.000000-1 2.521005+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.136345+6 7.000000-1 3.742127+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.142586+6 7.000000-1 3.697277+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.142930+6 7.000000-1 9.441586+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.145446+6 7.000000-1 1.890100+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.153957+6 7.000000-1 6.507070+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.165404+6 7.000000-1 2.397530+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.169624+6 7.000000-1 2.521247+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.172316+6 7.000000-1 4.367571+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.178299+6 7.000000-1 2.281827+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.186150+6 7.000000-1 2.175313+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.196735+6 7.000000-1 1.551876+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.197633+6 7.000000-1 1.472436+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.200900+6 7.000000-1 1.622288+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.203200+6 7.000000-1 2.045252+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.203900+6 7.000000-1 2.976736+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.211100+6 7.000000-1 8.314176+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.214200+6 7.000000-1 7.927457+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.228800+6 7.000000-1 1.436400+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.235900+6 7.000000-1 3.275355+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.246100+6 7.000000-1 3.361403+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.250500+6 7.000000-1 2.404723+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.258500+6 7.000000-1 8.707712+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.269500+6 7.000000-1 3.799778+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.288500+6 7.000000-1 4.659112+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.291400+6 7.000000-1 1.021775+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.293500+6 7.000000-1 7.503388+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.294400+6 7.000000-1 1.665153+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.295600+6 7.000000-1 4.178903+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.299000+6 7.000000-1 4.747361+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.299900+6 7.000000-1 3.440912+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.321001+6 7.000000-1 7.322842+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.325927+6 7.000000-1 1.240116+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.341047+6 7.000000-1 7.350127+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.365848+6 7.000000-1 5.123859-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.376402+6 7.000000-1 1.255995+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.384840+6 7.000000-1 3.267945+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.389363+6 7.000000-1 2.240603+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.390975+6 7.000000-1 4.612233+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.394987+6 7.000000-1 5.526015+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.415000+6 7.000000-1 3.262509+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.433400+6 7.000000-1 8.731223-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.437300+6 7.000000-1 2.112829+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.440300+6 7.000000-1 1.097058+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.458700+6 7.000000-1 1.608010+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.462300+6 7.000000-1 2.553174+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.470900+6 7.000000-1 1.946542+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    "-1.500000+0 0.000000+0          0          0         12          22625 2151     \n"
    " 1.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.000000+0 1.000000+0 5.000000-1 0.000000+0 5.437300-1 5.437300-12625 2151     \n"
    " 0.000000+0 0.000000+0          0        196       1176        1962625 2151     \n"
    " 3.099000+3 5.900000-1 1.400000-3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.358100+4 5.900000-1 1.750000-2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.927800+4 5.900000-1 2.750000-2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.303000+4 5.900000-1 2.675000-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.822000+4 5.900000-1 9.250000-2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.526000+4 5.900000-1 1.400000-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.844000+4 5.900000-1 1.400000+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.120000+4 5.900000-1 1.200000-2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.014000+4 5.900000-1 2.900000-2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.920000+4 5.900000-1 2.700000-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.876000+4 5.900000-1 1.500000-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.580000+4 4.000000-1 9.166800-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.720000+4 4.000000-1 2.743400+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.128000+4 5.900000-1 1.900000-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.315500+4 4.000000-1 3.876100+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.733000+4 5.900000-1 2.900000-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.135500+4 5.900000-1 1.300000-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.777000+4 5.900000-1 6.100000-2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.984000+4 5.900000-1 9.600000-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.017100+5 5.900000-1 1.900000-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.125800+5 5.900000-1 7.000000-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.129600+5 5.900000-1 3.300000-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.158500+5 5.500000-1 1.350000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.197800+5 5.300000-1 1.300000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.265300+5 8.000000-1 2.112100+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.379300+5 5.900000-1 6.300000-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.410100+5 5.900000-1 5.900000-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.455100+5 5.900000-1 4.100000-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.531900+5 4.000000-1 2.270600+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.572000+5 5.600000-1 1.085000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.645100+5 1.400000+0 5.649700+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.890900+5 5.900000-1 2.500000-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.944800+5 4.000000-1 4.153200+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.039700+5 4.000000-1 1.609000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.076000+5 5.900000-1 1.000000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.100100+5 5.900000-1 3.870000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.257100+5 5.900000-1 5.700000+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.312200+5 4.000000-1 1.590600+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.337100+5 5.900000-1 1.295000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.420600+5 5.900000-1 8.800000+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.452900+5 4.000000-1 2.665600+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.548900+5 5.900000-1 3.775000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.630200+5 4.000000-1 1.642900+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.829400+5 4.000000-1 8.750400+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.091400+5 5.900000-1 1.310000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.220200+5 4.000000-1 1.284700+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.326900+5 5.900000-1 1.455000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.452700+5 4.000000-1 1.177500+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.570900+5 4.000000-1 1.296200+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.610900+5 4.000000-1 2.322500+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.848400+5 4.000000-1 2.358900+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.897200+5 4.000000-1 5.620500+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.899100+5 4.000000-1 1.699400+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.979100+5 4.300000-1 2.239900+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.988000+5 4.000000-1 2.994400+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.123900+5 4.000000-1 9.687500+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.275300+5 5.900000-1 1.115000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.277000+5 5.900000-1 1.410000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.469700+5 4.000000-1 6.387000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.589700+5 5.900000-1 1.390000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.614600+5 5.900000-1 3.205000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.654500+5 5.900000-1 1.810000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.703300+5 4.000000-1 1.404700+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.786900+5 4.000000-1 8.111300+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.901500+5 4.000000-1 4.212400+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.980600+5 4.000000-1 1.854200+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.113799+5 4.000000-1 2.661896+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.242200+5 5.900000-1 7.000000+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.305283+5 4.000000-1 2.377581+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.525716+5 4.000000-1 1.806954+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.598908+5 4.000000-1 5.505216+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.651993+5 4.000000-1 7.960016+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.687365+5 4.000000-1 9.021971+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.809309+5 4.000000-1 4.985751+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.952300+5 5.900000-1 5.050000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.042492+5 4.000000-1 8.322616+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.044843+5 4.000000-1 4.213639+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.212139+5 4.000000-1 7.294993+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.385500+5 5.900000-1 2.755000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.470683+5 4.000000-1 5.606430+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.567200+5 5.900000-1 2.000000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.699057+5 4.000000-1 2.493408+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.768422+5 4.000000-1 3.262786+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.800487+5 4.000000-1 4.561021+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.805754+5 4.000000-1 3.858567+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.818110+5 4.000000-1 5.694376+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.835464+5 4.000000-1 4.575830+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.840877+5 4.000000-1 2.825082-4 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.876284+5 4.000000-1 4.352441+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.894095+5 4.000000-1 4.969923+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.978326+5 4.000000-1 6.647285+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.087790+5 4.000000-1 5.916424+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.226067+5 4.000000-1 1.205360+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.265519+5 4.000000-1 4.449364+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.387100+5 5.900000-1 2.675000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.494766+5 4.000000-1 8.185866+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.735575+5 4.000000-1 8.873871+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.909786+5 4.000000-1 6.627769+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.954090+5 4.000000-1 9.319217+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.070882+5 4.000000-1 3.330302+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.129077+5 4.000000-1 8.293181+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.186005+5 4.000000-1 9.917532+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.254062+5 4.000000-1 2.333022+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.290519+5 4.000000-1 6.492107+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.357682+5 4.000000-1 4.500824+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.412780+5 4.000000-1 1.254818+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.434657+5 4.000000-1 6.169918+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.474795+5 4.000000-1 2.071445+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.522371+5 4.000000-1 2.188346+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.591291+5 4.000000-1 4.349779+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.763763+5 4.000000-1 2.369700+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.881116+5 4.000000-1 1.688519+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.960101+5 4.000000-1 1.219338+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.977190+5 4.000000-1 1.028863+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.099892+5 4.000000-1 1.294740+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.144508+5 4.000000-1 2.595135+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.164981+5 4.000000-1 1.931471+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.187953+5 4.000000-1 1.313353+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.227527+5 4.000000-1 1.658686+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.273890+5 4.000000-1 2.096364+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.448935+5 4.000000-1 7.524011+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.636213+5 4.000000-1 1.307320+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.680020+5 4.000000-1 3.850489+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.784231+5 4.000000-1 8.367731+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.838438+5 4.000000-1 4.147539+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.900145+5 4.000000-1 5.543467-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.002368+6 7.000000-1 1.518759+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.005913+6 7.000000-1 7.659975+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.010327+6 7.000000-1 7.439483+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.012920+6 7.000000-1 7.671746+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.022815+6 7.000000-1 2.289961-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.026431+6 7.000000-1 1.268982+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.028429+6 7.000000-1 1.547053+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.031001+6 7.000000-1 3.223212+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.031465+6 7.000000-1 3.598004+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.032329+6 7.000000-1 4.183987+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.051001+6 7.000000-1 6.110995+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.056050+6 7.000000-1 4.184111+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.058038+6 7.000000-1 6.158646+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.062925+6 7.000000-1 2.761361+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.072620+6 7.000000-1 9.024670+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.083488+6 7.000000-1 1.088390+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.085559+6 7.000000-1 1.483430-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.091677+6 7.000000-1 1.194235+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.095462+6 7.000000-1 1.894608+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.099070+6 7.000000-1 5.528152+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.101233+6 7.000000-1 1.596818+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.104403+6 7.000000-1 1.760413+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.106613+6 7.000000-1 1.989903+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.107490+6 7.000000-1 1.274670+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.109121+6 7.000000-1 3.617540+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.110144+6 7.000000-1 2.370016+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.117903+6 7.000000-1 8.918474+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.118632+6 7.000000-1 7.202237+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.119608+6 7.000000-1 8.576966+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.125225+6 7.000000-1 3.283345+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.128490+6 7.000000-1 1.075546+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.135245+6 7.000000-1 6.675872-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.141235+6 7.000000-1 6.903190+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.144014+6 7.000000-1 8.474489+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.160850+6 7.000000-1 1.739419+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.166727+6 7.000000-1 3.483328+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.172094+6 7.000000-1 5.614840-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.179700+6 7.000000-1 1.049276+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.190647+6 7.000000-1 1.988563+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.191495+6 7.000000-1 9.753071-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.195039+6 7.000000-1 1.337952+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.197252+6 7.000000-1 2.097729+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.200500+6 7.000000-1 1.226610+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.214900+6 7.000000-1 6.750564+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.220700+6 7.000000-1 6.827222-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.274800+6 7.000000-1 5.396743+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.288100+6 7.000000-1 2.739994+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.297300+6 7.000000-1 1.083002+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.300387+6 7.000000-1 5.493265+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.306103+6 7.000000-1 1.778696+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.307309+6 7.000000-1 1.608504-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.312399+6 7.000000-1 6.039068+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.318355+6 7.000000-1 1.022869+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.324614+6 7.000000-1 7.662769+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.362008+6 7.000000-1 9.721002+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.370898+6 7.000000-1 1.373088+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.382499+6 7.000000-1 3.325686+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.396057+6 7.000000-1 1.046229+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.403500+6 7.000000-1 3.838244+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.407800+6 7.000000-1 8.793682+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.421700+6 7.000000-1 9.757445+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.430000+6 7.000000-1 9.268406-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.432400+6 7.000000-1 7.682903-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.439300+6 7.000000-1 1.981232+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.441300+6 7.000000-1 5.373680-3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.465800+6 7.000000-1 4.801503+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.467900+6 7.000000-1 8.062626+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.468900+6 7.000000-1 2.746135+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.475500+6 7.000000-1 2.636608+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.496700+6 7.000000-1 5.591693+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.500000+0 0.000000+0          0          0         12          22625 2151     \n"
    " 1.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.000000+0 2.000000+0 5.000000-1 0.000000+0 5.437300-1 5.437300-12625 2151     \n"
    " 0.000000+0 0.000000+0          0        174       1044        1742625 2151     \n"
    " 9.480000+3 2.700000-1 1.200000+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.118000+4 3.500000-1 3.850100+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.445000+4 3.500000-1 7.000200-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.064000+4 4.800000-1 5.000100+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.914000+4 4.100000-1 8.500200-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.354000+4 3.000000-1 8.502400+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.535000+4 3.400000-1 1.600400+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.206900+5 4.000000-1 9.119800-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.207700+5 4.500000-1 1.613900+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.428200+5 9.600000-1 1.200000+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.646500+5 4.000000-1 2.848500+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.824200+5 5.010000-1 4.015000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.921200+5 9.600000-1 1.200000+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.977000+5 9.600000-1 7.600000-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.099900+5 4.000000-1 2.473800+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.257700+5 4.000000-1 1.283500+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.258000+5 4.000000-1 1.312500+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.280900+5 4.000000-1 3.677300+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.293000+5 4.000000-1 5.195600+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.336800+5 4.000000-1 1.611700+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.549000+5 4.000000-1 2.828300+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.629000+5 4.000000-1 1.811800+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.642500+5 4.000000-1 8.732600+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.773200+5 4.000000-1 1.202800+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.059800+5 4.000000-1 1.287800+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.241600+5 4.000000-1 6.035600+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.677300+5 4.000000-1 2.369800+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.830000+5 4.000000-1 2.006800+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.907600+5 4.000000-1 1.392500+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.911900+5 4.000000-1 3.733800+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 3.989500+5 4.000000-1 2.217800+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.029600+5 4.000000-1 1.621900+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.095000+5 4.000000-1 3.039900+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.136700+5 4.000000-1 1.566300+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.167000+5 4.000000-1 4.239100+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.279900+5 4.000000-1 4.708300+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.321100+5 4.000000-1 1.917200+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.428100+5 4.000000-1 1.435500+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.460300+5 4.000000-1 2.966900+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.465700+5 4.000000-1 5.843200+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.590400+5 4.000000-1 1.394400+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.655400+5 4.000000-1 2.570000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.760200+5 4.000000-1 8.034900+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.762900+5 4.000000-1 2.767100+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.820800+5 4.000000-1 1.243900+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 4.857700+5 4.000000-1 5.517600+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.073475+5 4.000000-1 1.379720+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.102716+5 4.000000-1 6.736567+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.104292+5 4.000000-1 1.103354+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.154764+5 4.000000-1 2.801142+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.190815+5 4.000000-1 1.666742+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.216600+5 4.000000-1 4.820773+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.341374+5 4.000000-1 5.966375+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.385437+5 4.000000-1 1.679649+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.496540+5 4.000000-1 1.547438+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.561857+5 4.000000-1 4.058745+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.694581+5 4.000000-1 1.073030+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.716764+5 4.000000-1 2.133730+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.761440+5 4.000000-1 2.209928+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.768834+5 4.000000-1 3.343362+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.791168+5 4.000000-1 5.318386+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.852562+5 4.000000-1 3.838181+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.856907+5 4.000000-1 4.050193+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.882995+5 4.000000-1 4.559138+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.953607+5 4.000000-1 4.823196+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.972367+5 4.000000-1 5.953381+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 5.990237+5 4.000000-1 3.782001+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.167820+5 4.000000-1 2.131430+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.186659+5 4.000000-1 4.343671+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.348638+5 4.000000-1 1.106123+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.389130+5 4.000000-1 1.426816+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.428880+5 4.000000-1 5.749459+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.435311+5 4.000000-1 1.272912-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.485029+5 4.000000-1 1.072359+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.529115+5 4.000000-1 5.154546+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.547732+5 4.000000-1 3.807383+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.631077+5 4.000000-1 2.946856+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.648455+5 4.000000-1 7.303448+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.749092+5 4.000000-1 5.093406+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.986631+5 4.000000-1 1.384974+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.048440+5 4.000000-1 1.236817+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.172003+5 4.000000-1 1.853490+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.310510+5 4.000000-1 9.744808+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.415731+5 4.000000-1 5.859396+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.446768+5 4.000000-1 1.510914+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.537304+5 4.000000-1 3.635781+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.606821+5 4.000000-1 7.738128+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.695500+5 9.600000-1 2.750000+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.755800+5 9.600000-1 2.330000+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.839034+5 4.000000-1 6.394942+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.883167+5 4.000000-1 9.858968-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.070300+5 9.600000-1 7.950000+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.199535+5 4.000000-1 5.422480+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.622708+5 4.000000-1 2.318639+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.723056+5 4.000000-1 7.052260-2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.238092+5 4.000000-1 2.697234-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.329524+5 4.000000-1 1.072767+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.365564+5 4.000000-1 1.358183+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.616729+5 4.000000-1 1.123392+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.911457+5 4.000000-1 9.450324+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.020395+6 7.000000-1 1.286129+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.025416+6 7.000000-1 2.147429+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.032867+6 7.000000-1 8.407657-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.034654+6 7.000000-1 1.485668+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.036976+6 7.000000-1 2.008344+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.037442+6 7.000000-1 1.620103+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.044659+6 7.000000-1 2.324646-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.056393+6 7.000000-1 1.582706+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.063226+6 7.000000-1 1.769481+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.068528+6 7.000000-1 2.381345+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.069529+6 7.000000-1 2.270916+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.071085+6 7.000000-1 5.606077+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.083799+6 7.000000-1 6.748964+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.085802+6 7.000000-1 1.336727+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.087257+6 7.000000-1 4.625587+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.089521+6 7.000000-1 4.606651-2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.090365+6 7.000000-1 2.133190+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.097836+6 7.000000-1 4.604364-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.107098+6 7.000000-1 1.165450+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.110852+6 7.000000-1 3.184407-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.138229+6 7.000000-1 1.098227+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.138864+6 7.000000-1 2.995481+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.140031+6 7.000000-1 1.256156+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.144947+6 7.000000-1 5.943933-3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.160565+6 7.000000-1 1.941732+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.164436+6 7.000000-1 4.292274+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.172409+6 7.000000-1 2.776853+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.182174+6 7.000000-1 2.143896+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.183615+6 7.000000-1 2.970873+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.188969+6 7.000000-1 2.718275+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.193005+6 7.000000-1 1.868313+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.200100+6 7.000000-1 1.243752+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.206900+6 7.000000-1 1.720366+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.218000+6 7.000000-1 1.068909+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.226500+6 7.000000-1 1.388803+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.226900+6 7.000000-1 6.263633+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.231900+6 7.000000-1 4.681152-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.234700+6 7.000000-1 2.314244+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.242600+6 7.000000-1 2.928929+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.243000+6 7.000000-1 4.564524+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.245700+6 7.000000-1 3.598557-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.246500+6 7.000000-1 6.446429+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.247300+6 7.000000-1 1.382621+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.249300+6 7.000000-1 8.079770+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.250100+6 7.000000-1 8.699172+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.252900+6 7.000000-1 1.340302+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.255300+6 7.000000-1 3.059559-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.255700+6 7.000000-1 3.350943+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.258100+6 7.000000-1 6.970356+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.262200+6 7.000000-1 3.769258+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.264600+6 7.000000-1 1.110896+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.271100+6 7.000000-1 3.214584+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.272400+6 7.000000-1 2.032512+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.280200+6 7.000000-1 6.998520+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.285600+6 7.000000-1 8.680330+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.291900+6 7.000000-1 2.057113+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.296100+6 7.000000-1 3.140329+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.310522+6 7.000000-1 1.626881+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.342597+6 7.000000-1 2.866731+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.353687+6 7.000000-1 6.316532+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.373551+6 7.000000-1 1.034589-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.386529+6 7.000000-1 2.898360+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.392287+6 7.000000-1 2.060996+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.400200+6 7.000000-1 1.115285+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.404400+6 7.000000-1 1.298045+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.419800+6 7.000000-1 1.121747+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.432900+6 7.000000-1 5.891910-2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.435400+6 7.000000-1 5.562327-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.438800+6 7.000000-1 2.492189+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.455700+6 7.000000-1 2.276831+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.459800+6 7.000000-1 2.752649+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.465300+6 7.000000-1 1.902441+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.469400+6 7.000000-1 2.018815+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.480100+6 7.000000-1 1.766491+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.500000+0 0.000000+0          0          0         12          22625 2151     \n"
    " 1.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.000000+0 2.000000+0 5.000000-1 0.000000+0 5.437300-1 5.437300-12625 2151     \n"
    " 0.000000+0 0.000000+0          0         43        258         432625 2151     \n"
    " 1.264000+5 1.100000+0 2.900000+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.504700+5 9.600000-1 2.600000+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.779400+5 9.600000-1 1.400000+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.823600+5 4.000000-1 5.609500+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.945200+5 9.600000-1 5.300000-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 2.635000+5 4.000000-1 3.584800+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 6.241610+5 4.000000-1 1.280927+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.151255+5 4.000000-1 3.506706+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.196564+5 4.000000-1 4.798990+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.282568+5 4.000000-1 6.731346+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.406592+5 4.000000-1 1.917303+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.506336+5 4.000000-1 3.794667+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.555251+5 4.000000-1 2.078499+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.643228+5 4.000000-1 1.046145+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.671751+5 4.000000-1 1.961934+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.698952+5 4.000000-1 1.742457+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.759082+5 4.000000-1 1.815106+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.823973+5 4.000000-1 3.856682+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 7.965074+5 4.000000-1 6.345531+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.122098+5 4.000000-1 2.518495-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.307875+5 4.000000-1 7.280372+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.465205+5 4.000000-1 1.165599+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.564689+5 4.000000-1 1.551788+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.630488+5 4.000000-1 1.511368+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.699976+5 4.000000-1 1.737354+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.757905+5 4.000000-1 9.092705+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.814201+5 4.000000-1 1.593016+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 8.864225+5 4.000000-1 3.976289+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.012799+5 4.000000-1 1.208135+3 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.123647+5 4.000000-1 3.051593+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.382464+5 4.000000-1 5.450844+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.393989+5 4.000000-1 8.281391+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.413543+5 4.000000-1 1.625061+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.422685+5 4.000000-1 1.445374+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.430822+5 4.000000-1 1.347471+0 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.527008+5 4.000000-1 4.198989-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.555271+5 4.000000-1 1.052502+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.565799+5 4.000000-1 6.338464+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.584438+5 4.000000-1 1.034887+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.769656+5 4.000000-1 2.873313+2 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.794482+5 4.000000-1 8.738206+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 9.805248+5 4.000000-1 9.788140-1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 1.001152+6 4.000000-1 6.940802+1 0.000000+0 0.000000+0 0.000000+02625 2151     \n"
    " 0.000000+0 0.000000+0          0          0          0          02625 2  0     \n";
}

std::string Ca40() {

  // Ca40 ENDF/B-VIII.0 LRF=7 resonance evaluation

  return
    " 2.004000+4 3.961930+1          0          0          1          02025 2151     \n"
    " 2.004000+4 1.000000+0          0          0          1          02025 2151     \n"
    " 1.000000-5 1.500000+6          1          7          0          12025 2151     \n"
    " 0.000000+0 0.000000+0          0          3          5          02025 2151     \n"
    " 0.000000+0 0.000000+0          4          0         48          82025 2151     \n"
    " 0.000000+0 4.061929+1 0.000000+0 2.000000+1 1.000000+0 0.000000+02025 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 1.020000+2 0.000000+0 0.000000+02025 2151     \n"
    " 1.000000+0 3.961929+1 0.000000+0 2.000000+1 5.000000-1 0.000000+02025 2151     \n"
    " 0.000000+0 1.000000+0 0.000000+0 2.000000+0 0.000000+0 1.000000+02025 2151     \n"
    " 3.967131+0 3.664921+1 2.000000+0 1.800000+1 0.000000+0 1.500000+02025 2151     \n"
    " 1.747660+6 1.000000+0 0.000000+0 8.000000+2 1.000000+0 0.000000+02025 2151     \n"
    " 9.986235-1 3.962069+1 1.000000+0 1.900000+1 5.000000-1-4.000000+02025 2151     \n"
    "-5.285469+5 1.000000+0 0.000000+0 6.000000+2 0.000000+0 0.000000+02025 2151     \n"
    " 5.000000-1 0.000000+0          0          0         18          32025 2151     \n"
    " 1.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+02025 2151     \n"
    " 2.000000+0 0.000000+0 5.000000-1 0.000000+0 3.828792-1 4.993153-12025 2151     \n"
    " 3.000000+0 2.000000+0 1.500000+0 0.000000+0 4.885639-1 4.885639-12025 2151     \n"
    " 0.000000+0 0.000000+0          0         39        234         392025 2151     \n"
    "-4.586687+5 1.000091+0 9.809761+2 2.199782-3 0.000000+0 0.000000+02025 2151     \n"
    "-2.394214+5 1.000309+0 7.041383+2 1.219849-3 0.000000+0 0.000000+02025 2151     \n"
    "-1.929859+5 1.000388+0 5.076764+2 1.163776-3 0.000000+0 0.000000+02025 2151     \n"
    "-1.613693+5 1.000550+0 4.653400+2 1.150728-3 0.000000+0 0.000000+02025 2151     \n"
    "-3.033254+4 2.245017+0 6.554955+3 3.904361-3 0.000000+0 0.000000+02025 2151     \n"
    " 8.889281+4 8.150475-1 1.592992+2 9.999180-4 0.000000+0 0.000000+02025 2151     \n"
    " 1.326770+5 2.638134+0 2.510286+3 1.001433-3 0.000000+0 0.000000+02025 2151     \n"
    " 1.694194+5 2.435083+0 1.965314+3 1.001717-3 0.000000+0 0.000000+02025 2151     \n"
    " 2.176964+5 1.527959+0 4.751993+3 1.002430-3 0.000000+0 0.000000+02025 2151     \n"
    " 2.459636+5 1.520527+0 1.349507+4 1.005421-3 0.000000+0 0.000000+02025 2151     \n"
    " 2.918707+5 1.501569+0 1.368518+3 1.001707-3 0.000000+0 0.000000+02025 2151     \n"
    " 3.302023+5 1.007862+0 9.773980+3 1.005546-3 0.000000+0 0.000000+02025 2151     \n"
    " 3.557678+5 1.020744+0 1.549639+3 1.003722-3 0.000000+0 0.000000+02025 2151     \n"
    " 3.854519+5 1.016840+0 1.542018+2 1.003254-3 0.000000+0 0.000000+02025 2151     \n"
    " 4.394042+5 1.621389+0 6.101693+3 1.014057-3 0.000000+0 0.000000+02025 2151     \n"
    " 5.063640+5 1.624534+0 5.903202+3 1.011205-3 0.000000+0 0.000000+02025 2151     \n"
    " 5.915618+5 1.002107+0 3.368277+4 1.023269-3 0.000000+0 0.000000+02025 2151     \n"
    " 6.353157+5 1.006814+0 1.565760+3 1.008622-3 0.000000+0 0.000000+02025 2151     \n"
    " 6.746035+5 1.002731+0 1.969717+3 1.007068-3 0.000000+0 0.000000+02025 2151     \n"
    " 7.379045+5 1.566334+0 2.138456+3 9.881453-4 0.000000+0 0.000000+02025 2151     \n"
    " 7.430460+5 3.810711-1 2.869368+3 9.866825-4 0.000000+0 0.000000+02025 2151     \n"
    " 7.711377+5 9.999194-1 1.014552+4 1.012918-3 0.000000+0 0.000000+02025 2151     \n"
    " 7.924438+5 1.002873+0 1.842673+3 1.009080-3 0.000000+0 0.000000+02025 2151     \n"
    " 8.225723+5 1.000811+0 8.013577+1 1.015026-3 0.000000+0 0.000000+02025 2151     \n"
    " 8.620862+5 1.189165+0 2.037008+4 1.028543-3 0.000000+0 0.000000+02025 2151     \n"
    " 8.802359+5 1.000346+0 2.734472+4 1.023587-3 0.000000+0 0.000000+02025 2151     \n"
    " 9.700106+5 1.003842+0 5.871393+3 1.000075-3 0.000000+0 0.000000+02025 2151     \n"
    " 1.004426+6 7.164432-1 9.493242+3 9.811743-4 0.000000+0 0.000000+02025 2151     \n"
    " 1.093026+6 1.040785+0 1.817525+4 9.870769-4 0.000000+0 0.000000+02025 2151     \n"
    " 1.189083+6 1.051919+0 2.836113+3 1.000522-3 0.000000+0 0.000000+02025 2151     \n"
    " 1.283403+6 1.005859+0 2.510394+3 9.964546-4 0.000000+0 0.000000+02025 2151     \n"
    " 1.315693+6 1.000000+0 3.659372+3 9.866041-4 0.000000+0 0.000000+02025 2151     \n"
    " 1.373039+6 1.000000+0 3.518292+3 9.889280-4 0.000000+0 0.000000+02025 2151     \n"
    " 1.656895+6 1.000000+0 2.562519+5 1.319899-4 0.000000+0 0.000000+02025 2151     \n"
    " 1.767830+6 1.000000+0 1.601594+4 1.000914-3 0.000000+0 0.000000+02025 2151     \n"
    " 1.839939+6 1.000000+0 2.834579+4 1.012114-3 0.000000+0 0.000000+02025 2151     \n"
    " 1.850829+6 1.000000+0 7.005697+4 1.109617-3 0.000000+0 0.000000+02025 2151     \n"
    " 1.857524+6 1.000000+0 3.863731+5 1.440615-3 0.000000+0 0.000000+02025 2151     \n"
    " 1.913996+6 1.000000+0 1.878975+5 2.784131-4 0.000000+0 0.000000+02025 2151     \n"
    "-5.000000-1 0.000000+0          0          0         18          32025 2151     \n"
    " 1.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+02025 2151     \n"
    " 2.000000+0 1.000000+0 5.000000-1 0.000000+0 5.574746-1 4.993153-12025 2151     \n"
    " 3.000000+0 1.000000+0 1.500000+0 0.000000+0 4.885639-1 4.885639-12025 2151     \n"
    " 0.000000+0 0.000000+0          0         51        306         512025 2151     \n"
    " 2.880782+3 3.499989-1 1.810972-3 1.000029-3 0.000000+0 0.000000+02025 2151     \n"
    " 6.957171+3 3.500785-1 1.122320-2 1.000184-3 0.000000+0 0.000000+02025 2151     \n"
    " 1.003819+5 1.382810+0 1.388124+0 9.984372-4 0.000000+0 0.000000+02025 2151     \n"
    " 1.012155+5 6.278765-1 2.915016+1 9.973021-4 0.000000+0 0.000000+02025 2151     \n"
    " 1.295595+5 6.086651-1 4.388002-1 9.982139-4 0.000000+0 0.000000+02025 2151     \n"
    " 1.454277+5 6.016671-1 1.883684+2 1.000038-3 0.000000+0 0.000000+02025 2151     \n"
    " 1.580259+5 1.285404+0 3.974389-1 9.982235-4 0.000000+0 0.000000+02025 2151     \n"
    " 1.608229+5 1.197518+0 4.609762+1 1.002500-3 0.000000+0 0.000000+02025 2151     \n"
    " 2.114910+5 9.639369-1 3.426937+2 1.000291-3 0.000000+0 0.000000+02025 2151     \n"
    " 2.269925+5 1.915996+0 3.960588-1 1.009046-3 0.000000+0 0.000000+02025 2151     \n"
    " 2.295723+5 4.876817-1 3.569023+2 1.000273-3 0.000000+0 0.000000+02025 2151     \n"
    " 2.810462+5 1.084696+0 6.031804+0 9.988448-4 0.000000+0 0.000000+02025 2151     \n"
    " 2.938563+5 5.369699-1 1.612936+2 9.975840-4 0.000000+0 0.000000+02025 2151     \n"
    " 3.032079+5 6.505238-1 5.833486+1 9.999339-4 0.000000+0 0.000000+02025 2151     \n"
    " 3.059856+5 3.456795-1 6.955574+0 9.995323-4 0.000000+0 0.000000+02025 2151     \n"
    " 3.152028+5 3.444033-1 2.412008+0 9.951710-4 0.000000+0 0.000000+02025 2151     \n"
    " 3.322884+5 6.917381-1 1.990955+0 9.998567-4 0.000000+0 0.000000+02025 2151     \n"
    " 3.641640+5 1.000000+0 3.824830+0 1.000524-3 0.000000+0 0.000000+02025 2151     \n"
    " 4.194401+5 1.004392+0 1.049871+2 1.000059-3 0.000000+0 0.000000+02025 2151     \n"
    " 4.253701+5 1.000510+0 6.785937+0 9.972687-4 0.000000+0 0.000000+02025 2151     \n"
    " 4.371793+5 7.291594-1 5.576908+2 1.000090-3 0.000000+0 0.000000+02025 2151     \n"
    " 4.693849+5 3.813829-1 2.295565+2 1.000592-3 0.000000+0 0.000000+02025 2151     \n"
    " 4.931970+5 1.000000+0 5.898443+1 1.000429-3 0.000000+0 0.000000+02025 2151     \n"
    " 5.372131+5 5.111008-1 2.104851+2 1.002137-3 0.000000+0 0.000000+02025 2151     \n"
    " 6.026894+5 1.000000+0 1.229004+1 1.000982-3 0.000000+0 0.000000+02025 2151     \n"
    " 6.198572+5 1.000000+0 4.588795+1 1.000164-3 0.000000+0 0.000000+02025 2151     \n"
    " 6.307024+5 1.000000+0 4.875367+1 9.986493-4 0.000000+0 0.000000+02025 2151     \n"
    " 6.517651+5 3.495588-1 2.457345+1 9.942275-4 0.000000+0 0.000000+02025 2151     \n"
    " 6.585853+5 1.000000+0 6.445188+0 9.992971-4 0.000000+0 0.000000+02025 2151     \n"
    " 6.678979+5 3.495588-1 1.259035+2 1.001379-3 0.000000+0 0.000000+02025 2151     \n"
    " 6.942533+5 3.502237-1 8.926642+2 1.000665-3 0.000000+0 0.000000+02025 2151     \n"
    " 7.127810+5 3.800376-1 2.214914+2 9.035939-4 0.000000+0 0.000000+02025 2151     \n"
    " 8.234357+5 1.000496+0 2.024149+3 9.944059-4 0.000000+0 0.000000+02025 2151     \n"
    " 8.481097+5 9.006851-1 6.526597+1 1.001352-3 0.000000+0 0.000000+02025 2151     \n"
    " 8.674401+5 3.212123-1 4.117567+2 9.945676-4 0.000000+0 0.000000+02025 2151     \n"
    " 9.410248+5 3.500940-1 1.013363+3 1.003742-3 0.000000+0 0.000000+02025 2151     \n"
    " 9.650374+5 1.025332+0 1.609099+2 1.001175-3 0.000000+0 0.000000+02025 2151     \n"
    " 1.018558+6 9.932257-1 1.302725+3 9.572095-4 0.000000+0 0.000000+02025 2151     \n"
    " 1.024279+6 1.008507+0 2.021524+2 9.951288-4 0.000000+0 0.000000+02025 2151     \n"
    " 1.037285+6 1.006532+0 1.728796+3 1.001917-3 0.000000+0 0.000000+02025 2151     \n"
    " 1.057398+6 1.002368+0 8.359717+2 9.960064-4 0.000000+0 0.000000+02025 2151     \n"
    " 1.060311+6 1.003334+0 4.763447+2 9.888550-4 0.000000+0 0.000000+02025 2151     \n"
    " 1.100866+6 9.937983-1 4.382755+2 1.000022-3 0.000000+0 0.000000+02025 2151     \n"
    " 1.141864+6 9.202789-1 3.770592+2 9.967338-4 0.000000+0 0.000000+02025 2151     \n"
    " 1.242643+6 1.000077+0 1.730181+3 9.978615-4 0.000000+0 0.000000+02025 2151     \n"
    " 1.432040+6 1.000000+0 7.285957+2 9.796377-4 0.000000+0 0.000000+02025 2151     \n"
    " 1.439972+6 9.761681-1 3.451445+2 1.000397-3 0.000000+0 0.000000+02025 2151     \n"
    " 1.468021+6 9.389508-1 1.271009+2 9.842338-4 0.000000+0 0.000000+02025 2151     \n"
    " 1.470698+6 9.779341-1 1.445353+2 4.809440-4 0.000000+0 0.000000+02025 2151     \n"
    " 1.489149+6 1.109174+0 3.636708+2 1.000787-3 0.000000+0 0.000000+02025 2151     \n"
    " 1.496383+6 1.050878+0 1.448484+2 1.012630-3 0.000000+0 0.000000+02025 2151     \n"
    "-1.500000+0 0.000000+0          0          0         24          42025 2151     \n"
    " 1.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+02025 2151     \n"
    " 2.000000+0 1.000000+0 5.000000-1 0.000000+0 5.574746-1 4.993153-12025 2151     \n"
    " 4.000000+0 2.000000+0 3.500000+0 0.000000+0 4.993202-1 4.993202-12025 2151     \n"
    " 3.000000+0 1.000000+0 1.500000+0 0.000000+0 4.885639-1 4.885639-12025 2151     \n"
    " 0.000000+0 0.000000+0          0         52        312         522025 2151     \n"
    " 1.084167+4 6.385209-1 5.834810-1 0.000000+0 9.999433-4 0.000000+02025 2151     \n"
    " 1.933991+4 4.693232-1 7.927589-2 0.000000+0 1.000035-3 0.000000+02025 2151     \n"
    " 5.258763+4 3.948809-1 2.712932-1 0.000000+0 9.998658-4 0.000000+02025 2151     \n"
    " 8.505723+4 2.258309-1 1.717796+0 0.000000+0 9.938120-4 0.000000+02025 2151     \n"
    " 1.325016+5 3.097500-1 2.669534-1 0.000000+0 1.001213-3 0.000000+02025 2151     \n"
    " 1.515291+5 2.710230-1 1.357344+1 0.000000+0 9.999954-4 0.000000+02025 2151     \n"
    " 1.785771+5 5.279479-1 7.206823+0 0.000000+0 1.001054-3 0.000000+02025 2151     \n"
    " 2.215573+5 5.962590-1 2.379516+1 0.000000+0 1.007390-3 0.000000+02025 2151     \n"
    " 2.693832+5 3.962471-1 4.093579+1 0.000000+0 9.999634-4 0.000000+02025 2151     \n"
    " 2.884495+5 7.259105-1 7.754888-1 0.000000+0 9.994487-4 0.000000+02025 2151     \n"
    " 2.905101+5 2.733106-1 4.927727+1 0.000000+0 1.001839-3 0.000000+02025 2151     \n"
    " 3.210301+5 5.988639-1 2.819137-1 0.000000+0 9.976746-4 0.000000+02025 2151     \n"
    " 3.423883+5 3.902199-1 7.632160+1 0.000000+0 1.001149-3 0.000000+02025 2151     \n"
    " 3.457514+5 3.954839-1 1.226557+2 0.000000+0 1.001829-3 0.000000+02025 2151     \n"
    " 3.538809+5 1.011334+0 9.533776+1 0.000000+0 1.007196-3 0.000000+02025 2151     \n"
    " 3.715723+5 3.817660-1 1.081556+2 0.000000+0 1.001452-3 0.000000+02025 2151     \n"
    " 3.767322+5 1.000000+0 1.943685+1 0.000000+0 1.000646-3 0.000000+02025 2151     \n"
    " 3.973762+5 1.000000+0 1.337209+1 0.000000+0 9.997173-4 0.000000+02025 2151     \n"
    " 4.778264+5 3.811228-1 1.902890+2 0.000000+0 1.000948-3 0.000000+02025 2151     \n"
    " 4.871760+5 1.000000+0 3.910849+0 0.000000+0 1.007483-3 0.000000+02025 2151     \n"
    " 6.099266+5 1.000000+0 2.978815+1 1.00000-28 1.001369-3 0.000000+02025 2151     \n"
    " 6.406201+5 3.495588-1 6.761204+2 1.00000-23 9.996853-4 0.000000+02025 2151     \n"
    " 6.526481+5 1.000000+0 3.834609+1 1.00000-21 9.996850-4 0.000000+02025 2151     \n"
    " 6.886465+5 5.162762-1 1.024890+1 1.00000-20 9.999319-4 0.000000+02025 2151     \n"
    " 7.277809+5 4.139505-1 1.669951+2 5.000000-6 1.001089-3 0.000000+02025 2151     \n"
    " 7.643093+5 3.809662-1 1.351587+2 5.000000-6 9.988016-4 0.000000+02025 2151     \n"
    " 7.728142+5 1.000000+0 8.088075+1 1.500000-4 9.997971-4 0.000000+02025 2151     \n"
    " 8.256661+5 1.000496+0 3.727073+2 1.500000-4 1.000465-3 0.000000+02025 2151     \n"
    " 8.418854+5 8.794051-1 4.176216+2 1.500000-4 1.016792-3 0.000000+02025 2151     \n"
    " 8.523793+5 9.002441-1 2.022774+1 1.500000-4 9.984209-4 0.000000+02025 2151     \n"
    " 8.840479+5 4.158173-1 3.456310+2 1.500000-4 9.892646-4 0.000000+02025 2151     \n"
    " 9.075664+5 3.800084-1 6.364961+2 1.500000-4 9.961361-4 0.000000+02025 2151     \n"
    " 9.387972+5 1.002308+0 2.823054+1 1.500000-4 1.016902-3 0.000000+02025 2151     \n"
    " 9.585494+5 3.809018-1 4.532000+2 1.500000-4 9.982374-4 0.000000+02025 2151     \n"
    " 1.040712+6 1.000967+0 3.945334+2 1.500000-4 1.000298-3 0.000000+02025 2151     \n"
    " 1.082875+6 1.000676+0 7.546027+2 1.500000-4 9.991833-4 0.000000+02025 2151     \n"
    " 1.145045+6 1.060216+0 4.948621+2 1.500000-4 1.006641-3 0.000000+02025 2151     \n"
    " 1.159711+6 1.090563+0 4.579632+2 1.500000-4 9.981890-4 0.000000+02025 2151     \n"
    " 1.166595+6 1.083161+0 2.160592+2 1.500000-4 9.850694-4 0.000000+02025 2151     \n"
    " 1.185042+6 1.002685+0 1.879646+2 1.500000-4 1.018880-3 0.000000+02025 2151     \n"
    " 1.200423+6 1.159214+0 5.743287+3 1.500000-4 9.544719-4 0.000000+02025 2151     \n"
    " 1.214226+6 9.518523-1 6.166954+3 1.500000-4 9.991078-4 0.000000+02025 2151     \n"
    " 1.298746+6 1.000000+0 3.409365+2 1.500000-4 9.476255-4 0.000000+02025 2151     \n"
    " 1.301649+6 1.045509+0 1.917287+2 1.500000-4 9.971868-4 0.000000+02025 2151     \n"
    " 1.327471+6 1.000000+0 1.306775+2 1.000000-4 9.950707-4 0.000000+02025 2151     \n"
    " 1.347516+6 1.000000+0 1.159204+3 1.500000-4 9.972576-4 0.000000+02025 2151     \n"
    " 1.404028+6 1.000000+0 4.138882+1 1.500000-4 1.002915-3 0.000000+02025 2151     \n"
    " 1.406689+6 1.000000+0 4.318980+2 1.500000-4 9.989632-4 0.000000+02025 2151     \n"
    " 1.423854+6 1.000005+0 7.303676+2 1.500000-4 9.993239-4 0.000000+02025 2151     \n"
    " 1.445634+6 1.191034+0 1.491566+3 1.500000-4 1.000200-3 0.000000+02025 2151     \n"
    " 1.464709+6 1.016184+0 1.436454+2 1.500000-4 1.048600-3 0.000000+02025 2151     \n"
    " 1.477528+6 9.389508-1 1.603824+2 1.500000-4 1.008688-3 0.000000+02025 2151     \n"
    " 1.500000+0 0.000000+0          0          0         24          42025 2151     \n"
    " 1.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+02025 2151     \n"
    " 2.000000+0 2.000000+0 5.000000-1 0.000000+0 3.833842-1 4.993153-12025 2151     \n"
    " 3.000000+0 0.000000+0 1.500000+0 0.000000+0 4.885639-1 4.885639-12025 2151     \n"
    " 3.000000+0 2.000000+0 1.500000+0 0.000000+0 4.885639-1 4.885639-12025 2151     \n"
    " 0.000000+0 0.000000+0          0         54        324         542025 2151     \n"
    " 3.026620+3 5.099998-1 1.501024-4 1.000013-3 1.000013-3 0.000000+02025 2151     \n"
    " 4.657137+3 5.099760-1 4.947117-4 1.000004-3 1.000004-3 0.000000+02025 2151     \n"
    " 2.043708+4 7.210495-1 4.938822+0 9.999247-4 1.000164-3 0.000000+02025 2151     \n"
    " 4.623603+4 5.099926-1 2.026671-2 9.998318-4 9.998306-4 0.000000+02025 2151     \n"
    " 6.264166+4 1.121466+0 6.335089-2 9.998755-4 9.998746-4 0.000000+02025 2151     \n"
    " 6.557878+4 1.111692+0 5.383179-2 1.000216-3 1.000218-3 0.000000+02025 2151     \n"
    " 1.065398+5 4.707627-1 9.616531-2 1.000511-3 1.000668-3 0.000000+02025 2151     \n"
    " 1.385358+5 2.806143-1 9.549145-2 1.001714-3 1.000409-3 0.000000+02025 2151     \n"
    " 1.432626+5 3.596077-1 9.937466-1 1.000070-3 1.000070-3 0.000000+02025 2151     \n"
    " 1.500945+5 3.528211-1 6.032865-1 9.973723-4 9.968441-4 0.000000+02025 2151     \n"
    " 1.647324+5 5.324688-1 2.791209+0 1.000915-3 9.998078-4 0.000000+02025 2151     \n"
    " 1.722230+5 3.526227+0 8.406423-1 9.946027-4 9.967819-4 0.000000+02025 2151     \n"
    " 1.856688+5 1.239514+0 9.646669+0 9.975083-4 9.990511-4 0.000000+02025 2151     \n"
    " 1.924824+5 9.910902-1 3.775050-1 9.965835-4 9.999850-4 0.000000+02025 2151     \n"
    " 2.141463+5 6.883913-1 6.154837+0 1.000658-3 1.000664-3 0.000000+02025 2151     \n"
    " 2.227226+5 6.559927-1 9.108186+0 9.997763-4 9.980250-4 0.000000+02025 2151     \n"
    " 2.416358+5 1.093544+0 2.129680+1 1.001655-3 1.001670-3 0.000000+02025 2151     \n"
    " 2.805810+5 1.182492+0 7.251026+0 1.000530-3 9.941674-4 0.000000+02025 2151     \n"
    " 3.001983+5 6.199414-1 2.733042+1 1.000516-3 1.000521-3 0.000000+02025 2151     \n"
    " 3.798147+5 1.000000+0 1.287701+1 9.994548-4 9.994491-4 0.000000+02025 2151     \n"
    " 4.167255+5 5.128116-1 4.541446+0 1.004093-3 1.003108-3 0.000000+02025 2151     \n"
    " 4.700237+5 1.000000+0 4.799204+1 9.996709-4 9.996673-4 0.000000+02025 2151     \n"
    " 4.999691+5 1.000000+0 1.112486+1 9.980574-4 9.980354-4 0.000000+02025 2151     \n"
    " 5.035246+5 1.000000+0 4.336471+1 1.000728-3 1.000736-3 0.000000+02025 2151     \n"
    " 5.200488+5 1.000000+0 9.004734+0 1.000806-3 1.000815-3 0.000000+02025 2151     \n"
    " 5.250463+5 1.000000+0 5.904118+1 1.000435-3 1.000440-3 0.000000+02025 2151     \n"
    " 5.586053+5 5.100690-1 2.771598+2 1.004890-3 9.999780-4 0.000000+02025 2151     \n"
    " 6.293988+5 5.134487-1 1.499947+1 1.003408-3 1.001827-3 0.000000+02025 2151     \n"
    " 6.430281+5 1.000000+0 5.261832+1 1.000549-3 1.000556-3 0.000000+02025 2151     \n"
    " 7.472968+5 3.499870-1 4.616090+2 9.938156-4 9.984545-4 0.000000+02025 2151     \n"
    " 8.577094+5 5.106777-1 2.093856+2 1.012277-3 1.012457-3 0.000000+02025 2151     \n"
    " 8.736737+5 1.000000+0 9.569857+1 1.000208-3 1.000211-3 0.000000+02025 2151     \n"
    " 8.891634+5 9.002356-1 4.484038+1 9.926082-4 9.924987-4 0.000000+02025 2151     \n"
    " 8.914717+5 8.999572-1 4.317489+1 9.983370-4 9.983123-4 0.000000+02025 2151     \n"
    " 9.447745+5 5.119372-1 5.416408+2 9.977040-4 9.976686-4 0.000000+02025 2151     \n"
    " 9.688298+5 1.000000+0 2.363000-1 1.011136-3 1.011311-3 0.000000+02025 2151     \n"
    " 9.776050+5 1.012831+0 4.804160+1 9.963013-4 9.962432-4 0.000000+02025 2151     \n"
    " 9.933098+5 6.866854-1 9.094029+2 4.990968-5 1.000068-3 0.000000+02025 2151     \n"
    " 1.070749+6 1.003334+0 1.074788+2 5.031243-4 1.002502-3 0.000000+02025 2151     \n"
    " 1.094159+6 1.033640+0 3.819306+2 1.001033-3 1.001050-3 0.000000+02025 2151     \n"
    " 1.126744+6 1.013611+0 2.017282+2 1.004744-3 1.004826-3 0.000000+02025 2151     \n"
    " 1.128806+6 1.113102+0 8.801105+2 9.987654-4 9.987440-4 0.000000+02025 2151     \n"
    " 1.182505+6 1.001087+0 2.844307+2 1.020610-3 1.020983-3 0.000000+02025 2151     \n"
    " 1.231250+6 1.000157+0 8.685718+2 9.955052-4 9.954225-4 0.000000+02025 2151     \n"
    " 1.263048+6 1.000740+0 1.409171+3 1.008584-3 1.008746-3 0.000000+02025 2151     \n"
    " 1.281790+6 9.786015-1 2.563466+2 9.991459-4 9.991297-4 0.000000+02025 2151     \n"
    " 1.290230+6 1.000000+0 1.844469+3 9.909057-4 9.907326-4 0.000000+02025 2151     \n"
    " 1.333757+6 1.000000+0 5.978572+2 1.004338-3 9.995394-4 0.000000+02025 2151     \n"
    " 1.345129+6 1.000000+0 4.658871+2 9.996799-4 9.996736-4 0.000000+02025 2151     \n"
    " 1.380857+6 1.000000+0 5.188033+1 9.992091-4 9.991931-4 0.000000+02025 2151     \n"
    " 1.389886+6 1.000000+0 5.537658+2 5.014843-4 9.994047-4 0.000000+02025 2151     \n"
    " 1.432586+6 1.110417+0 2.342008+3 9.796801-4 9.984083-4 0.000000+02025 2151     \n"
    " 1.470814+6 1.061456+0 2.337906+2 1.015015-3 1.006809-3 0.000000+02025 2151     \n"
    " 1.515000+6 1.000000+0 2.500000+1 1.000000-3 1.000000-3 0.000000+02025 2151     \n"
    " 2.500000+0 0.000000+0          0          0         24          42025 2151     \n"
    " 1.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+02025 2151     \n"
    " 2.000000+0 2.000000+0 5.000000-1 0.000000+0 3.833842-1 4.993153-12025 2151     \n"
    " 4.000000+0 1.000000+0 3.500000+0 0.000000+0 4.993202-1 4.993202-12025 2151     \n"
    " 3.000000+0 2.000000+0 1.500000+0 0.000000+0 4.885639-1 4.885639-12025 2151     \n"
    " 0.000000+0 0.000000+0          0         50        300         502025 2151     \n"
    " 3.453854+4 1.306229-1 3.772489-1 0.000000+0 1.000300-3 0.000000+02025 2151     \n"
    " 4.209335+4 2.298787+0 5.142698-1 0.000000+0 1.000000-3 0.000000+02025 2151     \n"
    " 8.286527+4 4.473985-1 2.238498-2 0.000000+0 9.990652-4 0.000000+02025 2151     \n"
    " 1.102655+5 3.449248+0 8.887544-2 0.000000+0 1.001552-3 0.000000+02025 2151     \n"
    " 1.216381+5 9.203542-2 3.191748+0 0.000000+0 1.000441-3 0.000000+02025 2151     \n"
    " 1.591283+5 1.641489-2 9.477407-1 0.000000+0 1.001047-3 0.000000+02025 2151     \n"
    " 2.342025+5 1.246419+0 4.576535+0 0.000000+0 1.001027-3 0.000000+02025 2151     \n"
    " 2.584913+5 1.813201+0 1.722649+0 0.000000+0 9.874266-4 0.000000+02025 2151     \n"
    " 2.951489+5 1.021136+0 8.857884+1 0.000000+0 1.000587-3 0.000000+02025 2151     \n"
    " 3.085656+5 3.346849-1 1.108970+0 0.000000+0 1.000514-3 0.000000+02025 2151     \n"
    " 3.298763+5 6.653593-1 1.814720+0 0.000000+0 1.001090-3 0.000000+02025 2151     \n"
    " 3.397082+5 9.011701-1 2.576762-1 0.000000+0 9.968288-4 0.000000+02025 2151     \n"
    " 4.022092+5 5.160704-1 3.043052+0 0.000000+0 1.000203-3 0.000000+02025 2151     \n"
    " 4.751924+5 7.758394-1 4.342693+0 0.000000+0 1.000665-3 0.000000+02025 2151     \n"
    " 5.216675+5 1.000000+0 2.525939+1 0.000000+0 1.000546-3 0.000000+02025 2151     \n"
    " 5.699029+5 3.822176-1 1.924925+2 1.00000-44 9.991421-4 0.000000+02025 2151     \n"
    " 5.769344+5 1.000000+0 1.372460+1 1.00000-39 1.001567-3 0.000000+02025 2151     \n"
    " 5.941411+5 8.933261-1 1.082643+2 1.00000-31 1.003165-3 0.000000+02025 2151     \n"
    " 6.165436+5 1.000000+0 2.803231+1 1.00000-25 1.000673-3 0.000000+02025 2151     \n"
    " 6.230329+5 8.899142-1 7.007761+1 1.00000-24 1.003097-3 0.000000+02025 2151     \n"
    " 6.382422+5 8.906220-1 6.105679+1 1.00000-22 1.001349-3 0.000000+02025 2151     \n"
    " 7.411591+5 4.554076-1 7.018857+1 5.000000-6 1.002911-3 0.000000+02025 2151     \n"
    " 7.577087+5 1.169547+0 4.196178+2 5.000000-6 9.998576-4 0.000000+02025 2151     \n"
    " 7.898183+5 5.102981-1 7.601269+1 1.500000-4 9.979552-4 0.000000+02025 2151     \n"
    " 8.000581+5 8.905724-1 4.668667+1 1.500000-4 9.970821-4 0.000000+02025 2151     \n"
    " 8.297447+5 8.906578-1 3.507932+1 1.500000-4 1.002401-3 0.000000+02025 2151     \n"
    " 8.311155+5 8.901137-1 3.737703+1 1.500000-4 9.966001-4 0.000000+02025 2151     \n"
    " 8.417632+5 8.805770-1 4.232297+2 1.500000-4 1.000859-3 0.000000+02025 2151     \n"
    " 9.238486+5 1.008796+0 8.757183+0 1.500000-4 9.415124-4 0.000000+02025 2151     \n"
    " 9.244597+5 8.945398-1 2.227880+2 1.500000-4 9.989827-4 0.000000+02025 2151     \n"
    " 9.999858+5 9.224902-1 2.827924+1 1.500000-4 1.016953-3 0.000000+02025 2151     \n"
    " 1.045655+6 1.003072+0 2.247409+2 1.500000-4 9.989267-4 0.000000+02025 2151     \n"
    " 1.086756+6 1.005137+0 3.971382+1 1.500000-4 9.854369-4 0.000000+02025 2151     \n"
    " 1.098158+6 1.078827+0 8.173862+1 1.500000-4 1.001039-3 0.000000+02025 2151     \n"
    " 1.137219+6 1.222333+0 3.580701+2 1.500000-4 1.000237-3 0.000000+02025 2151     \n"
    " 1.169459+6 1.040866+0 4.714465+2 1.500000-4 9.993727-4 0.000000+02025 2151     \n"
    " 1.235112+6 1.000259+0 3.222590+2 1.500000-4 9.999556-4 0.000000+02025 2151     \n"
    " 1.249477+6 1.000056+0 9.000118+2 1.500000-4 9.968391-4 0.000000+02025 2151     \n"
    " 1.254854+6 9.938044-1 3.227746+1 1.500000-4 9.996581-4 0.000000+02025 2151     \n"
    " 1.257787+6 1.000058+0 1.629126+2 1.500000-4 1.001669-3 0.000000+02025 2151     \n"
    " 1.267411+6 9.997524-1 9.046679+2 1.500000-4 1.009044-3 0.000000+02025 2151     \n"
    " 1.306825+6 9.979807-1 4.509818+1 1.500000-4 1.007341-3 0.000000+02025 2151     \n"
    " 1.309848+6 1.000000+0 4.062436+2 1.500000-4 1.000226-3 0.000000+02025 2151     \n"
    " 1.338594+6 1.000000+0 2.752651+2 1.500000-4 1.000010-3 0.000000+02025 2151     \n"
    " 1.355534+6 1.000487+0 5.944995+0 1.500000-4 1.008256-3 0.000000+02025 2151     \n"
    " 1.356064+6 1.005737+0 1.275084+2 1.500000-4 1.002580-3 0.000000+02025 2151     \n"
    " 1.385422+6 9.992268-1 3.216466+1 1.500000-4 9.877915-4 0.000000+02025 2151     \n"
    " 1.458160+6 1.016184+0 1.566699+3 1.500000-4 1.033904-3 0.000000+02025 2151     \n"
    " 1.468385+6 1.061456+0 1.180005+2 1.500000-4 1.004753-3 0.000000+02025 2151     \n"
    " 1.530000+6 1.000000+0 4.500000+1 1.500000-4 1.000000-3 0.000000+02025 2151     \n"
    " 0.000000+0 0.000000+0          0          0          0          02025 2  0     \n";
}


std::string Cl35() {

  // Cl35 ENDF/B-VIII.0 LRF=7 resonance evaluation
  // particular features: - Z for n,gamma is set to 0.0 (should be 17.0)
  //                      - spin groups are not unitque (1- and 2- occur
  //                        multiple times)

  return
    " 1.703500+4 3.466845+1          0          0          1          01725 2151     \n"
    " 1.703500+4 1.000000+0          0          0          1          01725 2151     \n"
    " 1.000000-5 1.200000+6          1          7          0          01725 2151     \n"
    " 0.000000+0 0.000000+0          0          3          8          01725 2151     \n"
    " 0.000000+0 0.000000+0          3          0         36          61725 2151     \n"
    " 0.000000+0 3.565932+1 0.000000+0 0.000000+0 1.000000+0 0.000000+01725 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 1.020000+2 0.000000+0 0.000000+01725 2151     \n"
    " 1.000000+0 3.466845+1 0.000000+0 1.700000+1 5.000000-1 1.500000+01725 2151     \n"
    " 0.000000+0 1.000000+0 0.000000+0 2.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 9.986235-1 3.466863+1 1.000000+0 1.600000+1 5.000000-1 1.500000+01725 2151     \n"
    " 6.152200+5 1.000000+0 0.000000+0 6.000000+2 0.000000+0 0.000000+01725 2151     \n"
    " 1.000000+0 0.000000+0          0          0         18          31725 2151     \n"
    " 1.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.000000+0 0.000000+0 1.000000+0 0.000000+0 3.667980-1 4.822220-11725 2151     \n"
    " 3.000000+0 0.000000+0 1.000000+0 0.000000+0 3.667980-1 4.822220-11725 2151     \n"
    " 0.000000+0 0.000000+0          0         23        138         231725 2151     \n"
    " 5.493200+4 3.672600-1 4.644240+1 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 6.823616+4 3.933600-1 2.179040+2 1.000000-5 0.000000+0 0.000000+01725 2151     \n"
    " 1.150980+5 7.390000-1 4.307780+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.825230+5 7.451500-1 1.759740+3 4.000000-1 0.000000+0 0.000000+01725 2151     \n"
    " 2.397427+5 6.871600-1 2.685470+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 3.351287+5 3.583800-1 5.525660+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 3.991469+5 7.409700-1 1.093810+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 4.156650+5 3.286800-1 1.146260+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 4.506303+5 3.929900-1 4.613320+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 4.997896+5 6.704600-1 2.312410+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 5.422214+5 6.060000-1 5.220340+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 6.542074+5 6.060000-1 5.186030+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 6.652896+5 6.060000-1 1.406720+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 6.723573+5 6.060000-1 5.754080+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 6.780000+5 6.060000-1 1.331580+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 6.949358+5 6.060000-1 1.824950+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 7.251773+5 6.060000-1 1.672450+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 7.546912+5 6.060000-1 2.690680+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 7.814646+5 6.060000-1 1.294130+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 8.315787+5 8.600000-1 4.453040+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 8.829776+5 6.060000-1 9.780610+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.109188+6 6.060000-1 2.742350+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.205687+6 6.060000-1 6.425840+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.000000+0 0.000000+0          0          0         18          31725 2151     \n"
    " 1.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.000000+0 0.000000+0 2.000000+0 0.000000+0 3.667980-1 4.822220-11725 2151     \n"
    " 3.000000+0 0.000000+0 2.000000+0 0.000000+0 3.667980-1 4.822220-11725 2151     \n"
    " 0.000000+0 0.000000+0          0         32        192         321725 2151     \n"
    "-1.806500+2 5.301500-1 1.327700+1 5.992300-3 0.000000+0 0.000000+01725 2151     \n"
    " 1.480195+4 3.456800-1 3.259950+1 2.800020-2 0.000000+0 0.000000+01725 2151     \n"
    " 2.661579+4 3.041500-1 1.154980+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 7.549467+4 8.059200-1 7.991540-2 1.000000-5 0.000000+0 0.000000+01725 2151     \n"
    " 1.304435+5 7.593500-1 7.665700-1 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.149232+5 3.485100-1 6.528390+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.425998+5 9.018300-1 3.440090+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.846570+5 4.994200-1 1.194020+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 3.137503+5 4.763300-1 1.475680+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 3.310710+5 3.982400-1 3.276200+1 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 3.389756+5 2.383300+0 8.777960+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 3.802625+5 4.630100-1 1.178590+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 3.864177+5 4.180700-1 1.240150+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 4.593615+5 6.156500-1 3.859970+1 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 4.886274+5 4.081400-1 7.806800+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 5.143982+5 6.060000-1 5.284320+3 4.000000-1 0.000000+0 0.000000+01725 2151     \n"
    " 5.279391+5 6.060000-1 2.812390+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 5.590263+5 6.060000-1 1.813310+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 5.908121+5 6.060000-1 3.235630+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 6.216780+5 6.060000-1 6.758690+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 6.310363+5 6.060000-1 1.636760+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 6.408054+5 6.060000-1 8.597480+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 6.578143+5 6.060000-1 2.069420+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 6.810575+5 6.060000-1 6.654900+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 6.902708+5 6.060000-1 4.632110+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 7.799749+5 6.060000-1 5.311180+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 8.249523+5 6.060000-1 2.763840+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 8.609610+5 6.060000-1 7.202410+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.053393+6 6.060000-1 2.038130+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.137885+6 6.060000-1 1.835530+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.257680+6 6.060000-1 1.749830+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 7.563145+6 3.839800-1 6.219050+5 1.000000+3 0.000000+0 0.000000+01725 2151     \n"
    " 0.000000+0-1.000000+0          0          0         18          31725 2151     \n"
    " 1.000000+0 1.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.000000+0 1.000000+0 1.000000+0 0.000000+0 4.888750-1 4.822220-11725 2151     \n"
    " 3.000000+0 1.000000+0 1.000000+0 0.000000+0 4.888750-1 4.822220-11725 2151     \n"
    " 0.000000+0 0.000000+0          0          9         54          91725 2151     \n"
    " 2.239640+4 1.724800+0 9.663670-1 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 9.660415+4 1.565000+0 2.761380+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.300740+5 3.240100-1 8.118180+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.432153+5 7.032000-1 2.174300+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.454785+5 8.319700-1 6.556420+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 3.413389+5 6.549700-1 5.720320+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 4.812809+5 7.821300-1 1.774270+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 5.048567+5 8.600000-1 6.790940+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 5.478545+5 8.600000-1 7.640130+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    "-1.000000+0 0.000000+0          0          0         18          31725 2151     \n"
    " 1.000000+0 1.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.000000+0 1.000000+0 1.000000+0 0.000000+0 4.888750-1 4.822220-11725 2151     \n"
    " 3.000000+0 1.000000+0 1.000000+0 0.000000+0 4.888750-1 4.822220-11725 2151     \n"
    " 0.000000+0 0.000000+0          0         38        228         381725 2151     \n"
    " 4.250762+3 4.720000-1 6.280000-1 2.300000-1 0.000000+0 0.000000+01725 2151     \n"
    " 5.491020+3 9.702100-1 3.863540-3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 3.776792+4 1.912000-1 4.407970-1 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 4.416635+4 1.042900+0 3.054660+1 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 6.277945+4 6.209000-1 1.345600+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.035154+5 3.880500-1 3.819530+2 1.972980+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.351161+5 3.410900-1 1.875460+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.625608+5 6.471800-1 5.634890+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.654823+5 1.050100+0 2.072840+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.199922+5 4.001800-1 3.847750+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.241121+5 4.033900-1 7.567220+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.288860+5 5.941900-1 1.768540+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.615514+5 8.458200-1 1.063750+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.791149+5 3.767200-1 1.254490+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.845043+5 9.053700-1 4.161730+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 3.060336+5 5.628500-1 7.718360+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 3.074111+5 6.403400-1 1.172260+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 3.406957+5 5.654300-1 3.875200+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 3.819655+5 8.183600+0 1.602770+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 3.878895+5 9.112300-1 3.238570+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 4.044696+5 1.492300+0 1.010090+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 4.526283+5 8.600500-1 3.430810+3 4.000000-1 0.000000+0 0.000000+01725 2151     \n"
    " 4.653384+5 3.666800-1 6.938600+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 5.526254+5 8.600000-1 2.791240+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 5.592470+5 8.600000-1 5.857230+2 4.000000-1 0.000000+0 0.000000+01725 2151     \n"
    " 5.913358+5 8.600000-1 1.674080+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 6.008407+5 8.600000-1 1.208510+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 6.295837+5 8.600000-1 6.225330+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 6.593525+5 8.600000-1 1.356080+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 7.183788+5 8.600000-1 9.740400+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 7.328902+5 8.600000-1 1.955950+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 8.103378+5 8.600000-1 5.088680+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 8.715664+5 8.600000-1 1.987240+3 4.000000-1 0.000000+0 0.000000+01725 2151     \n"
    " 9.439464+5 8.600000-1 5.556340+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 9.845132+5 8.600000-1 1.801810+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.138720+6 8.600000-1 9.544250+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.353551+6 8.600000-1 3.879680+4 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.435502+6 8.600000-1 5.365630+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    "-2.000000+0 0.000000+0          0          0         18          31725 2151     \n"
    " 1.000000+0 1.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.000000+0 1.000000+0 1.000000+0 0.000000+0 4.888750-1 4.822220-11725 2151     \n"
    " 3.000000+0 1.000000+0 1.000000+0 0.000000+0 4.888750-1 4.822220-11725 2151     \n"
    " 0.000000+0 0.000000+0          0         67        402         671725 2151     \n"
    "-3.369334+5 5.340100-1 3.820180+4 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.734636+4 4.579210-1 6.027800+0 1.472210-1 0.000000+0 0.000000+01725 2151     \n"
    " 5.297412+4 5.624100-1 8.159470-1 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 5.781159+4 5.380650-1 1.073890+2 9.980470-1 0.000000+0 0.000000+01725 2151     \n"
    " 9.042023+4 7.164960-1 2.178860+1 2.740270-1 0.000000+0 0.000000+01725 2151     \n"
    " 9.052563+4 1.275920-1 4.235150+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.134033+5 3.369900-1 1.422880+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.408311+5 5.458500-1 9.871460+1 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.498295+5 7.555300-1 1.132080+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.926794+5 2.287600-1 3.381170+1 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.010949+5 2.905500-1 3.650400+1 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.066160+5 5.943300-1 5.594230-1 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.145480+5 2.323800-1 4.238200+1 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.153514+5 7.735500-1 4.558520+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.170986+5 6.185100-1 5.772230+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.213851+5 1.592700+0 4.073280+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.458488+5 7.646700-1 5.612490+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.501986+5 4.046000-1 4.345790+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.856091+5 8.420300-1 1.568700+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.900777+5 1.615200+0 1.520800+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 3.014487+5 4.017200-1 8.855620+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 3.038357+5 1.697600+0 6.957360+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 3.455490+5 3.479100-1 6.319440+1 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 3.543783+5 2.329000-1 8.058180+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 3.844697+5 5.785800+0 1.985550+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 4.078202+5 5.513100-1 2.040430+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 4.386472+5 6.448700-1 1.273690+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 4.444759+5 2.050000-1 1.810450+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 4.513762+5 1.759900+0 5.759080+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 4.578911+5 1.147400+0 7.322610+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 4.655024+5 5.580100-1 4.040480+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 4.691015+5 9.520300-1 3.761160+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 4.771962+5 5.403700-1 1.658970+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 4.854618+5 1.273200+0 4.799550+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 5.437824+5 8.600000-1 6.917280+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 5.645835+5 8.600000-1 5.359800+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 5.738796+5 8.600000-1 1.469740+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 6.330964+5 8.600000-1 1.445360+2 4.000000-1 0.000000+0 0.000000+01725 2151     \n"
    " 6.427422+5 8.600000-1 9.499240+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 7.125272+5 8.600000-1 3.300900+3 4.000000-1 0.000000+0 0.000000+01725 2151     \n"
    " 7.295798+5 8.600000-1 6.465640+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 7.398106+5 8.600000-1 3.042130+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 7.702915+5 8.600000-1 3.185460+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 7.985633+5 8.600000-1 2.910980+3 4.000000-1 0.000000+0 0.000000+01725 2151     \n"
    " 8.273814+5 6.060000-1 4.360650+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 8.322588+5 8.600000-1 1.382530+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 8.484064+5 8.600000-1 2.391310+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 8.524894+5 8.600000-1 1.726720+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 8.950719+5 8.600000-1 1.505930+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 9.156492+5 8.600000-1 1.033970+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 9.331357+5 8.600000-1 9.386810+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 9.533511+5 8.600000-1 1.000950+3 4.000000-1 0.000000+0 0.000000+01725 2151     \n"
    " 9.915402+5 8.600000-1 4.348100+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.006461+6 8.600000-1 3.591290+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.028543+6 8.600000-1 2.388860+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.055512+6 8.600000-1 3.888630+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.074009+6 8.600000-1 1.039250+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.080187+6 8.600000-1 7.511970+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.115662+6 8.600000-1 3.444140+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.131997+6 8.600000-1 1.523380+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.165343+6 8.600000-1 7.309520+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.198501+6 8.600000-1 3.243460+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.218130+6 8.600000-1 3.246660+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.243534+6 8.600000-1 3.096960+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.311787+6 8.600000-1 7.946430+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.354056+6 8.600000-1 1.103350+4 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.434336+6 8.600000-1 5.422940+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    "-1.000000+0 0.000000+0          0          0         18          31725 2151     \n"
    " 1.000000+0 1.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.000000+0 1.000000+0 2.000000+0 0.000000+0 4.888750-1 4.822220-11725 2151     \n"
    " 3.000000+0 1.000000+0 2.000000+0 0.000000+0 4.888750-1 4.822220-11725 2151     \n"
    " 0.000000+0 0.000000+0          0         18        108         181725 2151     \n"
    " 1.339884+5 2.314200+0 6.600970+2 4.000000-1 0.000000+0 0.000000+01725 2151     \n"
    " 2.251416+5 1.346200+0 5.685270+2 4.000000-1 0.000000+0 0.000000+01725 2151     \n"
    " 3.728997+5 2.107900+0 1.796360+3 4.000000-1 0.000000+0 0.000000+01725 2151     \n"
    " 4.019221+5 3.329500-1 1.212200+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 4.221047+5 2.868200+0 8.577140+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 4.836869+5 1.907300+0 2.735270+3 4.000000-1 0.000000+0 0.000000+01725 2151     \n"
    " 5.811770+5 8.600000-1 2.087720+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 6.142552+5 8.600000-1 1.153480+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 7.219908+5 8.600000-1 8.268190+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 7.484834+5 8.600000-1 1.023870+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 7.654137+5 8.600000-1 1.519860+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 7.927810+5 8.600000-1 1.040840+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 8.451494+5 8.600000-1 4.176760+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 9.366575+5 8.600000-1 6.061270+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 9.747505+5 8.600000-1 3.259090+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.091954+6 8.600000-1 1.974980+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.177028+6 8.600000-1 1.488160+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.225254+6 8.600000-1 1.807190+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    "-2.000000+0 0.000000+0          0          0         18          31725 2151     \n"
    " 1.000000+0 1.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.000000+0 1.000000+0 2.000000+0 0.000000+0 4.888750-1 4.822220-11725 2151     \n"
    " 3.000000+0 1.000000+0 2.000000+0 0.000000+0 4.888750-1 4.822220-11725 2151     \n"
    " 0.000000+0 0.000000+0          0         28        168         281725 2151     \n"
    " 3.978154+2 6.650000-1 5.050000-2 3.220000-1 0.000000+0 0.000000+01725 2151     \n"
    " 1.136064+5 2.947700-1 3.972030+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.430242+5 4.921800-1 3.147430+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.059671+5 4.506700-1 5.887780-1 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.837547+5 6.466700-1 4.237100+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.870089+5 5.536890-1 2.126030+1 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.949461+5 7.387600-1 4.747470+2 4.000000-1 0.000000+0 0.000000+01725 2151     \n"
    " 5.903574+5 8.600000-1 7.550480+2 4.000000-1 0.000000+0 0.000000+01725 2151     \n"
    " 6.181025+5 8.600000-1 1.075440+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 7.036239+5 8.600000-1 9.311240+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 7.350094+5 8.600000-1 3.777530+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 8.013281+5 8.600000-1 6.807230+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 8.358532+5 8.600000-1 1.181190+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 8.769947+5 8.600000-1 5.757280+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 9.225158+5 8.600000-1 1.571600+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 9.461013+5 8.600000-1 4.806310+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 9.810389+5 8.600000-1 4.229000+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 9.998479+5 8.600000-1 3.022270+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.016475+6 8.600000-1 4.019670+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.050581+6 8.600000-1 1.304140+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.088712+6 8.600000-1 5.239370+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.116025+6 8.600000-1 4.363740+4 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.172051+6 8.600000-1 1.494510+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.283751+6 8.600000-1 5.041400+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.337132+6 8.600000-1 4.850280+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.366148+6 8.600000-1 8.350230+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.403885+6 8.600000-1 8.067070+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.441365+6 8.600000-1 1.608740+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    "-3.000000+0 0.000000+0          0          0         18          31725 2151     \n"
    " 1.000000+0 1.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.000000+0 1.000000+0 2.000000+0 0.000000+0 4.888750-1 4.822220-11725 2151     \n"
    " 3.000000+0 1.000000+0 2.000000+0 0.000000+0 4.888750-1 4.822220-11725 2151     \n"
    " 0.000000+0 0.000000+0          0         57        342         571725 2151     \n"
    " 1.635612+4 3.865000-1 5.981800+0 1.640190-1 0.000000+0 0.000000+01725 2151     \n"
    " 1.713387+4 8.022800-1 1.409590+1 3.199990-2 0.000000+0 0.000000+01725 2151     \n"
    " 4.027028+4 5.774900-1 1.773460-1 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 5.160803+4 4.471700-2 2.417150+0 9.599660-2 0.000000+0 0.000000+01725 2151     \n"
    " 9.520675+4 4.531600-1 1.549000-1 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 9.944099+4 2.322700-1 2.393210+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.120451+5 3.240700-1 2.642330-1 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.400910+5 3.664800-1 3.778030+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.416415+5 3.131300-1 3.954240+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.529224+5 2.989300-1 3.820560-1 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.835405+5 3.333300-1 4.615080+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.852814+5 4.665000-1 4.482330+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.881523+5 4.975300-1 4.222480+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.901806+5 2.935300-1 1.024950+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.929420+5 7.559100-1 1.635970+1 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.991683+5 2.550000-1 2.775550+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 2.933823+5 1.802700+0 5.941830+0 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 3.367631+5 4.162100-1 2.890370+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 4.753035+5 8.300600-1 2.791800+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 5.097591+5 8.600000-1 3.870270+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 5.298958+5 8.600000-1 1.353670+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 5.352483+5 8.600000-1 4.481020+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 6.080558+5 8.600000-1 7.747700+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 6.657961+5 8.600000-1 2.099450+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 6.739830+5 6.060000-1 4.896460+1 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 6.816687+5 8.600000-1 6.028770+1 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 6.850220+5 8.600000-1 1.172880+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 7.010023+5 8.600000-1 1.251920+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 7.574622+5 8.600000-1 3.780250+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 7.611067+5 8.600000-1 2.238910+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 7.749692+5 8.600000-1 3.689300+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 8.068823+5 8.600000-1 6.298260+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 8.389452+5 8.600000-1 3.394540+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 8.626053+5 8.600000-1 6.945570+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 8.865767+5 8.600000-1 3.197900+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 9.058611+5 8.600000-1 4.411560+2 4.000000-1 0.000000+0 0.000000+01725 2151     \n"
    " 9.109009+5 8.600000-1 3.339470+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 9.504659+5 8.600000-1 6.451080+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.010401+6 8.600000-1 7.910460+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.033683+6 8.600000-1 3.549570+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.062148+6 8.600000-1 4.207760+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.071198+6 8.600000-1 1.944480+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.085425+6 8.600000-1 1.339760+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.103724+6 8.600000-1 1.609290+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.120932+6 8.600000-1 3.890220+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.126914+6 8.600000-1 6.155220+2 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.144865+6 8.600000-1 1.586600+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.155330+6 8.600000-1 1.577740+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.189685+6 8.600000-1 7.290220+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.209020+6 8.600000-1 3.485080+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.237194+6 8.600000-1 5.888380+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.268642+6 8.600000-1 2.391480+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.277668+6 8.600000-1 2.984480+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.315129+6 8.600000-1 1.210580+4 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.390745+6 8.600000-1 5.502450+3 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.425393+6 8.600000-1 1.381970+4 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    " 1.485128+6 8.600000-1 1.054090+4 0.000000+0 0.000000+0 0.000000+01725 2151     \n"
    "                                                                  1725 2  0     \n";
}
