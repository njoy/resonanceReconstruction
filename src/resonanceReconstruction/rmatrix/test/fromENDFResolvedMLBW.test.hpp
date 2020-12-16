std::string mlbwRh105();
std::string resolvedAg107();
std::string Tm168();
std::string Dy160();
std::string Dy162();

SCENARIO( "fromENDF - LRF2" ) {

  GIVEN( "valid ENDF data for Rh105" ) {

    std::string string = mlbwRh105();
    auto begin = string.begin();
    auto end = string.end();
    long lineNumber = 1;

    njoy::ENDFtk::HeadRecord head( begin, end, lineNumber );
    njoy::ENDFtk::section::Type< 2, 151 > endf( head, begin, end, lineNumber, 4531 );
    ResonanceRange endfResonanceRange = endf.isotopes().front().resonanceRanges().front();

    auto resonances = fromENDF( endfResonanceRange, neutronMass, elementaryCharge, ParticleID( "n" ), ParticleID( "Rh105" ) );

    double a = 0.123 * std::pow( 104.005 * 1.008664, 1. / 3. ) + 0.08;

    THEN( "the appropriate CompoundSystem is returned" ) {

      CHECK( true == resonances.isResolved() );
      CHECK( false == resonances.isUnresolved() );
      CHECK( 1e-5 == Approx( resonances.lowerEnergy().value ) );
      CHECK( 7.5 == Approx( resonances.upperEnergy().value ) );
      CHECK( false == bool( resonances.interpolation() ) );

      auto compoundsystem = std::get< legacy::resolved::CompoundSystem< MultiLevelBreitWigner > >( resonances.compoundSystem() );

      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // content verification
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      // spin groups
      auto spingroups = compoundsystem.spinGroups();
      CHECK( 1 == spingroups.size() );

      // check the minimal energy grid
      auto grid = compoundsystem.grid();
      CHECK( 3 == grid.size() );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0 -  lJ = 0,1+
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup0 = spingroups[0];

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0, elastic channel
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel00 = spingroup0.incidentChannel();
      CHECK( "n,Rh105->n,Rh105" == channel00.reactionID().symbol() );

      // incident particle pair
      const auto incident00 = channel00.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident00.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident00.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident00.particle().spin() ) );
      CHECK( +1 == incident00.particle().parity() );
      CHECK( 104.005 * 1.008664 == Approx( incident00.residual().mass().value ) );
      CHECK( 45.0 * 1.602e-19 == Approx( incident00.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident00.residual().spin() ) );
      CHECK( +1 == incident00.residual().parity() );
      CHECK( "n,Rh105" == incident00.pairID().symbol() );

      // particle pair
      const auto pair00 = channel00.particlePair();
      CHECK( 1.008664 == Approx( pair00.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair00.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair00.particle().spin() ) );
      CHECK( +1 == pair00.particle().parity() );
      CHECK( 104.005 * 1.008664 == Approx( pair00.residual().mass().value ) );
      CHECK( 45.0 * 1.602e-19 == Approx( pair00.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair00.residual().spin() ) );
      CHECK( +1 == pair00.residual().parity() );
      CHECK( "n,Rh105" == pair00.pairID().symbol() );

      // quantum numbers
      const auto numbers00 = channel00.quantumNumbers();
      CHECK( 0 == numbers00.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers00.spin() );
      CHECK( 1.0 == numbers00.totalAngularMomentum() );
      CHECK( +1 == numbers00.parity() );
      CHECK( "{0,1,1+}" == numbers00.toString() );

      // radii
      const auto radii00 = channel00.radii();
      CHECK( a == Approx( radii00.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( a == Approx( radii00.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .62 == Approx( radii00.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel00.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel00.Q().value ) );

      // resonance table
      auto table0 = spingroup0.resonanceTable();

      CHECK( 2 == table0.numberResonances() );
      CHECK( 2 == table0.resonances().size() );
      CHECK( 2 == table0.energies().size() );
      CHECK( -5. == Approx( table0.energies().front().value ) );
      CHECK( 5. == Approx( table0.energies().back().value ) );

      auto rfront0 = table0.resonances().front();
      CHECK( -5. == Approx( rfront0.energy().value ) );
      CHECK( 1.45 == Approx( rfront0.elastic().value ) );
      CHECK( 0.16 == Approx( rfront0.capture().value ) );
      CHECK( 0 == Approx( rfront0.fission().value ) );
      CHECK( 0 == Approx( rfront0.competition().value ) );

      auto rback0 = table0.resonances().back();
      CHECK( 5. == Approx( rback0.energy().value ) );
      CHECK( 0.33 == Approx( rback0.elastic().value ) );
      CHECK( 0.16 == Approx( rback0.capture().value ) );
      CHECK( 0 == Approx( rback0.fission().value ) );
      CHECK( 0 == Approx( rback0.competition().value ) );

      // check the minimal energy grid
      auto grid0 = spingroup0.grid();
      CHECK( 3 == grid0.size() );
      CHECK( 4.755 == Approx( grid0[0].value ) );
      CHECK( 5. == Approx( grid0[1].value ) );
      CHECK( 5.245 == Approx( grid0[2].value ) );
    } // THEN

    THEN( "cross sections can be reconstructed" ) {

      // values taken from NJOY2016 for ENDF/B-VII.1 Rh105

      ReactionID elas( "n,Rh105->n,Rh105" );
      ReactionID capt( "n,Rh105->capture" );
      Map< ReactionID, CrossSection > xs;

      xs = resonances( 1e-5 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 5268.5966369331500 == Approx( xs[ elas ].value ) );
      CHECK( 801565.16324338294  == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 5268.2765170745824 == Approx( xs[ elas ].value ) );
      CHECK( 253468.26154893031 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 5265.2500341956293 == Approx( xs[ elas ].value ) );
      CHECK( 80132.232264999941 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-2 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 5235.6252833826566 == Approx( xs[ elas ].value ) );
      CHECK( 25279.102640176850 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-1 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 4949.4375301958371 == Approx( xs[ elas ].value ) );
      CHECK( 7816.7206941266395 == Approx( xs[ capt ].value ) );

      xs = resonances( 1. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 2723.0907713948382 == Approx( xs[ elas ].value ) );
      CHECK( 2161.1909561504717 == Approx( xs[ capt ].value ) );

      xs = resonances( 4.755 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 74816.910028973463 == Approx( xs[ elas ].value ) );
      CHECK( 45893.526706445802 == Approx( xs[ capt ].value ) );

      xs = resonances( 5. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 185939.51714628452 == Approx( xs[ elas ].value ) );
      CHECK( 87782.287793509589 == Approx( xs[ capt ].value ) );

      xs = resonances( 5.245 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 111838.57454107259 == Approx( xs[ elas ].value ) );
      CHECK( 42264.081005233311 == Approx( xs[ capt ].value ) );
    } // THEN
  } // GIVEN

  GIVEN( "valid ENDF data for Ag107" ) {

    std::string string = resolvedAg107();
    auto begin = string.begin();
    auto end = string.end();
    long lineNumber = 1;

    njoy::ENDFtk::HeadRecord head( begin, end, lineNumber );
    njoy::ENDFtk::section::Type< 2, 151 > endf( head, begin, end, lineNumber, 4725 );
    ResonanceRange endfResonanceRange = endf.isotopes().front().resonanceRanges().front();

    auto resonances = fromENDF( endfResonanceRange, neutronMass, elementaryCharge, ParticleID( "n" ), ParticleID( "Ag107" ) );

    double a = 0.123 * std::pow( 105.987 * 1.008664, 1. / 3. ) + 0.08;

    THEN( "the appropriate CompoundSystem is returned" ) {

      CHECK( true == resonances.isResolved() );
      CHECK( false == resonances.isUnresolved() );
      CHECK( 1e-5 == Approx( resonances.lowerEnergy().value ) );
      CHECK( 6500. == Approx( resonances.upperEnergy().value ) );
      CHECK( false == bool( resonances.interpolation() ) );

      auto compoundsystem = std::get< legacy::resolved::CompoundSystem< MultiLevelBreitWigner > >( resonances.compoundSystem() );

      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // content verification
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      // spin groups
      auto spingroups = compoundsystem.spinGroups();
      CHECK( 5 == spingroups.size() );

      // check the minimal energy grid
      auto grid = compoundsystem.grid();
      CHECK( 1196 == grid.size() ); // ( 400 resonances - 1 negative ) * 3 - 1 duplicate

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0 -  lJ = 0,0+
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup0 = spingroups[0];

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0, elastic channel
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel00 = spingroup0.incidentChannel();
      CHECK( "n,Ag107->n,Ag107" == channel00.reactionID().symbol() );

      // incident particle pair
      const auto incident00 = channel00.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident00.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident00.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident00.particle().spin() ) );
      CHECK( +1 == incident00.particle().parity() );
      CHECK( 105.987 * 1.008664 == Approx( incident00.residual().mass().value ) );
      CHECK( 47.0 * 1.602e-19 == Approx( incident00.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident00.residual().spin() ) );
      CHECK( +1 == incident00.residual().parity() );
      CHECK( "n,Ag107" == incident00.pairID().symbol() );

      // particle pair
      const auto pair00 = channel00.particlePair();
      CHECK( 1.008664 == Approx( pair00.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair00.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair00.particle().spin() ) );
      CHECK( +1 == pair00.particle().parity() );
      CHECK( 105.987 * 1.008664 == Approx( pair00.residual().mass().value ) );
      CHECK( 47.0 * 1.602e-19 == Approx( pair00.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair00.residual().spin() ) );
      CHECK( +1 == pair00.residual().parity() );
      CHECK( "n,Ag107" == pair00.pairID().symbol() );

      // quantum numbers
      const auto numbers00 = channel00.quantumNumbers();
      CHECK( 0 == numbers00.orbitalAngularMomentum() );
      CHECK( 0 == numbers00.spin() );
      CHECK( 0. == numbers00.totalAngularMomentum() );
      CHECK( +1 == numbers00.parity() );
      CHECK( "{0,0,0+}" == numbers00.toString() );

      // radii - NRO=1 NAPS=0
      const auto radii00 = channel00.radii();
      CHECK( a == Approx( radii00.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( a == Approx( radii00.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .66 == Approx( radii00.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel00.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel00.Q().value ) );

      // resonance table
      auto table0 = spingroup0.resonanceTable();

      CHECK( 46 == table0.numberResonances() );
      CHECK( 46 == table0.resonances().size() );
      CHECK( 46 == table0.energies().size() );
      CHECK( 16.3 == Approx( table0.energies().front().value ) );
      CHECK( 6387.0 == Approx( table0.energies().back().value ) );

      auto rfront0 = table0.resonances().front();
      CHECK( 16.3 == Approx( rfront0.energy().value ) );
      CHECK( 1.160000e-2 == Approx( rfront0.elastic().value ) );
      CHECK( 1.350000e-1 == Approx( rfront0.capture().value ) );
      CHECK( 0 == Approx( rfront0.fission().value ) );
      CHECK( 0 == Approx( rfront0.competition().value ) );

      auto rback0 = table0.resonances().back();
      CHECK( 6387. == Approx( rback0.energy().value ) );
      CHECK( 3.000000e-1 == Approx( rback0.elastic().value ) );
      CHECK( 1.430000e-1 == Approx( rback0.capture().value ) );
      CHECK( 0 == Approx( rback0.fission().value ) );
      CHECK( 0 == Approx( rback0.competition().value ) );

      // check the minimal energy grid
      auto grid0 = spingroup0.grid();
      CHECK( 138 == grid0.size() ); // ( 46 resonances ) * 3

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1 -  lJ = 0,1+
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup1 = spingroups[1];

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1, elastic channel
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel10 = spingroup1.incidentChannel();
      CHECK( "n,Ag107->n,Ag107" == channel10.reactionID().symbol() );

      // incident particle pair
      const auto incident10 = channel10.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident10.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident10.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident10.particle().spin() ) );
      CHECK( +1 == incident10.particle().parity() );
      CHECK( 105.987 * 1.008664 == Approx( incident10.residual().mass().value ) );
      CHECK( 47.0 * 1.602e-19 == Approx( incident10.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident10.residual().spin() ) );
      CHECK( +1 == incident10.residual().parity() );
      CHECK( "n,Ag107" == incident10.pairID().symbol() );

      // particle pair
      const auto pair10 = channel10.particlePair();
      CHECK( 1.008664 == Approx( pair10.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair10.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair10.particle().spin() ) );
      CHECK( +1 == pair10.particle().parity() );
      CHECK( 105.987 * 1.008664 == Approx( pair10.residual().mass().value ) );
      CHECK( 47.0 * 1.602e-19 == Approx( pair10.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair10.residual().spin() ) );
      CHECK( +1 == pair10.residual().parity() );
      CHECK( "n,Ag107" == pair10.pairID().symbol() );

      // quantum numbers
      const auto numbers10 = channel10.quantumNumbers();
      CHECK( 0 == numbers10.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers10.spin() );
      CHECK( 1.0 == numbers10.totalAngularMomentum() );
      CHECK( +1 == numbers10.parity() );
      CHECK( "{0,1,1+}" == numbers10.toString() );

      // radii - NRO=1 NAPS=0
      const auto radii10 = channel10.radii();
      CHECK( a == Approx( radii10.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( a == Approx( radii10.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .66 == Approx( radii10.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel10.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel10.Q().value ) );

      // resonance table
      auto table1 = spingroup1.resonanceTable();

      CHECK( 158 == table1.numberResonances() );
      CHECK( 158 == table1.resonances().size() );
      CHECK( 158 == table1.energies().size() );
      CHECK( -12.6 == Approx( table1.energies().front().value ) );
      CHECK( 6499. == Approx( table1.energies().back().value ) );

      auto rfront1 = table1.resonances().front();
      CHECK( -12.6 == Approx( rfront1.energy().value ) );
      CHECK( 4.427000e-2  == Approx( rfront1.elastic().value ) );
      CHECK( 1.430000e-1 == Approx( rfront1.capture().value ) );
      CHECK( 0 == Approx( rfront1.fission().value ) );
      CHECK( 0 == Approx( rfront1.competition().value ) );

      auto rback1 = table1.resonances().back();
      CHECK( 6499. == Approx( rback1.energy().value ) );
      CHECK( 6.400000e-2 == Approx( rback1.elastic().value ) );
      CHECK( 1.430000e-1 == Approx( rback1.capture().value ) );
      CHECK( 0 == Approx( rback1.fission().value ) );
      CHECK( 0 == Approx( rback1.competition().value ) );

      // check the minimal energy grid
      auto grid1 = spingroup1.grid();
      CHECK( 471 == grid1.size() ); // ( 158 resonances - 1 negative ) * 3

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 2 -  lJ = 1,0-
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup2 = spingroups[2];

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 2, elastic channel
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel20 = spingroup2.incidentChannel();
      CHECK( "n,Ag107->n,Ag107" == channel20.reactionID().symbol() );

      // incident particle pair
      const auto incident20 = channel20.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident20.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident20.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident20.particle().spin() ) );
      CHECK( +1 == incident20.particle().parity() );
      CHECK( 105.987 * 1.008664 == Approx( incident20.residual().mass().value ) );
      CHECK( 47.0 * 1.602e-19 == Approx( incident20.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident20.residual().spin() ) );
      CHECK( +1 == incident20.residual().parity() );
      CHECK( "n,Ag107" == incident20.pairID().symbol() );

      // particle pair
      const auto pair20 = channel20.particlePair();
      CHECK( 1.008664 == Approx( pair20.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair20.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair20.particle().spin() ) );
      CHECK( +1 == pair20.particle().parity() );
      CHECK( 105.987 * 1.008664 == Approx( pair20.residual().mass().value ) );
      CHECK( 47.0 * 1.602e-19 == Approx( pair20.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair20.residual().spin() ) );
      CHECK( +1 == pair20.residual().parity() );
      CHECK( "n,Ag107" == pair20.pairID().symbol() );

      // quantum numbers
      const auto numbers20 = channel20.quantumNumbers();
      CHECK( 1 == numbers20.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers20.spin() );
      CHECK( 0.0 == numbers20.totalAngularMomentum() );
      CHECK( -1 == numbers20.parity() );
      CHECK( "{1,1,0-}" == numbers20.toString() );

      // radii - NRO=0 NAPS=0
      const auto radii20 = channel20.radii();
      CHECK( a == Approx( radii20.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( a == Approx( radii20.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .66 == Approx( radii20.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel20.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel20.Q().value ) );

      // resonance table
      auto table2 = spingroup2.resonanceTable();

      CHECK( 28 == table2.numberResonances() );
      CHECK( 28 == table2.resonances().size() );
      CHECK( 28 == table2.energies().size() );
      CHECK( 42.81 == Approx( table2.energies().front().value ) );
      CHECK( 6256. == Approx( table2.energies().back().value ) );

      auto rfront2 = table2.resonances().front();
      CHECK( 42.81 == Approx( rfront2.energy().value ) );
      CHECK( 1.960000e-5  == Approx( rfront2.elastic().value ) );
      CHECK( 0.143 == Approx( rfront2.capture().value ) );
      CHECK( 0 == Approx( rfront2.fission().value ) );
      CHECK( 0 == Approx( rfront2.competition().value ) );

      auto rback2 = table2.resonances().back();
      CHECK( 6256. == Approx( rback2.energy().value ) );
      CHECK( 0.142 == Approx( rback2.elastic().value ) );
      CHECK( 0.143 == Approx( rback2.capture().value ) );
      CHECK( 0 == Approx( rback2.fission().value ) );
      CHECK( 0 == Approx( rback2.competition().value ) );

      // check the minimal energy grid
      auto grid2 = spingroup2.grid();
      CHECK( 84 == grid2.size() ); // ( 27 resonances ) * 3

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 3 -  lJ = 1,1-
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup3 = spingroups[3];

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 3, elastic channel
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel30 = spingroup3.incidentChannel();
      CHECK( "n,Ag107->n,Ag107" == channel30.reactionID().symbol() );

      // incident particle pair
      const auto incident30 = channel30.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident30.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident30.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident30.particle().spin() ) );
      CHECK( +1 == incident30.particle().parity() );
      CHECK( 105.987 * 1.008664 == Approx( incident30.residual().mass().value ) );
      CHECK( 47.0 * 1.602e-19 == Approx( incident30.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident30.residual().spin() ) );
      CHECK( +1 == incident30.residual().parity() );
      CHECK( "n,Ag107" == incident30.pairID().symbol() );

      // particle pair
      const auto pair30 = channel30.particlePair();
      CHECK( 1.008664 == Approx( pair30.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair30.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair30.particle().spin() ) );
      CHECK( +1 == pair30.particle().parity() );
      CHECK( 105.987 * 1.008664 == Approx( pair30.residual().mass().value ) );
      CHECK( 47.0 * 1.602e-19 == Approx( pair30.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair30.residual().spin() ) );
      CHECK( +1 == pair30.residual().parity() );
      CHECK( "n,Ag107" == pair30.pairID().symbol() );

      // quantum numbers
      const auto numbers30 = channel30.quantumNumbers();
      CHECK( 1 == numbers30.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers30.spin() );
      CHECK( 1.0 == numbers30.totalAngularMomentum() );
      CHECK( -1 == numbers30.parity() );
      CHECK( "{1,1,1-}" == numbers30.toString() );

      // radii - NRO=0 NAPS=0
      const auto radii30 = channel30.radii();
      CHECK( a == Approx( radii30.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( a == Approx( radii30.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .66 == Approx( radii30.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel30.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel30.Q().value ) );

      // resonance table
      auto table3 = spingroup3.resonanceTable();

      CHECK( 60 == table3.numberResonances() );
      CHECK( 60 == table3.resonances().size() );
      CHECK( 60 == table3.energies().size() );
      CHECK( 20.3 == Approx( table3.energies().front().value ) );
      CHECK( 6140. == Approx( table3.energies().back().value ) );

      auto rfront3 = table3.resonances().front();
      CHECK( 20.3 == Approx( rfront3.energy().value ) );
      CHECK( 1.600000e-7  == Approx( rfront3.elastic().value ) );
      CHECK( 0.143 == Approx( rfront3.capture().value ) );
      CHECK( 0 == Approx( rfront3.fission().value ) );
      CHECK( 0 == Approx( rfront3.competition().value ) );

      auto rback3 = table3.resonances().back();
      CHECK( 6140. == Approx( rback3.energy().value ) );
      CHECK( 4.533334e-2 == Approx( rback3.elastic().value ) );
      CHECK( 0.143 == Approx( rback3.capture().value ) );
      CHECK( 0 == Approx( rback3.fission().value ) );
      CHECK( 0 == Approx( rback3.competition().value ) );

      // check the minimal energy grid
      auto grid3 = spingroup3.grid();
      CHECK( 180 == grid3.size() ); // ( 60 resonances ) * 3

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 4 -  lJ = 1,2-
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup4 = spingroups[4];

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 4, elastic channel
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel40 = spingroup4.incidentChannel();
      CHECK( "n,Ag107->n,Ag107" == channel40.reactionID().symbol() );

      // incident particle pair
      const auto incident40 = channel40.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident40.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident40.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident40.particle().spin() ) );
      CHECK( +1 == incident40.particle().parity() );
      CHECK( 105.987 * 1.008664 == Approx( incident40.residual().mass().value ) );
      CHECK( 47.0 * 1.602e-19 == Approx( incident40.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident40.residual().spin() ) );
      CHECK( +1 == incident40.residual().parity() );
      CHECK( "n,Ag107" == incident40.pairID().symbol() );

      // particle pair
      const auto pair40 = channel40.particlePair();
      CHECK( 1.008664 == Approx( pair40.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair40.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair40.particle().spin() ) );
      CHECK( +1 == pair40.particle().parity() );
      CHECK( 105.987 * 1.008664 == Approx( pair40.residual().mass().value ) );
      CHECK( 47.0 * 1.602e-19 == Approx( pair40.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair40.residual().spin() ) );
      CHECK( +1 == pair40.residual().parity() );
      CHECK( "n,Ag107" == pair40.pairID().symbol() );

      // quantum numbers
      const auto numbers40 = channel40.quantumNumbers();
      CHECK( 1 == numbers40.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers40.spin() );
      CHECK( 2.0 == numbers40.totalAngularMomentum() );
      CHECK( -1 == numbers40.parity() );
      CHECK( "{1,1,2-}" == numbers40.toString() );

      // radii - NRO=0 NAPS=0
      const auto radii40 = channel40.radii();
      CHECK( a == Approx( radii40.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( a == Approx( radii40.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .66 == Approx( radii40.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel40.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel40.Q().value ) );

      // resonance table
      auto table4 = spingroup4.resonanceTable();

      CHECK( 108 == table4.numberResonances() );
      CHECK( 108 == table4.resonances().size() );
      CHECK( 108 == table4.energies().size() );
      CHECK( 18.9 == Approx( table4.energies().front().value ) );
      CHECK( 6481. == Approx( table4.energies().back().value ) );

      auto rfront4 = table4.resonances().front();
      CHECK( 18.9 == Approx( rfront4.energy().value ) );
      CHECK( 8.800000e-8  == Approx( rfront4.elastic().value ) );
      CHECK( 0.143 == Approx( rfront4.capture().value ) );
      CHECK( 0 == Approx( rfront4.fission().value ) );
      CHECK( 0 == Approx( rfront4.competition().value ) );

      auto rback4 = table4.resonances().back();
      CHECK( 6481. == Approx( rback4.energy().value ) );
      CHECK( 1.280000e-2 == Approx( rback4.elastic().value ) );
      CHECK( 0.143 == Approx( rback4.capture().value ) );
      CHECK( 0 == Approx( rback4.fission().value ) );
      CHECK( 0 == Approx( rback4.competition().value ) );

      // check the minimal energy grid
      auto grid4 = spingroup4.grid();
      CHECK( 324 == grid4.size() ); // ( 108 resonances ) * 3
    } // THEN

    THEN( "cross sections can be reconstructed" ) {

      // values taken from NJOY2016 for ENDF/B-VIII.0 Ag107

      ReactionID elas( "n,Ag107->n,Ag107" );
      ReactionID capt( "n,Ag107->capture" );
      Map< ReactionID, CrossSection > xs;

      xs = resonances( 1e-5 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.3456830882900439 == Approx( xs[ elas ].value ) );
      CHECK( 1897.9525626737627 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.3456562608562566 == Approx( xs[ elas ].value ) );
      CHECK( 600.17761234155591 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.3453899443502610 == Approx( xs[ elas ].value ) );
      CHECK( 189.76856827723222 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-2 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.3427350936775975 == Approx( xs[ elas ].value ) );
      CHECK( 59.933525440428681 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-1 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.3164228743141617 == Approx( xs[ elas ].value ) );
      CHECK( 18.713595595710316 == Approx( xs[ capt ].value ) );

      xs = resonances( 1. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.0731912072392573 == Approx( xs[ elas ].value ) );
      CHECK( 5.2487444766558369 == Approx( xs[ capt ].value ) );

      xs = resonances( 10. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 5.5814855572949771 == Approx( xs[ elas ].value ) );
      CHECK( 1.1411368063239522 == Approx( xs[ capt ].value ) );

      xs = resonances( 100. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 5.0003302056688348 == Approx( xs[ elas ].value ) );
      CHECK( 3.8219152520813690E-002 == Approx( xs[ capt ].value ) );

      xs = resonances( 1000. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 5.2683665794806993 == Approx( xs[ elas ].value ) );
      CHECK( 1.2281057839360888E-002 == Approx( xs[ capt ].value ) );

      xs = resonances( 6499. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 33.065935889321828 == Approx( xs[ elas ].value ) );
      CHECK( 65.401834286193107 == Approx( xs[ capt ].value ) );
    } // THEN
  } // GIVEN

  GIVEN( "valid ENDF data for Tm168" ) {

    std::string string = Tm168();
    auto begin = string.begin();
    auto end = string.end();
    long lineNumber = 1;

    njoy::ENDFtk::HeadRecord head( begin, end, lineNumber );
    njoy::ENDFtk::section::Type< 2, 151 > endf( head, begin, end, lineNumber, 6922 );
    ResonanceRange endfResonanceRange = endf.isotopes().front().resonanceRanges().front();

    auto resonances = fromENDF( endfResonanceRange, neutronMass, elementaryCharge, ParticleID( "n" ), ParticleID( "Tm168" ) );

    double a = 0.123 * std::pow( 166.492 * 1.008664, 1. / 3. ) + 0.08;

    THEN( "the appropriate CompoundSystem is returned" ) {

      CHECK( true == resonances.isResolved() );
      CHECK( false == resonances.isUnresolved() );
      CHECK( 1e-5 == Approx( resonances.lowerEnergy().value ) );
      CHECK( 3.2 == Approx( resonances.upperEnergy().value ) );
      CHECK( false == bool( resonances.interpolation() ) );

      auto compoundsystem = std::get< legacy::resolved::CompoundSystem< MultiLevelBreitWigner > >( resonances.compoundSystem() );

      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // content verification
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      // spin groups
      auto spingroups = compoundsystem.spinGroups();
      CHECK( 1 == spingroups.size() );

      // check the minimal energy grid
      auto grid = compoundsystem.grid();
      CHECK( 6 == grid.size() );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0 -  lJ = 0,5/2+
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup0 = spingroups[0];

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0, elastic channel
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel00 = spingroup0.incidentChannel();
      CHECK( "n,Tm168->n,Tm168" == channel00.reactionID().symbol() );

      // incident particle pair
      const auto incident00 = channel00.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident00.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident00.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident00.particle().spin() ) );
      CHECK( +1 == incident00.particle().parity() );
      CHECK( 166.492 * 1.008664 == Approx( incident00.residual().mass().value ) );
      CHECK( 69.0 * 1.602e-19 == Approx( incident00.residual().charge().value ) );
      CHECK( 3. == Approx( incident00.residual().spin() ) );
      CHECK( +1 == incident00.residual().parity() );
      CHECK( "n,Tm168" == incident00.pairID().symbol() );

      // particle pair
      const auto pair00 = channel00.particlePair();
      CHECK( 1.008664 == Approx( pair00.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair00.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair00.particle().spin() ) );
      CHECK( +1 == pair00.particle().parity() );
      CHECK( 166.492 * 1.008664 == Approx( pair00.residual().mass().value ) );
      CHECK( 69.0 * 1.602e-19 == Approx( pair00.residual().charge().value ) );
      CHECK( 3. == Approx( pair00.residual().spin() ) );
      CHECK( +1 == pair00.residual().parity() );
      CHECK( "n,Tm168" == pair00.pairID().symbol() );

      // quantum numbers
      const auto numbers00 = channel00.quantumNumbers();
      CHECK( 0 == numbers00.orbitalAngularMomentum() );
      CHECK( 2.5 == numbers00.spin() );
      CHECK( 2.5 == numbers00.totalAngularMomentum() );
      CHECK( +1 == numbers00.parity() );
      CHECK( "{0,5/2,5/2+}" == numbers00.toString() );

      // radii - NRO=1 NAPS=0
      const auto radii00 = channel00.radii();
      CHECK( a == Approx( radii00.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( a == Approx( radii00.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( 1.2381 == Approx( radii00.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel00.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel00.Q().value ) );

      // resonance table
      auto table0 = spingroup0.resonanceTable();

      CHECK( 4 == table0.numberResonances() );
      CHECK( 4 == table0.resonances().size() );
      CHECK( 4 == table0.energies().size() );
      CHECK( -2.974700 == Approx( table0.energies().front().value ) );
      CHECK( 3.0253 == Approx( table0.energies().back().value ) );

      auto rfront0 = table0.resonances().front();
      CHECK( -2.974700 == Approx( rfront0.energy().value ) );
      CHECK( 4.616000e-4 == Approx( rfront0.elastic().value ) );
      CHECK( 7.800000e-2 == Approx( rfront0.capture().value ) );
      CHECK( 0 == Approx( rfront0.fission().value ) );
      CHECK( 0 == Approx( rfront0.competition().value ) );

      auto rback0 = table0.resonances().back();
      CHECK( 3.0253 == Approx( rback0.energy().value ) );
      CHECK( 4.616000e-4 == Approx( rback0.elastic().value ) );
      CHECK( 7.800000e-2 == Approx( rback0.capture().value ) );
      CHECK( 0 == Approx( rback0.fission().value ) );
      CHECK( 0 == Approx( rback0.competition().value ) );

      // check the minimal energy grid
      auto grid0 = spingroup0.grid();
      CHECK( 6 == grid0.size() ); // ( 4 resonances - 2 negative ) * 3
    } // THEN

    THEN( "cross sections can be reconstructed" ) {

      // values taken from NJOY2016 for ENDF/B-VII.1 Tm168

      ReactionID elas( "n,Tm168->n,Tm168" );
      ReactionID capt( "n,Tm168->capture" );
      Map< ReactionID, CrossSection > xs;

      xs = resonances( 1e-5 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 19.377220540239161 == Approx( xs[ elas ].value ) );
      CHECK( 6850.2432291685291  == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 19.376943488145745 == Approx( xs[ elas ].value ) );
      CHECK( 2166.2043524575993 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 19.374188263752348 == Approx( xs[ elas ].value ) );
      CHECK( 684.91264751523045 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-2 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 19.346727159961205 == Approx( xs[ elas ].value ) );
      CHECK( 216.32334287851921 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-1 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 19.074173987917121 == Approx( xs[ elas ].value ) );
      CHECK( 69.302277424063405 == Approx( xs[ capt ].value ) );

      xs = resonances( 1. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 29.995256452023909 == Approx( xs[ elas ].value ) );
      CHECK( 4612.1335908893370 == Approx( xs[ capt ].value ) );

      xs = resonances( 1.025300 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 57.823498652539691 == Approx( xs[ elas ].value ) );
      CHECK( 6446.7154333248773 == Approx( xs[ capt ].value ) );

      xs = resonances( 3.025300 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 33.221063754048949 == Approx( xs[ elas ].value ) );
      CHECK( 2185.2800476397319 == Approx( xs[ capt ].value ) );
    } // THEN
  } // GIVEN

  GIVEN( "valid ENDF data for Dy160" ) {

    std::string string = Dy160();
    auto begin = string.begin();
    auto end = string.end();
    long lineNumber = 1;

    njoy::ENDFtk::HeadRecord head( begin, end, lineNumber );
    njoy::ENDFtk::section::Type< 2, 151 > endf( head, begin, end, lineNumber, 6637 );
    ResonanceRange endfResonanceRange = endf.isotopes().front().resonanceRanges().front();

    auto resonances = fromENDF( endfResonanceRange, neutronMass, elementaryCharge, ParticleID( "n" ), ParticleID( "Dy160" ) );

    double a = 0.123 * std::pow( 158.5512 * 1.008664, 1. / 3. ) + 0.08;

    THEN( "the appropriate CompoundSystem is returned" ) {

      CHECK( true == resonances.isResolved() );
      CHECK( false == resonances.isUnresolved() );
      CHECK( 1e-5 == Approx( resonances.lowerEnergy().value ) );
      CHECK( 2008. == Approx( resonances.upperEnergy().value ) );
      CHECK( false == bool( resonances.interpolation() ) );

      auto compoundsystem = std::get< legacy::resolved::CompoundSystem< MultiLevelBreitWigner > >( resonances.compoundSystem() );

      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // content verification
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      // spin groups
      auto spingroups = compoundsystem.spinGroups();
      CHECK( 1 == spingroups.size() );

      // check the minimal energy grid
      auto grid = compoundsystem.grid();
      CHECK( 195 == grid.size() );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0 -  lJ = 0,1/2+
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup0 = spingroups[0];

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0, elastic channel
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel00 = spingroup0.incidentChannel();
      CHECK( "n,Dy160->n,Dy160" == channel00.reactionID().symbol() );

      // incident particle pair
      const auto incident00 = channel00.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident00.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident00.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident00.particle().spin() ) );
      CHECK( +1 == incident00.particle().parity() );
      CHECK( 158.5512 * 1.008664 == Approx( incident00.residual().mass().value ) );
      CHECK( 66.0 * 1.602e-19 == Approx( incident00.residual().charge().value ) );
      CHECK( 0. == Approx( incident00.residual().spin() ) );
      CHECK( +1 == incident00.residual().parity() );
      CHECK( "n,Dy160" == incident00.pairID().symbol() );

      // particle pair
      const auto pair00 = channel00.particlePair();
      CHECK( 1.008664 == Approx( pair00.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair00.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair00.particle().spin() ) );
      CHECK( +1 == pair00.particle().parity() );
      CHECK( 158.5512 * 1.008664 == Approx( pair00.residual().mass().value ) );
      CHECK( 66.0 * 1.602e-19 == Approx( pair00.residual().charge().value ) );
      CHECK( 0. == Approx( pair00.residual().spin() ) );
      CHECK( +1 == pair00.residual().parity() );
      CHECK( "n,Dy160" == pair00.pairID().symbol() );

      // quantum numbers
      const auto numbers00 = channel00.quantumNumbers();
      CHECK( 0 == numbers00.orbitalAngularMomentum() );
      CHECK( 0.5 == numbers00.spin() );
      CHECK( 0.5 == numbers00.totalAngularMomentum() );
      CHECK( +1 == numbers00.parity() );
      CHECK( "{0,1/2,1/2+}" == numbers00.toString() );

      // radii - NRO=1 NAPS=0
      const auto radii00 = channel00.radii();
      CHECK( a == Approx( radii00.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( a == Approx( radii00.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .7456127 == Approx( radii00.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel00.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel00.Q().value ) );

      // resonance table
      auto table0 = spingroup0.resonanceTable();

      CHECK( 66 == table0.numberResonances() );
      CHECK( 66 == table0.resonances().size() );
      CHECK( 66 == table0.energies().size() );
      CHECK( -5.819000e+1 == Approx( table0.energies().front().value ) );
      CHECK( 1994.3 == Approx( table0.energies().back().value ) );

      auto rfront0 = table0.resonances().front();
      CHECK( -5.819000e+1 == Approx( rfront0.energy().value ) );
      CHECK( 5.553000e-1 == Approx( rfront0.elastic().value ) );
      CHECK( 1.058000e-1 == Approx( rfront0.capture().value ) );
      CHECK( 0 == Approx( rfront0.fission().value ) );
      CHECK( 0 == Approx( rfront0.competition().value ) );

      auto rback0 = table0.resonances().back();
      CHECK( 1994.3 == Approx( rback0.energy().value ) );
      CHECK( 2.500000e-1 == Approx( rback0.elastic().value ) );
      CHECK( 1.058000e-1 == Approx( rback0.capture().value ) );
      CHECK( 0 == Approx( rback0.fission().value ) );
      CHECK( 0 == Approx( rback0.competition().value ) );

      // check the minimal energy grid
      auto grid0 = spingroup0.grid();
      CHECK( 195 == grid0.size() ); // ( 66 resonances - 1 negative ) * 3
    } // THEN

    THEN( "cross sections can be reconstructed" ) {

      // values taken from NJOY2016 for ENDF/B-VIII.0 Dy160

      ReactionID elas( "n,Dy160->n,Dy160" );
      ReactionID capt( "n,Dy160->capture" );
      Map< ReactionID, CrossSection > xs;

      xs = resonances( 1e-5 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 5.7145158652626975 == Approx( xs[ elas ].value ) );
      CHECK( 2785.3465249964465  == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 5.7144706453067666 == Approx( xs[ elas ].value ) );
      CHECK( 880.83729758558798  == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 5.7140201084068787 == Approx( xs[ elas ].value ) );
      CHECK( 278.65088951079861  == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-2 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 5.7095119231359446 == Approx( xs[ elas ].value ) );
      CHECK( 88.453629762421997  == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-1 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 5.6635745542418237 == Approx( xs[ elas ].value ) );
      CHECK( 29.112695931538980  == Approx( xs[ capt ].value ) );

      xs = resonances( 1. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 5.0349439486732308 == Approx( xs[ elas ].value ) );
      CHECK( 19.782038357730553  == Approx( xs[ capt ].value ) );

      xs = resonances( 10. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 45.290708367401329 == Approx( xs[ elas ].value ) );
      CHECK( 518.69376992948958  == Approx( xs[ capt ].value ) );

      xs = resonances( 100. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 9.4484330253560902 == Approx( xs[ elas ].value ) );
      CHECK( 0.43448307071166309  == Approx( xs[ capt ].value ) );

      xs = resonances( 1000. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 9.0062825973186911 == Approx( xs[ elas ].value ) );
      CHECK( 6.7745175723848902E-002  == Approx( xs[ capt ].value ) );

      xs = resonances( 1994.3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 655.54851943250458 == Approx( xs[ elas ].value ) );
      CHECK( 276.27064588721004  == Approx( xs[ capt ].value ) );
    } // THEN
  } // GIVEN

  GIVEN( "valid ENDF data for Dy162" ) {

    std::string string = Dy162();
    auto begin = string.begin();
    auto end = string.end();
    long lineNumber = 1;

    njoy::ENDFtk::HeadRecord head( begin, end, lineNumber );
    njoy::ENDFtk::section::Type< 2, 151 > endf( head, begin, end, lineNumber, 6643 );
    ResonanceRange endfResonanceRange = endf.isotopes().front().resonanceRanges().front();

    auto resonances = fromENDF( endfResonanceRange, neutronMass, elementaryCharge, ParticleID( "n" ), ParticleID( "Dy162" ) );

    double a = 0.123 * std::pow( 160.536 * 1.008664, 1. / 3. ) + 0.08;

    THEN( "the appropriate CompoundSystem is returned" ) {

      CHECK( true == resonances.isResolved() );
      CHECK( false == resonances.isUnresolved() );
      CHECK( 1e-5 == Approx( resonances.lowerEnergy().value ) );
      CHECK( 4845. == Approx( resonances.upperEnergy().value ) );
      CHECK( false == bool( resonances.interpolation() ) );

      auto compoundsystem = std::get< legacy::resolved::CompoundSystem< MultiLevelBreitWigner > >( resonances.compoundSystem() );

      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // content verification
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      // spin groups
      auto spingroups = compoundsystem.spinGroups();
      CHECK( 3 == spingroups.size() );

      // check the minimal energy grid
      auto grid = compoundsystem.grid();
      CHECK( 225 == grid.size() );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0 -  lJ = 0,1/2+
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup0 = spingroups[0];

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0, elastic channel
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel00 = spingroup0.incidentChannel();
      CHECK( "n,Dy162->n,Dy162" == channel00.reactionID().symbol() );

      // incident particle pair
      const auto incident00 = channel00.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident00.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident00.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident00.particle().spin() ) );
      CHECK( +1 == incident00.particle().parity() );
      CHECK( 160.536 * 1.008664 == Approx( incident00.residual().mass().value ) );
      CHECK( 66.0 * 1.602e-19 == Approx( incident00.residual().charge().value ) );
      CHECK( 0. == Approx( incident00.residual().spin() ) );
      CHECK( +1 == incident00.residual().parity() );
      CHECK( "n,Dy162" == incident00.pairID().symbol() );

      // particle pair
      const auto pair00 = channel00.particlePair();
      CHECK( 1.008664 == Approx( pair00.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair00.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair00.particle().spin() ) );
      CHECK( +1 == pair00.particle().parity() );
      CHECK( 160.536 * 1.008664 == Approx( pair00.residual().mass().value ) );
      CHECK( 66.0 * 1.602e-19 == Approx( pair00.residual().charge().value ) );
      CHECK( 0. == Approx( pair00.residual().spin() ) );
      CHECK( +1 == pair00.residual().parity() );
      CHECK( "n,Dy162" == pair00.pairID().symbol() );

      // quantum numbers
      const auto numbers00 = channel00.quantumNumbers();
      CHECK( 0 == numbers00.orbitalAngularMomentum() );
      CHECK( 0.5 == numbers00.spin() );
      CHECK( 0.5 == numbers00.totalAngularMomentum() );
      CHECK( +1 == numbers00.parity() );
      CHECK( "{0,1/2,1/2+}" == numbers00.toString() );

      // radii - NRO=1 NAPS=0
      const auto radii00 = channel00.radii();
      CHECK( a == Approx( radii00.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( a == Approx( radii00.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .59 == Approx( radii00.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel00.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel00.Q().value ) );

      // resonance table
      auto table0 = spingroup0.resonanceTable();

      CHECK( 68 == table0.numberResonances() );
      CHECK( 68 == table0.resonances().size() );
      CHECK( 68 == table0.energies().size() );
      CHECK( 5.44 == Approx( table0.energies().front().value ) );
      CHECK( 4844.6 == Approx( table0.energies().back().value ) );

      auto rfront0 = table0.resonances().front();
      CHECK( 5.44 == Approx( rfront0.energy().value ) );
      CHECK( 2.110000e-2  == Approx( rfront0.elastic().value ) );
      CHECK( 1.474000e-1 == Approx( rfront0.capture().value ) );
      CHECK( 0 == Approx( rfront0.fission().value ) );
      CHECK( 0 == Approx( rfront0.competition().value ) );

      auto rback0 = table0.resonances().back();
      CHECK( 4844.6 == Approx( rback0.energy().value ) );
      CHECK( 1.300000e-1 == Approx( rback0.elastic().value ) );
      CHECK( 1.168000e-1 == Approx( rback0.capture().value ) );
      CHECK( 0 == Approx( rback0.fission().value ) );
      CHECK( 0 == Approx( rback0.competition().value ) );

      // check the minimal energy grid
      auto grid0 = spingroup0.grid();
      CHECK( 204 == grid0.size() ); // ( 68 resonances ) * 3

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1 -  lJ = 1,1/2-
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup1 = spingroups[1];

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1, elastic channel
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel10 = spingroup1.incidentChannel();
      CHECK( "n,Dy162->n,Dy162" == channel10.reactionID().symbol() );

      // incident particle pair
      const auto incident10 = channel10.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident10.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident10.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident10.particle().spin() ) );
      CHECK( +1 == incident10.particle().parity() );
      CHECK( 160.536 * 1.008664 == Approx( incident10.residual().mass().value ) );
      CHECK( 66.0 * 1.602e-19 == Approx( incident10.residual().charge().value ) );
      CHECK( 0. == Approx( incident10.residual().spin() ) );
      CHECK( +1 == incident10.residual().parity() );
      CHECK( "n,Dy162" == incident10.pairID().symbol() );

      // particle pair
      const auto pair10 = channel10.particlePair();
      CHECK( 1.008664 == Approx( pair10.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair10.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair10.particle().spin() ) );
      CHECK( +1 == pair10.particle().parity() );
      CHECK( 160.536 * 1.008664 == Approx( pair10.residual().mass().value ) );
      CHECK( 66.0 * 1.602e-19 == Approx( pair10.residual().charge().value ) );
      CHECK( 0. == Approx( pair10.residual().spin() ) );
      CHECK( +1 == pair10.residual().parity() );
      CHECK( "n,Dy162" == pair10.pairID().symbol() );

      // quantum numbers
      const auto numbers10 = channel10.quantumNumbers();
      CHECK( 1 == numbers10.orbitalAngularMomentum() );
      CHECK( 0.5 == numbers10.spin() );
      CHECK( 0.5 == numbers10.totalAngularMomentum() );
      CHECK( -1 == numbers10.parity() );
      CHECK( "{1,1/2,1/2-}" == numbers10.toString() );

      // radii - NRO=1 NAPS=0
      const auto radii10 = channel10.radii();
      CHECK( a == Approx( radii10.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( a == Approx( radii10.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .59 == Approx( radii10.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel10.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel10.Q().value ) );

      // resonance table
      auto table1 = spingroup1.resonanceTable();

      CHECK( 1 == table1.numberResonances() );
      CHECK( 1 == table1.resonances().size() );
      CHECK( 1 == table1.energies().size() );
      CHECK( 3270.9 == Approx( table1.energies().front().value ) );
      CHECK( 3270.9 == Approx( table1.energies().back().value ) );

      auto rfront1 = table1.resonances().front();
      CHECK( 3270.9 == Approx( rfront1.energy().value ) );
      CHECK( 0.018  == Approx( rfront1.elastic().value ) );
      CHECK( 0.1168 == Approx( rfront1.capture().value ) );
      CHECK( 0 == Approx( rfront1.fission().value ) );
      CHECK( 0 == Approx( rfront1.competition().value ) );

      auto rback1 = table1.resonances().back();
      CHECK( 3270.9 == Approx( rback1.energy().value ) );
      CHECK( 0.018 == Approx( rback1.elastic().value ) );
      CHECK( 0.1168 == Approx( rback1.capture().value ) );
      CHECK( 0 == Approx( rback1.fission().value ) );
      CHECK( 0 == Approx( rback1.competition().value ) );

      // check the minimal energy grid
      auto grid1 = spingroup1.grid();
      CHECK( 3 == grid1.size() ); // ( 1 resonances ) * 3

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 2 -  lJ = 1,3/2-
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup2 = spingroups[2];

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 2, elastic channel
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel20 = spingroup2.incidentChannel();
      CHECK( "n,Dy162->n,Dy162" == channel20.reactionID().symbol() );

      // incident particle pair
      const auto incident20 = channel20.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident20.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident20.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident20.particle().spin() ) );
      CHECK( +1 == incident20.particle().parity() );
      CHECK( 160.536 * 1.008664 == Approx( incident20.residual().mass().value ) );
      CHECK( 66.0 * 1.602e-19 == Approx( incident20.residual().charge().value ) );
      CHECK( 0. == Approx( incident20.residual().spin() ) );
      CHECK( +1 == incident20.residual().parity() );
      CHECK( "n,Dy162" == incident20.pairID().symbol() );

      // particle pair
      const auto pair20 = channel20.particlePair();
      CHECK( 1.008664 == Approx( pair20.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair20.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair20.particle().spin() ) );
      CHECK( +1 == pair20.particle().parity() );
      CHECK( 160.536 * 1.008664 == Approx( pair20.residual().mass().value ) );
      CHECK( 66.0 * 1.602e-19 == Approx( pair20.residual().charge().value ) );
      CHECK( 0. == Approx( pair20.residual().spin() ) );
      CHECK( +1 == pair20.residual().parity() );
      CHECK( "n,Dy162" == pair20.pairID().symbol() );

      // quantum numbers
      const auto numbers20 = channel20.quantumNumbers();
      CHECK( 1 == numbers20.orbitalAngularMomentum() );
      CHECK( 0.5 == numbers20.spin() );
      CHECK( 1.5 == numbers20.totalAngularMomentum() );
      CHECK( -1 == numbers20.parity() );
      CHECK( "{1,1/2,3/2-}" == numbers20.toString() );

      // radii - NRO=0 NAPS=0
      const auto radii20 = channel20.radii();
      CHECK( a == Approx( radii20.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( a == Approx( radii20.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .59 == Approx( radii20.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel20.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel20.Q().value ) );

      // resonance table
      auto table2 = spingroup2.resonanceTable();

      CHECK( 6 == table2.numberResonances() );
      CHECK( 6 == table2.resonances().size() );
      CHECK( 6 == table2.energies().size() );
      CHECK( 1387.5 == Approx( table2.energies().front().value ) );
      CHECK( 3966.9 == Approx( table2.energies().back().value ) );

      auto rfront2 = table2.resonances().front();
      CHECK( 1387.5 == Approx( rfront2.energy().value ) );
      CHECK( 0.0024  == Approx( rfront2.elastic().value ) );
      CHECK( 1.168000e-1 == Approx( rfront2.capture().value ) );
      CHECK( 0 == Approx( rfront2.fission().value ) );
      CHECK( 0 == Approx( rfront2.competition().value ) );

      auto rback2 = table2.resonances().back();
      CHECK( 3966.9 == Approx( rback2.energy().value ) );
      CHECK( 1.300000e-2 == Approx( rback2.elastic().value ) );
      CHECK( 1.168000e-1 == Approx( rback2.capture().value ) );
      CHECK( 0 == Approx( rback2.fission().value ) );
      CHECK( 0 == Approx( rback2.competition().value ) );

      // check the minimal energy grid
      auto grid2 = spingroup2.grid();
      CHECK( 18 == grid2.size() ); // ( 6 resonances ) * 3
    } // THEN

    THEN( "cross sections can be reconstructed" ) {

      // values taken from NJOY2016 for ENDF/B-VIII.0 Dy162

      ReactionID elas( "n,Dy162->n,Dy162" );
      ReactionID capt( "n,Dy162->capture" );
      Map< ReactionID, CrossSection > xs;

      xs = resonances( 1e-5 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 0.16389960904605858 == Approx( xs[ elas ].value ) );
      CHECK( 9665.5552533233167 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 0.16391779597504672 == Approx( xs[ elas ].value ) );
      CHECK( 3056.6149238041858 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 0.16410333516673514 == Approx( xs[ elas ].value ) );
      CHECK( 966.89741896145790 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-2 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 0.16597922650957414 == Approx( xs[ elas ].value ) );
      CHECK( 306.74666951513154 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-1 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 0.18572129541119137 == Approx( xs[ elas ].value ) );
      CHECK( 100.21064013699589 == Approx( xs[ capt ].value ) );

      xs = resonances( 1. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 0.51481436908961931 == Approx( xs[ elas ].value ) );
      CHECK( 45.463932876992189 == Approx( xs[ capt ].value ) );

      xs = resonances( 10. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 6.0327932881822504 == Approx( xs[ elas ].value ) );
      CHECK( 13.725382293882630 == Approx( xs[ capt ].value ) );

      xs = resonances( 100. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.4167045411978370 == Approx( xs[ elas ].value ) );
      CHECK( 0.51983222202793633 == Approx( xs[ capt ].value ) );

      xs = resonances( 1000. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 6.1049261794790972 == Approx( xs[ elas ].value ) );
      CHECK( 0.16474127115670925 == Approx( xs[ capt ].value ) );

      xs = resonances( 3966.9 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 16.437892514295868 == Approx( xs[ elas ].value ) );
      CHECK( 119.81501088384495 == Approx( xs[ capt ].value ) );
    } // THEN
  } // GIVEN
} // SCENARIO

std::string mlbwRh105() {

  // Rh105 ENDF/B-VII LRF=1 resonance evaluation

  return
    " 4.510500+4 1.040000+2          0          0          1          04531 2151     \n"
    " 4.510500+4 1.000000+0          0          0          1          04531 2151     \n"
    " 1.000000-5 7.500000+0          1          2          0          04531 2151     \n"
    " 5.000000-1 6.200000-1          0          0          1          04531 2151     \n"
    " 1.040050+2 0.000000+0          0          0         12          24531 2151     \n"
    "-5.000000+0 1.000000+0 1.610000+0 1.450000+0 1.600000-1 0.000000+04531 2151     \n"
    " 5.000000+0 1.000000+0 4.900000-1 3.300000-1 1.600000-1 0.000000+04531 2151     \n"
    "                                                                  4531 2  0     \n";
}

std::string resolvedAg107() {

  // Cf252 ENDF/B-VII.0 LRF=1 resonance evaluation

  return
    " 4.710700+4 1.059870+2          0          0          1          04725 2151     \n"
    " 4.710700+4 1.000000+0          0          0          1          04725 2151     \n"
    " 1.000000-5 6.500000+3          1          2          0          04725 2151     \n"
    " 5.000000-1 6.600000-1          0          0          2          04725 2151     \n"
    " 1.059870+2 0.000000+0          0          0       1224        2044725 2151     \n"
    "-1.260000+1 1.000000+0 1.872700-1 4.427000-2 1.430000-1 0.000000+04725 2151     \n"
    " 1.630000+1 0.000000+0 1.466000-1 1.160000-2 1.350000-1 0.000000+04725 2151     \n"
    " 4.157000+1 1.000000+0 1.423667-1 5.366667-3 1.370000-1 0.000000+04725 2151     \n"
    " 4.490000+1 1.000000+0 1.482000-1 1.200000-3 1.470000-1 0.000000+04725 2151     \n"
    " 5.156000+1 1.000000+0 1.658667-1 2.386667-2 1.420000-1 0.000000+04725 2151     \n"
    " 1.442000+2 0.000000+0 1.470000-1 1.700000-2 1.300000-1 0.000000+04725 2151     \n"
    " 1.620000+2 1.000000+0 1.433733-1 3.733333-4 1.430000-1 0.000000+04725 2151     \n"
    " 1.737000+2 1.000000+0 1.493333-1 7.333333-3 1.420000-1 0.000000+04725 2151     \n"
    " 2.025000+2 1.000000+0 3.333333-1 1.793333-1 1.540000-1 0.000000+04725 2151     \n"
    " 2.513000+2 1.000000+0 1.643333-1 2.133333-2 1.430000-1 0.000000+04725 2151     \n"
    " 2.645000+2 1.000000+0 1.533333-1 3.333333-3 1.500000-1 0.000000+04725 2151     \n"
    " 3.108000+2 1.000000+0 2.436667-1 1.066667-1 1.370000-1 0.000000+04725 2151     \n"
    " 3.612000+2 1.000000+0 1.960000-1 2.100000-2 1.750000-1 0.000000+04725 2151     \n"
    " 4.447000+2 0.000000+0 2.184000-1 7.540000-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.614000+2 1.000000+0 1.584667-1 1.546667-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.667000+2 1.000000+0 2.293333-1 9.133333-2 1.380000-1 0.000000+04725 2151     \n"
    " 4.722000+2 0.000000+0 1.964000-1 5.340000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.154000+2 1.000000+0 1.990000-1 5.600000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.545000+2 0.000000+0 6.130000-1 4.700000-1 1.430000-1 0.000000+04725 2151     \n"
    " 5.769000+2 1.000000+0 1.696000-1 2.660000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.877000+2 1.000000+0 1.806667-1 6.266666-2 1.180000-1 0.000000+04725 2151     \n"
    " 6.259000+2 1.000000+0 1.621333-1 1.913333-2 1.430000-1 0.000000+04725 2151     \n"
    " 6.534000+2 1.000000+0 1.503600-1 7.360000-3 1.430000-1 0.000000+04725 2151     \n"
    " 6.746000+2 1.000000+0 1.733333-1 4.333333-2 1.300000-1 0.000000+04725 2151     \n"
    " 6.959000+2 1.000000+0 1.588000-1 1.580000-2 1.430000-1 0.000000+04725 2151     \n"
    " 7.221000+2 0.000000+0 1.504000-1 7.400000-3 1.430000-1 0.000000+04725 2151     \n"
    " 7.529000+2 1.000000+0 1.960000-1 5.300000-2 1.430000-1 0.000000+04725 2151     \n"
    " 7.802000+2 1.000000+0 1.526000-1 9.600000-3 1.430000-1 0.000000+04725 2151     \n"
    " 8.130000+2 0.000000+0 1.566000-1 1.360000-2 1.430000-1 0.000000+04725 2151     \n"
    " 8.440000+2 0.000000+0 1.560000-1 1.300000-2 1.430000-1 0.000000+04725 2151     \n"
    " 8.489000+2 1.000000+0 1.512667-1 8.266667-3 1.430000-1 0.000000+04725 2151     \n"
    " 8.823000+2 1.000000+0 2.223333-1 7.933334-2 1.430000-1 0.000000+04725 2151     \n"
    " 8.866000+2 0.000000+0 1.750000-1 3.200000-2 1.430000-1 0.000000+04725 2151     \n"
    " 9.094000+2 1.000000+0 1.530667-1 1.006667-2 1.430000-1 0.000000+04725 2151     \n"
    " 9.154000+2 1.000000+0 1.486000-1 5.600000-3 1.430000-1 0.000000+04725 2151     \n"
    " 9.239000+2 0.000000+0 1.634000-1 2.040000-2 1.430000-1 0.000000+04725 2151     \n"
    " 9.341000+2 1.000000+0 3.496667-1 2.066667-1 1.430000-1 0.000000+04725 2151     \n"
    " 1.017000+3 1.000000+0 1.552667-1 1.226667-2 1.430000-1 0.000000+04725 2151     \n"
    " 1.048000+3 0.000000+0 3.770000-1 2.340000-1 1.430000-1 0.000000+04725 2151     \n"
    " 1.068000+3 1.000000+0 1.526667-1 9.666666-3 1.430000-1 0.000000+04725 2151     \n"
    " 1.116000+3 1.000000+0 1.570667-1 1.406667-2 1.430000-1 0.000000+04725 2151     \n"
    " 1.135000+3 0.000000+0 3.490000-1 2.060000-1 1.430000-1 0.000000+04725 2151     \n"
    " 1.149000+3 0.000000+0 2.106000-1 6.760000-2 1.430000-1 0.000000+04725 2151     \n"
    " 1.179000+3 0.000000+0 5.630000-1 4.200000-1 1.430000-1 0.000000+04725 2151     \n"
    " 1.206000+3 1.000000+0 1.636000-1 2.060000-2 1.430000-1 0.000000+04725 2151     \n"
    " 1.228000+3 1.000000+0 1.548667-1 1.186667-2 1.430000-1 0.000000+04725 2151     \n"
    " 1.284000+3 1.000000+0 1.518000-1 8.800000-3 1.430000-1 0.000000+04725 2151     \n"
    " 1.290000+3 1.000000+0 1.523333-1 9.333334-3 1.430000-1 0.000000+04725 2151     \n"
    " 1.373000+3 1.000000+0 1.683333-1 2.533333-2 1.430000-1 0.000000+04725 2151     \n"
    " 1.399000+3 1.000000+0 1.568667-1 1.386667-2 1.430000-1 0.000000+04725 2151     \n"
    " 1.417000+3 1.000000+0 2.736667-1 1.306667-1 1.430000-1 0.000000+04725 2151     \n"
    " 1.439000+3 0.000000+0 2.570000-1 1.140000-1 1.430000-1 0.000000+04725 2151     \n"
    " 1.447000+3 1.000000+0 6.896667-1 5.466667-1 1.430000-1 0.000000+04725 2151     \n"
    " 1.497000+3 1.000000+0 1.583333-1 1.533333-2 1.430000-1 0.000000+04725 2151     \n"
    " 1.520000+3 1.000000+0 1.583333-1 1.533333-2 1.430000-1 0.000000+04725 2151     \n"
    " 1.566000+3 1.000000+0 1.530000-1 1.000000-2 1.430000-1 0.000000+04725 2151     \n"
    " 1.589000+3 1.000000+0 1.943333-1 5.133333-2 1.430000-1 0.000000+04725 2151     \n"
    " 1.628000+3 0.000000+0 4.110000-1 2.680000-1 1.430000-1 0.000000+04725 2151     \n"
    " 1.643000+3 1.000000+0 2.863334-1 1.433333-1 1.430000-1 0.000000+04725 2151     \n"
    " 1.706000+3 1.000000+0 1.723333-1 2.933333-2 1.430000-1 0.000000+04725 2151     \n"
    " 1.720000+3 1.000000+0 2.056667-1 6.266666-2 1.430000-1 0.000000+04725 2151     \n"
    " 1.834000+3 1.000000+0 2.070000-1 6.400000-2 1.430000-1 0.000000+04725 2151     \n"
    " 1.880000+3 0.000000+0 1.850000-1 4.200000-2 1.430000-1 0.000000+04725 2151     \n"
    " 1.906000+3 1.000000+0 3.710000-1 2.280000-1 1.430000-1 0.000000+04725 2151     \n"
    " 1.938000+3 0.000000+0 1.703000+0 1.560000+0 1.430000-1 0.000000+04725 2151     \n"
    " 1.958000+3 1.000000+0 1.896667-1 4.666667-2 1.430000-1 0.000000+04725 2151     \n"
    " 1.988000+3 1.000000+0 2.150000-1 7.200000-2 1.430000-1 0.000000+04725 2151     \n"
    " 2.027000+3 1.000000+0 1.610000-1 1.800000-2 1.430000-1 0.000000+04725 2151     \n"
    " 2.047000+3 1.000000+0 3.590000-1 2.160000-1 1.430000-1 0.000000+04725 2151     \n"
    " 2.074000+3 1.000000+0 1.636667-1 2.066667-2 1.430000-1 0.000000+04725 2151     \n"
    " 2.096000+3 1.000000+0 4.163333-1 2.733333-1 1.430000-1 0.000000+04725 2151     \n"
    " 2.105000+3 0.000000+0 2.010000-1 5.800000-2 1.430000-1 0.000000+04725 2151     \n"
    " 2.113000+3 1.000000+0 4.096667-1 2.666667-1 1.430000-1 0.000000+04725 2151     \n"
    " 2.130000+3 1.000000+0 1.563333-1 1.333333-2 1.430000-1 0.000000+04725 2151     \n"
    " 2.152000+3 1.000000+0 1.843333-1 4.133333-2 1.430000-1 0.000000+04725 2151     \n"
    " 2.177000+3 1.000000+0 4.003333-1 2.573333-1 1.430000-1 0.000000+04725 2151     \n"
    " 2.254000+3 1.000000+0 1.776667-1 3.466667-2 1.430000-1 0.000000+04725 2151     \n"
    " 2.314000+3 0.000000+0 4.690000-1 3.260000-1 1.430000-1 0.000000+04725 2151     \n"
    " 2.395000+3 1.000000+0 7.296667-1 5.866666-1 1.430000-1 0.000000+04725 2151     \n"
    " 2.434000+3 0.000000+0 3.930000-1 2.500000-1 1.430000-1 0.000000+04725 2151     \n"
    " 2.497000+3 1.000000+0 1.910000-1 4.800000-2 1.430000-1 0.000000+04725 2151     \n"
    " 2.510000+3 1.000000+0 2.243333-1 8.133333-2 1.430000-1 0.000000+04725 2151     \n"
    " 2.526000+3 1.000000+0 2.243333-1 8.133333-2 1.430000-1 0.000000+04725 2151     \n"
    " 2.555000+3 1.000000+0 1.750000-1 3.200000-2 1.430000-1 0.000000+04725 2151     \n"
    " 2.611000+3 0.000000+0 2.130000-1 7.000000-2 1.430000-1 0.000000+04725 2151     \n"
    " 2.662000+3 1.000000+0 5.195000-1 4.000000-1 1.195000-1 0.000000+04725 2151     \n"
    " 2.680000+3 1.000000+0 5.363333-1 3.933333-1 1.430000-1 0.000000+04725 2151     \n"
    " 2.714000+3 1.000000+0 1.856667-1 4.266667-2 1.430000-1 0.000000+04725 2151     \n"
    " 2.756000+3 1.000000+0 1.690000-1 2.600000-2 1.430000-1 0.000000+04725 2151     \n"
    " 2.769000+3 1.000000+0 2.296667-1 8.666667-2 1.430000-1 0.000000+04725 2151     \n"
    " 2.784000+3 0.000000+0 2.510000-1 1.080000-1 1.430000-1 0.000000+04725 2151     \n"
    " 2.835000+3 1.000000+0 2.943333-1 1.513333-1 1.430000-1 0.000000+04725 2151     \n"
    " 2.873000+3 1.000000+0 2.103333-1 6.733333-2 1.430000-1 0.000000+04725 2151     \n"
    " 2.905000+3 1.000000+0 1.930000-1 5.000000-2 1.430000-1 0.000000+04725 2151     \n"
    " 2.962000+3 1.000000+0 4.609334-1 3.053333-1 1.556000-1 0.000000+04725 2151     \n"
    " 2.983000+3 1.000000+0 1.650000-1 2.200000-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.020000+3 1.000000+0 2.010000-1 5.800000-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.103000+3 1.000000+0 2.496667-1 1.066667-1 1.430000-1 0.000000+04725 2151     \n"
    " 3.136000+3 1.000000+0 2.363333-1 9.333333-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.157000+3 1.000000+0 1.976667-1 5.466667-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.164000+3 1.000000+0 1.936667-1 5.066666-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.202000+3 1.000000+0 1.736667-1 3.066667-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.269000+3 1.000000+0 1.716667-1 2.866667-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.276000+3 1.000000+0 2.703333-1 1.273333-1 1.430000-1 0.000000+04725 2151     \n"
    " 3.311000+3 1.000000+0 1.736667-1 3.066667-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.328000+3 0.000000+0 2.630000-1 1.200000-1 1.430000-1 0.000000+04725 2151     \n"
    " 3.369000+3 1.000000+0 3.078333-1 1.273333-1 1.805000-1 0.000000+04725 2151     \n"
    " 3.411000+3 1.000000+0 1.803333-1 3.733334-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.428000+3 1.000000+0 1.863333-1 4.333333-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.452000+3 1.000000+0 1.843333-1 4.133333-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.485000+3 1.000000+0 2.063333-1 6.333333-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.491000+3 1.000000+0 3.410000-1 1.980000-1 1.430000-1 0.000000+04725 2151     \n"
    " 3.526000+3 1.000000+0 3.461333-1 1.933333-1 1.528000-1 0.000000+04725 2151     \n"
    " 3.616000+3 0.000000+0 8.290000-1 6.860000-1 1.430000-1 0.000000+04725 2151     \n"
    " 3.635000+3 0.000000+0 2.430000-1 1.000000-1 1.430000-1 0.000000+04725 2151     \n"
    " 3.649000+3 1.000000+0 2.403333-1 9.733333-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.670000+3 1.000000+0 3.803334-1 2.373333-1 1.430000-1 0.000000+04725 2151     \n"
    " 3.696000+3 1.000000+0 3.101333-1 1.653333-1 1.448000-1 0.000000+04725 2151     \n"
    " 3.702000+3 0.000000+0 6.060000-1 4.620000-1 1.440000-1 0.000000+04725 2151     \n"
    " 3.737000+3 1.000000+0 2.978000-1 1.880000-1 1.098000-1 0.000000+04725 2151     \n"
    " 3.785000+3 1.000000+0 1.970000-1 5.400000-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.796000+3 0.000000+0 4.730000-1 3.300000-1 1.430000-1 0.000000+04725 2151     \n"
    " 3.824000+3 1.000000+0 1.763333-1 3.333334-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.838000+3 1.000000+0 3.546667-1 1.706667-1 1.840000-1 0.000000+04725 2151     \n"
    " 3.869000+3 1.000000+0 1.890000-1 4.600000-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.890000+3 0.000000+0 2.430000-1 1.000000-1 1.430000-1 0.000000+04725 2151     \n"
    " 3.901000+3 1.000000+0 1.870000-1 4.400000-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.993000+3 1.000000+0 2.050000-1 6.200000-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.025000+3 1.000000+0 1.790000-1 3.600000-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.043000+3 0.000000+0 2.890000-1 1.460000-1 1.430000-1 0.000000+04725 2151     \n"
    " 4.056000+3 1.000000+0 4.979333-1 3.533333-1 1.446000-1 0.000000+04725 2151     \n"
    " 4.092000+3 1.000000+0 2.123334-1 6.933334-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.157000+3 1.000000+0 2.663333-1 1.233333-1 1.430000-1 0.000000+04725 2151     \n"
    " 4.197000+3 1.000000+0 1.810000-1 3.800000-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.311000+3 0.000000+0 7.530000-1 6.100000-1 1.430000-1 0.000000+04725 2151     \n"
    " 4.330000+3 1.000000+0 2.170000-1 7.400000-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.355000+3 1.000000+0 2.003333-1 5.733334-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.393000+3 1.000000+0 2.456667-1 1.026667-1 1.430000-1 0.000000+04725 2151     \n"
    " 4.404000+3 1.000000+0 2.043333-1 6.133333-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.429000+3 1.000000+0 3.877000-1 2.180000-1 1.697000-1 0.000000+04725 2151     \n"
    " 4.446000+3 0.000000+0 3.270000-1 1.840000-1 1.430000-1 0.000000+04725 2151     \n"
    " 4.493000+3 1.000000+0 2.250000-1 8.200000-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.510000+3 0.000000+0 6.810000-1 5.380000-1 1.430000-1 0.000000+04725 2151     \n"
    " 4.517000+3 1.000000+0 1.843333-1 4.133333-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.551000+3 1.000000+0 3.070000-1 1.640000-1 1.430000-1 0.000000+04725 2151     \n"
    " 4.573000+3 1.000000+0 1.823333-1 3.933333-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.604000+3 1.000000+0 2.030000-1 6.000000-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.637000+3 1.000000+0 1.963333-1 5.333333-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.659000+3 1.000000+0 2.116667-1 6.866667-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.686000+3 0.000000+0 3.090000-1 1.660000-1 1.430000-1 0.000000+04725 2151     \n"
    " 4.694000+3 1.000000+0 2.656667-1 1.226667-1 1.430000-1 0.000000+04725 2151     \n"
    " 4.707000+3 1.000000+0 3.350000-1 1.920000-1 1.430000-1 0.000000+04725 2151     \n"
    " 4.719000+3 0.000000+0 3.050000-1 1.620000-1 1.430000-1 0.000000+04725 2151     \n"
    " 4.756000+3 0.000000+0 3.330000-1 1.900000-1 1.430000-1 0.000000+04725 2151     \n"
    " 4.768000+3 0.000000+0 2.890000-1 1.460000-1 1.430000-1 0.000000+04725 2151     \n"
    " 4.774000+3 1.000000+0 2.363333-1 9.333333-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.788000+3 0.000000+0 3.070000-1 1.640000-1 1.430000-1 0.000000+04725 2151     \n"
    " 4.798000+3 0.000000+0 4.170000-1 2.740000-1 1.430000-1 0.000000+04725 2151     \n"
    " 4.827000+3 1.000000+0 1.896667-1 4.666667-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.840000+3 1.000000+0 2.230000-1 8.000000-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.863000+3 1.000000+0 4.163333-1 2.733333-1 1.430000-1 0.000000+04725 2151     \n"
    " 4.877000+3 1.000000+0 2.323333-1 8.933333-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.948000+3 1.000000+0 2.983333-1 1.553333-1 1.430000-1 0.000000+04725 2151     \n"
    " 4.968000+3 1.000000+0 2.296667-1 8.666667-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.003000+3 1.000000+0 3.987333-1 2.233333-1 1.754000-1 0.000000+04725 2151     \n"
    " 5.030000+3 1.000000+0 3.143333-1 1.713333-1 1.430000-1 0.000000+04725 2151     \n"
    " 5.055000+3 1.000000+0 5.096667-1 3.666667-1 1.430000-1 0.000000+04725 2151     \n"
    " 5.118000+3 1.000000+0 2.496667-1 1.066667-1 1.430000-1 0.000000+04725 2151     \n"
    " 5.183000+3 1.000000+0 1.936667-1 5.066666-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.191000+3 1.000000+0 2.003333-1 5.733334-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.216000+3 1.000000+0 2.856667-1 1.426667-1 1.430000-1 0.000000+04725 2151     \n"
    " 5.245000+3 1.000000+0 2.096667-1 6.666667-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.270000+3 1.000000+0 4.828000-1 3.400000-1 1.428000-1 0.000000+04725 2151     \n"
    " 5.370000+3 1.000000+0 4.582000-1 3.000000-1 1.582000-1 0.000000+04725 2151     \n"
    " 5.392000+3 1.000000+0 2.583333-1 1.153333-1 1.430000-1 0.000000+04725 2151     \n"
    " 5.433000+3 1.000000+0 5.637333-1 4.133333-1 1.504000-1 0.000000+04725 2151     \n"
    " 5.503000+3 1.000000+0 3.250000-1 1.820000-1 1.430000-1 0.000000+04725 2151     \n"
    " 5.610000+3 1.000000+0 5.970667-1 4.466667-1 1.504000-1 0.000000+04725 2151     \n"
    " 5.644000+3 1.000000+0 2.963333-1 1.533333-1 1.430000-1 0.000000+04725 2151     \n"
    " 5.689000+3 0.000000+0 3.090000-1 1.660000-1 1.430000-1 0.000000+04725 2151     \n"
    " 5.714000+3 0.000000+0 4.910000-1 3.480000-1 1.430000-1 0.000000+04725 2151     \n"
    " 5.731000+3 1.000000+0 3.696667-1 2.266667-1 1.430000-1 0.000000+04725 2151     \n"
    " 5.750000+3 1.000000+0 2.390000-1 9.599999-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.782000+3 1.000000+0 2.796667-1 1.366667-1 1.430000-1 0.000000+04725 2151     \n"
    " 5.844000+3 0.000000+0 3.710000-1 2.280000-1 1.430000-1 0.000000+04725 2151     \n"
    " 5.892000+3 1.000000+0 4.496667-1 3.066667-1 1.430000-1 0.000000+04725 2151     \n"
    " 5.912000+3 1.000000+0 2.663333-1 1.233333-1 1.430000-1 0.000000+04725 2151     \n"
    " 5.974000+3 1.000000+0 5.799000-1 4.200000-1 1.599000-1 0.000000+04725 2151     \n"
    " 6.026000+3 1.000000+0 2.496667-1 1.066667-1 1.430000-1 0.000000+04725 2151     \n"
    " 6.042000+3 1.000000+0 4.187000-1 2.400000-1 1.787000-1 0.000000+04725 2151     \n"
    " 6.054000+3 1.000000+0 5.912333-1 4.133333-1 1.779000-1 0.000000+04725 2151     \n"
    " 6.076000+3 1.000000+0 3.163334-1 1.733333-1 1.430000-1 0.000000+04725 2151     \n"
    " 6.090000+3 1.000000+0 2.363333-1 9.333333-2 1.430000-1 0.000000+04725 2151     \n"
    " 6.099000+3 1.000000+0 2.163333-1 7.333333-2 1.430000-1 0.000000+04725 2151     \n"
    " 6.112000+3 1.000000+0 2.496667-1 1.066667-1 1.430000-1 0.000000+04725 2151     \n"
    " 6.128000+3 0.000000+0 8.630000-1 7.200000-1 1.430000-1 0.000000+04725 2151     \n"
    " 6.224000+3 1.000000+0 2.963333-1 1.533333-1 1.430000-1 0.000000+04725 2151     \n"
    " 6.241000+3 1.000000+0 6.496667-1 5.066667-1 1.430000-1 0.000000+04725 2151     \n"
    " 6.314000+3 0.000000+0 6.830000-1 5.400000-1 1.430000-1 0.000000+04725 2151     \n"
    " 6.315000+3 1.000000+0 2.763334-1 1.333333-1 1.430000-1 0.000000+04725 2151     \n"
    " 6.347000+3 1.000000+0 2.896667-1 1.466667-1 1.430000-1 0.000000+04725 2151     \n"
    " 6.374000+3 0.000000+0 6.130000-1 4.700000-1 1.430000-1 0.000000+04725 2151     \n"
    " 6.387000+3 0.000000+0 4.430000-1 3.000000-1 1.430000-1 0.000000+04725 2151     \n"
    " 6.499000+3 1.000000+0 2.070000-1 6.400000-2 1.430000-1 0.000000+04725 2151     \n"
    " 1.059870+2 0.000000+0          1          0       1176        1964725 2151     \n"
    " 1.890000+1 2.000000+0 1.430001-1 8.800000-8 1.430000-1 0.000000+04725 2151     \n"
    " 2.030000+1 1.000000+0 1.430002-1 1.600000-7 1.430000-1 0.000000+04725 2151     \n"
    " 3.584000+1 2.000000+0 1.430003-1 2.720000-7 1.430000-1 0.000000+04725 2151     \n"
    " 4.281000+1 0.000000+0 1.430196-1 1.960000-5 1.430000-1 0.000000+04725 2151     \n"
    " 6.424000+1 1.000000+0 1.430240-1 2.400000-5 1.430000-1 0.000000+04725 2151     \n"
    " 6.474000+1 2.000000+0 1.430104-1 1.040000-5 1.430000-1 0.000000+04725 2151     \n"
    " 7.321000+1 1.000000+0 1.430360-1 3.600000-5 1.430000-1 0.000000+04725 2151     \n"
    " 8.355000+1 2.000000+0 1.430120-1 1.200000-5 1.430000-1 0.000000+04725 2151     \n"
    " 1.012000+2 1.000000+0 1.430053-1 5.333333-6 1.430000-1 0.000000+04725 2151     \n"
    " 1.076000+2 1.000000+0 1.430187-1 1.866667-5 1.430000-1 0.000000+04725 2151     \n"
    " 1.108000+2 2.000000+0 1.430640-1 6.400000-5 1.430000-1 0.000000+04725 2151     \n"
    " 1.251000+2 0.000000+0 1.430400-1 4.000000-5 1.430000-1 0.000000+04725 2151     \n"
    " 1.261000+2 1.000000+0 1.430240-1 2.400000-5 1.430000-1 0.000000+04725 2151     \n"
    " 1.285000+2 2.000000+0 1.430720-1 7.200000-5 1.430000-1 0.000000+04725 2151     \n"
    " 1.367000+2 2.000000+0 1.430224-1 2.240000-5 1.430000-1 0.000000+04725 2151     \n"
    " 1.415000+2 2.000000+0 1.430080-1 8.000000-6 1.430000-1 0.000000+04725 2151     \n"
    " 1.548000+2 1.000000+0 1.430333-1 3.333333-5 1.430000-1 0.000000+04725 2151     \n"
    " 1.669000+2 0.000000+0 1.437600-1 7.600000-4 1.430000-1 0.000000+04725 2151     \n"
    " 1.835000+2 1.000000+0 1.431733-1 1.733333-4 1.430000-1 0.000000+04725 2151     \n"
    " 2.010000+2 1.000000+0 1.433600-1 3.600000-4 1.430000-1 0.000000+04725 2151     \n"
    " 2.189000+2 1.000000+0 1.431133-1 1.133333-4 1.430000-1 0.000000+04725 2151     \n"
    " 2.283000+2 2.000000+0 1.430320-1 3.200000-5 1.430000-1 0.000000+04725 2151     \n"
    " 2.310000+2 2.000000+0 1.430416-1 4.160000-5 1.430000-1 0.000000+04725 2151     \n"
    " 2.355000+2 0.000000+0 1.431160-1 1.160000-4 1.430000-1 0.000000+04725 2151     \n"
    " 2.599000+2 1.000000+0 1.114333-1 3.333334-4 1.111000-1 0.000000+04725 2151     \n"
    " 2.699000+2 1.000000+0 1.432667-1 2.666667-4 1.430000-1 0.000000+04725 2151     \n"
    " 3.282000+2 2.000000+0 1.434800-1 4.800000-4 1.430000-1 0.000000+04725 2151     \n"
    " 3.468000+2 1.000000+0 1.435333-1 5.333333-4 1.430000-1 0.000000+04725 2151     \n"
    " 3.597000+2 1.000000+0 1.433467-1 3.466667-4 1.430000-1 0.000000+04725 2151     \n"
    " 3.725000+2 0.000000+0 1.437600-1 7.600000-4 1.430000-1 0.000000+04725 2151     \n"
    " 3.818000+2 1.000000+0 1.433867-1 3.866667-4 1.430000-1 0.000000+04725 2151     \n"
    " 3.849000+2 2.000000+0 1.430800-1 8.000000-5 1.430000-1 0.000000+04725 2151     \n"
    " 4.039000+2 1.000000+0 1.434000-1 4.000000-4 1.430000-1 0.000000+04725 2151     \n"
    " 4.100000+2 2.000000+0 1.432880-1 2.880000-4 1.430000-1 0.000000+04725 2151     \n"
    " 4.225000+2 0.000000+0 1.437200-1 7.200000-4 1.430000-1 0.000000+04725 2151     \n"
    " 4.793000+2 2.000000+0 1.433640-1 3.640000-4 1.430000-1 0.000000+04725 2151     \n"
    " 4.949000+2 2.000000+0 1.433840-1 3.840000-4 1.430000-1 0.000000+04725 2151     \n"
    " 5.038000+2 1.000000+0 1.433600-1 3.600000-4 1.430000-1 0.000000+04725 2151     \n"
    " 5.217000+2 1.000000+0 1.439133-1 9.133333-4 1.430000-1 0.000000+04725 2151     \n"
    " 5.244000+2 2.000000+0 1.438080-1 8.080000-4 1.430000-1 0.000000+04725 2151     \n"
    " 5.312000+2 1.000000+0 1.435000-1 5.000000-4 1.430000-1 0.000000+04725 2151     \n"
    " 5.924000+2 2.000000+0 1.435200-1 5.200000-4 1.430000-1 0.000000+04725 2151     \n"
    " 6.008000+2 2.000000+0 1.433640-1 3.640000-4 1.430000-1 0.000000+04725 2151     \n"
    " 6.082000+2 2.000000+0 1.455600-1 2.560000-3 1.430000-1 0.000000+04725 2151     \n"
    " 6.341000+2 2.000000+0 1.431680-1 1.680000-4 1.430000-1 0.000000+04725 2151     \n"
    " 6.479000+2 1.000000+0 1.445067-1 1.506667-3 1.430000-1 0.000000+04725 2151     \n"
    " 6.613000+2 1.000000+0 1.432000-1 2.000000-4 1.430000-1 0.000000+04725 2151     \n"
    " 7.031000+2 2.000000+0 1.453200-1 2.320000-3 1.430000-1 0.000000+04725 2151     \n"
    " 7.379000+2 1.000000+0 1.436067-1 6.066667-4 1.430000-1 0.000000+04725 2151     \n"
    " 1.010000+3 2.000000+0 1.434000-1 4.000000-4 1.430000-1 0.000000+04725 2151     \n"
    " 1.028000+3 2.000000+0 1.450000-1 2.000000-3 1.430000-1 0.000000+04725 2151     \n"
    " 1.199000+3 1.000000+0 1.482667-1 5.266666-3 1.430000-1 0.000000+04725 2151     \n"
    " 1.220000+3 2.000000+0 1.451200-1 2.120000-3 1.430000-1 0.000000+04725 2151     \n"
    " 1.238000+3 1.000000+0 1.482000-1 5.200000-3 1.430000-1 0.000000+04725 2151     \n"
    " 1.342000+3 2.000000+0 1.456000-1 2.600000-3 1.430000-1 0.000000+04725 2151     \n"
    " 1.346000+3 2.000000+0 1.463600-1 3.360000-3 1.430000-1 0.000000+04725 2151     \n"
    " 1.353000+3 2.000000+0 1.451200-1 2.120000-3 1.430000-1 0.000000+04725 2151     \n"
    " 1.461000+3 1.000000+0 1.496667-1 6.666666-3 1.430000-1 0.000000+04725 2151     \n"
    " 2.199000+3 1.000000+0 1.536667-1 1.066667-2 1.430000-1 0.000000+04725 2151     \n"
    " 2.203000+3 2.000000+0 1.470000-1 4.000000-3 1.430000-1 0.000000+04725 2151     \n"
    " 2.457000+3 2.000000+0 1.522000-1 9.200000-3 1.430000-1 0.000000+04725 2151     \n"
    " 2.469000+3 2.000000+0 1.526000-1 9.600000-3 1.430000-1 0.000000+04725 2151     \n"
    " 2.474000+3 2.000000+0 1.510000-1 7.999999-3 1.430000-1 0.000000+04725 2151     \n"
    " 2.570000+3 2.000000+0 1.514000-1 8.400000-3 1.430000-1 0.000000+04725 2151     \n"
    " 2.583000+3 2.000000+0 1.470000-1 4.000000-3 1.430000-1 0.000000+04725 2151     \n"
    " 2.592000+3 2.000000+0 1.486000-1 5.600000-3 1.430000-1 0.000000+04725 2151     \n"
    " 2.633000+3 1.000000+0 1.530000-1 1.000000-2 1.430000-1 0.000000+04725 2151     \n"
    " 2.673000+3 1.000000+0 1.467333-1 3.733333-3 1.430000-1 0.000000+04725 2151     \n"
    " 2.680000+3 1.000000+0 1.450000-1 2.000000-3 1.430000-1 0.000000+04725 2151     \n"
    " 2.699000+3 2.000000+0 1.468400-1 3.840000-3 1.430000-1 0.000000+04725 2151     \n"
    " 2.721000+3 0.000000+0 1.710000-1 2.800000-2 1.430000-1 0.000000+04725 2151     \n"
    " 2.737000+3 2.000000+0 1.530000-1 1.000000-2 1.430000-1 0.000000+04725 2151     \n"
    " 2.821000+3 2.000000+0 1.452000-1 2.200000-3 1.430000-1 0.000000+04725 2151     \n"
    " 2.849000+3 2.000000+0 1.542000-1 1.120000-2 1.430000-1 0.000000+04725 2151     \n"
    " 2.885000+3 2.000000+0 1.454000-1 2.400000-3 1.430000-1 0.000000+04725 2151     \n"
    " 2.915000+3 2.000000+0 1.550000-1 1.200000-2 1.430000-1 0.000000+04725 2151     \n"
    " 2.929000+3 1.000000+0 1.488000-1 5.800000-3 1.430000-1 0.000000+04725 2151     \n"
    " 2.944000+3 2.000000+0 1.494000-1 6.400000-3 1.430000-1 0.000000+04725 2151     \n"
    " 3.050000+3 1.000000+0 1.503333-1 7.333333-3 1.430000-1 0.000000+04725 2151     \n"
    " 3.062000+3 2.000000+0 1.470000-1 4.000000-3 1.430000-1 0.000000+04725 2151     \n"
    " 3.071000+3 1.000000+0 1.479333-1 4.933333-3 1.430000-1 0.000000+04725 2151     \n"
    " 3.084000+3 0.000000+0 1.650000-1 2.200000-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.096000+3 1.000000+0 1.478000-1 4.800000-3 1.430000-1 0.000000+04725 2151     \n"
    " 3.118000+3 0.000000+0 1.930000-1 5.000000-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.141000+3 2.000000+0 1.482000-1 5.200000-3 1.430000-1 0.000000+04725 2151     \n"
    " 3.150000+3 2.000000+0 1.506000-1 7.600000-3 1.430000-1 0.000000+04725 2151     \n"
    " 3.191000+3 2.000000+0 1.570000-1 1.400000-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.221000+3 0.000000+0 1.910000-1 4.800000-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.239000+3 2.000000+0 1.463200-1 3.320000-3 1.430000-1 0.000000+04725 2151     \n"
    " 3.252000+3 2.000000+0 1.562000-1 1.320000-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.256000+3 0.000000+0 2.050000-1 6.200000-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.288000+3 0.000000+0 1.558000-1 1.280000-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.317000+3 2.000000+0 1.518000-1 8.800000-3 1.430000-1 0.000000+04725 2151     \n"
    " 3.337000+3 0.000000+0 1.610000-1 1.800000-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.385000+3 1.000000+0 1.530000-1 1.000000-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.474000+3 0.000000+0 1.604000-1 1.740000-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.512000+3 0.000000+0 1.596000-1 1.660000-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.543000+3 2.000000+0 1.522000-1 9.200000-3 1.430000-1 0.000000+04725 2151     \n"
    " 3.551000+3 2.000000+0 1.554000-1 1.240000-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.567000+3 0.000000+0 1.536000-1 1.060000-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.585000+3 2.000000+0 1.570000-1 1.400000-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.595000+3 2.000000+0 1.448000-1 1.800000-3 1.430000-1 0.000000+04725 2151     \n"
    " 3.610000+3 1.000000+0 1.576667-1 1.466667-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.657000+3 1.000000+0 1.503333-1 7.333333-3 1.430000-1 0.000000+04725 2151     \n"
    " 3.711000+3 2.000000+0 1.430000-1 1.20000-28 1.430000-1 0.000000+04725 2151     \n"
    " 3.722000+3 0.000000+0 1.450000-1 2.000000-3 1.430000-1 0.000000+04725 2151     \n"
    " 3.759000+3 2.000000+0 1.474000-1 4.400000-3 1.430000-1 0.000000+04725 2151     \n"
    " 3.911000+3 1.000000+0 1.530000-1 1.000000-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.933000+3 0.000000+0 1.930000-1 5.000000-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.940000+3 0.000000+0 1.650000-1 2.200000-2 1.430000-1 0.000000+04725 2151     \n"
    " 3.971000+3 1.000000+0 1.690000-1 2.600000-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.005000+3 2.000000+0 1.482000-1 5.200000-3 1.430000-1 0.000000+04725 2151     \n"
    " 4.078000+3 2.000000+0 1.610000-1 1.800000-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.107000+3 1.000000+0 1.656667-1 2.266667-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.120000+3 2.000000+0 1.460400-1 3.040000-3 1.430000-1 0.000000+04725 2151     \n"
    " 4.167000+3 2.000000+0 1.486000-1 5.600000-3 1.430000-1 0.000000+04725 2151     \n"
    " 4.180000+3 2.000000+0 1.530000-1 1.000000-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.205000+3 1.000000+0 1.630000-1 2.000000-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.219000+3 1.000000+0 1.447333-1 1.733333-3 1.430000-1 0.000000+04725 2151     \n"
    " 4.226000+3 2.000000+0 1.606000-1 1.760000-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.236000+3 1.000000+0 1.458000-1 2.800000-3 1.430000-1 0.000000+04725 2151     \n"
    " 4.250000+3 2.000000+0 1.466800-1 3.680000-3 1.430000-1 0.000000+04725 2151     \n"
    " 4.264000+3 0.000000+0 1.690000-1 2.600000-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.292000+3 2.000000+0 1.463200-1 3.320000-3 1.430000-1 0.000000+04725 2151     \n"
    " 4.339000+3 2.000000+0 1.622000-1 1.920000-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.423000+3 2.000000+0 1.558000-1 1.280000-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.462000+3 2.000000+0 1.474000-1 4.400000-3 1.430000-1 0.000000+04725 2151     \n"
    " 4.477000+3 1.000000+0 1.630000-1 2.000000-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.532000+3 0.000000+0 2.090000-1 6.600000-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.559000+3 1.000000+0 1.510000-1 8.000000-3 1.430000-1 0.000000+04725 2151     \n"
    " 4.593000+3 1.000000+0 1.578667-1 1.486667-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.630000+3 2.000000+0 1.498000-1 6.800001-3 1.430000-1 0.000000+04725 2151     \n"
    " 4.669000+3 2.000000+0 1.650000-1 2.200000-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.733000+3 2.000000+0 1.442000-1 1.200000-3 1.430000-1 0.000000+04725 2151     \n"
    " 4.740000+3 2.000000+0 1.562000-1 1.320000-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.744000+3 2.000000+0 1.602000-1 1.720000-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.819000+3 0.000000+0 1.730000-1 3.000000-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.908000+3 2.000000+0 1.518000-1 8.800000-3 1.430000-1 0.000000+04725 2151     \n"
    " 4.924000+3 0.000000+0 2.170000-1 7.400000-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.934000+3 2.000000+0 1.574000-1 1.440000-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.941000+3 1.000000+0 1.623333-1 1.933333-2 1.430000-1 0.000000+04725 2151     \n"
    " 4.988000+3 2.000000+0 1.674000-1 2.440000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.020000+3 2.000000+0 1.482000-1 5.200000-3 1.430000-1 0.000000+04725 2151     \n"
    " 5.043000+3 2.000000+0 1.498000-1 6.800001-3 1.430000-1 0.000000+04725 2151     \n"
    " 5.065000+3 2.000000+0 1.514000-1 8.400000-3 1.430000-1 0.000000+04725 2151     \n"
    " 5.074000+3 1.000000+0 1.723333-1 2.933333-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.108000+3 1.000000+0 1.696667-1 2.666667-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.150000+3 0.000000+0 1.850000-1 4.200000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.161000+3 2.000000+0 1.464000-1 3.400000-3 1.430000-1 0.000000+04725 2151     \n"
    " 5.203000+3 2.000000+0 1.554000-1 1.240000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.229000+3 2.000000+0 1.478000-1 4.800000-3 1.430000-1 0.000000+04725 2151     \n"
    " 5.260000+3 2.000000+0 1.550000-1 1.200000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.292000+3 2.000000+0 1.534000-1 1.040000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.314000+3 2.000000+0 1.494000-1 6.400000-3 1.430000-1 0.000000+04725 2151     \n"
    " 5.327000+3 2.000000+0 1.686000-1 2.560000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.344000+3 2.000000+0 1.706000-1 2.760000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.359000+3 0.000000+0 1.710000-1 2.800000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.417000+3 1.000000+0 1.696667-1 2.666667-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.426000+3 2.000000+0 1.546000-1 1.160000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.474000+3 2.000000+0 1.554000-1 1.240000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.486000+3 2.000000+0 1.542000-1 1.120000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.495000+3 2.000000+0 1.570000-1 1.400000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.513000+3 2.000000+0 1.566000-1 1.360000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.524000+3 1.000000+0 1.876667-1 4.466667-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.556000+3 2.000000+0 1.630000-1 2.000000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.572000+3 1.000000+0 1.543333-1 1.133333-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.584000+3 2.000000+0 1.618000-1 1.880000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.598000+3 1.000000+0 1.616667-1 1.866667-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.632000+3 2.000000+0 1.506000-1 7.600000-3 1.430000-1 0.000000+04725 2151     \n"
    " 5.657000+3 2.000000+0 1.614000-1 1.840000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.669000+3 1.000000+0 1.716667-1 2.866667-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.680000+3 2.000000+0 1.666000-1 2.360000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.701000+3 2.000000+0 1.558000-1 1.280000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.794000+3 1.000000+0 1.603333-1 1.733333-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.810000+3 1.000000+0 1.510000-1 8.000000-3 1.430000-1 0.000000+04725 2151     \n"
    " 5.833000+3 2.000000+0 1.718000-1 2.880000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.865000+3 1.000000+0 1.530000-1 1.000000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.879000+3 1.000000+0 1.730000-1 3.000000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.933000+3 1.000000+0 1.643333-1 2.133333-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.948000+3 0.000000+0 2.350000-1 9.200000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.953000+3 2.000000+0 1.586000-1 1.560000-2 1.430000-1 0.000000+04725 2151     \n"
    " 5.996000+3 1.000000+0 1.570000-1 1.400000-2 1.430000-1 0.000000+04725 2151     \n"
    " 6.014000+3 2.000000+0 1.638000-1 2.080000-2 1.430000-1 0.000000+04725 2151     \n"
    " 6.140000+3 1.000000+0 1.883333-1 4.533334-2 1.430000-1 0.000000+04725 2151     \n"
    " 6.164000+3 2.000000+0 1.502000-1 7.200000-3 1.430000-1 0.000000+04725 2151     \n"
    " 6.185000+3 2.000000+0 1.706000-1 2.760000-2 1.430000-1 0.000000+04725 2151     \n"
    " 6.206000+3 0.000000+0 3.150000-1 1.720000-1 1.430000-1 0.000000+04725 2151     \n"
    " 6.256000+3 0.000000+0 2.850000-1 1.420000-1 1.430000-1 0.000000+04725 2151     \n"
    " 6.271000+3 2.000000+0 1.622000-1 1.920000-2 1.430000-1 0.000000+04725 2151     \n"
    " 6.334000+3 2.000000+0 1.658000-1 2.280000-2 1.430000-1 0.000000+04725 2151     \n"
    " 6.418000+3 2.000000+0 1.554000-1 1.240000-2 1.430000-1 0.000000+04725 2151     \n"
    " 6.435000+3 2.000000+0 1.469200-1 3.920000-3 1.430000-1 0.000000+04725 2151     \n"
    " 6.448000+3 2.000000+0 1.450000-1 2.000000-3 1.430000-1 0.000000+04725 2151     \n"
    " 6.461000+3 2.000000+0 1.650000-1 2.200000-2 1.430000-1 0.000000+04725 2151     \n"
    " 6.473000+3 2.000000+0 1.714000-1 2.840000-2 1.430000-1 0.000000+04725 2151     \n"
    " 6.481000+3 2.000000+0 1.558000-1 1.280000-2 1.430000-1 0.000000+04725 2151     \n"
    " 0.000000+0 0.000000+0          0          0          0          04725 2  0     \n";
}

std::string Tm168() {

  // Tm168 ENDF/B-VII.1 LRF=1 resonance evaluation with energy dependent radius

  return
    " 6.916800+4 1.664920+2          0          0          1          06922 2151     \n"
    " 6.916800+4 1.000000+0          0          0          1          06922 2151     \n"
    " 1.000000-5 3.200000+0          1          2          1          06922 2151     \n"
    " 0.000000+0 0.000000+0          0          0          1         506922 2151     \n"
    "         50          2                                            6922 2151     \n"
    " 1.000000-5 1.238100+0 4.000000+1 1.188400+0 5.000000+1 1.153200+06922 2151     \n"
    " 6.000000+1 1.126500+0 7.000000+1 1.105300+0 8.000000+1 1.087800+06922 2151     \n"
    " 9.000000+1 1.073100+0 1.000000+2 1.060500+0 2.000000+2 9.888000-16922 2151     \n"
    " 3.000000+2 9.547000-1 4.000000+2 9.334000-1 5.000000+2 9.184000-16922 2151     \n"
    " 6.000000+2 9.069000-1 7.000000+2 8.978000-1 8.000000+2 8.903000-16922 2151     \n"
    " 9.000000+2 8.839000-1 1.000000+3 8.783000-1 2.000000+3 8.456000-16922 2151     \n"
    " 3.000000+3 8.286000-1 4.000000+3 8.170000-1 5.000000+3 8.081000-16922 2151     \n"
    " 6.000000+3 8.008000-1 7.000000+3 7.946000-1 8.000000+3 7.892000-16922 2151     \n"
    " 9.000000+3 7.844000-1 1.000000+4 7.800000-1 1.200000+4 7.721000-16922 2151     \n"
    " 1.400000+4 7.653000-1 1.600000+4 7.592000-1 1.800000+4 7.536000-16922 2151     \n"
    " 2.000000+4 7.484000-1 2.500000+4 7.369000-1 3.000000+4 7.269000-16922 2151     \n"
    " 3.500000+4 7.180000-1 4.000000+4 7.098000-1 4.200000+4 7.067000-16922 2151     \n"
    " 4.400000+4 7.038000-1 4.600000+4 7.009000-1 4.800000+4 6.980000-16922 2151     \n"
    " 5.000000+4 6.953000-1 5.500000+4 6.888000-1 6.000000+4 6.826000-16922 2151     \n"
    " 6.500000+4 6.767000-1 7.000000+4 6.712000-1 7.500000+4 6.659000-16922 2151     \n"
    " 8.000000+4 6.608000-1 8.500000+4 6.560000-1 9.000000+4 6.513000-16922 2151     \n"
    " 9.500000+4 6.469000-1 2.000000+5 5.803000-1                      6922 2151     \n"
    " 3.000000+0 0.000000+0          0          0          1          06922 2151     \n"
    " 1.664920+2 0.000000+0          0          0         24          46922 2151     \n"
    "-2.974700+0 2.500000+0 7.846160-2 4.616000-4 7.800000-2 0.000000+06922 2151     \n"
    "-9.747000-1 2.500000+0 7.846160-2 4.616000-4 7.800000-2 0.000000+06922 2151     \n"
    " 1.025300+0 2.500000+0 7.846160-2 4.616000-4 7.800000-2 0.000000+06922 2151     \n"
    " 3.025300+0 2.500000+0 7.846160-2 4.616000-4 7.800000-2 0.000000+06922 2151     \n"
    " 0.000000+0 0.000000+0          0          0          0          06922 2  0     \n";
}

std::string Dy160() {

  return
    " 6.616000+4 1.585510+2          0          0          1          06637 2151     \n"
    " 6.616000+4 1.000000+0          0          0          1          06637 2151     \n"
    " 1.000000-5 2.008000+3          1          2          0          06637 2151     \n"
    " 0.000000+0 7.456127-1          0          0          1          06637 2151     \n"
    " 1.585512+2 0.000000+0          0          0        396         666637 2151     \n"
    "-5.819000+1 5.000000-1 6.611000-1 5.553000-1 1.058000-1 0.000000+06637 2151     \n"
    " 1.880000+0 5.000000-1 1.060000-1 2.000000-4 1.058000-1 0.000000+06637 2151     \n"
    " 1.045000+1 5.000000-1 1.165000-1 1.650000-2 1.000000-1 0.000000+06637 2151     \n"
    " 2.047000+1 5.000000-1 1.368000-1 3.100000-2 1.058000-1 0.000000+06637 2151     \n"
    " 3.490000+1 5.000000-1 1.069800-1 1.180000-3 1.058000-1 0.000000+06637 2151     \n"
    " 7.314000+1 5.000000-1 1.123000-1 6.500000-3 1.058000-1 0.000000+06637 2151     \n"
    " 8.560000+1 5.000000-1 2.030000-1 9.100000-2 1.120000-1 0.000000+06637 2151     \n"
    " 1.156000+2 5.000000-1 1.078400-1 2.040000-3 1.058000-1 0.000000+06637 2151     \n"
    " 1.364000+2 5.000000-1 1.220000-1 2.600000-2 9.600000-2 0.000000+06637 2151     \n"
    " 1.557000+2 5.000000-1 1.900000-1 8.500000-2 1.050000-1 0.000000+06637 2151     \n"
    " 1.783000+2 5.000000-1 1.507000-1 3.070000-2 1.200000-1 0.000000+06637 2151     \n"
    " 2.022000+2 5.000000-1 1.373000-1 3.130000-2 1.060000-1 0.000000+06637 2151     \n"
    " 2.463000+2 5.000000-1 1.146000-1 8.800000-3 1.058000-1 0.000000+06637 2151     \n"
    " 2.719000+2 5.000000-1 1.101000-1 4.300000-3 1.058000-1 0.000000+06637 2151     \n"
    " 3.207000+2 5.000000-1 1.189000-1 1.310000-2 1.058000-1 0.000000+06637 2151     \n"
    " 3.412000+2 5.000000-1 4.958000-1 3.900000-1 1.058000-1 0.000000+06637 2151     \n"
    " 3.793000+2 5.000000-1 3.978000-1 2.920000-1 1.058000-1 0.000000+06637 2151     \n"
    " 3.943000+2 5.000000-1 1.658000-1 6.000000-2 1.058000-1 0.000000+06637 2151     \n"
    " 3.988000+2 5.000000-1 3.158000-1 2.100000-1 1.058000-1 0.000000+06637 2151     \n"
    " 4.301000+2 5.000000-1 2.240000-1 1.140000-1 1.100000-1 0.000000+06637 2151     \n"
    " 5.222000+2 5.000000-1 2.400000-1 1.330000-1 1.070000-1 0.000000+06637 2151     \n"
    " 5.774000+2 5.000000-1 1.360000-1 3.600000-2 1.000000-1 0.000000+06637 2151     \n"
    " 5.900000+2 5.000000-1 1.143000-1 8.500000-3 1.058000-1 0.000000+06637 2151     \n"
    " 6.256000+2 5.000000-1 1.138000-1 8.000000-3 1.058000-1 0.000000+06637 2151     \n"
    " 6.536000+2 5.000000-1 1.114000-1 5.600000-3 1.058000-1 0.000000+06637 2151     \n"
    " 6.795000+2 5.000000-1 8.358001-1 7.300000-1 1.058000-1 0.000000+06637 2151     \n"
    " 7.036000+2 5.000000-1 9.808000-1 8.750000-1 1.058000-1 0.000000+06637 2151     \n"
    " 7.296000+2 5.000000-1 2.330000-1 1.130000-1 1.200000-1 0.000000+06637 2151     \n"
    " 7.364000+2 5.000000-1 1.355800+0 1.250000+0 1.058000-1 0.000000+06637 2151     \n"
    " 8.172000+2 5.000000-1 1.378000-1 3.200000-2 1.058000-1 0.000000+06637 2151     \n"
    " 8.452000+2 5.000000-1 2.158000-1 1.100000-1 1.058000-1 0.000000+06637 2151     \n"
    " 8.686000+2 5.000000-1 2.618000-1 1.560000-1 1.058000-1 0.000000+06637 2151     \n"
    " 8.891000+2 5.000000-1 1.237000-1 1.790000-2 1.058000-1 0.000000+06637 2151     \n"
    " 9.720000+2 5.000000-1 1.127000-1 6.900000-3 1.058000-1 0.000000+06637 2151     \n"
    " 9.893000+2 5.000000-1 1.939000-1 8.810000-2 1.058000-1 0.000000+06637 2151     \n"
    " 1.053500+3 5.000000-1 1.188000-1 1.300000-2 1.058000-1 0.000000+06637 2151     \n"
    " 1.085900+3 5.000000-1 7.958000-1 6.900000-1 1.058000-1 0.000000+06637 2151     \n"
    " 1.103000+3 5.000000-1 1.307000-1 2.490000-2 1.058000-1 0.000000+06637 2151     \n"
    " 1.133600+3 5.000000-1 2.708000-1 1.650000-1 1.058000-1 0.000000+06637 2151     \n"
    " 1.139200+3 5.000000-1 1.139000-1 8.100000-3 1.058000-1 0.000000+06637 2151     \n"
    " 1.161500+3 5.000000-1 1.278000-1 2.200000-2 1.058000-1 0.000000+06637 2151     \n"
    " 1.201300+3 5.000000-1 1.368000-1 3.100000-2 1.058000-1 0.000000+06637 2151     \n"
    " 1.244300+3 5.000000-1 1.508000-1 4.500000-2 1.058000-1 0.000000+06637 2151     \n"
    " 1.279400+3 5.000000-1 1.228000-1 1.700000-2 1.058000-1 0.000000+06637 2151     \n"
    " 1.313100+3 5.000000-1 3.858000-1 2.800000-1 1.058000-1 0.000000+06637 2151     \n"
    " 1.367400+3 5.000000-1 1.238000-1 1.800000-2 1.058000-1 0.000000+06637 2151     \n"
    " 1.401400+3 5.000000-1 2.408000-1 1.350000-1 1.058000-1 0.000000+06637 2151     \n"
    " 1.411700+3 5.000000-1 2.148000-1 1.090000-1 1.058000-1 0.000000+06637 2151     \n"
    " 1.437200+3 5.000000-1 2.838000-1 1.780000-1 1.058000-1 0.000000+06637 2151     \n"
    " 1.495300+3 5.000000-1 3.058000-1 2.000000-1 1.058000-1 0.000000+06637 2151     \n"
    " 1.541700+3 5.000000-1 2.358000-1 1.300000-1 1.058000-1 0.000000+06637 2151     \n"
    " 1.592100+3 5.000000-1 4.848000-1 3.790000-1 1.058000-1 0.000000+06637 2151     \n"
    " 1.625300+3 5.000000-1 1.188000-1 1.300000-2 1.058000-1 0.000000+06637 2151     \n"
    " 1.648300+3 5.000000-1 4.148000-1 3.090000-1 1.058000-1 0.000000+06637 2151     \n"
    " 1.710300+3 5.000000-1 4.158000-1 3.100000-1 1.058000-1 0.000000+06637 2151     \n"
    " 1.721000+3 5.000000-1 3.658000-1 2.600000-1 1.058000-1 0.000000+06637 2151     \n"
    " 1.746200+3 5.000000-1 1.478000-1 4.200000-2 1.058000-1 0.000000+06637 2151     \n"
    " 1.763600+3 5.000000-1 1.318000-1 2.600000-2 1.058000-1 0.000000+06637 2151     \n"
    " 1.785600+3 5.000000-1 1.608000-1 5.500000-2 1.058000-1 0.000000+06637 2151     \n"
    " 1.829700+3 5.000000-1 2.558000-1 1.500000-1 1.058000-1 0.000000+06637 2151     \n"
    " 1.861600+3 5.000000-1 1.395800+0 1.290000+0 1.058000-1 0.000000+06637 2151     \n"
    " 1.874900+3 5.000000-1 2.658000-1 1.600000-1 1.058000-1 0.000000+06637 2151     \n"
    " 1.899600+3 5.000000-1 3.858000-1 2.800000-1 1.058000-1 0.000000+06637 2151     \n"
    " 1.937800+3 5.000000-1 7.258000-1 6.200000-1 1.058000-1 0.000000+06637 2151     \n"
    " 1.944300+3 5.000000-1 3.458000-1 2.400000-1 1.058000-1 0.000000+06637 2151     \n"
    " 1.994300+3 5.000000-1 3.558000-1 2.500000-1 1.058000-1 0.000000+06637 2151     \n"
    "                                                                  6637 2  0     \n";
  }

std::string Dy162() {

  return
    " 6.616200+4 1.605360+2          0          0          1          06643 2151     \n"
    " 6.616200+4 1.000000+0          0          0          1          06643 2151     \n"
    " 1.000000-5 4.845000+3          1          2          0          06643 2151     \n"
    " 0.000000+0 5.900000-1          0          0          2          06643 2151     \n"
    " 1.605360+2 0.000000+0          0          0        408         686643 2151     \n"
    " 5.440000+0 5.000000-1 1.685000-1 2.110000-2 1.474000-1 0.000000+06643 2151     \n"
    " 7.110000+1 5.000000-1 5.250000-1 4.000000-1 1.250000-1 0.000000+06643 2151     \n"
    " 1.172200+2 5.000000-1 1.295000-1 9.500000-3 1.200000-1 0.000000+06643 2151     \n"
    " 2.079700+2 5.000000-1 1.290000-1 2.400000-2 1.050000-1 0.000000+06643 2151     \n"
    " 2.232900+2 5.000000-1 1.520000-1 3.700000-2 1.150000-1 0.000000+06643 2151     \n"
    " 2.693600+2 5.000000-1 7.368000-1 6.200000-1 1.168000-1 0.000000+06643 2151     \n"
    " 3.570200+2 5.000000-1 1.313000-1 2.630000-2 1.050000-1 0.000000+06643 2151     \n"
    " 4.127500+2 5.000000-1 2.500000-1 1.500000-1 1.000000-1 0.000000+06643 2151     \n"
    " 4.703200+2 5.000000-1 1.248000-1 8.000000-3 1.168000-1 0.000000+06643 2151     \n"
    " 5.298300+2 5.000000-1 4.000000-1 2.800000-1 1.200000-1 0.000000+06643 2151     \n"
    " 6.328000+2 5.000000-1 1.626800+0 1.510000+0 1.168000-1 0.000000+06643 2151     \n"
    " 6.859600+2 5.000000-1 5.550000-1 4.450000-1 1.100000-1 0.000000+06643 2151     \n"
    " 7.164700+2 5.000000-1 7.200000-1 5.900000-1 1.300000-1 0.000000+06643 2151     \n"
    " 7.664200+2 5.000000-1 9.768000-1 8.600000-1 1.168000-1 0.000000+06643 2151     \n"
    " 8.660000+2 5.000000-1 2.616800+0 2.500000+0 1.168000-1 0.000000+06643 2151     \n"
    " 9.520800+2 5.000000-1 3.040000-1 1.940000-1 1.100000-1 0.000000+06643 2151     \n"
    " 1.005200+3 5.000000-1 1.630000-1 5.100000-2 1.120000-1 0.000000+06643 2151     \n"
    " 1.066500+3 5.000000-1 1.478000-1 3.100000-2 1.168000-1 0.000000+06643 2151     \n"
    " 1.110600+3 5.000000-1 3.680000-1 2.600000-1 1.080000-1 0.000000+06643 2151     \n"
    " 1.261400+3 5.000000-1 9.668000-1 8.500000-1 1.168000-1 0.000000+06643 2151     \n"
    " 1.360200+3 5.000000-1 2.400000-1 1.300000-1 1.100000-1 0.000000+06643 2151     \n"
    " 1.431600+3 5.000000-1 1.616800+0 1.500000+0 1.168000-1 0.000000+06643 2151     \n"
    " 1.526200+3 5.000000-1 3.250000-1 2.200000-1 1.050000-1 0.000000+06643 2151     \n"
    " 1.567300+3 5.000000-1 7.468000-1 6.300000-1 1.168000-1 0.000000+06643 2151     \n"
    " 1.680700+3 5.000000-1 1.298000-1 1.300000-2 1.168000-1 0.000000+06643 2151     \n"
    " 1.724500+3 5.000000-1 2.230000-1 1.080000-1 1.150000-1 0.000000+06643 2151     \n"
    " 1.794900+3 5.000000-1 1.348000-1 1.800000-2 1.168000-1 0.000000+06643 2151     \n"
    " 1.835300+3 5.000000-1 3.290000-1 2.190000-1 1.100000-1 0.000000+06643 2151     \n"
    " 1.935900+3 5.000000-1 1.036800+0 9.200000-1 1.168000-1 0.000000+06643 2151     \n"
    " 1.950200+3 5.000000-1 1.878000-1 7.100000-2 1.168000-1 0.000000+06643 2151     \n"
    " 2.003200+3 5.000000-1 3.868000-1 2.700000-1 1.168000-1 0.000000+06643 2151     \n"
    " 2.081500+3 5.000000-1 2.968000-1 1.800000-1 1.168000-1 0.000000+06643 2151     \n"
    " 2.181400+3 5.000000-1 1.918000-1 7.500000-2 1.168000-1 0.000000+06643 2151     \n"
    " 2.240900+3 5.000000-1 1.066800+0 9.500000-1 1.168000-1 0.000000+06643 2151     \n"
    " 2.275600+3 5.000000-1 9.268000-1 8.100000-1 1.168000-1 0.000000+06643 2151     \n"
    " 2.341500+3 5.000000-1 8.868000-1 7.700000-1 1.168000-1 0.000000+06643 2151     \n"
    " 2.363400+3 5.000000-1 5.368000-1 4.200000-1 1.168000-1 0.000000+06643 2151     \n"
    " 2.492600+3 5.000000-1 5.368000-1 4.200000-1 1.168000-1 0.000000+06643 2151     \n"
    " 2.526800+3 5.000000-1 3.368000-1 2.200000-1 1.168000-1 0.000000+06643 2151     \n"
    " 2.576200+3 5.000000-1 2.868000-1 1.700000-1 1.168000-1 0.000000+06643 2151     \n"
    " 2.673400+3 5.000000-1 8.368000-1 7.200000-1 1.168000-1 0.000000+06643 2151     \n"
    " 2.796800+3 5.000000-1 2.816800+0 2.700000+0 1.168000-1 0.000000+06643 2151     \n"
    " 2.847600+3 5.000000-1 9.168000-1 8.000000-1 1.168000-1 0.000000+06643 2151     \n"
    " 2.885900+3 5.000000-1 1.186800+0 1.070000+0 1.168000-1 0.000000+06643 2151     \n"
    " 2.944200+3 5.000000-1 2.668000-1 1.500000-1 1.168000-1 0.000000+06643 2151     \n"
    " 2.957200+3 5.000000-1 2.668000-1 1.500000-1 1.168000-1 0.000000+06643 2151     \n"
    " 3.080800+3 5.000000-1 1.588000-1 4.200000-2 1.168000-1 0.000000+06643 2151     \n"
    " 3.183400+3 5.000000-1 2.368000-1 1.200000-1 1.168000-1 0.000000+06643 2151     \n"
    " 3.243800+3 5.000000-1 8.568000-1 7.400000-1 1.168000-1 0.000000+06643 2151     \n"
    " 3.327300+3 5.000000-1 3.768000-1 2.600000-1 1.168000-1 0.000000+06643 2151     \n"
    " 3.363900+3 5.000000-1 3.468000-1 2.300000-1 1.168000-1 0.000000+06643 2151     \n"
    " 3.554400+3 5.000000-1 2.516800+0 2.400000+0 1.168000-1 0.000000+06643 2151     \n"
    " 3.624000+3 5.000000-1 1.518000-1 3.500000-2 1.168000-1 0.000000+06643 2151     \n"
    " 3.666000+3 5.000000-1 6.268000-1 5.100000-1 1.168000-1 0.000000+06643 2151     \n"
    " 3.798900+3 5.000000-1 3.016800+0 2.900000+0 1.168000-1 0.000000+06643 2151     \n"
    " 3.894600+3 5.000000-1 1.558000-1 3.900000-2 1.168000-1 0.000000+06643 2151     \n"
    " 4.012700+3 5.000000-1 1.516800+0 1.400000+0 1.168000-1 0.000000+06643 2151     \n"
    " 4.081800+3 5.000000-1 1.516800+0 1.400000+0 1.168000-1 0.000000+06643 2151     \n"
    " 4.096100+3 5.000000-1 5.068000-1 3.900000-1 1.168000-1 0.000000+06643 2151     \n"
    " 4.114900+3 5.000000-1 2.716800+0 2.600000+0 1.168000-1 0.000000+06643 2151     \n"
    " 4.239000+3 5.000000-1 9.668000-1 8.500000-1 1.168000-1 0.000000+06643 2151     \n"
    " 4.273300+3 5.000000-1 6.368000-1 5.200000-1 1.168000-1 0.000000+06643 2151     \n"
    " 4.393500+3 5.000000-1 3.268000-1 2.100000-1 1.168000-1 0.000000+06643 2151     \n"
    " 4.452300+3 5.000000-1 2.216800+0 2.100000+0 1.168000-1 0.000000+06643 2151     \n"
    " 4.554900+3 5.000000-1 1.848000-1 6.800000-2 1.168000-1 0.000000+06643 2151     \n"
    " 4.657000+3 5.000000-1 2.128000-1 9.600000-2 1.168000-1 0.000000+06643 2151     \n"
    " 4.785600+3 5.000000-1 5.968000-1 4.800000-1 1.168000-1 0.000000+06643 2151     \n"
    " 4.844600+3 5.000000-1 2.468000-1 1.300000-1 1.168000-1 0.000000+06643 2151     \n"
    " 1.605360+2 0.000000+0          1          0         42          76643 2151     \n"
    " 1.387500+3 1.500000+0 1.192000-1 2.400000-3 1.168000-1 0.000000+06643 2151     \n"
    " 1.483300+3 1.500000+0 1.203000-1 3.500000-3 1.168000-1 0.000000+06643 2151     \n"
    " 2.047600+3 1.500000+0 1.213000-1 4.500000-3 1.168000-1 0.000000+06643 2151     \n"
    " 3.028900+3 1.500000+0 1.258000-1 9.000000-3 1.168000-1 0.000000+06643 2151     \n"
    " 3.270900+3 5.000000-1 1.348000-1 1.800000-2 1.168000-1 0.000000+06643 2151     \n"
    " 3.513800+3 1.500000+0 1.298000-1 1.300000-2 1.168000-1 0.000000+06643 2151     \n"
    " 3.966900+3 1.500000+0 1.298000-1 1.300000-2 1.168000-1 0.000000+06643 2151     \n"
    "                                                                  6643 2  0     \n";
  }
