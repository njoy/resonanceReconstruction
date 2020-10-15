std::string Rh105();
std::string Cf252();

SCENARIO( "fromENDF - LRF1" ) {

  GIVEN( "valid ENDF data for Rh105" ) {

    std::string string = Rh105();
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

      auto compoundsystem = std::get< legacy::resolved::CompoundSystem< SingleLevelBreitWigner > >( resonances.compoundSystem() );

      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // content verification
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      // spin groups
      auto spingroups = compoundsystem.spinGroups();
      CHECK( 1 == spingroups.size() );

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

      // values taken from NJOY2016 for ENDF/B-VII.1 Rh105 (set to SLBW)

      ReactionID elas( "n,Rh105->n,Rh105" );
      ReactionID capt( "n,Rh105->capture" );
      std::map< ReactionID, CrossSection > xs;

      xs = resonances( 1e-5 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 9075.762 == Approx( xs[ elas ].value ) );
      CHECK( 801565.16324338294  == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 9075.340 == Approx( xs[ elas ].value ) );
      CHECK( 253468.26154893031 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 9071.981 == Approx( xs[ elas ].value ) );
      CHECK( 80132.232264999941 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-2 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 9041.196 == Approx( xs[ elas ].value ) );
      CHECK( 25279.102640176850 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-1 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 8751.658 == Approx( xs[ elas ].value ) );
      CHECK( 7816.7206941266395 == Approx( xs[ capt ].value ) );

      xs = resonances( 1. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 6657.999 == Approx( xs[ elas ].value ) );
      CHECK( 2161.1909561504717 == Approx( xs[ capt ].value ) );

      xs = resonances( 4.755 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 9.335218e+4 == Approx( xs[ elas ].value ) );
      CHECK( 45893.526706445802 == Approx( xs[ capt ].value ) );

      xs = resonances( 5. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 1.828299e+5 == Approx( xs[ elas ].value ) );
      CHECK( 87782.287793509589 == Approx( xs[ capt ].value ) );

      xs = resonances( 5.245 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 9.177689e+4 == Approx( xs[ elas ].value ) );
      CHECK( 42264.081005233311 == Approx( xs[ capt ].value ) );
    } // THEN
  } // GIVEN

  GIVEN( "valid ENDF data for Cf252" ) {

    std::string string = Cf252();
    auto begin = string.begin();
    auto end = string.end();
    long lineNumber = 1;

    njoy::ENDFtk::HeadRecord head( begin, end, lineNumber );
    njoy::ENDFtk::section::Type< 2, 151 > endf( head, begin, end, lineNumber, 9861 );
    ResonanceRange endfResonanceRange = endf.isotopes().front().resonanceRanges().front();

    auto resonances = fromENDF( endfResonanceRange, neutronMass, elementaryCharge, ParticleID( "n" ), ParticleID( "Cf252" ) );

    double a = 0.123 * std::pow( 249.9160 * 1.008664, 1. / 3. ) + 0.08;

    THEN( "the appropriate CompoundSystem is returned" ) {

      CHECK( true == resonances.isResolved() );
      CHECK( false == resonances.isUnresolved() );
      CHECK( 1e-5 == Approx( resonances.lowerEnergy().value ) );
      CHECK( 366.5 == Approx( resonances.upperEnergy().value ) );

      auto compoundsystem = std::get< legacy::resolved::CompoundSystem< SingleLevelBreitWigner > >( resonances.compoundSystem() );

      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // content verification
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      // spin groups
      auto spingroups = compoundsystem.spinGroups();
      CHECK( 1 == spingroups.size() );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0 -  lJ = 0,1+
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup0 = spingroups[0];

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0, elastic channel
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel00 = spingroup0.incidentChannel();
      CHECK( "n,Cf252->n,Cf252" == channel00.reactionID().symbol() );

      // incident particle pair
      const auto incident00 = channel00.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident00.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident00.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident00.particle().spin() ) );
      CHECK( +1 == incident00.particle().parity() );
      CHECK( 249.9160 * 1.008664 == Approx( incident00.residual().mass().value ) );
      CHECK( 98.0 * 1.602e-19 == Approx( incident00.residual().charge().value ) );
      CHECK( 0. == Approx( incident00.residual().spin() ) );
      CHECK( +1 == incident00.residual().parity() );
      CHECK( "n,Cf252" == incident00.pairID().symbol() );

      // particle pair
      const auto pair00 = channel00.particlePair();
      CHECK( 1.008664 == Approx( pair00.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair00.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair00.particle().spin() ) );
      CHECK( +1 == pair00.particle().parity() );
      CHECK( 249.9160 * 1.008664 == Approx( pair00.residual().mass().value ) );
      CHECK( 98.0 * 1.602e-19 == Approx( pair00.residual().charge().value ) );
      CHECK( 0. == Approx( pair00.residual().spin() ) );
      CHECK( +1 == pair00.residual().parity() );
      CHECK( "n,Cf252" == pair00.pairID().symbol() );

      // quantum numbers
      const auto numbers00 = channel00.quantumNumbers();
      CHECK( 0 == numbers00.orbitalAngularMomentum() );
      CHECK( 0.5 == numbers00.spin() );
      CHECK( 0.5 == numbers00.totalAngularMomentum() );
      CHECK( +1 == numbers00.parity() );
      CHECK( "{0,1/2,1/2+}" == numbers00.toString() );

      // radii
      const auto radii00 = channel00.radii();
      CHECK( a == Approx( radii00.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( a == Approx( radii00.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .883096 == Approx( radii00.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel00.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel00.Q().value ) );

      // resonance table
      auto table0 = spingroup0.resonanceTable();

      CHECK( 21 == table0.numberResonances() );
      CHECK( 21 == table0.resonances().size() );
      CHECK( 21 == table0.energies().size() );
      CHECK( -2.399000 == Approx( table0.energies().front().value ) );
      CHECK( 349. == Approx( table0.energies().back().value ) );

      auto rfront0 = table0.resonances().front();
      CHECK( -2.399000 == Approx( rfront0.energy().value ) );
      CHECK( 1.862000e-3 == Approx( rfront0.elastic().value ) );
      CHECK( 2.354000e-2 == Approx( rfront0.capture().value ) );
      CHECK( 3.624000e-2 == Approx( rfront0.fission().value ) );
      CHECK( 0 == Approx( rfront0.competition().value ) );

      auto rback0 = table0.resonances().back();
      CHECK( 349. == Approx( rback0.energy().value ) );
      CHECK( 2.246000e-2 == Approx( rback0.elastic().value ) );
      CHECK( 2.354000e-2 == Approx( rback0.capture().value ) );
      CHECK( 5.953000e-2 == Approx( rback0.fission().value ) );
      CHECK( 0 == Approx( rback0.competition().value ) );

      // check the minimal energy grid
      auto grid0 = spingroup0.grid();

      CHECK( 60 == grid0.size() ); // ( 21 resonances - 1 negative ) * 3
    } // THEN

    THEN( "cross sections can be reconstructed" ) {

      // values taken from NJOY2016 for ENDF/B-VII.0 Cf252

      ReactionID elas( "n,Cf252->n,Cf252" );
      ReactionID fiss( "n,Cf252->fission" );
      ReactionID capt( "n,Cf252->capture" );
      std::map< ReactionID, CrossSection > xs;

      xs = resonances( 1e-5 * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 11.243977400001121 == Approx( xs[ elas ].value ) );
      CHECK( 1650.6912248981080 == Approx( xs[ fiss ].value ) );
      CHECK( 1051.8579209354205 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-4 * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 11.243866300001121 == Approx( xs[ elas ].value ) );
      CHECK( 521.95732192271782 == Approx( xs[ fiss ].value ) );
      CHECK( 332.60254596822426 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-3 * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 11.242759200001121 == Approx( xs[ elas ].value ) );
      CHECK( 164.94027048892630 == Approx( xs[ fiss ].value ) );
      CHECK( 105.10191714896118 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-2 * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 11.231747200001122 == Approx( xs[ elas ].value ) );
      CHECK( 51.790643523388539 == Approx( xs[ fiss ].value ) );
      CHECK( 32.996562401239352 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-1 * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 11.126130000001110 == Approx( xs[ elas ].value ) );
      CHECK( 15.282510117379221 == Approx( xs[ fiss ].value ) );
      CHECK( 9.7214460025030842 == Approx( xs[ capt ].value ) );

      xs = resonances( 1. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 10.379848340001038 == Approx( xs[ elas ].value ) );
      CHECK( 2.7519626584302150 == Approx( xs[ fiss ].value ) );
      CHECK( 1.7169953446951749 == Approx( xs[ capt ].value ) );

      xs = resonances( 10. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 8.3995451000008394 == Approx( xs[ elas ].value ) );
      CHECK( 0.39961771351488529 == Approx( xs[ fiss ].value ) );
      CHECK( 0.17296504320680961 == Approx( xs[ capt ].value ) );

      xs = resonances( 17.06514 * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 179.13129340001791 == Approx( xs[ elas ].value ) );
      CHECK( 2877.6510322349923 == Approx( xs[ fiss ].value ) );
      CHECK( 1137.9167082430763 == Approx( xs[ capt ].value ) );

      xs = resonances( 17.11 * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 497.02052280004966 == Approx( xs[ elas ].value ) );
      CHECK( 5857.8372713158751 == Approx( xs[ fiss ].value ) );
      CHECK( 2316.37433360366289 == Approx( xs[ capt ].value ) );

      xs = resonances( 17.154 * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 322.28217230003219 == Approx( xs[ elas ].value ) );
      CHECK( 2926.3819811416593 == Approx( xs[ fiss ].value ) );
      CHECK( 1157.1863759917796 == Approx( xs[ capt ].value ) );

      xs = resonances( 100. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 8.5078650000008498 == Approx( xs[ elas ].value ) );
      CHECK( 0.28588768344614179 == Approx( xs[ fiss ].value ) );
      CHECK( 0.11311810319326356 == Approx( xs[ capt ].value ) );

      xs = resonances( 349. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 347.55912570003477 == Approx( xs[ elas ].value ) );
      CHECK( 902.96112615519689 == Approx( xs[ fiss ].value ) );
      CHECK( 357.05871152655800 == Approx( xs[ capt ].value ) );
    } // THEN
  } // GIVEN
} // SCENARIO

std::string Rh105() {

  // Rh105 ENDF/B-VII LRF=1 resonance evaluation

  return
    " 4.510500+4 1.040000+2          0          0          1          04531 2151     \n"
    " 4.510500+4 1.000000+0          0          0          1          04531 2151     \n"
    " 1.000000-5 7.500000+0          1          1          0          04531 2151     \n"
    " 5.000000-1 6.200000-1          0          0          1          04531 2151     \n"
    " 1.040050+2 0.000000+0          0          0         12          24531 2151     \n"
    "-5.000000+0 1.000000+0 1.610000+0 1.450000+0 1.600000-1 0.000000+04531 2151     \n"
    " 5.000000+0 1.000000+0 4.900000-1 3.300000-1 1.600000-1 0.000000+04531 2151     \n"
    "                                                                  4531 2  0     \n";
}

std::string Cf252() {

  // Cf252 ENDF/B-VII.0 LRF=1 resonance evaluation

  return
    " 9.825200+4 2.499160+2          0          0          1          09861 2151     \n"
    " 9.825200+4 1.000000+0          0          1          1          09861 2151     \n"
    " 1.000000-5 3.665000+2          1          1          0          09861 2151     \n"
    " 0.000000+0 8.830960-1          0          0          1          09861 2151     \n"
    " 2.499160+2 0.000000+0          0          0        126         219861 2151     \n"
    "-2.399000+0 5.000000-1 6.164000-2 1.862000-3 2.354000-2 3.624000-29861 2151     \n"
    " 1.711000+1 5.000000-1 8.805000-2 4.972000-3 2.354000-2 5.953000-29861 2151     \n"
    " 3.458000+1 5.000000-1 9.014000-2 7.068000-3 2.354000-2 5.953000-29861 2151     \n"
    " 5.205000+1 5.000000-1 9.174000-2 8.671000-3 2.354000-2 5.953000-29861 2151     \n"
    " 6.952000+1 5.000000-1 9.310000-2 1.002000-2 2.354000-2 5.953000-29861 2151     \n"
    " 8.699000+1 5.000000-1 9.428000-2 1.121000-2 2.354000-2 5.953000-29861 2151     \n"
    " 1.045000+2 5.000000-1 9.536000-2 1.228000-2 2.354000-2 5.953000-29861 2151     \n"
    " 1.219000+2 5.000000-1 9.635000-2 1.327000-2 2.354000-2 5.953000-29861 2151     \n"
    " 1.394000+2 5.000000-1 9.726000-2 1.419000-2 2.354000-2 5.953000-29861 2151     \n"
    " 1.569000+2 5.000000-1 9.813000-2 1.505000-2 2.354000-2 5.953000-29861 2151     \n"
    " 1.743000+2 5.000000-1 9.894000-2 1.587000-2 2.354000-2 5.953000-29861 2151     \n"
    " 1.918000+2 5.000000-1 9.972000-2 1.665000-2 2.354000-2 5.953000-29861 2151     \n"
    " 2.093000+2 5.000000-1 1.005000-1 1.739000-2 2.354000-2 5.953000-29861 2151     \n"
    " 2.267000+2 5.000000-1 1.012000-1 1.810000-2 2.354000-2 5.953000-29861 2151     \n"
    " 2.442000+2 5.000000-1 1.019000-1 1.878000-2 2.354000-2 5.953000-29861 2151     \n"
    " 2.617000+2 5.000000-1 1.025000-1 1.944000-2 2.354000-2 5.953000-29861 2151     \n"
    " 2.792000+2 5.000000-1 1.032000-1 2.008000-2 2.354000-2 5.953000-29861 2151     \n"
    " 2.966000+2 5.000000-1 1.038000-1 2.070000-2 2.354000-2 5.953000-29861 2151     \n"
    " 3.141000+2 5.000000-1 1.044000-1 2.130000-2 2.354000-2 5.953000-29861 2151     \n"
    " 3.316000+2 5.000000-1 1.050000-1 2.189000-2 2.354000-2 5.953000-29861 2151     \n"
    " 3.490000+2 5.000000-1 1.055000-1 2.246000-2 2.354000-2 5.953000-29861 2151     \n"
    "                                                                  9861 2  0     \n";
}
