SCENARIO( "SpinGroup" ) {

  GIVEN( "valid data for a SpinGroup" ) {

    // test based on Rh105 ENDF/B-VII.1 resolved resonance evaluation

    // scattering radius from equation D.14
    double a = 0.123 * std::pow( 104.005 * 1.008664, 1. / 3. ) + 0.08;

    // particles
    Particle neutron( ParticleID( "n" ), neutronMass, 0.0 * coulombs, 0.5, +1);
    Particle rh105( ParticleID( "Rh105" ), 104.005 * neutronMass,
                    45.0 * elementary, 0.5, +1);

    // particle pairs
    ParticlePair in( neutron, rh105 );

    // channels
    Channel< Neutron > elastic( in, in, 0. * electronVolt,
                                { 0, 0.5, 1.0, +1 },
                                { a * rootBarn, 0.62 * rootBarn } );

    // resonance table
    ResonanceTable table(
      { Resonance( -5. * electronVolt,
                   1.45 * electronVolt, 0.16 * electronVolt,
                   0. * electronVolt, 0. * electronVolt,
                   elastic.penetrability( -5. * electronVolt ),
                   elastic.penetrability( -5. * electronVolt ),
                   elastic.shiftFactor( -5. * electronVolt ) ),
        Resonance( 5. * electronVolt,
                   0.33 * electronVolt, 0.16 * electronVolt,
                   0. * electronVolt, 0. * electronVolt,
                   elastic.penetrability( 5. * electronVolt ),
                   elastic.penetrability( 5. * electronVolt ),
                   elastic.shiftFactor( 5. * electronVolt ) ) } );

    THEN( "a SpinGroup can be constructed" ) {

      Energy energy = 1e-5 * electronVolts;

      SpinGroup< SingleLevelBreitWigner > group( std::move( elastic ), std::move( table ), 0. * electronVolt );

      auto channel = group.incidentChannel();

      CHECK( "n,Rh105{0,1/2,1+}" == channel.channelID() );
      CHECK( "n,Rh105->n,Rh105" == channel.reactionID().symbol() );

      CHECK( true == channel.isIncidentChannel() );

      auto incident = channel.incidentParticlePair();
      CHECK( 1.00866491582 == Approx( incident.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident.particle().spin() ) );
      CHECK( +1 == incident.particle().parity() );
      CHECK( 104.005 * 1.008664 == Approx( incident.residual().mass().value ) );
      CHECK( 45. * e == Approx( incident.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident.residual().spin() ) );
      CHECK( +1 == incident.residual().parity() );
      CHECK( "n,Rh105" == incident.pairID().symbol() );

      auto pair = channel.particlePair();
      CHECK( 1.00866491582 == Approx( pair.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair.particle().spin() ) );
      CHECK( +1 == pair.particle().parity() );
      CHECK( 104.005 * 1.008664 == Approx( pair.residual().mass().value ) );
      CHECK( 45. * e == Approx( pair.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair.residual().spin() ) );
      CHECK( +1 == pair.residual().parity() );
      CHECK( "n,Rh105" == pair.pairID().symbol() );

      auto numbers = channel.quantumNumbers();
      CHECK( 0 == numbers.orbitalAngularMomentum() );
      CHECK( 0.5 == numbers.spin() );
      CHECK( 1.0 == numbers.totalAngularMomentum() );
      CHECK( +1 == numbers.parity() );

      auto radii = channel.radii();
      CHECK( a == Approx( radii.penetrabilityRadius( energy ).value ) );
      CHECK( a == Approx( radii.shiftFactorRadius( energy ).value ) );
      CHECK( .62 == Approx( radii.phaseShiftRadius( energy ).value ) );

      CHECK( 0.0 == channel.boundaryCondition() );

      CHECK( 0.0 == Approx( channel.Q().value ) );

      // resonance data
      auto table = group.resonanceTable();

      CHECK( 2 == table.numberResonances() );
      CHECK( 2 == table.resonances().size() );
      CHECK( 2 == table.energies().size() );
      CHECK( -5. == Approx( table.energies()[0].value ) );
      CHECK( 5. == Approx( table.energies()[1].value ) );

      auto resonance = table.resonances()[0];
      CHECK( -5. == Approx( resonance.energy().value ) );
      CHECK( 1.45 == Approx( resonance.elastic().value ) );
      CHECK( 0.16 == Approx( resonance.capture().value ) );
      CHECK( 0 == Approx( resonance.fission().value ) );
      CHECK( 0 == Approx( resonance.competition().value ) );

      resonance = table.resonances()[1];
      CHECK( 5. == Approx( resonance.energy().value ) );
      CHECK( 0.33 == Approx( resonance.elastic().value ) );
      CHECK( 0.16 == Approx( resonance.capture().value ) );
      CHECK( 0 == Approx( resonance.fission().value ) );
      CHECK( 0 == Approx( resonance.competition().value ) );

      // check the reaction identifiers
      auto reactions = group.reactionIDs();

      CHECK( 2 == reactions.size() );
      CHECK( "n,Rh105->n,Rh105" == reactions[0].symbol() );
      CHECK( "n,Rh105->capture" == reactions[1].symbol() );

      CHECK( false == group.hasFission() );

      // check the minimal energy grid
      auto grid = group.grid();

      CHECK( 3 == grid.size() );
      CHECK( 4.755 == Approx( grid[0].value ) );
      CHECK( 5. == Approx( grid[1].value ) );
      CHECK( 5.245 == Approx( grid[2].value ) );
    } // THEN
  } // GIVEN
} // SCENARIO
