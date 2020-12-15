SCENARIO( "CompoundSystem" ) {

  GIVEN( "valid data for a CompoundSystem" ) {

    // particles
    Particle neutron( ParticleID( "n" ), neutronMass, 0.0 * coulombs, 0.5, +1);
    Particle pu239( ParticleID( "Pu239" ), 239. * daltons,
                    94.0 * elementary, 1.5, +1);

    // particle pairs
    ParticlePair in( neutron, pu239 );

    // channels
    Channel< Neutron > elastic1( in, in, 0. * electronVolt,
                                 { 0, 0.5, 0.0, +1 },
                                 { 0.5 * rootBarn, 0.946 * rootBarn } );
    Channel< Neutron > elastic2( in, in, 0. * electronVolt,
                                 { 1, 0.5, 1.0, +1 },
                                 { 0.5 * rootBarn, 0.946 * rootBarn } );

    // resonance table
    ResonanceTable table1(
      { Resonance( 2500. * electronVolt, 1. * electronVolt,
                   2. * rootElectronVolt, 3. * electronVolt,
                   4. * electronVolt, 5. * electronVolt ),
        Resonance( 30000. * electronVolt, 5. * electronVolt,
                   4. * rootElectronVolt, 3 * electronVolt,
                   2 * electronVolt, 1. * electronVolt ) },
      { 1, 0, 2, 0 } );
    ResonanceTable table2(
      { Resonance( 2500. * electronVolt, 10. * electronVolt,
                   20. * rootElectronVolt, 30. * electronVolt,
                   40. * electronVolt, 50. * electronVolt ),
        Resonance( 30000. * electronVolt, 50. * electronVolt,
                   40. * rootElectronVolt, 30 * electronVolt,
                   20 * electronVolt, 10. * electronVolt ) },
      { 2, 0, 1, 0 } );

    SpinGroup group1( std::move( elastic1 ), std::move( table1 ) );
    SpinGroup group2( std::move( elastic2 ), std::move( table2 ) );

    THEN( "a CompoundSystem can be constructed" ) {

      Energy energy = 1e-5 * electronVolts;

      CompoundSystem system( { group1, group2 }, 5 );

      CHECK( 2 == system.spinGroups().size() );
      CHECK( 5 == system.interpolation() );

      // group 1 - l,J = 0,0
      auto group = system.spinGroups()[0];

      auto channel = group.incidentChannel();

      CHECK( "n,Pu239{0,1/2,0+}" == channel.channelID() );
      CHECK( "n,Pu239->n,Pu239" == channel.reactionID().symbol() );

      CHECK( true == channel.isIncidentChannel() );

      auto incident = channel.incidentParticlePair();
      CHECK( 1.00866491582 == Approx( incident.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident.particle().spin() ) );
      CHECK( +1 == incident.particle().parity() );
      CHECK( 239. == Approx( incident.residual().mass().value ) );
      CHECK( 94. * e == Approx( incident.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident.residual().spin() ) );
      CHECK( +1 == incident.residual().parity() );
      CHECK( "n,Pu239" == incident.pairID().symbol() );

      auto pair = channel.particlePair();
      CHECK( 1.00866491582 == Approx( pair.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair.particle().spin() ) );
      CHECK( +1 == pair.particle().parity() );
      CHECK( 239. == Approx( pair.residual().mass().value ) );
      CHECK( 94. * e == Approx( pair.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair.residual().spin() ) );
      CHECK( +1 == pair.residual().parity() );
      CHECK( "n,Pu239" == pair.pairID().symbol() );

      auto numbers = channel.quantumNumbers();
      CHECK( 0 == numbers.orbitalAngularMomentum() );
      CHECK( 0.5 == numbers.spin() );
      CHECK( 0.0 == numbers.totalAngularMomentum() );
      CHECK( +1 == numbers.parity() );

      auto radii = channel.radii();
      CHECK( .5 == Approx( radii.penetrabilityRadius( energy ).value ) );
      CHECK( .5 == Approx( radii.shiftFactorRadius( energy ).value ) );
      CHECK( .946 == Approx( radii.phaseShiftRadius( energy ).value ) );

      CHECK( 0.0 == channel.boundaryCondition() );

      CHECK( 0.0 == Approx( channel.Q().value ) );

      // resonance data
      auto table = group.resonanceTable();
      CHECK( 1 == table.degreesOfFreedom().elastic );
      CHECK( 0 == table.degreesOfFreedom().capture );
      CHECK( 2 == table.degreesOfFreedom().fission );
      CHECK( 0 == table.degreesOfFreedom().competition );

      CHECK( 2 == table.numberResonances() );
      CHECK( 2 == table.resonances().size() );
      CHECK( 2 == table.energies().size() );
      CHECK( 2500. == Approx( table.energies()[0].value ) );
      CHECK( 30000. == Approx( table.energies()[1].value ) );

      auto resonance = table.resonances()[0];
      CHECK( 2500. == Approx( resonance.energy().value ) );
      CHECK( 1. == Approx( resonance.levelSpacing().value ) );
      CHECK( 2. == Approx( resonance.elastic().value ) );
      CHECK( 3. == Approx( resonance.capture().value ) );
      CHECK( 4. == Approx( resonance.fission().value ) );
      CHECK( 5. == Approx( resonance.competition().value ) );

      resonance = table.resonances()[1];
      CHECK( 30000. == Approx( resonance.energy().value ) );
      CHECK( 5. == Approx( resonance.levelSpacing().value ) );
      CHECK( 4. == Approx( resonance.elastic().value ) );
      CHECK( 3. == Approx( resonance.capture().value ) );
      CHECK( 2. == Approx( resonance.fission().value ) );
      CHECK( 1. == Approx( resonance.competition().value ) );

      // group 2 - l,J = 1,1
      group = system.spinGroups()[1];

      channel = group.incidentChannel();

      CHECK( "n,Pu239{1,1/2,1+}" == channel.channelID() );
      CHECK( "n,Pu239->n,Pu239" == channel.reactionID().symbol() );

      CHECK( true == channel.isIncidentChannel() );

      incident = channel.incidentParticlePair();
      CHECK( 1.00866491582 == Approx( incident.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident.particle().spin() ) );
      CHECK( +1 == incident.particle().parity() );
      CHECK( 239. == Approx( incident.residual().mass().value ) );
      CHECK( 94. * e == Approx( incident.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident.residual().spin() ) );
      CHECK( +1 == incident.residual().parity() );
      CHECK( "n,Pu239" == incident.pairID().symbol() );

      pair = channel.particlePair();
      CHECK( 1.00866491582 == Approx( pair.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair.particle().spin() ) );
      CHECK( +1 == pair.particle().parity() );
      CHECK( 239. == Approx( pair.residual().mass().value ) );
      CHECK( 94. * e == Approx( pair.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair.residual().spin() ) );
      CHECK( +1 == pair.residual().parity() );
      CHECK( "n,Pu239" == pair.pairID().symbol() );

      numbers = channel.quantumNumbers();
      CHECK( 1 == numbers.orbitalAngularMomentum() );
      CHECK( 0.5 == numbers.spin() );
      CHECK( 1.0 == numbers.totalAngularMomentum() );
      CHECK( +1 == numbers.parity() );

      radii = channel.radii();
      CHECK( .5 == Approx( radii.penetrabilityRadius( energy ).value ) );
      CHECK( .5 == Approx( radii.shiftFactorRadius( energy ).value ) );
      CHECK( .946 == Approx( radii.phaseShiftRadius( energy ).value ) );

      CHECK( 0.0 == channel.boundaryCondition() );

      CHECK( 0.0 == Approx( channel.Q().value ) );

      // resonance data
      table = group.resonanceTable();
      CHECK( 2 == table.degreesOfFreedom().elastic );
      CHECK( 0 == table.degreesOfFreedom().capture );
      CHECK( 1 == table.degreesOfFreedom().fission );
      CHECK( 0 == table.degreesOfFreedom().competition );

      CHECK( 2 == table.numberResonances() );
      CHECK( 2 == table.resonances().size() );
      CHECK( 2 == table.energies().size() );
      CHECK( 2500. == Approx( table.energies()[0].value ) );
      CHECK( 30000. == Approx( table.energies()[1].value ) );

      resonance = table.resonances()[0];
      CHECK( 2500. == Approx( resonance.energy().value ) );
      CHECK( 10. == Approx( resonance.levelSpacing().value ) );
      CHECK( 20. == Approx( resonance.elastic().value ) );
      CHECK( 30. == Approx( resonance.capture().value ) );
      CHECK( 40. == Approx( resonance.fission().value ) );
      CHECK( 50. == Approx( resonance.competition().value ) );

      resonance = table.resonances()[1];
      CHECK( 30000. == Approx( resonance.energy().value ) );
      CHECK( 50. == Approx( resonance.levelSpacing().value ) );
      CHECK( 40. == Approx( resonance.elastic().value ) );
      CHECK( 30. == Approx( resonance.capture().value ) );
      CHECK( 20. == Approx( resonance.fission().value ) );
      CHECK( 10. == Approx( resonance.competition().value ) );

      // check the minimal energy grid
      auto grid = system.grid();

      CHECK( 15 == grid.size() );
      CHECK( 2500. == Approx( grid[0].value ) );
      CHECK( 3000. == Approx( grid[1].value ) );
      CHECK( 3500. == Approx( grid[2].value ) );
      CHECK( 4000. == Approx( grid[3].value ) );
      CHECK( 5000. == Approx( grid[4].value ) );
      CHECK( 6000. == Approx( grid[5].value ) );
      CHECK( 7200. == Approx( grid[6].value ) );
      CHECK( 8500. == Approx( grid[7].value ) );
      CHECK( 10000. == Approx( grid[8].value ) );
      CHECK( 12500. == Approx( grid[9].value ) );
      CHECK( 15000. == Approx( grid[10].value ) );
      CHECK( 17000. == Approx( grid[11].value ) );
      CHECK( 20000. == Approx( grid[12].value ) );
      CHECK( 25000. == Approx( grid[13].value ) );
      CHECK( 30000. == Approx( grid[14].value ) );
    } // THEN
  } // GIVEN

  GIVEN( "data for a CompoundSystem with errors" ) {

    // particles
    Particle neutron( ParticleID( "n" ), neutronMass, 0.0 * coulombs, 0.5, +1);
    Particle pu239( ParticleID( "Pu239" ), 239. * daltons,
                    94.0 * elementary, 1.5, +1);

    // particle pairs
    ParticlePair in( neutron, pu239 );

    // channels
    Channel< Neutron > elastic( in, in, 0. * electronVolt,
                                { 0, 0.5, 0.0, +1 },
                                { 0.5 * rootBarn, 0.946 * rootBarn } );

    // resonance table
    ResonanceTable table(
      { Resonance( 2500. * electronVolt, 1. * electronVolt,
                   2. * rootElectronVolt, 3. * electronVolt,
                   4. * electronVolt, 5. * electronVolt ),
        Resonance( 30000. * electronVolt, 5. * electronVolt,
                   4. * rootElectronVolt, 3 * electronVolt,
                   2 * electronVolt, 1. * electronVolt ) },
      { 1, 0, 2, 0 } );

    SpinGroup group( std::move( elastic ), std::move( table ) );

    THEN( "an exception is thrown at construction when there are no spin "
          "groups" ) {

      CHECK_THROWS( CompoundSystem( std::vector< SpinGroup >{}, 2 ) );
    } // THEN

    THEN( "an exception is thrown at construction when the spin groups "
          "are not unique" ) {

      CHECK_THROWS( CompoundSystem( { group, group }, 2 ) );
    } // THEN
  } // GIVEN
} // SCENARIO
