SCENARIO( "SpinGroup" ) {

  GIVEN( "valid data for a SpinGroup" ) {

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
      { Resonance( 2.500000e+3 * electronVolt, 8.917200e+0 * electronVolt,
                   9.508100e-4 * rootElectronVolt, 4.070000e-2 * electronVolt,
                   2.842000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 3.000000e+4 * electronVolt, 8.465900e+0 * electronVolt,
                   9.465000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt,
                   2.693000e+0 * electronVolt, 0. * electronVolt ) },
      { 1, 0, 2, 0 } );

    THEN( "a SpinGroup can be constructed" ) {

      Energy energy = 1e-5 * electronVolts;

      SpinGroup group( std::move( elastic ), std::move( table ) );

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
      CHECK( 8.917200e+0 == Approx( resonance.levelSpacing().value ) );
      CHECK( 9.508100e-4 == Approx( resonance.elastic().value ) );
      CHECK( 4.070000e-2 == Approx( resonance.capture().value ) );
      CHECK( 2.842000e+0 == Approx( resonance.fission().value ) );
      CHECK( 0 == Approx( resonance.competition().value ) );

      resonance = table.resonances()[1];
      CHECK( 30000. == Approx( resonance.energy().value ) );
      CHECK( 8.465900e+0 == Approx( resonance.levelSpacing().value ) );
      CHECK( 9.465000e-4 == Approx( resonance.elastic().value ) );
      CHECK( 4.070000e-2 == Approx( resonance.capture().value ) );
      CHECK( 2.693000e+0 == Approx( resonance.fission().value ) );
      CHECK( 0. == Approx( resonance.competition().value ) );
    } // THEN
  } // GIVEN
} // SCENARIO
