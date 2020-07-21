SCENARIO( "SpinGroup" ) {

  GIVEN( "valid data for a SpinGroup" ) {

    // particles
    Particle photon( ParticleID( "g" ), 0.0 * daltons, 0.0 * coulombs, 1., +1);
    Particle neutron( ParticleID( "n" ), 1.00866491582 * daltons,
                      0.0 * coulombs, 0.5, +1);
    Particle proton( ParticleID( "p" ), 1.00727647 * daltons,
                     elementary, 0.5, +1);
    Particle cl36( ParticleID( "Cl36" ), 35.968306822 * daltons,
                   17.0 * elementary, 0., +1);
    Particle cl35( ParticleID( "Cl35" ), 34.968852694 * daltons,
                   17.0 * elementary, 1.5, +1);
    Particle cl35_e1( ParticleID( "Cl35_e1" ), 34.968852694 * daltons,
                      17.0 * elementary, 1.5, +1);
    Particle s36( ParticleID( "S36" ), 35.967080699 * daltons,
                  16.0 * elementary, 1.5, +1);

    // particle pairs
    ParticlePair elasticPair( neutron, cl35 );
    ParticlePair inelasticPair( neutron, cl35_e1 );
    ParticlePair capturePair( photon, cl36 );
    ParticlePair protonEmissionPair( proton, s36 );

    // Q values
    QValue elasticQ = 0.0 * electronVolt;
    QValue inelasticQ = -1.219440e+6 * electronVolt;
    QValue captureQ = 0.0 * electronVolt;
    QValue protonEmissionQ = 6.152200e+5 * electronVolt;

    // quantum numbers
    ChannelQuantumNumbers elasticNumbers( 0, 1.0, 1.0, +1 );
    ChannelQuantumNumbers inelasticNumbers( 0, 1.0, 1.0, +1 );
    ChannelQuantumNumbers captureNumbers( 0, 0.0, 1.0, +1 );
    ChannelQuantumNumbers protonEmissionNumbers( 0, 1.0, 1.0, +1 );

    // channel radii
    ChannelRadii elasticRadii( 4.822220e-1 * rootBarn,
                               3.667980e-1 * rootBarn );
    ChannelRadii inelasticRadii( 4.822220e-1 * rootBarn,
                               3.667980e-1 * rootBarn );
    ChannelRadii captureRadii( 0.0 * rootBarn );
    ChannelRadii protonEmissionRadii( 4.822220e-1 * rootBarn,
                                      3.667980e-1 * rootBarn );

    // channels
    Channel< Photon > capture( elasticPair, capturePair, captureQ,
                               captureNumbers, captureRadii );
    Channel< Neutron > elastic( elasticPair, elasticPair, elasticQ,
                                elasticNumbers, elasticRadii );
    Channel< Neutron > inelastic( elasticPair, inelasticPair, inelasticQ,
                                 inelasticNumbers, inelasticRadii );
    Channel< ChargedParticle > protonEmission( elasticPair,
                                               protonEmissionPair,
                                               protonEmissionQ,
                                               protonEmissionNumbers,
                                               protonEmissionRadii );

    // single resonance table
    ResonanceTable single(
      { elastic.channelID(), inelastic.channelID(), protonEmission.channelID() },
      { Resonance( 0.25 * electronVolt,
                   { 1.0 * rootElectronVolt, 2.0 * rootElectronVolt,
                     3.0 * rootElectronVolt },
                   0.5 * rootElectronVolt ) } );

    // usefull lambdas
    auto channelID = [] ( const auto& channel )
                        { return channel.channelID(); };

    THEN( "a SpinGroup can be constructed" ) {

      Energy energy = 1e-5 * electronVolts;

      SpinGroup< ReichMoore, ShiftFactor >
          group( { elastic, inelastic, protonEmission }, std::move( single ) );

      CHECK( 1 == group.incidentChannels().size() );
      CHECK( elastic.channelID()
             == std::visit( channelID, group.incidentChannels().front() ) );

      CHECK( elasticPair.pairID() == group.incidentPair().pairID() );

      CHECK( 3 == group.channels().size() );

      // channel 1 - elastic
      auto channel1 = std::get< Channel< Neutron > >( group.channels()[0] );

      CHECK( "n,Cl35{0,1,1+}" == channel1.channelID() );
      CHECK( "n,Cl35->n,Cl35" == channel1.reactionID().symbol() );

      CHECK( true == channel1.isIncidentChannel() );

      auto incident = channel1.incidentParticlePair();
      CHECK( 1.00866491582 == Approx( incident.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident.particle().spin() ) );
      CHECK( +1 == incident.particle().parity() );
      CHECK( 34.968852694 == Approx( incident.residual().mass().value ) );
      CHECK( 17. * e == Approx( incident.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident.residual().spin() ) );
      CHECK( +1 == incident.residual().parity() );
      CHECK( "n,Cl35" == incident.pairID().symbol() );

      auto pair = channel1.particlePair();
      CHECK( 1.00866491582 == Approx( pair.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair.particle().spin() ) );
      CHECK( +1 == pair.particle().parity() );
      CHECK( 34.968852694 == Approx( pair.residual().mass().value ) );
      CHECK( 17. * e == Approx( pair.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair.residual().spin() ) );
      CHECK( +1 == pair.residual().parity() );
      CHECK( "n,Cl35" == pair.pairID().symbol() );

      auto numbers = channel1.quantumNumbers();
      CHECK( 0 == numbers.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers.spin() );
      CHECK( 1.0 == numbers.totalAngularMomentum() );
      CHECK( +1 == numbers.parity() );

      auto radii = channel1.radii();
      CHECK( .4822220 == Approx( radii.penetrabilityRadius( energy ).value ) );
      CHECK( .4822220 == Approx( radii.shiftFactorRadius( energy ).value ) );
      CHECK( .3667980 == Approx( radii.phaseShiftRadius( energy ).value ) );

      CHECK( 0.0 == channel1.boundaryCondition() );

      CHECK( 0.0 == Approx( channel1.Q().value ) );

      // channel 2 - inelastic
      auto channel2 = std::get< Channel< Neutron > >( group.channels()[1] );

      CHECK( "n,Cl35_e1{0,1,1+}" == channel2.channelID() );
      CHECK( "n,Cl35->n,Cl35_e1" == channel2.reactionID().symbol() );

      CHECK( false == channel2.isIncidentChannel() );

      incident = channel2.incidentParticlePair();
      CHECK( 1.00866491582 == Approx( incident.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident.particle().spin() ) );
      CHECK( +1 == incident.particle().parity() );
      CHECK( 34.968852694 == Approx( incident.residual().mass().value ) );
      CHECK( 17. * e == Approx( incident.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident.residual().spin() ) );
      CHECK( +1 == incident.residual().parity() );
      CHECK( "n,Cl35" == incident.pairID().symbol() );

      pair = channel2.particlePair();
      CHECK( 1.00866491582 == Approx( pair.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair.particle().spin() ) );
      CHECK( +1 == pair.particle().parity() );
      CHECK( 34.968852694 == Approx( pair.residual().mass().value ) );
      CHECK( 17. * e == Approx( pair.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair.residual().spin() ) );
      CHECK( +1 == pair.residual().parity() );
      CHECK( "n,Cl35_e1" == pair.pairID().symbol() );

      numbers = channel2.quantumNumbers();
      CHECK( 0 == numbers.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers.spin() );
      CHECK( 1.0 == numbers.totalAngularMomentum() );
      CHECK( +1 == numbers.parity() );

      radii = channel2.radii();
      CHECK( .4822220 == Approx( radii.penetrabilityRadius( energy ).value ) );
      CHECK( .4822220 == Approx( radii.shiftFactorRadius( energy ).value ) );
      CHECK( .3667980 == Approx( radii.phaseShiftRadius( energy ).value ) );

      CHECK( 0.0 == channel2.boundaryCondition() );

      CHECK( -1.219440e+6 == Approx( channel2.Q().value ) );

      // channel 3 - proton emission
      auto channel3 = std::get< Channel< ChargedParticle > >( group.channels()[2] );

      CHECK( "p,S36{0,1,1+}" == channel3.channelID() );
      CHECK( "n,Cl35->p,S36" == channel3.reactionID().symbol() );

      CHECK( false == channel3.isIncidentChannel() );

      incident = channel3.incidentParticlePair();
      CHECK( 1.00866491582 == Approx( incident.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident.particle().spin() ) );
      CHECK( +1 == incident.particle().parity() );
      CHECK( 34.968852694 == Approx( incident.residual().mass().value ) );
      CHECK( 17. * e == Approx( incident.residual().charge().value ) );
      CHECK( 1.5 == Approx( incident.residual().spin() ) );
      CHECK( +1 == incident.residual().parity() );
      CHECK( "n,Cl35" == incident.pairID().symbol() );

      pair = channel3.particlePair();
      CHECK( 1.00727647 == Approx( pair.particle().mass().value ) );
      CHECK( e == Approx( pair.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair.particle().spin() ) );
      CHECK( +1 == pair.particle().parity() );
      CHECK( 35.967080699 == Approx( pair.residual().mass().value ) );
      CHECK( 16. * e == Approx( pair.residual().charge().value ) );
      CHECK( 1.5 == Approx( pair.residual().spin() ) );
      CHECK( +1 == pair.residual().parity() );
      CHECK( "p,S36" == pair.pairID().symbol() );

      numbers = channel3.quantumNumbers();
      CHECK( 0 == numbers.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers.spin() );
      CHECK( 1.0 == numbers.totalAngularMomentum() );
      CHECK( +1 == numbers.parity() );

      radii = channel3.radii();
      CHECK( .4822220 == Approx( radii.penetrabilityRadius( energy ).value ) );
      CHECK( .4822220 == Approx( radii.shiftFactorRadius( energy ).value ) );
      CHECK( .3667980 == Approx( radii.phaseShiftRadius( energy ).value ) );

      CHECK( 0.0 == channel3.boundaryCondition() );

      CHECK( 6.152200e+5 == Approx( channel3.Q().value ) );

      // reactions in the spin group
      auto reactions = group.reactionIDs();
      CHECK( 4 == reactions.size() );
      CHECK( "n,Cl35->n,Cl35" == reactions[0].symbol() );
      CHECK( "n,Cl35->n,Cl35_e1" == reactions[1].symbol() );
      CHECK( "n,Cl35->p,S36" == reactions[2].symbol() );
      CHECK( "n,Cl35->capture" == reactions[3].symbol() );

      // resonance data
      auto table = group.resonanceTable();
      CHECK( 3 == table.numberChannels() );
      CHECK( 3 == table.channels().size() );
      CHECK( "n,Cl35{0,1,1+}" == table.channels()[0] );
      CHECK( "n,Cl35_e1{0,1,1+}" == table.channels()[1] );
      CHECK( "p,S36{0,1,1+}" == table.channels()[2] );
      CHECK( 1 == table.numberResonances() );
      CHECK( 1 == table.resonances().size() );
      CHECK( 1 == table.energies().size() );
      CHECK( 0.25 == Approx( table.energies()[0].value ) );
      auto resonance = table.resonances()[0];
      CHECK( 0.25 == Approx( resonance.energy().value ) );
      CHECK( 0.5 == Approx( resonance.eliminatedWidth().value ) );
      CHECK( 3 == resonance.widths().size() );
      CHECK( 1.0 == Approx( resonance.widths()[0].value ) );
      CHECK( 2.0 == Approx( resonance.widths()[1].value ) );
      CHECK( 3.0 == Approx( resonance.widths()[2].value ) );
    }
  } // GIVEN

  GIVEN( "data for a SpinGroup with errors" ) {

    // particles
    Particle photon( ParticleID( "g" ), 0.0 * daltons, 0.0 * coulombs, 1., +1);
    Particle neutron( ParticleID( "n" ), 1.00866491582 * daltons,
                      0.0 * coulombs, 0.5, +1);
    Particle proton( ParticleID( "p" ), 1.00727647 * daltons,
                     elementary, 0.5, +1);
    Particle cl36( ParticleID( "Cl36" ), 35.968306822 * daltons,
                   17.0 * elementary, 0., +1);
    Particle cl35( ParticleID( "Cl35" ), 34.968852694 * daltons,
                   17.0 * elementary, 1.5, +1);
    Particle cl35_e1( ParticleID( "Cl35_e1" ), 34.968852694 * daltons,
                      17.0 * elementary, 1.5, +1);
    Particle s36( ParticleID( "S36" ), 35.967080699 * daltons,
                  16.0 * elementary, 1.5, +1);

    // particle pairs
    ParticlePair elasticPair( neutron, cl35 );
    ParticlePair inelasticPair( neutron, cl35_e1 );
    ParticlePair capturePair( photon, cl36 );
    ParticlePair protonEmissionPair( proton, s36 );

    // Q values
    QValue elasticQ = 0.0 * electronVolt;
    QValue inelasticQ = -1.219440e+6 * electronVolt;
    QValue captureQ = 0.0 * electronVolt;
    QValue protonEmissionQ = 6.152200e+5 * electronVolt;

    // quantum numbers
    ChannelQuantumNumbers elasticNumbers( 0, 1.0, 1.0, +1 );
    ChannelQuantumNumbers inelasticNumbers( 0, 1.0, 1.0, +1 );
    ChannelQuantumNumbers captureNumbers( 0, 0.0, 1.0, +1 );
    ChannelQuantumNumbers protonEmissionNumbers( 0, 1.0, 1.0, +1 );

    // channel radii
    ChannelRadii elasticRadii( 4.822220e-1 * rootBarn,
                               3.667980e-1 * rootBarn );
    ChannelRadii inelasticRadii( 4.822220e-1 * rootBarn,
                               3.667980e-1 * rootBarn );
    ChannelRadii captureRadii( 0.0 * rootBarn );
    ChannelRadii protonEmissionRadii( 4.822220e-1 * rootBarn,
                                      3.667980e-1 * rootBarn );

    // channels (the third channel has an incorrect J,pi with respect to the
    // other channels and the sixth has the wrong incident pair)
    Channel< Photon > capture( elasticPair, capturePair, captureQ,
                               captureNumbers, captureRadii );
    Channel< Neutron > elastic( elasticPair, elasticPair, elasticQ,
                                elasticNumbers, elasticRadii );
    Channel< Neutron > incorrect( elasticPair, elasticPair, elasticQ,
                                  { 0, 1.0, 0.0, +1 }, elasticRadii );
    Channel< Neutron > inelastic( elasticPair, inelasticPair, inelasticQ,
                                 inelasticNumbers, inelasticRadii );
    Channel< ChargedParticle > protonEmission( elasticPair,
                                               protonEmissionPair,
                                               protonEmissionQ,
                                               protonEmissionNumbers,
                                               protonEmissionRadii );
    Channel< ChargedParticle > wrong( inelasticPair, protonEmissionPair,
                                      protonEmissionQ, protonEmissionNumbers,
                                      protonEmissionRadii );

    // single resonance table
    ResonanceTable single(
      { elastic.channelID(), inelastic.channelID(), protonEmission.channelID() },
      { Resonance( 0.25 * electronVolt,
                   { 1.0 * rootElectronVolt, 2.0 * rootElectronVolt,
                     3.0 * rootElectronVolt },
                   0.5 * rootElectronVolt ) } );
    ResonanceTable wrongorder(
        { elastic.channelID(), protonEmission.channelID(), inelastic.channelID() },
        { Resonance( 0.25 * electronVolt,
                     { 1.0 * rootElectronVolt, 2.0 * rootElectronVolt,
                       3.0 * rootElectronVolt },
                     0.5 * rootElectronVolt ) } );
    ResonanceTable wrongsize(
      { elastic.channelID(), inelastic.channelID() },
      { Resonance( 0.25 * electronVolt,
                   { 1.0 * rootElectronVolt, 2.0 * rootElectronVolt },
                   0.5 * rootElectronVolt ) } );

    // no need to test with other template parameter: construction of SpinGroup
    // is independent of the template parameters
    using ReichMooreShiftFactorSpinGroup = SpinGroup< ReichMoore, ShiftFactor >;

    THEN( "an exception is thrown at construction for inconsistent Jpi" ) {

      ResonanceTable copy = single;
      CHECK_THROWS( ReichMooreShiftFactorSpinGroup
                          ( { incorrect, inelastic, protonEmission },
                            std::move( copy ) ) );
    }

    THEN( "an exception is thrown at construction for out of order channels "
          "between the SpinGroup's Channels and ResonanceTable" ) {

      CHECK_THROWS( ReichMooreShiftFactorSpinGroup
                          ( { incorrect, inelastic, protonEmission },
                            std::move( wrongorder ) ) );
    }

    THEN( "an exception is thrown at construction for inconsistent sizes  "
          "between the SpinGroup's Channels and ResonanceTable" ) {

      ResonanceTable copy = wrongsize;
      CHECK_THROWS( ReichMooreShiftFactorSpinGroup
                          ( { incorrect, inelastic, protonEmission },
                            std::move( copy ) ) );
    }

    THEN( "an exception is thrown at construction when the channels are not "
          "unique" ) {

      ResonanceTable copy = wrongsize;
      CHECK_THROWS( ReichMooreShiftFactorSpinGroup
                          ( { incorrect, inelastic, inelastic, protonEmission },
                            std::move( wrongsize ) ) );
    }

    THEN( "an exception is thrown at construction when the incident particle "
          "pair is not the same for each channel" ) {

      ResonanceTable copy = single;
      CHECK_THROWS( ReichMooreShiftFactorSpinGroup
                          ( { incorrect, inelastic, wrong },
                            std::move( copy ) ) );
    }
  } // GIVEN
} // SCENARIO
