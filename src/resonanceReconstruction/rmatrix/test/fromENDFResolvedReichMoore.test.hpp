std::string Pu239LRF3();
std::string Si29();

SCENARIO( "fromENDF - LRF3" ) {

  GIVEN( "valid ENDF data for Pu239" ) {

    std::string string = Pu239LRF3();
    auto begin = string.begin();
    auto end = string.end();
    long lineNumber = 1;

    njoy::ENDFtk::HeadRecord head( begin, end, lineNumber );
    njoy::ENDFtk::section::Type< 2, 151 > endf( head, begin, end, lineNumber, 9437 );
    ResonanceRange endfResonanceRange = endf.isotopes().front().resonanceRanges().front();

    auto resonances = fromENDF( endfResonanceRange, neutronMass, elementaryCharge, ParticleID( "n" ), ParticleID( "Pu239" ) );

    THEN( "the appropriate CompoundSystem is returned" ) {

      CHECK( true == resonances.isResolved() );
      CHECK( false == resonances.isUnresolved() );
      CHECK( 1e-5 == Approx( resonances.lowerEnergy().value ) );
      CHECK( 2500. == Approx( resonances.upperEnergy().value ) );

      auto compoundsystem = std::get< CompoundSystem< ReichMoore, ShiftFactor > >( resonances.compoundSystem() );

      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // content verification
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      // spin groups
      auto spingroups = compoundsystem.spinGroups();
      CHECK( 9 == spingroups.size() );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0 -  empty spin group
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup0 = spingroups[0];

      // channels
      auto channels0 = spingroup0.channels();

      CHECK( 1 == channels0.size() ); // 1 normal channel + 0 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel00 = std::get< Channel< Neutron > >( channels0[0] );
      CHECK( "n,Pu239->n,Pu239" == channel00.reactionID().symbol() );

      // incident particle pair
      const auto incident00 = channel00.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident00.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident00.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident00.particle().spin() ) );
      CHECK( +1 == incident00.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( incident00.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( incident00.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident00.residual().spin() ) );
      CHECK( +1 == incident00.residual().parity() );
      CHECK( "n,Pu239" == incident00.pairID().symbol() );

      // particle pair
      const auto pair00 = channel00.particlePair();
      CHECK( 1.008664 == Approx( pair00.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair00.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair00.particle().spin() ) );
      CHECK( +1 == pair00.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( pair00.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( pair00.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair00.residual().spin() ) );
      CHECK( +1 == pair00.residual().parity() );
      CHECK( "n,Pu239" == pair00.pairID().symbol() );

      // quantum numbers
      const auto numbers00 = channel00.quantumNumbers();
      CHECK( 1 == numbers00.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers00.spin() );
      CHECK( 0.0 == numbers00.totalAngularMomentum() );
      CHECK( -1 == numbers00.parity() );
      CHECK( "{1,1,0-}" == numbers00.toString() );

      // radii
      const auto radii00 = channel00.radii();
      CHECK( .9410000 == Approx( radii00.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii00.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii00.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel00.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel00.Q().value ) );

      // resonance table
      auto table0 = spingroup0.resonanceTable();

      CHECK( 1 == table0.numberChannels() ); // 1 normal channel + 0 eliminated
      CHECK( 0 == table0.numberResonances() );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup1 = spingroups[1];

      // channels
      auto channels1 = spingroup1.channels();

      CHECK( 3 == channels1.size() ); // 3 normal channel + 1 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel10 = std::get< Channel< Neutron > >( channels1[0] );
      CHECK( "n,Pu239->n,Pu239" == channel10.reactionID().symbol() );

      // incident particle pair
      const auto incident10 = channel10.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident10.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident10.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident10.particle().spin() ) );
      CHECK( +1 == incident10.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( incident10.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( incident10.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident10.residual().spin() ) );
      CHECK( +1 == incident10.residual().parity() );
      CHECK( "n,Pu239" == incident10.pairID().symbol() );

      // particle pair
      const auto pair10 = channel10.particlePair();
      CHECK( 1.008664 == Approx( pair10.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair10.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair10.particle().spin() ) );
      CHECK( +1 == pair10.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( pair10.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( pair10.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair10.residual().spin() ) );
      CHECK( +1 == pair10.residual().parity() );
      CHECK( "n,Pu239" == pair10.pairID().symbol() );

      // quantum numbers
      const auto numbers10 = channel10.quantumNumbers();
      CHECK( 0 == numbers10.orbitalAngularMomentum() );
      CHECK( 0.0 == numbers10.spin() );
      CHECK( 0.0 == numbers10.totalAngularMomentum() );
      CHECK( +1 == numbers10.parity() );
      CHECK( "{0,0,0+}" == numbers10.toString() );

      // radii
      const auto radii10 = channel10.radii();
      CHECK( .9410000 == Approx( radii10.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii10.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii10.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel10.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel10.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1, channel 1: fission1
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel11 = std::get< Channel< Fission > >( channels1[1] );
      CHECK( "n,Pu239->fission" == channel11.reactionID().symbol() );

      // incident particle pair
      const auto incident11 = channel11.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident11.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident11.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident11.particle().spin() ) );
      CHECK( +1 == incident11.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( incident11.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( incident11.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident11.residual().spin() ) );
      CHECK( +1 == incident11.residual().parity() );
      CHECK( "n,Pu239" == incident11.pairID().symbol() );

      // quantum numbers
      const auto numbers11 = channel11.quantumNumbers();
      CHECK( 0 == numbers11.orbitalAngularMomentum() );
      CHECK( 0.0 == numbers11.spin() );
      CHECK( 0.0 == numbers11.totalAngularMomentum() );
      CHECK( +1 == numbers11.parity() );
      CHECK( "{0,0,0+}" == numbers11.toString() );

      // radii
      const auto radii11 = channel11.radii();
      CHECK( .9410000 == Approx( radii11.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii11.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii11.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel11.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel11.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1, channel 2: fission2
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel12 = std::get< Channel< Fission > >( channels1[2] );
      CHECK( "n,Pu239->fission" == channel12.reactionID().symbol() );

      // incident particle pair
      const auto incident12 = channel12.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident12.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident12.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident12.particle().spin() ) );
      CHECK( +1 == incident12.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( incident12.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( incident12.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident12.residual().spin() ) );
      CHECK( +1 == incident12.residual().parity() );
      CHECK( "n,Pu239" == incident12.pairID().symbol() );

      // quantum numbers
      const auto numbers12 = channel12.quantumNumbers();
      CHECK( 0 == numbers12.orbitalAngularMomentum() );
      CHECK( 0.0 == numbers12.spin() );
      CHECK( 0.0 == numbers12.totalAngularMomentum() );
      CHECK( +1 == numbers12.parity() );
      CHECK( "{0,0,0+}" == numbers12.toString() );

      // radii
      const auto radii12 = channel12.radii();
      CHECK( .9410000 == Approx( radii12.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii12.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii12.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel12.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel12.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1, resonance table
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      auto table1 = spingroup1.resonanceTable();

      CHECK( 3 == table1.numberChannels() ); // 1 normal channel + 1 eliminated
      CHECK( 231 == table1.numberResonances() );

      auto energies1 = table1.energies();
      CHECK( -6.908700 == Approx( energies1.front().value ) );
      CHECK( 2.511660e+3 == Approx( energies1.back().value ) );

      auto resonances1 = table1.resonances();
      CHECK( -6.908700 == Approx( resonances1.front().energy().value ) );
      CHECK( 2.511660e+3 == Approx( resonances1.back().energy().value ) );
      CHECK( 3 == resonances1.front().widths().size() );
      CHECK( 3 == resonances1.back().widths().size() );
      CHECK( std::sqrt( 1.805854e-2 / 2. / channel10.penetrability( -6.908700 * electronVolt ) )
             == Approx( resonances1.front().widths()[0].value ) );
      CHECK( std::sqrt( 3.123558e-1 / 2. / channel10.penetrability( 2.511660e+3 * electronVolt ) )
             == Approx( resonances1.back().widths()[0].value ) );
      CHECK( -std::sqrt( 9.235517e-1 / 2. )
             == Approx( resonances1.front().widths()[1].value ) );
      CHECK( std::sqrt( 1.002582 / 2. )
             == Approx( resonances1.back().widths()[1].value ) );
      CHECK( std::sqrt( 3.464891e-1 / 2. )
             == Approx( resonances1.front().widths()[2].value ) );
      CHECK( std::sqrt( 1.971943e-2 / 2. )
             == Approx( resonances1.back().widths()[2].value ) );
      CHECK( std::sqrt( 6.012769e-2 / 2. ) == Approx( resonances1.front().eliminatedWidth().value ) );
      CHECK( std::sqrt( 3.926596e-2 / 2. ) == Approx( resonances1.back().eliminatedWidth().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 2 -  empty spin group, 2 elastic channels
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup2 = spingroups[2];

      // channels
      auto channels2 = spingroup2.channels();

      CHECK( 2 == channels2.size() ); // 1 normal channel + 0 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 2, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel20 = std::get< Channel< Neutron > >( channels2[0] );
      CHECK( "n,Pu239->n,Pu239" == channel20.reactionID().symbol() );

      // incident particle pair
      const auto incident20 = channel20.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident20.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident20.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident20.particle().spin() ) );
      CHECK( +1 == incident20.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( incident20.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( incident20.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident20.residual().spin() ) );
      CHECK( +1 == incident20.residual().parity() );
      CHECK( "n,Pu239" == incident20.pairID().symbol() );

      // particle pair
      const auto pair20 = channel20.particlePair();
      CHECK( 1.008664 == Approx( pair20.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair20.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair20.particle().spin() ) );
      CHECK( +1 == pair20.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( pair20.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( pair20.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair20.residual().spin() ) );
      CHECK( +1 == pair20.residual().parity() );
      CHECK( "n,Pu239" == pair20.pairID().symbol() );

      // quantum numbers
      const auto numbers20 = channel20.quantumNumbers();
      CHECK( 1 == numbers20.orbitalAngularMomentum() );
      CHECK( 0.0 == numbers20.spin() );
      CHECK( 1.0 == numbers20.totalAngularMomentum() );
      CHECK( -1 == numbers20.parity() );
      CHECK( "{1,0,1-}" == numbers20.toString() );

      // radii
      const auto radii20 = channel20.radii();
      CHECK( .9410000 == Approx( radii20.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii20.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii20.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel20.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel20.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 2, channel 1: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel21 = std::get< Channel< Neutron > >( channels2[1] );
      CHECK( "n,Pu239->n,Pu239" == channel21.reactionID().symbol() );

      // incident particle pair
      const auto incident21 = channel21.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident21.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident21.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident21.particle().spin() ) );
      CHECK( +1 == incident21.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( incident21.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( incident21.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident21.residual().spin() ) );
      CHECK( +1 == incident21.residual().parity() );
      CHECK( "n,Pu239" == incident21.pairID().symbol() );

      // particle pair
      const auto pair21 = channel21.particlePair();
      CHECK( 1.008664 == Approx( pair21.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair21.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair21.particle().spin() ) );
      CHECK( +1 == pair21.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( pair21.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( pair21.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair21.residual().spin() ) );
      CHECK( +1 == pair21.residual().parity() );
      CHECK( "n,Pu239" == pair21.pairID().symbol() );

      // quantum numbers
      const auto numbers21 = channel21.quantumNumbers();
      CHECK( 1 == numbers21.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers21.spin() );
      CHECK( 1.0 == numbers21.totalAngularMomentum() );
      CHECK( -1 == numbers21.parity() );
      CHECK( "{1,1,1-}" == numbers21.toString() );

      // radii
      const auto radii21 = channel21.radii();
      CHECK( .9410000 == Approx( radii21.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii21.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii21.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel21.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel21.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 2, resonance table
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      auto table2 = spingroup2.resonanceTable();

      CHECK( 2 == table2.numberChannels() ); // 1 normal channel + 0 eliminated
      CHECK( 0 == table2.numberResonances() );

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
      CHECK( "n,Pu239->n,Pu239" == channel30.reactionID().symbol() );

      // incident particle pair
      const auto incident30 = channel30.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident30.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident30.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident30.particle().spin() ) );
      CHECK( +1 == incident30.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( incident30.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( incident30.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident30.residual().spin() ) );
      CHECK( +1 == incident30.residual().parity() );
      CHECK( "n,Pu239" == incident30.pairID().symbol() );

      // particle pair
      const auto pair30 = channel30.particlePair();
      CHECK( 1.008664 == Approx( pair30.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair30.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair30.particle().spin() ) );
      CHECK( +1 == pair30.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( pair30.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( pair30.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair30.residual().spin() ) );
      CHECK( +1 == pair30.residual().parity() );
      CHECK( "n,Pu239" == pair30.pairID().symbol() );

      // quantum numbers
      const auto numbers30 = channel30.quantumNumbers();
      CHECK( 0 == numbers30.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers30.spin() );
      CHECK( 1.0 == numbers30.totalAngularMomentum() );
      CHECK( +1 == numbers30.parity() );
      CHECK( "{0,1,1+}" == numbers30.toString() );

      // radii
      const auto radii30 = channel30.radii();
      CHECK( .9410000 == Approx( radii30.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii30.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii30.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel30.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel30.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 3, channel 1: fission1
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel31 = std::get< Channel< Fission > >( channels3[1] );
      CHECK( "n,Pu239->fission" == channel31.reactionID().symbol() );

      // incident particle pair
      const auto incident31 = channel31.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident31.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident31.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident31.particle().spin() ) );
      CHECK( +1 == incident31.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( incident31.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( incident31.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident31.residual().spin() ) );
      CHECK( +1 == incident31.residual().parity() );
      CHECK( "n,Pu239" == incident31.pairID().symbol() );

      // quantum numbers
      const auto numbers31 = channel31.quantumNumbers();
      CHECK( 0 == numbers31.orbitalAngularMomentum() );
      CHECK( 0.0 == numbers31.spin() );
      CHECK( 1.0 == numbers31.totalAngularMomentum() );
      CHECK( +1 == numbers31.parity() );
      CHECK( "{0,0,1+}" == numbers31.toString() );

      // radii
      const auto radii31 = channel31.radii();
      CHECK( .9410000 == Approx( radii31.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii31.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii31.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel31.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel31.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 3, channel 2: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel32 = std::get< Channel< Neutron > >( channels3[2] );
      CHECK( "n,Pu239->n,Pu239" == channel32.reactionID().symbol() );

      // incident particle pair
      const auto incident32 = channel32.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident32.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident32.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident32.particle().spin() ) );
      CHECK( +1 == incident32.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( incident32.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( incident32.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident32.residual().spin() ) );
      CHECK( +1 == incident32.residual().parity() );
      CHECK( "n,Pu239" == incident32.pairID().symbol() );

      // particle pair
      const auto pair32 = channel32.particlePair();
      CHECK( 1.008664 == Approx( pair32.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair32.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair32.particle().spin() ) );
      CHECK( +1 == pair32.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( pair32.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( pair32.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair32.residual().spin() ) );
      CHECK( +1 == pair32.residual().parity() );
      CHECK( "n,Pu239" == pair32.pairID().symbol() );

      // quantum numbers
      const auto numbers32 = channel32.quantumNumbers();
      CHECK( 2 == numbers32.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers32.spin() );
      CHECK( 1.0 == numbers32.totalAngularMomentum() );
      CHECK( +1 == numbers32.parity() );
      CHECK( "{2,1,1+}" == numbers32.toString() );

      // radii
      const auto radii32 = channel32.radii();
      CHECK( .9410000 == Approx( radii32.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii32.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii32.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel32.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel32.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 3, resonance table
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      auto table3 = spingroup3.resonanceTable();

      CHECK( 3 == table3.numberChannels() ); // 1 normal channel + 1 eliminated
      CHECK( 812 == table3.numberResonances() );

      auto energies3 = table3.energies();
      CHECK( -1.500200e+2 == Approx( energies3.front().value ) );
      CHECK( 2.523800e+3 == Approx( energies3.back().value ) );

      auto resonances3 = table3.resonances();
      CHECK( -1.500200e+2 == Approx( resonances3.front().energy().value ) );
      CHECK( 2.523800e+3 == Approx( resonances3.back().energy().value ) );
      CHECK( 3 == resonances3.front().widths().size() );
      CHECK( 3 == resonances3.back().widths().size() );
      CHECK( std::sqrt( .3726588 / 2. / channel30.penetrability( -1.500200e+2 * electronVolt ) )
             == Approx( resonances3.front().widths()[0].value ) );
      CHECK( std::sqrt( .1829061 / 2. / channel30.penetrability( 2.523800e+3 * electronVolt ) )
             == Approx( resonances3.back().widths()[0].value ) );
      CHECK( std::sqrt( 1.352425e-1 / 2. )
             == Approx( resonances3.front().widths()[1].value ) );
      CHECK( std::sqrt( 4.074383e-1 / 2. )
             == Approx( resonances3.back().widths()[1].value ) );
      CHECK( 0. == Approx( resonances3.front().widths()[2].value ) );
      CHECK( 0. == Approx( resonances3.back().widths()[2].value ) );
      CHECK( std::sqrt( 4.685986e-2 / 2. ) == Approx( resonances3.front().eliminatedWidth().value ) );
      CHECK( std::sqrt( 3.926666e-2 / 2. ) == Approx( resonances3.back().eliminatedWidth().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 4 -  empty spin group, 2 elastic channels
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup4 = spingroups[4];

      // channels
      auto channels4 = spingroup4.channels();

      CHECK( 2 == channels4.size() ); // 1 normal channel + 0 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 4, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel40 = std::get< Channel< Neutron > >( channels4[0] );
      CHECK( "n,Pu239->n,Pu239" == channel40.reactionID().symbol() );

      // incident particle pair
      const auto incident40 = channel40.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident40.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident40.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident40.particle().spin() ) );
      CHECK( +1 == incident40.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( incident40.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( incident40.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident40.residual().spin() ) );
      CHECK( +1 == incident40.residual().parity() );
      CHECK( "n,Pu239" == incident40.pairID().symbol() );

      // particle pair
      const auto pair40 = channel40.particlePair();
      CHECK( 1.008664 == Approx( pair40.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair40.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair40.particle().spin() ) );
      CHECK( +1 == pair40.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( pair40.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( pair40.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair40.residual().spin() ) );
      CHECK( +1 == pair40.residual().parity() );
      CHECK( "n,Pu239" == pair40.pairID().symbol() );

      // quantum numbers
      const auto numbers40 = channel40.quantumNumbers();
      CHECK( 1 == numbers40.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers40.spin() );
      CHECK( 2.0 == numbers40.totalAngularMomentum() );
      CHECK( -1 == numbers40.parity() );
      CHECK( "{1,1,2-}" == numbers40.toString() );

      // radii
      const auto radii40 = channel40.radii();
      CHECK( .9410000 == Approx( radii40.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii40.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii40.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel40.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel40.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 4, channel 1: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel41 = std::get< Channel< Neutron > >( channels4[1] );
      CHECK( "n,Pu239->n,Pu239" == channel41.reactionID().symbol() );

      // incident particle pair
      const auto incident41 = channel41.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident41.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident41.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident41.particle().spin() ) );
      CHECK( +1 == incident41.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( incident41.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( incident41.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident41.residual().spin() ) );
      CHECK( +1 == incident41.residual().parity() );
      CHECK( "n,Pu239" == incident41.pairID().symbol() );

      // particle pair
      const auto pair41 = channel41.particlePair();
      CHECK( 1.008664 == Approx( pair41.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair41.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair41.particle().spin() ) );
      CHECK( +1 == pair41.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( pair41.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( pair41.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair41.residual().spin() ) );
      CHECK( +1 == pair41.residual().parity() );
      CHECK( "n,Pu239" == pair41.pairID().symbol() );

      // quantum numbers
      const auto numbers41 = channel41.quantumNumbers();
      CHECK( 3 == numbers41.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers41.spin() );
      CHECK( 2.0 == numbers41.totalAngularMomentum() );
      CHECK( -1 == numbers41.parity() );
      CHECK( "{3,1,2-}" == numbers41.toString() );

      // radii
      const auto radii41 = channel41.radii();
      CHECK( .9410000 == Approx( radii41.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii41.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii41.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel41.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel41.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 4, resonance table
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      auto table4 = spingroup4.resonanceTable();

      CHECK( 2 == table4.numberChannels() ); // 1 normal channel + 0 eliminated
      CHECK( 0 == table4.numberResonances() );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 5 -  empty spin group, 2 elastic channels
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup5 = spingroups[5];

      // channels
      auto channels5 = spingroup5.channels();

      CHECK( 2 == channels5.size() ); // 1 normal channel + 0 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 5, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel50 = std::get< Channel< Neutron > >( channels5[0] );
      CHECK( "n,Pu239->n,Pu239" == channel50.reactionID().symbol() );

      // incident particle pair
      const auto incident50 = channel50.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident50.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident50.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident50.particle().spin() ) );
      CHECK( +1 == incident50.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( incident50.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( incident50.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident50.residual().spin() ) );
      CHECK( +1 == incident50.residual().parity() );
      CHECK( "n,Pu239" == incident50.pairID().symbol() );

      // particle pair
      const auto pair50 = channel50.particlePair();
      CHECK( 1.008664 == Approx( pair50.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair50.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair50.particle().spin() ) );
      CHECK( +1 == pair50.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( pair50.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( pair50.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair50.residual().spin() ) );
      CHECK( +1 == pair50.residual().parity() );
      CHECK( "n,Pu239" == pair50.pairID().symbol() );

      // quantum numbers
      const auto numbers50 = channel50.quantumNumbers();
      CHECK( 2 == numbers50.orbitalAngularMomentum() );
      CHECK( 0.0 == numbers50.spin() );
      CHECK( 2.0 == numbers50.totalAngularMomentum() );
      CHECK( +1 == numbers50.parity() );
      CHECK( "{2,0,2+}" == numbers50.toString() );

      // radii
      const auto radii50 = channel50.radii();
      CHECK( .9410000 == Approx( radii50.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii50.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii50.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel50.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel50.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 5, channel 1: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel51 = std::get< Channel< Neutron > >( channels5[1] );
      CHECK( "n,Pu239->n,Pu239" == channel51.reactionID().symbol() );

      // incident particle pair
      const auto incident51 = channel51.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident51.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident51.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident51.particle().spin() ) );
      CHECK( +1 == incident51.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( incident51.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( incident51.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident51.residual().spin() ) );
      CHECK( +1 == incident51.residual().parity() );
      CHECK( "n,Pu239" == incident51.pairID().symbol() );

      // particle pair
      const auto pair51 = channel51.particlePair();
      CHECK( 1.008664 == Approx( pair51.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair51.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair51.particle().spin() ) );
      CHECK( +1 == pair51.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( pair51.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( pair51.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair51.residual().spin() ) );
      CHECK( +1 == pair51.residual().parity() );
      CHECK( "n,Pu239" == pair51.pairID().symbol() );

      // quantum numbers
      const auto numbers51 = channel51.quantumNumbers();
      CHECK( 2 == numbers51.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers51.spin() );
      CHECK( 2.0 == numbers51.totalAngularMomentum() );
      CHECK( +1 == numbers51.parity() );
      CHECK( "{2,1,2+}" == numbers51.toString() );

      // radii
      const auto radii51 = channel51.radii();
      CHECK( .9410000 == Approx( radii51.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii51.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii51.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel51.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel51.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 5, resonance table
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      auto table5 = spingroup5.resonanceTable();

      CHECK( 2 == table5.numberChannels() ); // 1 normal channel + 0 eliminated
      CHECK( 0 == table5.numberResonances() );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 6 -  empty spin group, 2 elastic channels
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup6 = spingroups[6];

      // channels
      auto channels6 = spingroup6.channels();

      CHECK( 2 == channels6.size() ); // 1 normal channel + 0 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 6, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel60 = std::get< Channel< Neutron > >( channels6[0] );
      CHECK( "n,Pu239->n,Pu239" == channel60.reactionID().symbol() );

      // incident particle pair
      const auto incident60 = channel60.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident60.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident60.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident60.particle().spin() ) );
      CHECK( +1 == incident60.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( incident60.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( incident60.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident60.residual().spin() ) );
      CHECK( +1 == incident60.residual().parity() );
      CHECK( "n,Pu239" == incident60.pairID().symbol() );

      // particle pair
      const auto pair60 = channel60.particlePair();
      CHECK( 1.008664 == Approx( pair60.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair60.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair60.particle().spin() ) );
      CHECK( +1 == pair60.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( pair60.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( pair60.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair60.residual().spin() ) );
      CHECK( +1 == pair60.residual().parity() );
      CHECK( "n,Pu239" == pair60.pairID().symbol() );

      // quantum numbers
      const auto numbers60 = channel60.quantumNumbers();
      CHECK( 3 == numbers60.orbitalAngularMomentum() );
      CHECK( 0.0 == numbers60.spin() );
      CHECK( 3.0 == numbers60.totalAngularMomentum() );
      CHECK( -1 == numbers60.parity() );
      CHECK( "{3,0,3-}" == numbers60.toString() );

      // radii
      const auto radii60 = channel60.radii();
      CHECK( .9410000 == Approx( radii60.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii60.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii60.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel60.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel60.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 6, channel 1: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel61 = std::get< Channel< Neutron > >( channels6[1] );
      CHECK( "n,Pu239->n,Pu239" == channel61.reactionID().symbol() );

      // incident particle pair
      const auto incident61 = channel61.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident61.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident61.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident61.particle().spin() ) );
      CHECK( +1 == incident61.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( incident61.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( incident61.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident61.residual().spin() ) );
      CHECK( +1 == incident61.residual().parity() );
      CHECK( "n,Pu239" == incident61.pairID().symbol() );

      // particle pair
      const auto pair61 = channel61.particlePair();
      CHECK( 1.008664 == Approx( pair61.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair61.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair61.particle().spin() ) );
      CHECK( +1 == pair61.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( pair61.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( pair61.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair61.residual().spin() ) );
      CHECK( +1 == pair61.residual().parity() );
      CHECK( "n,Pu239" == pair61.pairID().symbol() );

      // quantum numbers
      const auto numbers61 = channel61.quantumNumbers();
      CHECK( 3 == numbers61.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers61.spin() );
      CHECK( 3.0 == numbers61.totalAngularMomentum() );
      CHECK( -1 == numbers61.parity() );
      CHECK( "{3,1,3-}" == numbers61.toString() );

      // radii
      const auto radii61 = channel61.radii();
      CHECK( .9410000 == Approx( radii61.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii61.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii61.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel61.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel61.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 6, resonance table
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      auto table6 = spingroup6.resonanceTable();

      CHECK( 2 == table6.numberChannels() ); // 1 normal channel + 0 eliminated
      CHECK( 0 == table6.numberResonances() );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 7 -  empty spin group, 1 elastic channels
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup7 = spingroups[7];

      // channels
      auto channels7 = spingroup7.channels();

      CHECK( 1 == channels7.size() ); // 1 normal channel + 0 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 7, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel70 = std::get< Channel< Neutron > >( channels7[0] );
      CHECK( "n,Pu239->n,Pu239" == channel70.reactionID().symbol() );

      // incident particle pair
      const auto incident70 = channel70.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident70.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident70.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident70.particle().spin() ) );
      CHECK( +1 == incident70.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( incident70.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( incident70.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident70.residual().spin() ) );
      CHECK( +1 == incident70.residual().parity() );
      CHECK( "n,Pu239" == incident70.pairID().symbol() );

      // particle pair
      const auto pair70 = channel70.particlePair();
      CHECK( 1.008664 == Approx( pair70.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair70.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair70.particle().spin() ) );
      CHECK( +1 == pair70.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( pair70.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( pair70.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair70.residual().spin() ) );
      CHECK( +1 == pair70.residual().parity() );
      CHECK( "n,Pu239" == pair70.pairID().symbol() );

      // quantum numbers
      const auto numbers70 = channel70.quantumNumbers();
      CHECK( 2 == numbers70.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers70.spin() );
      CHECK( 3.0 == numbers70.totalAngularMomentum() );
      CHECK( +1 == numbers70.parity() );
      CHECK( "{2,1,3+}" == numbers70.toString() );

      // radii
      const auto radii70 = channel70.radii();
      CHECK( .9410000 == Approx( radii70.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii70.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii70.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel70.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel70.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 7, resonance table
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      auto table7 = spingroup7.resonanceTable();

      CHECK( 1 == table7.numberChannels() ); // 1 normal channel + 0 eliminated
      CHECK( 0 == table7.numberResonances() );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 8 -  empty spin group, 2 elastic channels
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup8 = spingroups[8];

      // channels
      auto channels8 = spingroup8.channels();

      CHECK( 1 == channels8.size() ); // 1 normal channel + 0 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 8, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel80 = std::get< Channel< Neutron > >( channels8[0] );
      CHECK( "n,Pu239->n,Pu239" == channel80.reactionID().symbol() );

      // incident particle pair
      const auto incident80 = channel80.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident80.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident80.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident80.particle().spin() ) );
      CHECK( +1 == incident80.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( incident80.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( incident80.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident80.residual().spin() ) );
      CHECK( +1 == incident80.residual().parity() );
      CHECK( "n,Pu239" == incident80.pairID().symbol() );

      // particle pair
      const auto pair80 = channel80.particlePair();
      CHECK( 1.008664 == Approx( pair80.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair80.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair80.particle().spin() ) );
      CHECK( +1 == pair80.particle().parity() );
      CHECK( 236.9986 * 1.008664 == Approx( pair80.residual().mass().value ) );
      CHECK( 94.0 * 1.602e-19 == Approx( pair80.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair80.residual().spin() ) );
      CHECK( +1 == pair80.residual().parity() );
      CHECK( "n,Pu239" == pair80.pairID().symbol() );

      // quantum numbers
      const auto numbers80 = channel80.quantumNumbers();
      CHECK( 3 == numbers80.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers80.spin() );
      CHECK( 4.0 == numbers80.totalAngularMomentum() );
      CHECK( -1 == numbers80.parity() );
      CHECK( "{3,1,4-}" == numbers80.toString() );

      // radii
      const auto radii80 = channel80.radii();
      CHECK( .9410000 == Approx( radii80.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii80.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .9410000 == Approx( radii80.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel80.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel80.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 8, resonance table
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      auto table8 = spingroup8.resonanceTable();

      CHECK( 1 == table8.numberChannels() ); // 1 normal channel + 0 eliminated
      CHECK( 0 == table8.numberResonances() );
    } // THEN

    THEN( "cross sections can be reconstructed" ) {

      // values taken from NJOY2016 PeNDF tape for ENDF/B-VIII.0 Pu239

      ReactionID elas( "n,Pu239->n,Pu239" );
      ReactionID fiss( "n,Pu239->fission" );
      ReactionID capt( "n,Pu239->capture" );
      std::map< ReactionID, CrossSection > xs;

      xs = resonances( 1e-5 * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 8.152130 == Approx( xs[ elas ].value ) );
      CHECK( 3.456462e+4 == Approx( xs[ fiss ].value ) );
      CHECK( 1.284211e+4 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-4 * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 8.151873 == Approx( xs[ elas ].value ) );
      CHECK( 1.093649e+4 == Approx( xs[ fiss ].value ) );
      CHECK( 4.060599e+3 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-3 * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 8.149338 == Approx( xs[ elas ].value ) );
      CHECK( 3.477441e+3 == Approx( xs[ fiss ].value ) );
      CHECK( 1.282863e+3 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-2 * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 8.120662 == Approx( xs[ elas ].value ) );
      CHECK( 1.144071e+3 == Approx( xs[ fiss ].value ) );
      CHECK( 4.072596e+2 == Approx( xs[ capt ].value ) );

      xs = resonances( 0.0253 * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 8.055773 == Approx( xs[ elas ].value ) );
      CHECK( 7.469471e+2 == Approx( xs[ fiss ].value ) );
      CHECK( 2.697482e+2 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-1 * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 7.542244 == Approx( xs[ elas ].value ) );
      CHECK( 4.771902e+2 == Approx( xs[ fiss ].value ) );
      CHECK( 2.320353e+2 == Approx( xs[ capt ].value ) );

      xs = resonances( 1. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 1.015843e+1 == Approx( xs[ elas ].value ) );
      CHECK( 3.948663e+1 == Approx( xs[ fiss ].value ) );
      CHECK( 7.661602 == Approx( xs[ capt ].value ) );

      xs = resonances( 7.81580 * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 3.078505e+1 == Approx( xs[ elas ].value ) );
      CHECK( 1.286856e+3 == Approx( xs[ fiss ].value ) );
      CHECK( 1.109709e+3 == Approx( xs[ capt ].value ) );

      xs = resonances( 10. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 6.422677 == Approx( xs[ elas ].value ) );
      CHECK( 1.080674e+1 == Approx( xs[ fiss ].value ) );
      CHECK( 4.830660 == Approx( xs[ capt ].value ) );

      xs = resonances( 10.928 * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 2.126328e+1 == Approx( xs[ elas ].value ) );
      CHECK( 1.353841e+3 == Approx( xs[ fiss ].value ) );
      CHECK( 3.245726e+2 == Approx( xs[ capt ].value ) );

      xs = resonances( 100. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 1.085331e+1 == Approx( xs[ elas ].value ) );
      CHECK( 1.786214e+1 == Approx( xs[ fiss ].value ) );
      CHECK( 3.802623e-1 == Approx( xs[ capt ].value ) );

      xs = resonances( 1000. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 1.813166e+1 == Approx( xs[ elas ].value ) );
      CHECK( 2.852688 == Approx( xs[ fiss ].value ) );
      CHECK( 6.884886e-1 == Approx( xs[ capt ].value ) );

      xs = resonances( 1.960508e+3 * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 3.428219e+1 == Approx( xs[ elas ].value ) );
      CHECK( 7.707020e+1 == Approx( xs[ fiss ].value ) );
      CHECK( 5.625787e+1 == Approx( xs[ capt ].value ) );
    } // THEN
  } // GIVEN

  GIVEN( "valid ENDF data for Si29" ) {

    std::string string = Si29();
    auto begin = string.begin();
    auto end = string.end();
    long lineNumber = 1;

    njoy::ENDFtk::HeadRecord head( begin, end, lineNumber );
    njoy::ENDFtk::section::Type< 2, 151 > endf( head, begin, end, lineNumber, 1428 );
    ResonanceRange endfResonanceRange = endf.isotopes().front().resonanceRanges().front();

    auto resonances = fromENDF( endfResonanceRange, neutronMass, elementaryCharge, ParticleID( "n" ), ParticleID( "Si29" ) );

    THEN( "the appropriate CompoundSystem is returned" ) {

      CHECK( true == resonances.isResolved() );
      CHECK( false == resonances.isUnresolved() );
      CHECK( 1e-5 == Approx( resonances.lowerEnergy().value ) );
      CHECK( 1.3e+6 == Approx( resonances.upperEnergy().value ) );

      auto compoundsystem = std::get< CompoundSystem< ReichMoore, ShiftFactor > >( resonances.compoundSystem() );

      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // content verification
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      // spin groups
      auto spingroups = compoundsystem.spinGroups();
      CHECK( 7 == spingroups.size() );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0 -  empty spin group, 1 elastic channels
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup0 = spingroups[0];

      // channels
      auto channels0 = spingroup0.channels();

      CHECK( 1 == channels0.size() ); // 1 normal channel + 0 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel00 = std::get< Channel< Neutron > >( channels0[0] );
      CHECK( "n,Si29->n,Si29" == channel00.reactionID().symbol() );

      // incident particle pair
      const auto incident00 = channel00.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident00.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident00.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident00.particle().spin() ) );
      CHECK( +1 == incident00.particle().parity() );
      CHECK( 28.72800 * 1.008664 == Approx( incident00.residual().mass().value ) );
      CHECK( 14.0 * 1.602e-19 == Approx( incident00.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident00.residual().spin() ) );
      CHECK( +1 == incident00.residual().parity() );
      CHECK( "n,Si29" == incident00.pairID().symbol() );

      // particle pair
      const auto pair00 = channel00.particlePair();
      CHECK( 1.008664 == Approx( pair00.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair00.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair00.particle().spin() ) );
      CHECK( +1 == pair00.particle().parity() );
      CHECK( 28.72800 * 1.008664 == Approx( pair00.residual().mass().value ) );
      CHECK( 14.0 * 1.602e-19 == Approx( pair00.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair00.residual().spin() ) );
      CHECK( +1 == pair00.residual().parity() );
      CHECK( "n,Si29" == pair00.pairID().symbol() );

      // quantum numbers
      const auto numbers00 = channel00.quantumNumbers();
      CHECK( 1 == numbers00.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers00.spin() );
      CHECK( 0.0 == numbers00.totalAngularMomentum() );
      CHECK( -1 == numbers00.parity() );
      CHECK( "{1,1,0-}" == numbers00.toString() );

      // radii
      const auto radii00 = channel00.radii();
      CHECK( .44 == Approx( radii00.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .44 == Approx( radii00.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .44 == Approx( radii00.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel00.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel00.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0, resonance table
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      auto table0 = spingroup0.resonanceTable();

      CHECK( 1 == table0.numberChannels() ); // 1 normal channel + 0 eliminated
      CHECK( 0 == table0.numberResonances() );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1 -  1 elastic channels + 1 eliminated
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup1 = spingroups[1];

      // channels
      auto channels1 = spingroup1.channels();

      CHECK( 1 == channels1.size() ); // 1 normal channel + 1 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1, channel 1: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel10 = std::get< Channel< Neutron > >( channels1[0] );
      CHECK( "n,Si29->n,Si29" == channel10.reactionID().symbol() );

      // incident particle pair
      const auto incident10 = channel10.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident10.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident10.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident10.particle().spin() ) );
      CHECK( +1 == incident10.particle().parity() );
      CHECK( 28.72800 * 1.008664 == Approx( incident10.residual().mass().value ) );
      CHECK( 14.0 * 1.602e-19 == Approx( incident10.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident10.residual().spin() ) );
      CHECK( +1 == incident10.residual().parity() );
      CHECK( "n,Si29" == incident10.pairID().symbol() );

      // particle pair
      const auto pair10 = channel10.particlePair();
      CHECK( 1.008664 == Approx( pair10.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair10.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair10.particle().spin() ) );
      CHECK( +1 == pair10.particle().parity() );
      CHECK( 28.72800 * 1.008664 == Approx( pair10.residual().mass().value ) );
      CHECK( 14.0 * 1.602e-19 == Approx( pair10.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair10.residual().spin() ) );
      CHECK( +1 == pair10.residual().parity() );
      CHECK( "n,Si29" == pair10.pairID().symbol() );

      // quantum numbers
      const auto numbers10 = channel10.quantumNumbers();
      CHECK( 0 == numbers10.orbitalAngularMomentum() );
      CHECK( 0.0 == numbers10.spin() );
      CHECK( 0.0 == numbers10.totalAngularMomentum() );
      CHECK( +1 == numbers10.parity() );
      CHECK( "{0,0,0+}" == numbers10.toString() );

      // radii
      const auto radii10 = channel10.radii();
      CHECK( .44 == Approx( radii10.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .44 == Approx( radii10.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .44 == Approx( radii10.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel10.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel10.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1, resonance table
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      auto table1 = spingroup1.resonanceTable();

      CHECK( 1 == table1.numberChannels() ); // 1 normal channel + 1 eliminated
      CHECK( 7 == table1.numberResonances() );

      auto energies1 = table1.energies();
      CHECK( -2.179600e+6 == Approx( energies1.front().value ) );
      CHECK( 2.248487e+6 == Approx( energies1.back().value ) );

      auto resonances1 = table1.resonances();
      CHECK( -2.179600e+6 == Approx( resonances1.front().energy().value ) );
      CHECK( 2.248487e+6 == Approx( resonances1.back().energy().value ) );
      CHECK( 1 == resonances1.front().widths().size() );
      CHECK( 1 == resonances1.back().widths().size() );
      CHECK( std::sqrt( 1.722200e+6 / 2. / channel10.penetrability( -2.179600e+6 * electronVolt ) )
             == Approx( resonances1.front().widths()[0].value ) );
      CHECK( std::sqrt( 169.32 / 2. / channel10.penetrability( 2.248487e+6 * electronVolt ) )
             == Approx( resonances1.back().widths()[0].value ) );
      CHECK( std::sqrt( 409.08 / 2. ) == Approx( resonances1.front().eliminatedWidth().value ) );
      CHECK( std::sqrt( 3. / 2. ) == Approx( resonances1.back().eliminatedWidth().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 2 -  2 elastic channels + 1 eliminated
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup2 = spingroups[2];

      // channels
      auto channels2 = spingroup2.channels();

      CHECK( 2 == channels2.size() ); // 2 normal channel + 1 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 2, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel20 = std::get< Channel< Neutron > >( channels2[0] );
      CHECK( "n,Si29->n,Si29" == channel20.reactionID().symbol() );

      // incident particle pair
      const auto incident20 = channel20.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident20.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident20.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident20.particle().spin() ) );
      CHECK( +1 == incident20.particle().parity() );
      CHECK( 28.72800 * 1.008664 == Approx( incident20.residual().mass().value ) );
      CHECK( 14.0 * 1.602e-19 == Approx( incident20.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident20.residual().spin() ) );
      CHECK( +1 == incident20.residual().parity() );
      CHECK( "n,Si29" == incident20.pairID().symbol() );

      // particle pair
      const auto pair20 = channel20.particlePair();
      CHECK( 1.008664 == Approx( pair20.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair20.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair20.particle().spin() ) );
      CHECK( +1 == pair20.particle().parity() );
      CHECK( 28.72800 * 1.008664 == Approx( pair20.residual().mass().value ) );
      CHECK( 14.0 * 1.602e-19 == Approx( pair20.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair20.residual().spin() ) );
      CHECK( +1 == pair20.residual().parity() );
      CHECK( "n,Si29" == pair20.pairID().symbol() );

      // quantum numbers
      const auto numbers20 = channel20.quantumNumbers();
      CHECK( 1 == numbers20.orbitalAngularMomentum() );
      CHECK( 0.0 == numbers20.spin() );
      CHECK( 1.0 == numbers20.totalAngularMomentum() );
      CHECK( -1 == numbers20.parity() );
      CHECK( "{1,0,1-}" == numbers20.toString() );

      // radii
      const auto radii20 = channel20.radii();
      CHECK( .44 == Approx( radii20.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .44 == Approx( radii20.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .44 == Approx( radii20.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel20.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel20.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 2, channel 1: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel21 = std::get< Channel< Neutron > >( channels2[1] );
      CHECK( "n,Si29->n,Si29" == channel21.reactionID().symbol() );

      // incident particle pair
      const auto incident21 = channel21.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident21.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident21.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident21.particle().spin() ) );
      CHECK( +1 == incident21.particle().parity() );
      CHECK( 28.72800 * 1.008664 == Approx( incident21.residual().mass().value ) );
      CHECK( 14.0 * 1.602e-19 == Approx( incident21.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident21.residual().spin() ) );
      CHECK( +1 == incident21.residual().parity() );
      CHECK( "n,Si29" == incident21.pairID().symbol() );

      // particle pair
      const auto pair21 = channel21.particlePair();
      CHECK( 1.008664 == Approx( pair21.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair21.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair21.particle().spin() ) );
      CHECK( +1 == pair21.particle().parity() );
      CHECK( 28.72800 * 1.008664 == Approx( pair21.residual().mass().value ) );
      CHECK( 14.0 * 1.602e-19 == Approx( pair21.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair21.residual().spin() ) );
      CHECK( +1 == pair21.residual().parity() );
      CHECK( "n,Si29" == pair21.pairID().symbol() );

      // quantum numbers
      const auto numbers21 = channel21.quantumNumbers();
      CHECK( 1 == numbers21.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers21.spin() );
      CHECK( 1.0 == numbers21.totalAngularMomentum() );
      CHECK( -1 == numbers21.parity() );
      CHECK( "{1,1,1-}" == numbers21.toString() );

      // radii
      const auto radii21 = channel21.radii();
      CHECK( .44 == Approx( radii21.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .44 == Approx( radii21.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .44 == Approx( radii21.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel21.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel21.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 2, resonance table
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      auto table2 = spingroup2.resonanceTable();

      CHECK( 2 == table2.numberChannels() ); // 2 normal channel + 1 eliminated
      CHECK( 10 == table2.numberResonances() );

      auto energies2 = table2.energies();
      CHECK( 1.528200e+4 == Approx( energies2.front().value ) );
      CHECK( 1.192268e+6 == Approx( energies2.back().value ) );

      auto resonances2 = table2.resonances();
      CHECK( 1.528200e+4 == Approx( resonances2.front().energy().value ) );
      CHECK( 1.192268e+6 == Approx( resonances2.back().energy().value ) );
      CHECK( 2 == resonances2.front().widths().size() );
      CHECK( 2 == resonances2.back().widths().size() );
      CHECK( std::sqrt( 10. / 2. / channel20.penetrability( 1.528200e+4 * electronVolt ) )
             == Approx( resonances2.front().widths()[0].value ) );
      CHECK( 0. == Approx( resonances2.front().widths()[1].value ) );
      CHECK( 0. == Approx( resonances2.back().widths()[0].value ) );
      CHECK( std::sqrt( 375.06 / 2. / channel21.penetrability( 1.192268e+6 * electronVolt ) )
             == Approx( resonances2.back().widths()[1].value ) );
      CHECK( std::sqrt( 1.646 / 2. ) == Approx( resonances2.front().eliminatedWidth().value ) );
      CHECK( std::sqrt( .3 / 2. ) == Approx( resonances2.back().eliminatedWidth().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 3 -  2 elastic channels + 1 eliminated
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup3 = spingroups[3];

      // channels
      auto channels3 = spingroup3.channels();

      CHECK( 2 == channels3.size() ); // 2 normal channel + 1 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 3, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel30 = std::get< Channel< Neutron > >( channels3[0] );
      CHECK( "n,Si29->n,Si29" == channel30.reactionID().symbol() );

      // incident particle pair
      const auto incident30 = channel30.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident30.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident30.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident30.particle().spin() ) );
      CHECK( +1 == incident30.particle().parity() );
      CHECK( 28.72800 * 1.008664 == Approx( incident30.residual().mass().value ) );
      CHECK( 14.0 * 1.602e-19 == Approx( incident30.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident30.residual().spin() ) );
      CHECK( +1 == incident30.residual().parity() );
      CHECK( "n,Si29" == incident30.pairID().symbol() );

      // particle pair
      const auto pair30 = channel30.particlePair();
      CHECK( 1.008664 == Approx( pair30.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair30.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair30.particle().spin() ) );
      CHECK( +1 == pair30.particle().parity() );
      CHECK( 28.72800 * 1.008664 == Approx( pair30.residual().mass().value ) );
      CHECK( 14.0 * 1.602e-19 == Approx( pair30.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair30.residual().spin() ) );
      CHECK( +1 == pair30.residual().parity() );
      CHECK( "n,Si29" == pair30.pairID().symbol() );

      // quantum numbers
      const auto numbers30 = channel30.quantumNumbers();
      CHECK( 0 == numbers30.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers30.spin() );
      CHECK( 1.0 == numbers30.totalAngularMomentum() );
      CHECK( +1 == numbers30.parity() );
      CHECK( "{0,1,1+}" == numbers30.toString() );

      // radii
      const auto radii30 = channel30.radii();
      CHECK( .44 == Approx( radii30.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .44 == Approx( radii30.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .44 == Approx( radii30.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel30.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel30.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 3, channel 1: elastic channel
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel31 = std::get< Channel< Neutron > >( channels3[1] );
      CHECK( "n,Si29->n,Si29" == channel31.reactionID().symbol() );

      // incident particle pair
      const auto incident31 = channel31.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident31.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident31.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident31.particle().spin() ) );
      CHECK( +1 == incident31.particle().parity() );
      CHECK( 28.72800 * 1.008664 == Approx( incident31.residual().mass().value ) );
      CHECK( 14.0 * 1.602e-19 == Approx( incident31.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident31.residual().spin() ) );
      CHECK( +1 == incident31.residual().parity() );
      CHECK( "n,Si29" == incident31.pairID().symbol() );

      // particle pair
      const auto pair31 = channel31.particlePair();
      CHECK( 1.008664 == Approx( pair31.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair31.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair31.particle().spin() ) );
      CHECK( +1 == pair31.particle().parity() );
      CHECK( 28.72800 * 1.008664 == Approx( pair31.residual().mass().value ) );
      CHECK( 14.0 * 1.602e-19 == Approx( pair31.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair31.residual().spin() ) );
      CHECK( +1 == pair31.residual().parity() );
      CHECK( "n,Si29" == pair31.pairID().symbol() );

      // quantum numbers
      const auto numbers31 = channel31.quantumNumbers();
      CHECK( 2 == numbers31.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers31.spin() );
      CHECK( 1.0 == numbers31.totalAngularMomentum() );
      CHECK( +1 == numbers31.parity() );
      CHECK( "{2,1,1+}" == numbers31.toString() );

      // radii
      const auto radii31 = channel31.radii();
      CHECK( .44 == Approx( radii31.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .44 == Approx( radii31.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .44 == Approx( radii31.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel31.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel31.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 3, resonance table
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      auto table3 = spingroup3.resonanceTable();

      CHECK( 2 == table3.numberChannels() ); // 2 normal channel + 1 eliminated
      CHECK( 5 == table3.numberResonances() );

      auto energies3 = table3.energies();
      CHECK( 3.857642e+5 == Approx( energies3.front().value ) );
      CHECK( 1.388859e+6 == Approx( energies3.back().value ) );

      auto resonances3 = table3.resonances();
      CHECK( 3.857642e+5 == Approx( resonances3.front().energy().value ) );
      CHECK( 1.388859e+6 == Approx( resonances3.back().energy().value ) );
      CHECK( 2 == resonances3.front().widths().size() );
      CHECK( 2 == resonances3.back().widths().size() );
      CHECK( std::sqrt( 24133. / 2. / channel30.penetrability( 3.857642e+5 * electronVolt ) )
             == Approx( resonances3.front().widths()[0].value ) );
      CHECK( 0. == Approx( resonances3.front().widths()[1].value ) );
      CHECK( 0. == Approx( resonances3.back().widths()[0].value ) );
      CHECK( std::sqrt( 4271.4 / 2. / channel31.penetrability( 1.388859e+6 * electronVolt ) )
             == Approx( resonances3.back().widths()[1].value ) );
      CHECK( std::sqrt( 4.67 / 2. ) == Approx( resonances3.front().eliminatedWidth().value ) );
      CHECK( std::sqrt( 3. / 2. ) == Approx( resonances3.back().eliminatedWidth().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 4 -  empty spin group, 1 elastic channels
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup4 = spingroups[4];

      // channels
      auto channels4 = spingroup4.channels();

      CHECK( 1 == channels4.size() ); // 1 normal channel + 0 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 4, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel40 = std::get< Channel< Neutron > >( channels4[0] );
      CHECK( "n,Si29->n,Si29" == channel40.reactionID().symbol() );

      // incident particle pair
      const auto incident40 = channel40.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident40.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident40.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident40.particle().spin() ) );
      CHECK( +1 == incident40.particle().parity() );
      CHECK( 28.72800 * 1.008664 == Approx( incident40.residual().mass().value ) );
      CHECK( 14.0 * 1.602e-19 == Approx( incident40.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident40.residual().spin() ) );
      CHECK( +1 == incident40.residual().parity() );
      CHECK( "n,Si29" == incident40.pairID().symbol() );

      // particle pair
      const auto pair40 = channel40.particlePair();
      CHECK( 1.008664 == Approx( pair40.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair40.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair40.particle().spin() ) );
      CHECK( +1 == pair40.particle().parity() );
      CHECK( 28.72800 * 1.008664 == Approx( pair40.residual().mass().value ) );
      CHECK( 14.0 * 1.602e-19 == Approx( pair40.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair40.residual().spin() ) );
      CHECK( +1 == pair40.residual().parity() );
      CHECK( "n,Si29" == pair40.pairID().symbol() );

      // quantum numbers
      const auto numbers40 = channel40.quantumNumbers();
      CHECK( 1 == numbers40.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers40.spin() );
      CHECK( 2.0 == numbers40.totalAngularMomentum() );
      CHECK( -1 == numbers40.parity() );
      CHECK( "{1,1,2-}" == numbers40.toString() );

      // radii
      const auto radii40 = channel40.radii();
      CHECK( .44 == Approx( radii40.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .44 == Approx( radii40.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .44 == Approx( radii40.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel40.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel40.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 5, resonance table
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      auto table4 = spingroup4.resonanceTable();

      CHECK( 1 == table4.numberChannels() ); // 1 normal channel + 0 eliminated
      CHECK( 7 == table4.numberResonances() );

      auto energies4 = table4.energies();
      CHECK( 3.881900e+4 == Approx( energies4.front().value ) );
      CHECK( 1.207629e+6 == Approx( energies4.back().value ) );

      auto resonances4 = table4.resonances();
      CHECK( 3.881900e+4 == Approx( resonances4.front().energy().value ) );
      CHECK( 1.207629e+6 == Approx( resonances4.back().energy().value ) );
      CHECK( 1 == resonances4.front().widths().size() );
      CHECK( 1 == resonances4.back().widths().size() );
      CHECK( std::sqrt( 75.926 / 2. / channel40.penetrability( 3.881900e+4 * electronVolt ) )
             == Approx( resonances4.front().widths()[0].value ) );
      CHECK( std::sqrt( 19795. / 2. / channel40.penetrability( 1.207629e+6 * electronVolt ) )
             == Approx( resonances4.back().widths()[0].value ) );
      CHECK( std::sqrt( 2.4 / 2. ) == Approx( resonances4.front().eliminatedWidth().value ) );
      CHECK( std::sqrt( .3 / 2. ) == Approx( resonances4.back().eliminatedWidth().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 5 -  empty spin group, 1 elastic channels
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup5 = spingroups[5];

      // channels
      auto channels5 = spingroup5.channels();

      CHECK( 2 == channels5.size() ); // 1 normal channel + 0 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 5, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel50 = std::get< Channel< Neutron > >( channels5[0] );
      CHECK( "n,Si29->n,Si29" == channel50.reactionID().symbol() );

      // incident particle pair
      const auto incident50 = channel50.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident50.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident50.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident50.particle().spin() ) );
      CHECK( +1 == incident50.particle().parity() );
      CHECK( 28.72800 * 1.008664 == Approx( incident50.residual().mass().value ) );
      CHECK( 14.0 * 1.602e-19 == Approx( incident50.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident50.residual().spin() ) );
      CHECK( +1 == incident50.residual().parity() );
      CHECK( "n,Si29" == incident50.pairID().symbol() );

      // particle pair
      const auto pair50 = channel50.particlePair();
      CHECK( 1.008664 == Approx( pair50.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair50.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair50.particle().spin() ) );
      CHECK( +1 == pair50.particle().parity() );
      CHECK( 28.72800 * 1.008664 == Approx( pair50.residual().mass().value ) );
      CHECK( 14.0 * 1.602e-19 == Approx( pair50.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair50.residual().spin() ) );
      CHECK( +1 == pair50.residual().parity() );
      CHECK( "n,Si29" == pair50.pairID().symbol() );

      // quantum numbers
      const auto numbers50 = channel50.quantumNumbers();
      CHECK( 2 == numbers50.orbitalAngularMomentum() );
      CHECK( 0.0 == numbers50.spin() );
      CHECK( 2.0 == numbers50.totalAngularMomentum() );
      CHECK( +1 == numbers50.parity() );
      CHECK( "{2,0,2+}" == numbers50.toString() );

      // radii
      const auto radii50 = channel50.radii();
      CHECK( .44 == Approx( radii50.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .44 == Approx( radii50.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .44 == Approx( radii50.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel50.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel50.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 5, channel 1: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel51 = std::get< Channel< Neutron > >( channels5[1] );
      CHECK( "n,Si29->n,Si29" == channel51.reactionID().symbol() );

      // incident particle pair
      const auto incident51 = channel51.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident51.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident51.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident51.particle().spin() ) );
      CHECK( +1 == incident51.particle().parity() );
      CHECK( 28.72800 * 1.008664 == Approx( incident51.residual().mass().value ) );
      CHECK( 14.0 * 1.602e-19 == Approx( incident51.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident51.residual().spin() ) );
      CHECK( +1 == incident51.residual().parity() );
      CHECK( "n,Si29" == incident51.pairID().symbol() );

      // particle pair
      const auto pair51 = channel51.particlePair();
      CHECK( 1.008664 == Approx( pair51.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair51.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair51.particle().spin() ) );
      CHECK( +1 == pair51.particle().parity() );
      CHECK( 28.72800 * 1.008664 == Approx( pair51.residual().mass().value ) );
      CHECK( 14.0 * 1.602e-19 == Approx( pair51.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair51.residual().spin() ) );
      CHECK( +1 == pair51.residual().parity() );
      CHECK( "n,Si29" == pair51.pairID().symbol() );

      // quantum numbers
      const auto numbers51 = channel51.quantumNumbers();
      CHECK( 2 == numbers51.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers51.spin() );
      CHECK( 2.0 == numbers51.totalAngularMomentum() );
      CHECK( +1 == numbers51.parity() );
      CHECK( "{2,1,2+}" == numbers51.toString() );

      // radii
      const auto radii51 = channel51.radii();
      CHECK( .44 == Approx( radii51.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .44 == Approx( radii51.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .44 == Approx( radii51.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel51.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel51.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 5, resonance table
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      auto table5 = spingroup5.resonanceTable();

      CHECK( 2 == table5.numberChannels() ); // 1 normal channel + 1 eliminated
      CHECK( 0 == table5.numberResonances() );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 6 -  empty spin group, 1 elastic channels
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup6 = spingroups[6];

      // channels
      auto channels6 = spingroup6.channels();

      CHECK( 1 == channels6.size() ); // 1 normal channel + 0 eliminated

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 6, channel 0: elastic
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const auto channel60 = std::get< Channel< Neutron > >( channels6[0] );
      CHECK( "n,Si29->n,Si29" == channel60.reactionID().symbol() );

      // incident particle pair
      const auto incident60 = channel60.incidentParticlePair();
      CHECK( 1.008664 == Approx( incident60.particle().mass().value ) );
      CHECK( 0.0 == Approx( incident60.particle().charge().value ) );
      CHECK( 0.5 == Approx( incident60.particle().spin() ) );
      CHECK( +1 == incident60.particle().parity() );
      CHECK( 28.72800 * 1.008664 == Approx( incident60.residual().mass().value ) );
      CHECK( 14.0 * 1.602e-19 == Approx( incident60.residual().charge().value ) );
      CHECK( 0.5 == Approx( incident60.residual().spin() ) );
      CHECK( +1 == incident60.residual().parity() );
      CHECK( "n,Si29" == incident60.pairID().symbol() );

      // particle pair
      const auto pair60 = channel60.particlePair();
      CHECK( 1.008664 == Approx( pair60.particle().mass().value ) );
      CHECK( 0.0 == Approx( pair60.particle().charge().value ) );
      CHECK( 0.5 == Approx( pair60.particle().spin() ) );
      CHECK( +1 == pair60.particle().parity() );
      CHECK( 28.72800 * 1.008664 == Approx( pair60.residual().mass().value ) );
      CHECK( 14.0 * 1.602e-19 == Approx( pair60.residual().charge().value ) );
      CHECK( 0.5 == Approx( pair60.residual().spin() ) );
      CHECK( +1 == pair60.residual().parity() );
      CHECK( "n,Si29" == pair60.pairID().symbol() );

      // quantum numbers
      const auto numbers60 = channel60.quantumNumbers();
      CHECK( 2 == numbers60.orbitalAngularMomentum() );
      CHECK( 1.0 == numbers60.spin() );
      CHECK( 3.0 == numbers60.totalAngularMomentum() );
      CHECK( +1 == numbers60.parity() );
      CHECK( "{2,1,3+}" == numbers60.toString() );

      // radii
      const auto radii60 = channel60.radii();
      CHECK( .44 == Approx( radii60.penetrabilityRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .44 == Approx( radii60.shiftFactorRadius( 1e-5 * electronVolt ).value ) );
      CHECK( .44 == Approx( radii60.phaseShiftRadius( 1e-5 * electronVolt ).value ) );

      // boundary conditions
      CHECK( 0. == channel60.boundaryCondition() );

      // Q value
      CHECK( 0.0 == Approx( channel60.Q().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 6, resonance table
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      auto table6 = spingroup6.resonanceTable();

      CHECK( 1 == table6.numberChannels() ); // 1 normal channel + 0 eliminated
      CHECK( 0 == table6.numberResonances() );
    } // THEN

    THEN( "cross sections can be reconstructed" ) {

      // values taken from NJOY2016 PENDF tape for ENDF/B-VIII.0 Si29

      ReactionID elas( "n,Si29->n,Si29" );
      ReactionID capt( "n,Si29->capture" );
      std::map< ReactionID, CrossSection > xs;

      xs = resonances( 1e-5 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 2.584133 == Approx( xs[ elas ].value ) );
      CHECK( 6.033624 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 2.584133 == Approx( xs[ elas ].value ) );
      CHECK( 1.907999 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 2.584133 == Approx( xs[ elas ].value ) );
      CHECK( 6.033624e-1 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-2 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 2.584133 == Approx( xs[ elas ].value ) );
      CHECK( 1.907999e-1 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e-1 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 2.584133 == Approx( xs[ elas ].value ) );
      CHECK( 6.033631e-2 == Approx( xs[ capt ].value ) );

      xs = resonances( 1. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 2.584129 == Approx( xs[ elas ].value ) );
      CHECK( 1.908002e-2 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e+1 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 2.584092 == Approx( xs[ elas ].value ) );
      CHECK( 6.033695e-3 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e+2 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 2.583726 == Approx( xs[ elas ].value ) );
      CHECK( 1.908231e-3 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 2.580076 == Approx( xs[ elas ].value ) );
      CHECK( 6.042010e-4 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 2.544204 == Approx( xs[ elas ].value ) );
      CHECK( 2.087144e-4 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e+5 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 2.224929 == Approx( xs[ elas ].value ) );
      CHECK( 6.188437e-5 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e+6 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 2.126174 == Approx( xs[ elas ].value ) );
      CHECK( 7.969512e-6 == Approx( xs[ capt ].value ) );

      xs = resonances( 1.299999e+6 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 2.886449 == Approx( xs[ elas ].value ) );
      CHECK( 6.120396e-6 == Approx( xs[ capt ].value ) );
    } // THEN
  } // GIVEN
} // SCENARIO

std::string Pu239LRF3() {

  // Pu239 ENDF/B-VIII.0 LRF=3 resonance evaluation

  return
    " 9.423900+4 2.369986+2          0          0          1          09437 2151     \n"
    " 9.423900+4 1.000000+0          0          1          1          09437 2151     \n"
    " 1.000000-5 2.500000+3          1          3          0          19437 2151     \n"
    " 5.000000-1 9.410000-1          1          0          1          49437 2151     \n"
    " 2.369986+2 9.410000-1          0          0       6258       10439437 2151     \n"
    "-1.500200+2 1.000000+0 3.726588-1 4.685986-2 1.352425-1 0.000000+09437 2151     \n"
    "-8.068561+0 1.000000+0 1.474790-4 4.963603-2-1.528359-3 0.000000+09437 2151     \n"
    "-6.908700+0 0.000000+0 1.805854-2 6.012769-2-9.235517-1 3.464891-19437 2151     \n"
    "-2.194400-1 0.000000+0 3.497563-5 2.013091-2 1.443355-2-7.743345-19437 2151     \n"
    "-2.000000-2 0.000000+0 6.59741-11 2.102923-2-2.353050-2 1.000010-79437 2151     \n"
    " 2.956243-1 1.000000+0 7.947046-5 3.982423-2 5.619673-2 0.000000+09437 2151     \n"
    " 7.815800+0 1.000000+0 8.058332-4 3.859953-2-4.473833-2 0.000000+09437 2151     \n"
    " 1.092800+1 1.000000+0 1.807943-3 3.652167-2-1.552867-1 0.000000+09437 2151     \n"
    " 1.189800+1 1.000000+0 9.794152-4 3.800918-2 2.066556-2 0.000000+09437 2151     \n"
    " 1.432900+1 1.000000+0 6.081492-4 2.896772-2 5.873738-2 0.000000+09437 2151     \n"
    " 1.467800+1 1.000000+0 1.922615-3 3.952068-2 3.095395-2 0.000000+09437 2151     \n"
    " 1.541700+1 0.000000+0 2.056203-3 4.054259-2-1.093928-6 7.550000-19437 2151     \n"
    " 1.765700+1 1.000000+0 1.806885-3 3.833472-2-3.623305-2 0.000000+09437 2151     \n"
    " 2.226600+1 1.000000+0 2.592015-3 4.161996-2-6.217418-2 0.000000+09437 2151     \n"
    " 2.393300+1 1.000000+0 8.549952-5 3.697286-2 2.693898-2 0.000000+09437 2151     \n"
    " 2.626900+1 1.000000+0 1.548603-3 3.868245-2 3.930623-2 0.000000+09437 2151     \n"
    " 2.728800+1 1.000000+0 1.498496-4 3.866958-2 2.618306-3 0.000000+09437 2151     \n"
    " 3.232700+1 0.000000+0 8.678823-4 4.182541-2 5.235058-3-1.279000-19437 2151     \n"
    " 3.548600+1 1.000000+0 2.770046-4 3.954051-2 3.231034-3 0.000000+09437 2151     \n"
    " 4.145700+1 1.000000+0 3.314222-3 4.905718-2-6.496140-3 0.000000+09437 2151     \n"
    " 4.173600+1 1.000000+0 9.388960-4 3.794632-2 6.238103-2 0.000000+09437 2151     \n"
    " 4.453100+1 1.000000+0 6.144083-3 4.052580-2-4.514990-3 0.000000+09437 2151     \n"
    " 4.753400+1 0.000000+0 5.171861-3 2.938826-2 5.548812-1-1.274000-79437 2151     \n"
    " 4.957600+1 0.000000+0 4.339861-3 4.448069-2-1.116076+0 4.983000-39437 2151     \n"
    " 5.014400+1 1.000000+0 3.184552-3 2.177446-2-4.395953-3 0.000000+09437 2151     \n"
    " 5.264800+1 1.000000+0 9.474846-3 4.462903-2-9.570580-3 0.000000+09437 2151     \n"
    " 5.570400+1 1.000000+0 1.497682-3 3.753685-2 2.490887-2 0.000000+09437 2151     \n"
    " 5.692400+1 0.000000+0 1.103568-2 5.672256-2-1.948257+0-5.051426-19437 2151     \n"
    " 5.929100+1 1.000000+0 4.578011-3 3.388975-2 1.055869-1 0.000000+09437 2151     \n"
    " 6.162100+1 0.000000+0 2.451508-2 3.898194-2 7.859031+0-3.320667-39437 2151     \n"
    " 6.317000+1 1.000000+0 6.250878-4 3.858910-2 7.765319-2 0.000000+09437 2151     \n"
    " 6.545400+1 0.000000+0 4.060900-3 3.720344-2 4.602538-1-2.139877-99437 2151     \n"
    " 6.579300+1 1.000000+0 9.292900-3 5.970370-2 9.664072-2 0.000000+09437 2151     \n"
    " 7.416700+1 1.000000+0 3.151923-3 3.175752-2-2.862053-2 0.000000+09437 2151     \n"
    " 7.488500+1 0.000000+0 1.662930-3 3.142321-2 5.836330-1 1.693882-19437 2151     \n"
    " 7.503400+1 1.000000+0 2.103244-2 4.103538-2-9.277256-2 0.000000+09437 2151     \n"
    " 7.908500+1 1.000000+0 5.045545-5 6.344780-2 2.312717-2 0.000000+09437 2151     \n"
    " 8.086779+1 0.000000+0 4.376256-3 4.200000-2 1.802467+0-8.000096-49437 2151     \n"
    " 8.277396+1 1.000000+0 3.313543-4 4.685000-2 7.649002-4 0.000000+09437 2151     \n"
    " 8.536073+1 0.000000+0 4.908061-2 5.100000-2-1.850222+0 1.557989-39437 2151     \n"
    " 8.562257+1 1.000000+0 1.146209-2 2.618000-2 4.951489-3 0.000000+09437 2151     \n"
    " 9.085484+1 1.000000+0 1.216419-2 3.399000-2 7.105036-3 0.000000+09437 2151     \n"
    " 9.307943+1 1.000000+0 6.859974-4 3.764000-2-2.866529-3 0.000000+09437 2151     \n"
    " 9.550517+1 1.000000+0 2.063776-3 6.316000-2-3.150635-2 0.000000+09437 2151     \n"
    " 9.623591+1 0.000000+0 6.640183-3 4.200000-2 1.499218+0-6.429247-19437 2151     \n"
    " 1.002558+2 0.000000+0 3.033674-2 4.200000-2 8.692401-4 1.157935+19437 2151     \n"
    " 1.031473+2 1.000000+0 1.613360-3 3.901790-2 7.138305-3 0.000000+09437 2151     \n"
    " 1.054592+2 1.000000+0 4.878843-3 3.695720-2 3.928122-3 0.000000+09437 2151     \n"
    " 1.068318+2 1.000000+0 1.049666-2 4.018590-2-2.037800-2 0.000000+09437 2151     \n"
    " 1.105574+2 1.000000+0 5.098524-4 3.900000-2 1.432129-2 0.000000+09437 2151     \n"
    " 1.140747+2 0.000000+0 2.472140-3 3.900000-2-2.806990+0 3.209759-19437 2151     \n"
    " 1.153194+2 1.000000+0 5.265794-4 3.900000-2 6.675931-1 0.000000+09437 2151     \n"
    " 1.161893+2 0.000000+0 1.021982-2 3.900000-2-1.276005-1-2.553139-29437 2151     \n"
    " 1.190095+2 1.000000+0 1.816584-2 4.050010-2-3.685477-2 0.000000+09437 2151     \n"
    " 1.211816+2 1.000000+0 2.668518-3 3.475467-2 2.517246-2 0.000000+09437 2151     \n"
    " 1.236246+2 1.000000+0 5.732016-4 5.448998-2-5.410530-2 0.000000+09437 2151     \n"
    " 1.263933+2 1.000000+0 1.733803-3 4.433934-2-9.743127-3 0.000000+09437 2151     \n"
    " 1.277214+2 1.000000+0 5.429897-4 3.807883-2 1.062898-2 0.000000+09437 2151     \n"
    " 1.315982+2 0.000000+0 3.731837-2 4.380129-2 4.207085+0-5.409000-39437 2151     \n"
    " 1.339861+2 1.000000+0 5.973919-3 3.706473-2-3.237859-3 0.000000+09437 2151     \n"
    " 1.369355+2 0.000000+0 1.152383-2 3.349058-2 1.400187-4-9.181000-29437 2151     \n"
    " 1.394440+2 1.000000+0 9.430501-5 4.046013-2-2.338238-2 0.000000+09437 2151     \n"
    " 1.431607+2 1.000000+0 3.341474-3 2.919127-2 5.135081-2 0.000000+09437 2151     \n"
    " 1.436778+2 1.000000+0 4.504148-3 3.703204-2 2.938467-2 0.000000+09437 2151     \n"
    " 1.464736+2 1.000000+0 7.532541-3 3.943747-2 8.051087-3 0.000000+09437 2151     \n"
    " 1.473577+2 0.000000+0 5.194162-3 4.010731-2 2.460544+0 1.586000-39437 2151     \n"
    " 1.484448+2 1.000000+0 3.229834-4 4.105206-2 2.605463-2 0.000000+09437 2151     \n"
    " 1.496479+2 1.000000+0 1.606244-3 3.747636-2 2.653004-2 0.000000+09437 2151     \n"
    " 1.571648+2 0.000000+0 3.060764-2 3.601941-2 2.070989-1 5.684000-19437 2151     \n"
    " 1.571674+2 1.000000+0 1.386961-3 3.918334-2 1.345284-2 0.000000+09437 2151     \n"
    " 1.621466+2 1.000000+0 1.120623-4 4.025781-2 1.946938-2 0.000000+09437 2151     \n"
    " 1.647110+2 1.000000+0 2.826146-2 4.456044-2-7.874520-3 0.000000+09437 2151     \n"
    " 1.672877+2 1.000000+0 6.252791-3 4.167247-2 6.804219-2 0.000000+09437 2151     \n"
    " 1.700732+2 1.000000+0 1.289689-4 3.764917-2 2.943757-3 0.000000+09437 2151     \n"
    " 1.706638+2 1.000000+0 4.443652-4 3.683407-2 4.745634-2 0.000000+09437 2151     \n"
    " 1.709384+2 0.000000+0 3.604519-3 4.350615-2-2.345839-2 1.600000+09437 2151     \n"
    " 1.713205+2 1.000000+0 2.345664-5 3.909383-2 2.931810-3 0.000000+09437 2151     \n"
    " 1.761759+2 1.000000+0 2.205127-3 4.739790-2 3.159702-2 0.000000+09437 2151     \n"
    " 1.774081+2 1.000000+0 3.704962-3 4.074615-2 6.314532-3 0.000000+09437 2151     \n"
    " 1.790932+2 1.000000+0 1.235103-3 3.635043-2 1.009297-2 0.000000+09437 2151     \n"
    " 1.838349+2 1.000000+0 1.602775-3 3.813775-2 1.505707-2 0.000000+09437 2151     \n"
    " 1.846230+2 0.000000+0 2.040746-2 4.277047-2-2.794653+0 6.720000-29437 2151     \n"
    " 1.884706+2 1.000000+0 6.479084-4 4.169582-2 1.365360-2 0.000000+09437 2151     \n"
    " 1.908416+2 1.000000+0 1.710176-3 4.384543-2-1.393745-2 0.000000+09437 2151     \n"
    " 1.924995+2 0.000000+0 1.168135-3 3.939159-2 2.795820-2 3.018000+09437 2151     \n"
    " 1.953805+2 0.000000+0 5.310713-2 2.229571-2-6.798810-1 3.389000-49437 2151     \n"
    " 1.969217+2 1.000000+0 3.997362-3 3.745528-2 1.560094-2 0.000000+09437 2151     \n"
    " 1.996114+2 1.000000+0 9.175631-3 4.432589-2 7.307720-2 0.000000+09437 2151     \n"
    " 2.036364+2 1.000000+0 1.487327-3 3.239595-2 1.529437-2 0.000000+09437 2151     \n"
    " 2.040760+2 0.000000+0 6.599902-2 5.250208-2 1.528452-1 3.372000-19437 2151     \n"
    " 2.075994+2 1.000000+0 7.342238-3 4.434600-2 6.509585-3 0.000000+09437 2151     \n"
    " 2.109839+2 0.000000+0 3.274738-3 3.907835-2-1.422234-1 1.304000+09437 2151     \n"
    " 2.135009+2 1.000000+0 3.550441-4 4.107989-2 3.640293-2 0.000000+09437 2151     \n"
    " 2.137371+2 0.000000+0 9.024432-3 3.943158-2 1.164560+1-2.859000-19437 2151     \n"
    " 2.154522+2 1.000000+0 3.028941-5 4.160829-2 3.112646-2 0.000000+09437 2151     \n"
    " 2.167739+2 1.000000+0 6.387411-3 3.567796-2 5.411154-3 0.000000+09437 2151     \n"
    " 2.197130+2 1.000000+0 3.612230-3 3.774929-2 2.188833-2 0.000000+09437 2151     \n"
    " 2.204710+2 1.000000+0 7.778440-3 3.827461-2-1.163443-2 0.000000+09437 2151     \n"
    " 2.234077+2 1.000000+0 3.367105-3 3.546319-2 5.132089-3 0.000000+09437 2151     \n"
    " 2.251327+2 1.000000+0 1.789707-3 4.765063-2 1.719275-2 0.000000+09437 2151     \n"
    " 2.262561+2 0.000000+0 1.247322-2 4.150060-2-2.338262+0-3.405000+09437 2151     \n"
    " 2.281417+2 1.000000+0 1.840713-3 4.090050-2-1.988443-2 0.000000+09437 2151     \n"
    " 2.302855+2 0.000000+0 3.799153-3 3.956492-2 3.000851+0-2.844000-29437 2151     \n"
    " 2.316503+2 1.000000+0 1.281097-2 3.779987-2 6.481896-3 0.000000+09437 2151     \n"
    " 2.328542+2 1.000000+0 3.840981-4 3.347397-2 2.318010-2 0.000000+09437 2151     \n"
    " 2.340550+2 0.000000+0 9.456099-3 3.964911-2 4.633585+0-1.852000-19437 2151     \n"
    " 2.345686+2 1.000000+0 1.031779-2 3.360264-2 7.018691-3 0.000000+09437 2151     \n"
    " 2.393066+2 1.000000+0 5.647113-3 3.857938-2-1.145553-2 0.000000+09437 2151     \n"
    " 2.407957+2 1.000000+0 3.870289-5 3.966599-2 4.061064-2 0.000000+09437 2151     \n"
    " 2.431433+2 1.000000+0 6.732343-3 3.593803-2-6.022179-2 0.000000+09437 2151     \n"
    " 2.477775+2 1.000000+0 7.878040-4 3.878265-2 1.280003-1 0.000000+09437 2151     \n"
    " 2.491250+2 1.000000+0 1.509038-2 3.627066-2-4.951815-3 0.000000+09437 2151     \n"
    " 2.514959+2 1.000000+0 2.957727-2 3.337172-2 1.016623-2 0.000000+09437 2151     \n"
    " 2.548502+2 1.000000+0 2.995218-3 3.966213-2 3.211628-2 0.000000+09437 2151     \n"
    " 2.562022+2 0.000000+0 1.133117-3 4.009078-2 1.173512+0 9.896000-49437 2151     \n"
    " 2.563798+2 1.000000+0 5.781383-3 4.063503-2-1.152601-2 0.000000+09437 2151     \n"
    " 2.593035+2 1.000000+0 2.598063-4 4.121986-2-4.070325-2 0.000000+09437 2151     \n"
    " 2.622701+2 0.000000+0 8.387607-2 4.019009-2-2.669389+0 2.877000+09437 2151     \n"
    " 2.630155+2 1.000000+0 2.330800-3 3.913850-2-1.063069-2 0.000000+09437 2151     \n"
    " 2.646131+2 1.000000+0 1.501117-4 3.556779-2 2.380099-2 0.000000+09437 2151     \n"
    " 2.693427+2 1.000000+0 1.089763-3 3.571531-2-5.384277-2 0.000000+09437 2151     \n"
    " 2.698319+2 1.000000+0 4.231678-3 3.143269-2 2.591366-2 0.000000+09437 2151     \n"
    " 2.729132+2 1.000000+0 2.641724-2 3.938930-2-4.291288-2 0.000000+09437 2151     \n"
    " 2.751073+2 0.000000+0 3.983328-2 4.340210-2 6.853877-2 9.469000-19437 2151     \n"
    " 2.758839+2 1.000000+0 2.190006-2 3.998626-2 5.387438-2 0.000000+09437 2151     \n"
    " 2.761631+2 0.000000+0 3.064318-2 3.931186-2 9.661600+0 1.754000-29437 2151     \n"
    " 2.798595+2 1.000000+0 7.949437-3 3.090389-2-2.566541-2 0.000000+09437 2151     \n"
    " 2.832207+2 1.000000+0 2.761073-2 3.860662-2 5.968918-3 0.000000+09437 2151     \n"
    " 2.860432+2 1.000000+0 1.821999-4 4.025295-2 8.731618-3 0.000000+09437 2151     \n"
    " 2.883364+2 1.000000+0 8.369419-5 3.952563-2-2.873943-2 0.000000+09437 2151     \n"
    " 2.908393+2 0.000000+0 2.465293-2 3.681574-2-5.622453+0 1.590000-29437 2151     \n"
    " 2.926630+2 1.000000+0 3.581342-3 3.176599-2 2.903242-2 0.000000+09437 2151     \n"
    " 2.927490+2 0.000000+0 4.961382-4 3.844300-2 1.065324-4-6.108000-19437 2151     \n"
    " 2.968010+2 1.000000+0 3.456354-3 4.347556-2-2.963552-2 0.000000+09437 2151     \n"
    " 2.989242+2 1.000000+0 1.057251-2 4.831328-2 2.653474-2 0.000000+09437 2151     \n"
    " 3.021599+2 1.000000+0 1.954505-2 4.360863-2-4.492502-2 0.000000+09437 2151     \n"
    " 3.085677+2 1.000000+0 3.212258-3 3.185492-2 1.100244-1 0.000000+09437 2151     \n"
    " 3.093524+2 1.000000+0 1.404507-2 4.918974-2 2.313565-2 0.000000+09437 2151     \n"
    " 3.105521+2 0.000000+0 2.760651-3 3.975185-2 4.002067+0 1.006000-39437 2151     \n"
    " 3.114854+2 1.000000+0 2.784909-4 3.254685-2 3.647669-3 0.000000+09437 2151     \n"
    " 3.139725+2 1.000000+0 1.485417-2 3.573590-2-9.039714-3 0.000000+09437 2151     \n"
    " 3.170193+2 1.000000+0 5.824822-3 3.782196-2 2.015555-2 0.000000+09437 2151     \n"
    " 3.237350+2 1.000000+0 1.931797-2 8.143225-2-2.637698-2 0.000000+09437 2151     \n"
    " 3.241357+2 0.000000+0 1.652905-2 4.039980-2-3.246838+0-9.607000-19437 2151     \n"
    " 3.256704+2 1.000000+0 8.390583-3 2.796051-2-1.665561-2 0.000000+09437 2151     \n"
    " 3.300431+2 0.000000+0 6.143988-3 4.257778-2 1.463603+0-6.337000-19437 2151     \n"
    " 3.343006+2 1.000000+0 6.144359-3 6.126662-2 1.232670-2 0.000000+09437 2151     \n"
    " 3.363112+2 1.000000+0 1.534526-2 4.309744-2-1.525378-2 0.000000+09437 2151     \n"
    " 3.383400+2 1.000000+0 8.548506-3 5.866363-2 1.234824-2 0.000000+09437 2151     \n"
    " 3.396264+2 1.000000+0 3.719531-3 3.644792-2-1.475337-2 0.000000+09437 2151     \n"
    " 3.435826+2 1.000000+0 1.667234-2 3.477121-2 1.981101-2 0.000000+09437 2151     \n"
    " 3.465818+2 0.000000+0 1.094014-2 3.751078-2-1.188357-3-1.204000+09437 2151     \n"
    " 3.507091+2 1.000000+0 2.217382-2 3.652749-2 2.926571-2 0.000000+09437 2151     \n"
    " 3.532109+2 1.000000+0 3.939216-3 3.803297-2 1.710778-2 0.000000+09437 2151     \n"
    " 3.553274+2 1.000000+0 3.665337-4 4.112377-2 1.559974-2 0.000000+09437 2151     \n"
    " 3.604106+2 1.000000+0 1.195691-3 5.284381-2-1.887661-2 0.000000+09437 2151     \n"
    " 3.616462+2 0.000000+0 3.287063-3 3.919859-2-7.974583-1 5.595000-19437 2151     \n"
    " 3.617950+2 1.000000+0 2.068719-4 3.967719-2 2.867348-3 0.000000+09437 2151     \n"
    " 3.658672+2 0.000000+0 6.246070-2 4.428747-2-2.294659+1 2.723000-19437 2151     \n"
    " 3.686572+2 1.000000+0 2.601594-4 4.040011-2 1.403694-2 0.000000+09437 2151     \n"
    " 3.702134+2 0.000000+0 8.159673-4 4.060435-2 1.116570-3 2.023000-19437 2151     \n"
    " 3.707489+2 1.000000+0 2.986568-3 5.416670-2 1.371487-2 0.000000+09437 2151     \n"
    " 3.737744+2 0.000000+0 2.644426-3 4.374821-2 2.416770+0 1.977000+09437 2151     \n"
    " 3.754594+2 1.000000+0 2.893906-3 3.746791-2 6.968103-3 0.000000+09437 2151     \n"
    " 3.775271+2 1.000000+0 1.829649-3 3.925996-2-2.214816-2 0.000000+09437 2151     \n"
    " 3.784878+2 1.000000+0 1.081091-3 5.074818-2 7.526413-3 0.000000+09437 2151     \n"
    " 3.828254+2 1.000000+0 3.556091-4 3.973070-2 7.862884-3 0.000000+09437 2151     \n"
    " 3.847135+2 1.000000+0 6.087171-3 2.802594-2-2.932436-2 0.000000+09437 2151     \n"
    " 3.870952+2 1.000000+0 2.052078-4 4.373636-2 1.420371-2 0.000000+09437 2151     \n"
    " 3.877882+2 0.000000+0 3.121119-3 4.177090-2 5.677094+0-8.237000-19437 2151     \n"
    " 3.899771+2 1.000000+0 1.486202-3 3.601872-2-1.300681-2 0.000000+09437 2151     \n"
    " 3.919691+2 1.000000+0 1.033213-3 3.862653-2 2.077667-2 0.000000+09437 2151     \n"
    " 3.948903+2 1.000000+0 6.979024-3 4.859272-2-4.153027-2 0.000000+09437 2151     \n"
    " 3.973586+2 1.000000+0 1.930823-3 3.395739-2 4.277199-2 0.000000+09437 2151     \n"
    " 4.020717+2 1.000000+0 1.961485-2 3.660109-2-1.582770-1 0.000000+09437 2151     \n"
    " 4.046630+2 1.000000+0 2.283186-2 3.707319-2-6.766713-2 0.000000+09437 2151     \n"
    " 4.065219+2 1.000000+0 1.228077-3 3.799614-2-6.881958-3 0.000000+09437 2151     \n"
    " 4.075087+2 1.000000+0 8.261139-4 3.587424-2 1.851917-1 0.000000+09437 2151     \n"
    " 4.092060+2 1.000000+0 1.168092-3 4.107570-2-1.612711-2 0.000000+09437 2151     \n"
    " 4.128026+2 1.000000+0 9.268847-3 3.540607-2 3.297267-2 0.000000+09437 2151     \n"
    " 4.161594+2 1.000000+0 3.000471-3 4.063137-2-4.564097-3 0.000000+09437 2151     \n"
    " 4.181058+2 1.000000+0 1.276656-3 4.513880-2-1.277358-2 0.000000+09437 2151     \n"
    " 4.203287+2 1.000000+0 5.601041-3 3.614157-2-2.583655-2 0.000000+09437 2151     \n"
    " 4.205241+2 0.000000+0 8.041953-3 3.846921-2 1.967152-1 2.204000+09437 2151     \n"
    " 4.261147+2 1.000000+0 2.545217-4 4.071665-2 1.543716-2 0.000000+09437 2151     \n"
    " 4.264019+2 0.000000+0 2.051742-2 4.147834-2 1.476389-1 7.585000+09437 2151     \n"
    " 4.301462+2 1.000000+0 3.830818-3 4.353171-2 6.306675-1 0.000000+09437 2151     \n"
    " 4.333524+2 1.000000+0 6.099460-4 3.994403-2-3.433602-3 0.000000+09437 2151     \n"
    " 4.335373+2 0.000000+0 1.709126-2 4.072026-2 2.919823+0-1.630000+09437 2151     \n"
    " 4.382523+2 1.000000+0 2.872592-3 3.756978-2-6.524479-3 0.000000+09437 2151     \n"
    " 4.393027+2 1.000000+0 3.211845-3 2.377536-2 2.830593-3 0.000000+09437 2151     \n"
    " 4.406952+2 1.000000+0 3.375108-4 3.531040-2-1.043885-1 0.000000+09437 2151     \n"
    " 4.429984+2 0.000000+0 2.246280-2 4.269684-2-1.668125-1 2.798000-19437 2151     \n"
    " 4.471214+2 1.000000+0 1.371618-4 3.566869-2 4.074136-2 0.000000+09437 2151     \n"
    " 4.503265+2 1.000000+0 9.739199-4 3.880860-2-3.673654-2 0.000000+09437 2151     \n"
    " 4.518799+2 1.000000+0 1.463614-2 3.327980-2 3.196286-3 0.000000+09437 2151     \n"
    " 4.552305+2 1.000000+0 8.430497-4 4.127246-2-2.170071-1 0.000000+09437 2151     \n"
    " 4.562900+2 0.000000+0 7.923493-2 3.973328-2 1.848849-1-1.557000-19437 2151     \n"
    " 4.578712+2 1.000000+0 7.158732-3 5.290556-2-1.511008-1 0.000000+09437 2151     \n"
    " 4.593781+2 1.000000+0 4.501516-3 3.363276-2 2.746769-2 0.000000+09437 2151     \n"
    " 4.618091+2 1.000000+0 2.100687-3 3.567143-2-7.658610-2 0.000000+09437 2151     \n"
    " 4.632010+2 1.000000+0 6.042520-4 3.542058-2 1.214446-2 0.000000+09437 2151     \n"
    " 4.663832+2 1.000000+0 1.839130-4 3.901289-2 1.835137-1 0.000000+09437 2151     \n"
    " 4.686568+2 0.000000+0 1.977734-2 4.352881-2-3.274579+0 6.616000-39437 2151     \n"
    " 4.736615+2 1.000000+0 4.440218-3 3.530062-2-2.914237-3 0.000000+09437 2151     \n"
    " 4.756315+2 0.000000+0 3.366957-2 4.301155-2 3.911347-3 3.619000+09437 2151     \n"
    " 4.769055+2 0.000000+0 5.775075-3 4.213734-2 3.664827+0 7.732000+09437 2151     \n"
    " 4.797173+2 1.000000+0 2.081310-4 4.044506-2 2.085881-2 0.000000+09437 2151     \n"
    " 4.847187+2 1.000000+0 2.620010-3 3.861902-2-6.742402-3 0.000000+09437 2151     \n"
    " 4.883689+2 0.000000+0 2.005321-2 4.331299-2-1.492059-2 7.037000-19437 2151     \n"
    " 4.918849+2 0.000000+0 4.235353-2 4.335213-2 1.457965+0 1.339000+09437 2151     \n"
    " 4.946964+2 1.000000+0 5.195154-3 4.506322-2-8.330242-2 0.000000+09437 2151     \n"
    " 4.962332+2 1.000000+0 6.833097-4 3.975851-2 8.023701-3 0.000000+09437 2151     \n"
    " 5.011171+2 1.000000+0 3.829380-3 3.390470-2 5.134930-2 0.000000+09437 2151     \n"
    " 5.034406+2 1.000000+0 1.316134-2 8.365237-2 7.676739-2 0.000000+09437 2151     \n"
    " 5.064420+2 1.000000+0 3.250601-4 4.713329-2-1.560367-1 0.000000+09437 2151     \n"
    " 5.101366+2 0.000000+0 4.772952-5 3.813398-2 3.622382-1 4.306000+09437 2151     \n"
    " 5.103755+2 0.000000+0 1.182814-1 4.308014-2 1.387333-1 1.989000+09437 2151     \n"
    " 5.103792+2 1.000000+0 3.356940-2 4.209886-2-2.902579-1 0.000000+09437 2151     \n"
    " 5.156943+2 1.000000+0 4.804840-4 4.036645-2-1.141041-2 0.000000+09437 2151     \n"
    " 5.173710+2 1.000000+0 6.516255-5 3.903992-2 5.897690-2 0.000000+09437 2151     \n"
    " 5.186573+2 1.000000+0 6.522346-4 4.149632-2-1.115439-1 0.000000+09437 2151     \n"
    " 5.208523+2 1.000000+0 1.571244-2 4.378416-2-2.644072-2 0.000000+09437 2151     \n"
    " 5.248497+2 1.000000+0 3.906948-2 2.738258-2-1.062748-2 0.000000+09437 2151     \n"
    " 5.251023+2 0.000000+0 2.624918-1 4.346788-2 1.206997+1-4.821000-39437 2151     \n"
    " 5.267096+2 1.000000+0 1.464488-3 3.904258-2 2.898874-3 0.000000+09437 2151     \n"
    " 5.280369+2 1.000000+0 1.622745-3 4.180917-2 2.474641-2 0.000000+09437 2151     \n"
    " 5.291459+2 1.000000+0 5.625398-4 4.069623-2 3.142131-3 0.000000+09437 2151     \n"
    " 5.309268+2 0.000000+0 2.816228-5 4.269371-2-1.672814+0 1.496000-19437 2151     \n"
    " 5.311855+2 1.000000+0 5.696320-2 4.038763-2-2.988686-2 0.000000+09437 2151     \n"
    " 5.398403+2 1.000000+0 1.283026-2 2.387914-2 2.930494-3 0.000000+09437 2151     \n"
    " 5.412448+2 1.000000+0 1.821610-3 3.710530-2-3.575465-2 0.000000+09437 2151     \n"
    " 5.422341+2 1.000000+0 4.727604-3 4.040394-2 5.233756-1 0.000000+09437 2151     \n"
    " 5.437541+2 1.000000+0 1.234252-2 7.566264-2-1.129908-2 0.000000+09437 2151     \n"
    " 5.464036+2 0.000000+0 3.749981-2 4.073720-2-1.563868-1 1.316000+09437 2151     \n"
    " 5.503441+2 1.000000+0 1.201780-2 5.236334-2 8.641241-3 0.000000+09437 2151     \n"
    " 5.542487+2 1.000000+0 1.264631-2 3.816429-2 2.763424-3 0.000000+09437 2151     \n"
    " 5.548716+2 0.000000+0 1.146108-1 4.283265-2 1.397478+0-8.340000-29437 2151     \n"
    " 5.562775+2 1.000000+0 1.534854-3 4.155845-2-3.421197-3 0.000000+09437 2151     \n"
    " 5.598431+2 1.000000+0 2.772224-2 5.408135-2 8.159544-3 0.000000+09437 2151     \n"
    " 5.631675+2 0.000000+0 4.588218-2 4.013043-2-3.010188-2-1.424000+09437 2151     \n"
    " 5.635616+2 1.000000+0 2.216000-2 4.407004-2-6.704767-2 0.000000+09437 2151     \n"
    " 5.646892+2 1.000000+0 8.081678-3 4.508040-2 2.750859-3 0.000000+09437 2151     \n"
    " 5.665009+2 1.000000+0 9.856418-3 4.280396-2-1.767463-3 0.000000+09437 2151     \n"
    " 5.717909+2 1.000000+0 8.870227-3 3.899157-2-3.710186-2 0.000000+09437 2151     \n"
    " 5.738044+2 0.000000+0 8.856669-3 4.079210-2 8.763332-4 1.221000+09437 2151     \n"
    " 5.746812+2 0.000000+0 1.522849-1 5.311546-2 1.786412-1 3.834000-49437 2151     \n"
    " 5.765348+2 1.000000+0 4.203656-2 4.097592-2 7.704668-3 0.000000+09437 2151     \n"
    " 5.787065+2 1.000000+0 1.710701-3 3.111873-2 4.321328-3 0.000000+09437 2151     \n"
    " 5.797690+2 1.000000+0 6.380274-3 3.523548-2 6.884277-3 0.000000+09437 2151     \n"
    " 5.855314+2 1.000000+0 6.935521-4 4.261583-2-3.175822-2 0.000000+09437 2151     \n"
    " 5.888298+2 1.000000+0 1.233003-2 6.964469-2 1.167238-2 0.000000+09437 2151     \n"
    " 5.905356+2 1.000000+0 2.358120-4 3.756981-2 1.591305-1 0.000000+09437 2151     \n"
    " 5.943140+2 1.000000+0 1.865545-3 3.204129-2 1.300247-2 0.000000+09437 2151     \n"
    " 5.980233+2 0.000000+0 2.714507-2 4.079652-2 1.364107+0 4.519000+09437 2151     \n"
    " 5.980618+2 1.000000+0 8.579505-3 3.115084-2-4.829549-3 0.000000+09437 2151     \n"
    " 6.047507+2 1.000000+0 2.440448-2 2.968405-2-2.946350-3 0.000000+09437 2151     \n"
    " 6.084163+2 1.000000+0 9.747669-3 3.135006-2-8.250941-3 0.000000+09437 2151     \n"
    " 6.100461+2 1.000000+0 1.627087-2 3.454292-2-5.398584-3 0.000000+09437 2151     \n"
    " 6.135804+2 1.000000+0 5.730294-3 4.875523-2-1.399719-2 0.000000+09437 2151     \n"
    " 6.216027+2 1.000000+0 1.165588-2 3.741849-2-5.683957-3 0.000000+09437 2151     \n"
    " 6.233525+2 1.000000+0 8.990367-3 6.051346-2 1.009166-2 0.000000+09437 2151     \n"
    " 6.259406+2 1.000000+0 7.857106-3 4.082782-2-1.083304-2 0.000000+09437 2151     \n"
    " 6.289611+2 1.000000+0 1.606048-3 4.114645-2 8.584463-3 0.000000+09437 2151     \n"
    " 6.332451+2 0.000000+0 7.142701-2 4.792362-2 2.921175-3-4.258000+09437 2151     \n"
    " 6.372534+2 1.000000+0 4.707769-3 5.216232-2 8.263370-3 0.000000+09437 2151     \n"
    " 6.400480+2 1.000000+0 8.507634-3 3.110222-2 8.872551-3 0.000000+09437 2151     \n"
    " 6.425418+2 1.000000+0 2.188865-4 4.241774-2 9.493947-3 0.000000+09437 2151     \n"
    " 6.431958+2 0.000000+0 6.206144-3 4.419358-2 5.446254-3 3.398000+09437 2151     \n"
    " 6.457371+2 1.000000+0 7.035179-3 6.397754-2 8.689677-3 0.000000+09437 2151     \n"
    " 6.470887+2 0.000000+0 3.032679-3 4.265360-2 1.866054+0 1.006000-19437 2151     \n"
    " 6.472857+2 1.000000+0 4.473188-4 4.368392-2 1.069871-2 0.000000+09437 2151     \n"
    " 6.512968+2 1.000000+0 2.586034-4 3.800996-2 9.747060-3 0.000000+09437 2151     \n"
    " 6.532892+2 1.000000+0 1.923674-4 3.894073-2 1.044627-2 0.000000+09437 2151     \n"
    " 6.587032+2 0.000000+0 4.653284-2 3.677751-2-3.732910-2-4.353000-69437 2151     \n"
    " 6.593726+2 1.000000+0 5.150331-2 7.703811-2-1.551313-2 0.000000+09437 2151     \n"
    " 6.626196+2 1.000000+0 3.246936-4 3.107767-2 1.793063-2 0.000000+09437 2151     \n"
    " 6.671783+2 1.000000+0 2.236421-3 4.731737-2-4.126501-3 0.000000+09437 2151     \n"
    " 6.699020+2 1.000000+0 3.456058-3 4.195664-2 3.651246-3 0.000000+09437 2151     \n"
    " 6.721552+2 1.000000+0 1.550818-2 2.828754-2 3.630133-3 0.000000+09437 2151     \n"
    " 6.730107+2 0.000000+0 4.700499-3 3.939835-2 6.901635-5-5.076000+09437 2151     \n"
    " 6.734819+2 1.000000+0 4.843144-4 3.700536-2 1.361054-2 0.000000+09437 2151     \n"
    " 6.751003+2 1.000000+0 2.157156-3 3.916991-2 1.930566-2 0.000000+09437 2151     \n"
    " 6.763804+2 1.000000+0 8.782993-4 4.146520-2 4.192667-2 0.000000+09437 2151     \n"
    " 6.821301+2 1.000000+0 1.536255-3 3.863383-2-1.091844-2 0.000000+09437 2151     \n"
    " 6.828396+2 0.000000+0 3.207899-2 4.731281-2-3.714050+0 1.000000-79437 2151     \n"
    " 6.845661+2 0.000000+0 2.015138-2 4.528294-2 2.775889-2-1.822000+09437 2151     \n"
    " 6.846920+2 1.000000+0 6.831326-3 4.037643-2 4.289992-4 0.000000+09437 2151     \n"
    " 6.866216+2 1.000000+0 7.557227-4 3.856379-2 3.395549-3 0.000000+09437 2151     \n"
    " 6.885271+2 1.000000+0 1.404854-2 2.496324-2 9.645329-3 0.000000+09437 2151     \n"
    " 6.927008+2 1.000000+0 4.919868-3 4.719673-2 8.829941-2 0.000000+09437 2151     \n"
    " 6.939349+2 0.000000+0 1.046527-2 6.109776-2 6.124281-1 1.094000-29437 2151     \n"
    " 6.983477+2 1.000000+0 2.301353-4 3.433936-2 3.134930-2 0.000000+09437 2151     \n"
    " 7.009578+2 0.000000+0 4.408097-3 4.293751-2 1.548784-2-2.320000+09437 2151     \n"
    " 7.070560+2 1.000000+0 2.586042-3 3.880190-2 8.531837-3 0.000000+09437 2151     \n"
    " 7.076714+2 1.000000+0 1.016568-2 3.818873-2 1.144961-2 0.000000+09437 2151     \n"
    " 7.091374+2 1.000000+0 2.039718-4 4.114392-2 7.748456-2 0.000000+09437 2151     \n"
    " 7.106205+2 1.000000+0 2.051329-4 4.407402-2 3.973504-2 0.000000+09437 2151     \n"
    " 7.132141+2 1.000000+0 1.028906-2 4.684606-2 1.862477-2 0.000000+09437 2151     \n"
    " 7.140326+2 1.000000+0 6.196447-3 5.309077-2 3.857828-2 0.000000+09437 2151     \n"
    " 7.174767+2 1.000000+0 1.168726-2 3.044653-2 6.736366-3 0.000000+09437 2151     \n"
    " 7.188522+2 1.000000+0 7.322007-3 3.560543-2 3.749761-3 0.000000+09437 2151     \n"
    " 7.206331+2 1.000000+0 1.512364-4 4.443304-2 2.689663-2 0.000000+09437 2151     \n"
    " 7.265005+2 0.000000+0 1.356282-3 4.309591-2 3.829210-1 1.000000-79437 2151     \n"
    " 7.278973+2 1.000000+0 2.910703-3 3.153174-2-1.214504-2 0.000000+09437 2151     \n"
    " 7.322661+2 0.000000+0 1.277716-2 4.485396-2-2.303915-1-3.278000+09437 2151     \n"
    " 7.332079+2 1.000000+0 8.462035-3 3.908557-2 5.597476-3 0.000000+09437 2151     \n"
    " 7.347009+2 1.000000+0 7.874644-4 3.581889-2 4.465288-2 0.000000+09437 2151     \n"
    " 7.384343+2 0.000000+0 2.924567-2 6.246509-2 1.059092+0 8.949000-19437 2151     \n"
    " 7.396529+2 1.000000+0 4.327295-3 6.556145-2 1.760294-1 0.000000+09437 2151     \n"
    " 7.446443+2 1.000000+0 1.270830-2 4.369908-2 2.405756-2 0.000000+09437 2151     \n"
    " 7.470574+2 1.000000+0 1.459901-2 4.119231-2 9.399227-3 0.000000+09437 2151     \n"
    " 7.493107+2 1.000000+0 1.189432-3 3.523531-2 4.683992-3 0.000000+09437 2151     \n"
    " 7.507972+2 1.000000+0 2.052612-3 3.295338-2 8.131421-3 0.000000+09437 2151     \n"
    " 7.534326+2 1.000000+0 3.974090-2 3.702278-2 8.585514-3 0.000000+09437 2151     \n"
    " 7.552556+2 1.000000+0 1.723552-3 4.486952-2 4.717091-2 0.000000+09437 2151     \n"
    " 7.577346+2 0.000000+0 7.720173-2 5.366048-2-9.166854-4 2.873000+09437 2151     \n"
    " 7.580851+2 1.000000+0 4.305152-3 4.778312-2 3.059650-2 0.000000+09437 2151     \n"
    " 7.647686+2 1.000000+0 8.850957-3 6.629032-2-1.107648-2 0.000000+09437 2151     \n"
    " 7.671131+2 1.000000+0 3.115478-3 3.707256-2 2.108656-2 0.000000+09437 2151     \n"
    " 7.738660+2 1.000000+0 6.915749-3 4.423599-2 1.233725-2 0.000000+09437 2151     \n"
    " 7.758683+2 0.000000+0 3.331844-3 3.612785-2 3.187576-1 1.760000-29437 2151     \n"
    " 7.776782+2 1.000000+0 1.052153-2 5.083407-2 3.723105-2 0.000000+09437 2151     \n"
    " 7.810239+2 1.000000+0 2.170843-3 3.461434-2 5.462915-2 0.000000+09437 2151     \n"
    " 7.818048+2 1.000000+0 1.365182-2 3.818024-2 3.893343-2 0.000000+09437 2151     \n"
    " 7.826366+2 1.000000+0 2.223440-3 4.374059-2 1.919071-2 0.000000+09437 2151     \n"
    " 7.846973+2 1.000000+0 2.495388-3 3.784303-2 4.507207-3 0.000000+09437 2151     \n"
    " 7.871323+2 1.000000+0 9.598206-3 3.852299-2 1.133270-2 0.000000+09437 2151     \n"
    " 7.894108+2 0.000000+0 6.156494-3 3.877782-2-6.660915-1 1.498000-39437 2151     \n"
    " 7.905304+2 1.000000+0 1.465146-3 3.970871-2 2.677170-2 0.000000+09437 2151     \n"
    " 7.957684+2 1.000000+0 7.531680-3 3.121138-2 3.414007-2 0.000000+09437 2151     \n"
    " 7.971194+2 0.000000+0 1.396401-1 4.084429-2-5.126747+0 3.903000+09437 2151     \n"
    " 8.009453+2 1.000000+0 3.772263-3 3.966945-2 4.591165-2 0.000000+09437 2151     \n"
    " 8.049771+2 0.000000+0 2.888690-2 4.227652-2-1.308175+0 6.736000+09437 2151     \n"
    " 8.053325+2 1.000000+0 5.893815-4 4.413116-2 3.214769-2 0.000000+09437 2151     \n"
    " 8.070266+2 1.000000+0 7.674716-3 4.286197-2 9.801186-3 0.000000+09437 2151     \n"
    " 8.109521+2 1.000000+0 7.321001-4 2.970891-2 2.485183-2 0.000000+09437 2151     \n"
    " 8.166403+2 1.000000+0 1.382405-2 2.854985-2-9.801570-3 0.000000+09437 2151     \n"
    " 8.206701+2 1.000000+0 5.667000-3 2.553094-2 9.381838-3 0.000000+09437 2151     \n"
    " 8.244933+2 0.000000+0 2.274675-2 4.522227-2-1.524007+0 1.668000+09437 2151     \n"
    " 8.252255+2 1.000000+0 5.050805-3 3.915242-2-3.954821-2 0.000000+09437 2151     \n"
    " 8.297894+2 1.000000+0 2.960282-2 4.969616-2 3.207590-2 0.000000+09437 2151     \n"
    " 8.315103+2 0.000000+0 5.724634-3 3.964252-2 1.257763+0 1.993000+09437 2151     \n"
    " 8.330281+2 1.000000+0 1.576122-3 2.699405-2 2.114227-1 0.000000+09437 2151     \n"
    " 8.361280+2 1.000000+0 2.609465-3 4.195789-2 2.099327-3 0.000000+09437 2151     \n"
    " 8.437938+2 1.000000+0 1.132195-3 4.192703-2 2.714001-2 0.000000+09437 2151     \n"
    " 8.438249+2 0.000000+0 1.277081-1 5.256548-2 1.743974-1 8.447000+09437 2151     \n"
    " 8.455773+2 1.000000+0 5.433371-3 4.840749-2-5.156047-2 0.000000+09437 2151     \n"
    " 8.482100+2 1.000000+0 6.802721-3 6.375749-2 7.287004-2 0.000000+09437 2151     \n"
    " 8.534863+2 0.000000+0 1.790161-2 4.242655-2 3.424360-1 1.315000+09437 2151     \n"
    " 8.545286+2 1.000000+0 6.003936-4 4.015874-2 8.956956-2 0.000000+09437 2151     \n"
    " 8.558689+2 1.000000+0 1.909437-3 3.450445-2 5.866353-2 0.000000+09437 2151     \n"
    " 8.581511+2 0.000000+0 2.641483-2 4.026735-2 8.438115+0 4.701000-29437 2151     \n"
    " 8.586356+2 1.000000+0 2.256582-3 3.368855-2 2.675680-2 0.000000+09437 2151     \n"
    " 8.615052+2 1.000000+0 2.031885-3 2.830008-2 4.636393-2 0.000000+09437 2151     \n"
    " 8.659304+2 1.000000+0 1.841591-3 4.236314-2 2.298692-2 0.000000+09437 2151     \n"
    " 8.674831+2 1.000000+0 3.858275-3 4.743730-2 1.453961-1 0.000000+09437 2151     \n"
    " 8.701566+2 1.000000+0 4.077732-3 3.849981-2 5.119784-2 0.000000+09437 2151     \n"
    " 8.744073+2 1.000000+0 1.919706-2 3.929730-2 2.636258-2 0.000000+09437 2151     \n"
    " 8.756697+2 1.000000+0 1.484189-2 3.786754-2 2.179769-2 0.000000+09437 2151     \n"
    " 8.785926+2 1.000000+0 4.334581-3 3.541187-2 8.854534-3 0.000000+09437 2151     \n"
    " 8.825889+2 0.000000+0 2.095104-3 3.959547-2 1.012134+1-6.319000-49437 2151     \n"
    " 8.855485+2 1.000000+0 1.308273-2 5.860365-2 1.158294-2 0.000000+09437 2151     \n"
    " 8.919243+2 1.000000+0 7.033367-2 2.672966-2 4.298756-3 0.000000+09437 2151     \n"
    " 8.958218+2 1.000000+0 7.181442-3 6.382839-2 2.924603-1 0.000000+09437 2151     \n"
    " 8.970871+2 1.000000+0 8.054007-3 4.686219-2 8.069735-2 0.000000+09437 2151     \n"
    " 8.980231+2 0.000000+0 1.248787-3 3.833600-2 2.103582-1-4.047000-29437 2151     \n"
    " 9.037608+2 1.000000+0 1.114268-2 3.528469-2-1.579406-2 0.000000+09437 2151     \n"
    " 9.059955+2 1.000000+0 4.054574-3 3.941024-2-1.007810-2 0.000000+09437 2151     \n"
    " 9.062446+2 0.000000+0 6.928705-2 4.046201-2-6.842530+0 1.203000+09437 2151     \n"
    " 9.085206+2 1.000000+0 4.246855-3 3.334666-2 3.214916-2 0.000000+09437 2151     \n"
    " 9.113147+2 0.000000+0 8.986288-4 4.045764-2 1.779802+0 2.499000-19437 2151     \n"
    " 9.128831+2 0.000000+0 4.379334-2 3.668998-2 1.271080-1 3.177000+09437 2151     \n"
    " 9.161950+2 1.000000+0 4.673094-3 3.238282-2 1.671639-1 0.000000+09437 2151     \n"
    " 9.199744+2 1.000000+0 1.157340-2 3.795798-2-3.287997-1 0.000000+09437 2151     \n"
    " 9.207157+2 1.000000+0 8.896472-3 4.031465-2 1.132565-2 0.000000+09437 2151     \n"
    " 9.231620+2 1.000000+0 2.539198-2 2.564613-2 5.598312-3 0.000000+09437 2151     \n"
    " 9.258159+2 0.000000+0 1.244373-2 4.262094-2 1.915237+0 3.428000-39437 2151     \n"
    " 9.279503+2 1.000000+0 3.569999-3 3.491557-2 9.383353-3 0.000000+09437 2151     \n"
    " 9.296086+2 1.000000+0 1.314013-3 2.947874-2 8.413547-1 0.000000+09437 2151     \n"
    " 9.329601+2 1.000000+0 1.207484-2 2.627164-2-3.388488-3 0.000000+09437 2151     \n"
    " 9.375769+2 1.000000+0 1.888264-3 3.807920-2 1.160439-1 0.000000+09437 2151     \n"
    " 9.376646+2 0.000000+0 5.437480-2 3.995892-2-6.999390+0-7.002000-49437 2151     \n"
    " 9.396823+2 1.000000+0 1.173529-2 5.022624-2 7.957792-2 0.000000+09437 2151     \n"
    " 9.412113+2 1.000000+0 8.135165-3 4.177474-2-1.637486-2 0.000000+09437 2151     \n"
    " 9.441814+2 1.000000+0 8.690877-3 3.836778-2-3.404893-3 0.000000+09437 2151     \n"
    " 9.461826+2 0.000000+0 2.387788-3 3.777355-2 4.105881+0-7.463000-19437 2151     \n"
    " 9.464395+2 1.000000+0 6.895291-3 4.064718-2 3.103748-3 0.000000+09437 2151     \n"
    " 9.517085+2 1.000000+0 6.324384-3 3.552037-2 1.720662-2 0.000000+09437 2151     \n"
    " 9.550253+2 1.000000+0 8.053171-3 4.713097-2 9.471669-2 0.000000+09437 2151     \n"
    " 9.589294+2 1.000000+0 4.134156-2 7.258728-2 4.561590-2 0.000000+09437 2151     \n"
    " 9.624734+2 1.000000+0 2.542085-4 4.767302-2 1.229616-1 0.000000+09437 2151     \n"
    " 9.664055+2 1.000000+0 2.771426-2 4.935312-2 1.285373-1 0.000000+09437 2151     \n"
    " 9.721599+2 0.000000+0 4.635552-2 2.559839-2 6.194969-1 7.061000-49437 2151     \n"
    " 9.744375+2 1.000000+0 6.122124-2 2.302036-2 4.354165-2 0.000000+09437 2151     \n"
    " 9.779830+2 1.000000+0 3.438589-3 3.356349-2 4.123361-2 0.000000+09437 2151     \n"
    " 9.804002+2 0.000000+0 9.077394-2 4.623096-2-4.066149-3 1.991000+09437 2151     \n"
    " 9.854019+2 1.000000+0 7.924796-2 2.959419-2-1.924065-2 0.000000+09437 2151     \n"
    " 9.895984+2 1.000000+0 4.336639-3 6.094465-2 5.594205-1 0.000000+09437 2151     \n"
    " 9.908873+2 1.000000+0 5.942128-3 7.377183-2 8.539720-2 0.000000+09437 2151     \n"
    " 9.942781+2 1.000000+0 2.002483-2 3.663357-2 3.323906-2 0.000000+09437 2151     \n"
    " 9.978914+2 1.000000+0 3.163647-2 2.676963-2-2.148713-1 0.000000+09437 2151     \n"
    " 9.986819+2 1.000000+0 3.741182-2 4.854021-2-7.244252-2 0.000000+09437 2151     \n"
    " 1.002792+3 0.000000+0 5.418385-2 4.416785-2 4.808016-1 1.595000-19437 2151     \n"
    " 1.006220+3 1.000000+0 2.419441-3 3.993902-2 2.795912-2 0.000000+09437 2151     \n"
    " 1.007128+3 1.000000+0 4.571709-3 3.197707-2 1.632066-2 0.000000+09437 2151     \n"
    " 1.009760+3 0.000000+0 5.157725-3 3.971568-2 3.246132+0 1.000000-79437 2151     \n"
    " 1.012270+3 1.000000+0 6.009370-3 5.060078-2 1.136921-2 0.000000+09437 2151     \n"
    " 1.015743+3 1.000000+0 5.806556-4 3.948526-2 3.128093-2 0.000000+09437 2151     \n"
    " 1.017048+3 1.000000+0 1.089299-2 4.660113-2 1.155269-2 0.000000+09437 2151     \n"
    " 1.020197+3 1.000000+0 1.937959-2 5.279650-2 1.325690-1 0.000000+09437 2151     \n"
    " 1.024750+3 1.000000+0 1.485195-3 2.408570-2 2.869970-2 0.000000+09437 2151     \n"
    " 1.025467+3 0.000000+0 2.908738-2 5.805025-2-2.622716-1 3.576000-19437 2151     \n"
    " 1.028311+3 1.000000+0 4.210154-2 5.952788-2-1.052268-1 0.000000+09437 2151     \n"
    " 1.031670+3 1.000000+0 2.321410-2 5.229206-2-6.607440-3 0.000000+09437 2151     \n"
    " 1.033527+3 1.000000+0 3.677093-3 2.865867-2-1.402218-2 0.000000+09437 2151     \n"
    " 1.034869+3 1.000000+0 7.486226-4 2.739936-2-7.840193-3 0.000000+09437 2151     \n"
    " 1.035783+3 1.000000+0 1.303857-3 4.258370-2 2.157418-2 0.000000+09437 2151     \n"
    " 1.036146+3 1.000000+0 7.933356-4 3.638979-2 2.774860-2 0.000000+09437 2151     \n"
    " 1.038911+3 1.000000+0 1.304352-3 4.183032-2 8.174552-2 0.000000+09437 2151     \n"
    " 1.039793+3 1.000000+0 2.038068-3 4.198675-2-1.224094-1 0.000000+09437 2151     \n"
    " 1.041855+3 1.000000+0 8.663990-4 5.249027-2 4.550311-3 0.000000+09437 2151     \n"
    " 1.044168+3 1.000000+0 2.669653-2 4.660496-2-6.298306-2 0.000000+09437 2151     \n"
    " 1.048095+3 0.000000+0 9.750596-2 4.554386-2-1.697648+0 8.337000-39437 2151     \n"
    " 1.048516+3 1.000000+0 2.791158-3 3.663145-2 4.596489-3 0.000000+09437 2151     \n"
    " 1.049584+3 1.000000+0 3.275866-3 4.083969-2 3.043420-3 0.000000+09437 2151     \n"
    " 1.050715+3 1.000000+0 2.717892-3 4.133065-2 1.761958-3 0.000000+09437 2151     \n"
    " 1.056301+3 0.000000+0 4.526807-2 4.241760-2 7.983416-3 1.258000-19437 2151     \n"
    " 1.057277+3 1.000000+0 1.394616-2 4.381212-2 3.227394-2 0.000000+09437 2151     \n"
    " 1.059196+3 1.000000+0 7.059741-4 3.793672-2 3.959612-1 0.000000+09437 2151     \n"
    " 1.061212+3 1.000000+0 1.182077-2 5.067575-2 5.494278-2 0.000000+09437 2151     \n"
    " 1.065221+3 1.000000+0 5.821621-3 5.089913-2-9.976552-3 0.000000+09437 2151     \n"
    " 1.068214+3 1.000000+0 4.205709-4 2.650773-2 2.801431-2 0.000000+09437 2151     \n"
    " 1.070474+3 1.000000+0 5.314930-2 2.488663-2 1.413428-2 0.000000+09437 2151     \n"
    " 1.071723+3 0.000000+0 4.034340-3 4.323569-2 8.442486-2 1.000000-79437 2151     \n"
    " 1.073224+3 0.000000+0 1.026392-2 3.841129-2-2.160738-2 2.606000+09437 2151     \n"
    " 1.074667+3 1.000000+0 7.312618-4 4.233884-2 7.421563-3 0.000000+09437 2151     \n"
    " 1.075984+3 1.000000+0 1.063645-2 2.984706-2-7.343807-3 0.000000+09437 2151     \n"
    " 1.077492+3 0.000000+0 2.274203-3 3.953940-2 9.995465-8-5.533000-29437 2151     \n"
    " 1.079573+3 1.000000+0 1.583330-3 1.962183-2 4.594267-3 0.000000+09437 2151     \n"
    " 1.082616+3 1.000000+0 4.763651-3 3.119178-2 8.718702-4 0.000000+09437 2151     \n"
    " 1.085933+3 1.000000+0 2.277898-2 3.369031-2 3.269608-3 0.000000+09437 2151     \n"
    " 1.089352+3 1.000000+0 1.188682-2 5.466009-2 2.601755-2 0.000000+09437 2151     \n"
    " 1.091579+3 0.000000+0 5.795018-2 4.321113-2-1.149973-2 8.929000+09437 2151     \n"
    " 1.093606+3 1.000000+0 1.384288-3 4.382469-2 4.494405-3 0.000000+09437 2151     \n"
    " 1.097424+3 1.000000+0 5.602531-2 3.856839-2 3.397026-2 0.000000+09437 2151     \n"
    " 1.104250+3 0.000000+0 3.592222-2 3.932480-2-4.537268-3 1.737000+09437 2151     \n"
    " 1.104546+3 1.000000+0 1.038196-2 4.366787-2 6.598540-3 0.000000+09437 2151     \n"
    " 1.104756+3 0.000000+0 6.134121-2 4.046482-2 3.513028+0-4.938000-29437 2151     \n"
    " 1.107374+3 1.000000+0 1.615755-2 5.323221-2-1.362418-2 0.000000+09437 2151     \n"
    " 1.108562+3 1.000000+0 6.392414-3 3.913625-2 6.490083-2 0.000000+09437 2151     \n"
    " 1.111672+3 1.000000+0 1.064808-3 3.849324-2 2.086300-2 0.000000+09437 2151     \n"
    " 1.113971+3 1.000000+0 8.654903-4 4.220656-2-4.371675-2 0.000000+09437 2151     \n"
    " 1.115474+3 0.000000+0 1.858563-2 4.207079-2-5.178448+0 5.796000-29437 2151     \n"
    " 1.117131+3 1.000000+0 1.983723-3 4.015283-2 6.779944-3 0.000000+09437 2151     \n"
    " 1.119293+3 1.000000+0 8.189147-3 3.737087-2-6.165791-3 0.000000+09437 2151     \n"
    " 1.120826+3 1.000000+0 3.447248-3 3.754282-2 6.844607-3 0.000000+09437 2151     \n"
    " 1.124825+3 1.000000+0 1.160017-2 4.117854-2 2.793102-2 0.000000+09437 2151     \n"
    " 1.125101+3 0.000000+0 2.345373-2 4.021372-2 3.769129+0-2.944000-19437 2151     \n"
    " 1.127678+3 1.000000+0 5.913835-4 3.497989-2 2.518119-2 0.000000+09437 2151     \n"
    " 1.134155+3 1.000000+0 6.038607-3 4.407209-2 8.607518-3 0.000000+09437 2151     \n"
    " 1.137336+3 0.000000+0 1.426220-2 3.901862-2-1.716324-1 4.593000+09437 2151     \n"
    " 1.137741+3 1.000000+0 2.303103-3 4.276173-2 9.614522-3 0.000000+09437 2151     \n"
    " 1.139774+3 0.000000+0 2.118515-1 5.014858-2-1.495298+0 9.752000-29437 2151     \n"
    " 1.139813+3 1.000000+0 4.923996-3 4.077958-2 2.907653-3 0.000000+09437 2151     \n"
    " 1.140504+3 1.000000+0 6.828612-3 4.271839-2-3.065406-2 0.000000+09437 2151     \n"
    " 1.142855+3 1.000000+0 1.667720-2 3.833275-2 5.507973-5 0.000000+09437 2151     \n"
    " 1.145319+3 1.000000+0 2.472863-2 2.588293-2-2.257052-3 0.000000+09437 2151     \n"
    " 1.149055+3 1.000000+0 1.987767-2 2.637251-2 3.396376-2 0.000000+09437 2151     \n"
    " 1.151305+3 1.000000+0 3.521314-3 4.249763-2 1.771267-2 0.000000+09437 2151     \n"
    " 1.153184+3 1.000000+0 1.207355-2 3.663405-2-2.729817-3 0.000000+09437 2151     \n"
    " 1.156107+3 1.000000+0 1.181229-3 3.764953-2 5.499649-2 0.000000+09437 2151     \n"
    " 1.158581+3 1.000000+0 2.477591-2 3.631595-2-4.337992-2 0.000000+09437 2151     \n"
    " 1.160890+3 1.000000+0 3.054862-3 2.899469-2-6.893028-3 0.000000+09437 2151     \n"
    " 1.166888+3 0.000000+0 1.107015-3 4.002132-2 1.421171-1 8.714000-19437 2151     \n"
    " 1.170153+3 1.000000+0 8.919303-3 4.132105-2 3.202669-2 0.000000+09437 2151     \n"
    " 1.171026+3 1.000000+0 2.293567-3 4.153504-2 2.289843-2 0.000000+09437 2151     \n"
    " 1.175274+3 0.000000+0 4.318852-2 4.343279-2 2.218665-1 5.874000-29437 2151     \n"
    " 1.178180+3 1.000000+0 1.757446-3 3.623931-2 9.754615-2 0.000000+09437 2151     \n"
    " 1.182049+3 1.000000+0 6.626862-4 3.937227-2-4.628978-2 0.000000+09437 2151     \n"
    " 1.184761+3 1.000000+0 1.013614-3 3.682441-2-6.832193-3 0.000000+09437 2151     \n"
    " 1.184801+3 0.000000+0 1.082958-1 5.229950-2 1.645405+0 8.022000-29437 2151     \n"
    " 1.188711+3 1.000000+0 1.430187-2 4.807617-2 3.770984-2 0.000000+09437 2151     \n"
    " 1.191487+3 1.000000+0 1.079773-3 5.059838-2 7.719867-3 0.000000+09437 2151     \n"
    " 1.194968+3 1.000000+0 3.080157-2 6.611447-2-1.180578-1 0.000000+09437 2151     \n"
    " 1.200378+3 1.000000+0 8.109422-3 3.571965-2 5.570775-3 0.000000+09437 2151     \n"
    " 1.201624+3 1.000000+0 2.795457-3 3.905805-2 1.756895-2 0.000000+09437 2151     \n"
    " 1.205908+3 1.000000+0 2.485851-3 3.967647-2 1.097108-2 0.000000+09437 2151     \n"
    " 1.206104+3 0.000000+0 1.748331-2 3.846141-2-1.602432+0 2.573000-19437 2151     \n"
    " 1.211658+3 0.000000+0 1.209425-1 4.018796-2 7.346228-3 9.532000+09437 2151     \n"
    " 1.212099+3 1.000000+0 5.810250-4 3.881586-2 2.991043-3 0.000000+09437 2151     \n"
    " 1.215993+3 1.000000+0 3.176405-3 4.004388-2-1.045799-2 0.000000+09437 2151     \n"
    " 1.218653+3 1.000000+0 4.131952-2 3.787240-2 1.259410-2 0.000000+09437 2151     \n"
    " 1.220664+3 1.000000+0 1.848112-3 4.021047-2-2.670147-2 0.000000+09437 2151     \n"
    " 1.225531+3 1.000000+0 1.153397-2 4.504651-2-1.258392-2 0.000000+09437 2151     \n"
    " 1.229282+3 1.000000+0 3.740226-3 3.734811-2-7.751557-2 0.000000+09437 2151     \n"
    " 1.231996+3 1.000000+0 4.913699-4 3.323057-2-1.065455-1 0.000000+09437 2151     \n"
    " 1.233772+3 0.000000+0 4.081459-2 4.210576-2 3.617883+0 3.392000-19437 2151     \n"
    " 1.236201+3 1.000000+0 1.215348-2 4.109389-2 6.749051-4 0.000000+09437 2151     \n"
    " 1.237716+3 0.000000+0 2.703245-2 4.232187-2-1.607494+0 4.916000-19437 2151     \n"
    " 1.241586+3 1.000000+0 8.537022-4 3.653658-2-9.793084-2 0.000000+09437 2151     \n"
    " 1.243525+3 1.000000+0 6.095570-4 3.133136-2 1.183038-1 0.000000+09437 2151     \n"
    " 1.245190+3 0.000000+0 7.560257-3 3.664903-2-1.959706-1 5.496000-19437 2151     \n"
    " 1.245232+3 1.000000+0 6.090600-4 3.848958-2 3.070407-3 0.000000+09437 2151     \n"
    " 1.247100+3 0.000000+0 4.172798-2 4.544878-2 3.361983+0 7.908000-19437 2151     \n"
    " 1.247702+3 1.000000+0 8.000590-4 4.155773-2 2.987402-3 0.000000+09437 2151     \n"
    " 1.253453+3 1.000000+0 1.316400-2 4.597079-2 1.494993-2 0.000000+09437 2151     \n"
    " 1.255371+3 1.000000+0 9.202503-3 3.088866-2-2.718695-2 0.000000+09437 2151     \n"
    " 1.258578+3 1.000000+0 3.999747-2 4.573369-2-7.274452-3 0.000000+09437 2151     \n"
    " 1.260566+3 1.000000+0 3.622800-2 5.627086-2-8.467946-3 0.000000+09437 2151     \n"
    " 1.262823+3 1.000000+0 9.839656-4 4.656442-2 1.932349-1 0.000000+09437 2151     \n"
    " 1.267221+3 1.000000+0 1.417106-2 4.271810-2 3.370146-2 0.000000+09437 2151     \n"
    " 1.267353+3 0.000000+0 8.253073-2 4.712075-2-2.645507+0 1.000000-79437 2151     \n"
    " 1.268561+3 1.000000+0 1.301689-2 4.107229-2 4.258010-3 0.000000+09437 2151     \n"
    " 1.269883+3 1.000000+0 3.016181-2 4.566185-2 3.624078-3 0.000000+09437 2151     \n"
    " 1.272911+3 1.000000+0 8.699510-4 4.149174-2 8.584960-2 0.000000+09437 2151     \n"
    " 1.274341+3 1.000000+0 5.542726-3 3.368623-2-3.145915-2 0.000000+09437 2151     \n"
    " 1.277209+3 1.000000+0 1.896647-2 5.787735-2 1.204599-2 0.000000+09437 2151     \n"
    " 1.283074+3 0.000000+0 8.352874-2 2.250198-2 4.294795-2 1.000000-79437 2151     \n"
    " 1.285110+3 1.000000+0 1.681953-3 3.752073-2-8.719130-2 0.000000+09437 2151     \n"
    " 1.287934+3 0.000000+0 3.406277-2 5.114911-2 7.497450-1 1.000000-79437 2151     \n"
    " 1.290356+3 1.000000+0 1.728290-2 4.285103-2-1.029925-2 0.000000+09437 2151     \n"
    " 1.291517+3 1.000000+0 8.189160-4 4.099910-2 9.729493-2 0.000000+09437 2151     \n"
    " 1.294628+3 1.000000+0 6.243523-3 4.276460-2 3.732114-2 0.000000+09437 2151     \n"
    " 1.297717+3 0.000000+0 4.614642-2 3.630446-2 1.228426-2-5.527000+09437 2151     \n"
    " 1.300049+3 1.000000+0 2.188750-2 3.736956-2 5.680086-2 0.000000+09437 2151     \n"
    " 1.300514+3 0.000000+0 9.104001-3 4.064681-2 5.443526-1 1.271000-29437 2151     \n"
    " 1.302662+3 1.000000+0 1.294304-2 3.861520-2 8.110314-3 0.000000+09437 2151     \n"
    " 1.304157+3 1.000000+0 3.217173-3 3.539263-2 2.021898-2 0.000000+09437 2151     \n"
    " 1.306310+3 1.000000+0 1.349048-3 3.909708-2 3.020218-3 0.000000+09437 2151     \n"
    " 1.308443+3 0.000000+0 3.906571-2 3.724239-2 2.139081-2 3.719000+09437 2151     \n"
    " 1.308672+3 1.000000+0 9.500518-4 3.641406-2 2.989611-3 0.000000+09437 2151     \n"
    " 1.311915+3 1.000000+0 4.689760-3 4.210387-2-2.196384-2 0.000000+09437 2151     \n"
    " 1.313732+3 1.000000+0 1.228129-2 4.188955-2 2.179573-2 0.000000+09437 2151     \n"
    " 1.316814+3 0.000000+0 1.001387-2 4.097227-2 8.670481-1 9.966000-39437 2151     \n"
    " 1.317996+3 1.000000+0 1.913387-2 5.208355-2 1.268465-1 0.000000+09437 2151     \n"
    " 1.320015+3 1.000000+0 2.620732-2 6.588413-2-1.625356-1 0.000000+09437 2151     \n"
    " 1.321554+3 0.000000+0 6.809721-3 4.360998-2 2.290857+0 8.977000-29437 2151     \n"
    " 1.327949+3 1.000000+0 1.329424-3 4.214751-2 2.801285-2 0.000000+09437 2151     \n"
    " 1.330488+3 1.000000+0 3.213795-3 3.954934-2 4.541337-2 0.000000+09437 2151     \n"
    " 1.332914+3 0.000000+0 5.536113-2 4.154205-2-2.257164-1 1.030000+19437 2151     \n"
    " 1.333232+3 1.000000+0 2.051351-3 3.655557-2 1.622965-2 0.000000+09437 2151     \n"
    " 1.335366+3 0.000000+0 8.484850-2 4.512400-2-1.691700+0 1.000000-79437 2151     \n"
    " 1.335956+3 1.000000+0 3.713072-2 4.886544-2 4.764897-2 0.000000+09437 2151     \n"
    " 1.337099+3 1.000000+0 1.761371-2 4.479144-2-2.131117-1 0.000000+09437 2151     \n"
    " 1.339669+3 1.000000+0 1.237069-2 3.953003-2 1.646022-2 0.000000+09437 2151     \n"
    " 1.342511+3 1.000000+0 8.454184-3 3.883858-2 6.592526-2 0.000000+09437 2151     \n"
    " 1.344335+3 0.000000+0 8.253742-3 4.269867-2 9.332751-1 1.000000-79437 2151     \n"
    " 1.345704+3 1.000000+0 2.007146-3 3.794868-2-7.841676-2 0.000000+09437 2151     \n"
    " 1.347948+3 0.000000+0 1.072089-2 3.593943-2 2.132422-1 8.768000-29437 2151     \n"
    " 1.348909+3 1.000000+0 3.419954-3 4.049284-2-1.905398-2 0.000000+09437 2151     \n"
    " 1.353692+3 1.000000+0 9.188736-3 4.260532-2 1.104630-1 0.000000+09437 2151     \n"
    " 1.356541+3 0.000000+0 1.010317-1 4.117607-2 3.495376+0-1.473000+09437 2151     \n"
    " 1.357213+3 1.000000+0 1.094646-2 4.096426-2 1.047258-1 0.000000+09437 2151     \n"
    " 1.358497+3 1.000000+0 2.535807-2 4.735169-2-1.771032-1 0.000000+09437 2151     \n"
    " 1.359994+3 1.000000+0 8.863324-3 3.916782-2 1.631789-1 0.000000+09437 2151     \n"
    " 1.361662+3 1.000000+0 6.966761-3 3.425137-2 4.110364-2 0.000000+09437 2151     \n"
    " 1.363090+3 1.000000+0 7.068760-4 3.790795-2-8.414727-2 0.000000+09437 2151     \n"
    " 1.365128+3 1.000000+0 1.577826-3 4.091052-2 3.045546-2 0.000000+09437 2151     \n"
    " 1.369429+3 1.000000+0 1.798773-3 4.037452-2-1.516176-1 0.000000+09437 2151     \n"
    " 1.371539+3 1.000000+0 4.191972-3 3.752292-2 4.974530-2 0.000000+09437 2151     \n"
    " 1.373620+3 1.000000+0 2.603815-3 3.824310-2-4.861209-2 0.000000+09437 2151     \n"
    " 1.376450+3 1.000000+0 8.548504-4 2.814377-2-7.838756-2 0.000000+09437 2151     \n"
    " 1.378450+3 0.000000+0 1.199506-2 3.375397-2 5.067559-1 1.000000-79437 2151     \n"
    " 1.382199+3 1.000000+0 1.054209-2 3.934033-2 6.259976-2 0.000000+09437 2151     \n"
    " 1.383526+3 0.000000+0 4.523258-3 3.966689-2 7.159907-3 5.596000-19437 2151     \n"
    " 1.383686+3 1.000000+0 2.912823-3 4.283493-2 2.616728-2 0.000000+09437 2151     \n"
    " 1.386164+3 1.000000+0 7.652707-3 4.160727-2-1.913886-1 0.000000+09437 2151     \n"
    " 1.390418+3 1.000000+0 9.103578-3 3.715078-2 3.571474-2 0.000000+09437 2151     \n"
    " 1.391250+3 0.000000+0 3.301191-2 4.102170-2-3.815832+0 1.000000-79437 2151     \n"
    " 1.393123+3 1.000000+0 3.695952-2 3.869859-2-5.670782-2 0.000000+09437 2151     \n"
    " 1.396075+3 1.000000+0 4.180059-2 3.984447-2 6.572094-2 0.000000+09437 2151     \n"
    " 1.399899+3 1.000000+0 1.133497-2 3.169490-2 1.127578-2 0.000000+09437 2151     \n"
    " 1.404410+3 1.000000+0 2.259442-2 4.026523-2-1.525842-2 0.000000+09437 2151     \n"
    " 1.406105+3 1.000000+0 2.692652-3 3.886409-2 1.238911-1 0.000000+09437 2151     \n"
    " 1.407640+3 1.000000+0 5.446369-4 3.658161-2 1.349788-2 0.000000+09437 2151     \n"
    " 1.409268+3 1.000000+0 1.295331-2 3.765340-2-5.800456-3 0.000000+09437 2151     \n"
    " 1.411233+3 1.000000+0 2.475757-3 3.832449-2-1.247774-1 0.000000+09437 2151     \n"
    " 1.414704+3 1.000000+0 3.539864-3 4.080207-2-9.379240-2 0.000000+09437 2151     \n"
    " 1.415283+3 0.000000+0 6.576623-2 3.990486-2 6.576892-1 6.214000+09437 2151     \n"
    " 1.416198+3 1.000000+0 9.159507-3 3.983855-2 2.501844-2 0.000000+09437 2151     \n"
    " 1.418055+3 1.000000+0 2.846565-3 3.922851-2 1.174080-2 0.000000+09437 2151     \n"
    " 1.421691+3 1.000000+0 3.701877-2 4.157403-2-8.610752-3 0.000000+09437 2151     \n"
    " 1.423327+3 1.000000+0 3.700302-3 4.029777-2 2.341964-1 0.000000+09437 2151     \n"
    " 1.425140+3 1.000000+0 9.676214-4 4.109447-2-2.832965-3 0.000000+09437 2151     \n"
    " 1.427256+3 1.000000+0 7.453983-3 3.389074-2 5.156030-3 0.000000+09437 2151     \n"
    " 1.428556+3 0.000000+0 3.757522-3 3.978262-2-1.863232+0 1.000000-79437 2151     \n"
    " 1.431145+3 1.000000+0 1.373239-3 4.150563-2 2.595045-2 0.000000+09437 2151     \n"
    " 1.434525+3 1.000000+0 2.403975-2 4.020248-2 2.839590-2 0.000000+09437 2151     \n"
    " 1.436426+3 1.000000+0 1.234656-3 4.909969-2 6.751061-2 0.000000+09437 2151     \n"
    " 1.439029+3 1.000000+0 2.784679-2 4.323870-2-9.826583-3 0.000000+09437 2151     \n"
    " 1.441248+3 1.000000+0 5.523182-4 3.847212-2 9.887014-2 0.000000+09437 2151     \n"
    " 1.443846+3 1.000000+0 3.030246-3 3.993088-2 4.846489-2 0.000000+09437 2151     \n"
    " 1.444806+3 0.000000+0 3.086597-2 4.180735-2 1.354942+0 1.632000-19437 2151     \n"
    " 1.446600+3 1.000000+0 1.130697-2 4.213830-2-6.485488-3 0.000000+09437 2151     \n"
    " 1.451997+3 0.000000+0 1.064282-1 3.934766-2 4.402744+0 2.891000-29437 2151     \n"
    " 1.456536+3 0.000000+0 3.230381-2 3.480047-2-2.529174-4-1.176000+09437 2151     \n"
    " 1.457850+3 1.000000+0 2.117082-2 4.049420-2 1.982139-2 0.000000+09437 2151     \n"
    " 1.464874+3 1.000000+0 5.393109-3 4.070394-2-2.160258-2 0.000000+09437 2151     \n"
    " 1.467679+3 1.000000+0 1.360271-2 4.052552-2 4.679536-2 0.000000+09437 2151     \n"
    " 1.471343+3 0.000000+0 2.660117-2 4.193597-2-3.851847-1 1.000000-79437 2151     \n"
    " 1.472932+3 1.000000+0 1.530723-2 4.143996-2 2.087962-2 0.000000+09437 2151     \n"
    " 1.474796+3 1.000000+0 1.867846-3 4.409829-2 2.703464-1 0.000000+09437 2151     \n"
    " 1.477136+3 1.000000+0 1.624181-2 4.054755-2 1.466782-2 0.000000+09437 2151     \n"
    " 1.481132+3 1.000000+0 3.186048-3 4.155791-2 1.587946-2 0.000000+09437 2151     \n"
    " 1.482196+3 0.000000+0 7.107138-2 4.051360-2-2.985736+0 1.152000-19437 2151     \n"
    " 1.482302+3 1.000000+0 2.798962-3 3.846897-2 1.034943-2 0.000000+09437 2151     \n"
    " 1.484294+3 1.000000+0 2.076087-2 4.152714-2-8.585384-3 0.000000+09437 2151     \n"
    " 1.487245+3 1.000000+0 5.947312-3 3.676911-2-7.586540-3 0.000000+09437 2151     \n"
    " 1.491041+3 1.000000+0 2.020616-3 3.882261-2 4.626620-3 0.000000+09437 2151     \n"
    " 1.492638+3 0.000000+0 1.621608-2 3.897466-2 2.288279+0 1.000000-79437 2151     \n"
    " 1.493682+3 1.000000+0 3.105727-3 4.113223-2 3.966791-1 0.000000+09437 2151     \n"
    " 1.496019+3 1.000000+0 1.299879-3 3.697397-2 3.419558-2 0.000000+09437 2151     \n"
    " 1.499237+3 1.000000+0 7.726135-3 3.983164-2-1.213939-2 0.000000+09437 2151     \n"
    " 1.500043+3 0.000000+0 3.450431-2 4.137739-2 3.687901+0 5.323000-19437 2151     \n"
    " 1.504031+3 0.000000+0 1.182715-2 4.866050-2 4.914865-3 1.649000+09437 2151     \n"
    " 1.506536+3 1.000000+0 4.566242-3 3.917675-2-2.764903-2 0.000000+09437 2151     \n"
    " 1.508870+3 0.000000+0 1.091200-2 4.218317-2 5.276705-1 1.000000-79437 2151     \n"
    " 1.511006+3 1.000000+0 7.125366-3 3.828455-2 2.771777-2 0.000000+09437 2151     \n"
    " 1.512586+3 1.000000+0 8.539281-3 3.900351-2 2.462037-2 0.000000+09437 2151     \n"
    " 1.514396+3 1.000000+0 1.509658-2 4.031653-2 1.108774-2 0.000000+09437 2151     \n"
    " 1.515859+3 1.000000+0 6.610973-3 3.573535-2 1.454309-1 0.000000+09437 2151     \n"
    " 1.519166+3 1.000000+0 2.435101-3 3.601195-2 1.746900-2 0.000000+09437 2151     \n"
    " 1.521717+3 1.000000+0 1.720391-2 4.199410-2 9.072181-3 0.000000+09437 2151     \n"
    " 1.523582+3 0.000000+0 8.347506-3 4.366472-2-2.290847-1 1.582000-19437 2151     \n"
    " 1.523853+3 0.000000+0 8.603553-2 4.044025-2 2.833562+0 6.558000+09437 2151     \n"
    " 1.526531+3 1.000000+0 1.118886-3 4.116832-2 2.354257-1 0.000000+09437 2151     \n"
    " 1.527757+3 1.000000+0 6.964398-3 4.122918-2 2.599506-3 0.000000+09437 2151     \n"
    " 1.532156+3 1.000000+0 4.765866-3 3.622506-2 1.759666-2 0.000000+09437 2151     \n"
    " 1.536685+3 1.000000+0 3.020981-3 4.056626-2 4.070932-3 0.000000+09437 2151     \n"
    " 1.537701+3 1.000000+0 1.564994-2 4.441470-2-5.400205-3 0.000000+09437 2151     \n"
    " 1.538378+3 0.000000+0 6.976159-3 4.064984-2 3.657506-1 1.000000-79437 2151     \n"
    " 1.540620+3 1.000000+0 1.876396-2 4.100180-2-2.872074-2 0.000000+09437 2151     \n"
    " 1.542270+3 1.000000+0 5.650441-4 4.359750-2 1.097142-1 0.000000+09437 2151     \n"
    " 1.546824+3 1.000000+0 2.778844-4 4.307969-2 8.938288-2 0.000000+09437 2151     \n"
    " 1.550835+3 1.000000+0 9.966095-3 3.637146-2-4.115382-2 0.000000+09437 2151     \n"
    " 1.552123+3 1.000000+0 1.172530-2 4.348547-2-1.610061-2 0.000000+09437 2151     \n"
    " 1.557744+3 1.000000+0 3.688738-4 5.567004-2 2.211260-1 0.000000+09437 2151     \n"
    " 1.559574+3 1.000000+0 3.155519-3 3.289076-2-3.232218-2 0.000000+09437 2151     \n"
    " 1.563102+3 0.000000+0 3.916074-3 4.360683-2-2.081452-1 1.000000-79437 2151     \n"
    " 1.564524+3 1.000000+0 5.849923-4 3.977061-2 1.040441-1 0.000000+09437 2151     \n"
    " 1.566241+3 1.000000+0 8.702600-3 4.092420-2 2.046184-2 0.000000+09437 2151     \n"
    " 1.568356+3 1.000000+0 3.631899-3 4.502741-2 5.297826-2 0.000000+09437 2151     \n"
    " 1.569575+3 1.000000+0 2.156036-2 4.211808-2 9.385665-3 0.000000+09437 2151     \n"
    " 1.570128+3 1.000000+0 1.934353-2 3.865555-2-3.255291-2 0.000000+09437 2151     \n"
    " 1.571690+3 1.000000+0 7.109076-3 4.584187-2-2.168033-2 0.000000+09437 2151     \n"
    " 1.574788+3 1.000000+0 1.995084-2 3.890381-2 1.427501-2 0.000000+09437 2151     \n"
    " 1.576279+3 1.000000+0 1.090782-3 4.686034-2 4.216976-2 0.000000+09437 2151     \n"
    " 1.578347+3 1.000000+0 3.839442-2 4.088444-2-1.660528-2 0.000000+09437 2151     \n"
    " 1.580890+3 0.000000+0 5.897379-3 3.896900-2 3.863468-1 1.000000-79437 2151     \n"
    " 1.583039+3 1.000000+0 2.255286-2 3.862757-2 1.821174-2 0.000000+09437 2151     \n"
    " 1.585052+3 0.000000+0 1.240116-2 3.971420-2 6.088996-1 1.000000-79437 2151     \n"
    " 1.587999+3 1.000000+0 7.124930-3 3.314613-2 7.417490-3 0.000000+09437 2151     \n"
    " 1.590350+3 1.000000+0 8.230963-3 2.897816-2-1.320837-2 0.000000+09437 2151     \n"
    " 1.592474+3 0.000000+0 1.276194-3 4.029342-2 7.139953-1 1.000000-79437 2151     \n"
    " 1.594976+3 1.000000+0 2.984195-3 2.760058-2 3.003551-2 0.000000+09437 2151     \n"
    " 1.598915+3 1.000000+0 3.340676-2 3.894831-2-6.901969-3 0.000000+09437 2151     \n"
    " 1.602916+3 0.000000+0 2.948916-2 4.580371-2-1.246428+0 1.000000-79437 2151     \n"
    " 1.605806+3 1.000000+0 2.014086-2 3.951579-2 7.408464-3 0.000000+09437 2151     \n"
    " 1.607970+3 1.000000+0 3.012886-3 3.971982-2-3.011533-2 0.000000+09437 2151     \n"
    " 1.612740+3 1.000000+0 7.474933-2 4.300882-2 2.254005-2 0.000000+09437 2151     \n"
    " 1.613364+3 1.000000+0 2.477356-2 4.228493-2 1.759119-2 0.000000+09437 2151     \n"
    " 1.616698+3 1.000000+0 3.150004-2 4.148632-2-4.907410-3 0.000000+09437 2151     \n"
    " 1.619706+3 0.000000+0 6.178951-2 4.391703-2 1.253031+0 1.000000-79437 2151     \n"
    " 1.620946+3 1.000000+0 1.678829-2 4.091634-2-2.153566-2 0.000000+09437 2151     \n"
    " 1.621931+3 0.000000+0 5.953826-2 3.924169-2 1.300149+0 1.745000+09437 2151     \n"
    " 1.624314+3 1.000000+0 2.059338-2 4.049823-2 3.038709-2 0.000000+09437 2151     \n"
    " 1.626173+3 1.000000+0 6.330705-3 3.894329-2-9.780800-2 0.000000+09437 2151     \n"
    " 1.627701+3 1.000000+0 3.232090-3 3.896631-2 2.568641-1 0.000000+09437 2151     \n"
    " 1.630860+3 0.000000+0 7.036835-2 3.685052-2 3.167514+0 1.000000-79437 2151     \n"
    " 1.634041+3 1.000000+0 2.049784-3 3.431361-2-8.111342-3 0.000000+09437 2151     \n"
    " 1.636987+3 1.000000+0 2.056415-2 3.931203-2-8.573754-3 0.000000+09437 2151     \n"
    " 1.639065+3 1.000000+0 1.188407-4 3.814640-2 1.151873-1 0.000000+09437 2151     \n"
    " 1.639870+3 0.000000+0 5.682098-5 4.167975-2-5.253923-1 1.000000-79437 2151     \n"
    " 1.646270+3 0.000000+0 2.956361-2 4.590595-2-2.923025-2 1.012000+09437 2151     \n"
    " 1.650389+3 1.000000+0 6.056144-3 2.888134-2-1.659557-2 0.000000+09437 2151     \n"
    " 1.652180+3 1.000000+0 1.911769-3 3.952538-2 9.032362-2 0.000000+09437 2151     \n"
    " 1.655306+3 0.000000+0 4.351512-3 4.342360-2-1.029571-3-1.188000+09437 2151     \n"
    " 1.657664+3 1.000000+0 5.728018-3 4.137134-2-1.300393-2 0.000000+09437 2151     \n"
    " 1.659273+3 1.000000+0 3.256254-2 4.276305-2 1.155499-2 0.000000+09437 2151     \n"
    " 1.660763+3 1.000000+0 1.689118-3 4.102038-2 3.309430-2 0.000000+09437 2151     \n"
    " 1.662708+3 1.000000+0 9.724620-4 4.205909-2 8.522015-2 0.000000+09437 2151     \n"
    " 1.664461+3 1.000000+0 2.720784-2 4.210845-2-1.140968-2 0.000000+09437 2151     \n"
    " 1.666295+3 1.000000+0 5.330123-4 4.047333-2 3.694936-1 0.000000+09437 2151     \n"
    " 1.669267+3 0.000000+0 7.730555-2 4.227635-2-1.271293+0 1.000000-79437 2151     \n"
    " 1.669394+3 1.000000+0 1.278070-3 3.916490-2 3.003652-3 0.000000+09437 2151     \n"
    " 1.670667+3 1.000000+0 3.474737-3 3.742640-2-1.846061-1 0.000000+09437 2151     \n"
    " 1.673367+3 1.000000+0 5.409977-2 4.144087-2 1.186652-2 0.000000+09437 2151     \n"
    " 1.674893+3 0.000000+0 1.388118-2 3.947483-2-5.032517-1 1.640000+09437 2151     \n"
    " 1.679674+3 1.000000+0 6.544559-3 3.676866-2 2.005412-2 0.000000+09437 2151     \n"
    " 1.683068+3 1.000000+0 5.349629-3 4.350298-2 3.460519-1 0.000000+09437 2151     \n"
    " 1.684552+3 1.000000+0 5.112396-2 4.048631-2-2.170392-2 0.000000+09437 2151     \n"
    " 1.688258+3 0.000000+0 3.709590-2 3.843721-2 7.044237+0 1.000000-79437 2151     \n"
    " 1.690072+3 1.000000+0 9.290946-4 3.468935-2 3.991099-3 0.000000+09437 2151     \n"
    " 1.693605+3 1.000000+0 2.453897-2 4.100139-2 4.586737-3 0.000000+09437 2151     \n"
    " 1.699208+3 1.000000+0 1.134523-3 5.635610-2-2.231435-1 0.000000+09437 2151     \n"
    " 1.702071+3 1.000000+0 4.152680-2 4.248617-2-1.037341-2 0.000000+09437 2151     \n"
    " 1.704027+3 1.000000+0 2.499437-2 4.156703-2-5.853833-2 0.000000+09437 2151     \n"
    " 1.706112+3 1.000000+0 2.594870-2 4.185049-2 4.547195-2 0.000000+09437 2151     \n"
    " 1.707785+3 1.000000+0 2.108760-3 4.966060-2 1.966808-1 0.000000+09437 2151     \n"
    " 1.710424+3 1.000000+0 4.150933-2 4.091020-2-1.238229-2 0.000000+09437 2151     \n"
    " 1.713803+3 1.000000+0 1.863292-3 3.613587-2 8.541425-2 0.000000+09437 2151     \n"
    " 1.717586+3 1.000000+0 9.607198-4 4.435042-2 6.340364-2 0.000000+09437 2151     \n"
    " 1.719238+3 1.000000+0 4.832266-3 3.794273-2-8.954762-3 0.000000+09437 2151     \n"
    " 1.720387+3 1.000000+0 1.552932-3 4.386218-2 1.486304-2 0.000000+09437 2151     \n"
    " 1.722976+3 1.000000+0 8.688860-3 4.018474-2-1.750815-2 0.000000+09437 2151     \n"
    " 1.724384+3 1.000000+0 2.655551-3 4.416814-2 3.378296-2 0.000000+09437 2151     \n"
    " 1.726803+3 1.000000+0 2.264317-2 4.341697-2-1.544182-2 0.000000+09437 2151     \n"
    " 1.727982+3 1.000000+0 7.505759-3 3.706324-2-1.811601-1 0.000000+09437 2151     \n"
    " 1.729387+3 1.000000+0 5.494941-3 3.800494-2 1.955786-1 0.000000+09437 2151     \n"
    " 1.730909+3 1.000000+0 9.553811-3 4.194650-2 5.204453-2 0.000000+09437 2151     \n"
    " 1.731450+3 1.000000+0 2.157623-3 3.843515-2-5.688960-2 0.000000+09437 2151     \n"
    " 1.732643+3 1.000000+0 1.630371-3 3.822084-2 2.772574-1 0.000000+09437 2151     \n"
    " 1.735278+3 0.000000+0 4.917505-3 3.921070-2 1.964753-1 5.991000+09437 2151     \n"
    " 1.735738+3 1.000000+0 6.345286-3 3.124483-2 3.498890-2 0.000000+09437 2151     \n"
    " 1.736132+3 1.000000+0 1.639569-2 4.585585-2-4.218490-3 0.000000+09437 2151     \n"
    " 1.739937+3 1.000000+0 2.362743-2 4.099582-2-1.050688-2 0.000000+09437 2151     \n"
    " 1.742836+3 1.000000+0 2.759219-4 3.760036-2 1.940767-2 0.000000+09437 2151     \n"
    " 1.745263+3 1.000000+0 3.733285-4 4.092034-2-1.250771-1 0.000000+09437 2151     \n"
    " 1.747604+3 1.000000+0 3.735222-2 4.010153-2 2.278720-2 0.000000+09437 2151     \n"
    " 1.751009+3 0.000000+0 8.153499-3 3.932911-2-3.727060+0 1.000000-79437 2151     \n"
    " 1.754350+3 1.000000+0 1.255327-2 3.955920-2-1.286135-2 0.000000+09437 2151     \n"
    " 1.757384+3 1.000000+0 1.097770-2 3.338594-2 1.928317-2 0.000000+09437 2151     \n"
    " 1.759685+3 1.000000+0 7.007912-4 4.635319-2-4.641222-3 0.000000+09437 2151     \n"
    " 1.761748+3 0.000000+0 7.580063-3 3.981958-2 1.253306+0 1.000000-79437 2151     \n"
    " 1.763383+3 1.000000+0 6.735378-3 3.610739-2-4.326686-3 0.000000+09437 2151     \n"
    " 1.764556+3 1.000000+0 2.818982-2 4.261843-2-9.595198-3 0.000000+09437 2151     \n"
    " 1.767725+3 1.000000+0 8.995019-4 4.151672-2 7.294351-2 0.000000+09437 2151     \n"
    " 1.770328+3 1.000000+0 3.853142-3 3.851691-2-1.086816-1 0.000000+09437 2151     \n"
    " 1.772942+3 1.000000+0 2.038642-2 3.940415-2 6.015464-3 0.000000+09437 2151     \n"
    " 1.775200+3 1.000000+0 1.646264-3 4.575969-2 5.187251-2 0.000000+09437 2151     \n"
    " 1.778522+3 1.000000+0 8.861833-3 4.449717-2-4.033587-2 0.000000+09437 2151     \n"
    " 1.779766+3 1.000000+0 3.833288-3 3.949910-2 4.296136-2 0.000000+09437 2151     \n"
    " 1.781136+3 1.000000+0 1.658331-2 3.837457-2-2.182673-1 0.000000+09437 2151     \n"
    " 1.783443+3 1.000000+0 2.117389-3 3.816408-2-2.313702-1 0.000000+09437 2151     \n"
    " 1.786370+3 0.000000+0 7.680299-3 4.047164-2 5.586091-2 7.609000-19437 2151     \n"
    " 1.788145+3 1.000000+0 2.589544-2 4.204621-2-1.908880-2 0.000000+09437 2151     \n"
    " 1.788651+3 0.000000+0 1.196252-1 3.945081-2-2.000092+1 1.000000-79437 2151     \n"
    " 1.790004+3 1.000000+0 3.733801-3 4.023356-2-6.123860-2 0.000000+09437 2151     \n"
    " 1.794294+3 1.000000+0 4.739712-3 4.173772-2-3.931534-2 0.000000+09437 2151     \n"
    " 1.796373+3 1.000000+0 1.180302-1 4.194016-2 3.567523-2 0.000000+09437 2151     \n"
    " 1.800653+3 1.000000+0 2.733420-2 4.118328-2 2.756310-3 0.000000+09437 2151     \n"
    " 1.803328+3 1.000000+0 4.560241-2 4.126434-2-5.360029-3 0.000000+09437 2151     \n"
    " 1.807995+3 1.000000+0 7.914902-2 4.298813-2-6.707160-2 0.000000+09437 2151     \n"
    " 1.809487+3 1.000000+0 5.622799-2 3.840444-2-9.354808-2 0.000000+09437 2151     \n"
    " 1.813266+3 0.000000+0 9.599051-3 3.821830-2 1.153139+0 1.000000-79437 2151     \n"
    " 1.815816+3 1.000000+0 1.330026-2 3.839700-2 3.430883-2 0.000000+09437 2151     \n"
    " 1.819480+3 1.000000+0 1.828068-3 4.093588-2-1.872442-1 0.000000+09437 2151     \n"
    " 1.825811+3 0.000000+0 1.103033-1 5.126171-2 6.761695-3 9.139000-19437 2151     \n"
    " 1.827046+3 1.000000+0 1.626799-2 4.290127-2-7.932833-2 0.000000+09437 2151     \n"
    " 1.827930+3 0.000000+0 2.639295-2 3.889639-2-2.144782+0 1.000000-79437 2151     \n"
    " 1.832219+3 1.000000+0 3.394646-2 4.146287-2 8.351866-3 0.000000+09437 2151     \n"
    " 1.834919+3 1.000000+0 1.447181-3 3.894226-2 3.197874-3 0.000000+09437 2151     \n"
    " 1.835118+3 1.000000+0 6.738376-2 4.463785-2-3.932375-2 0.000000+09437 2151     \n"
    " 1.837947+3 1.000000+0 7.508285-4 4.028248-2 2.659130-1 0.000000+09437 2151     \n"
    " 1.840528+3 1.000000+0 5.669787-3 3.774066-2-6.529168-2 0.000000+09437 2151     \n"
    " 1.844185+3 1.000000+0 6.634817-3 3.871627-2 2.582649-3 0.000000+09437 2151     \n"
    " 1.845372+3 0.000000+0 2.303310-2 3.923165-2 8.107456-2 3.068000+09437 2151     \n"
    " 1.846396+3 1.000000+0 2.271896-2 4.309301-2-2.038301-2 0.000000+09437 2151     \n"
    " 1.847579+3 0.000000+0 6.843698-3 3.981562-2-2.265982+0 1.000000-79437 2151     \n"
    " 1.849440+3 1.000000+0 4.417666-3 3.845386-2 3.248611-3 0.000000+09437 2151     \n"
    " 1.851833+3 1.000000+0 1.555685-2 3.973718-2 1.278448-2 0.000000+09437 2151     \n"
    " 1.853212+3 1.000000+0 1.608623-2 3.999510-2-9.825459-3 0.000000+09437 2151     \n"
    " 1.856685+3 0.000000+0 1.450978-1 4.391652-2 4.225984+0 1.000000-79437 2151     \n"
    " 1.858037+3 1.000000+0 1.352448-2 3.870988-2-4.021323-2 0.000000+09437 2151     \n"
    " 1.862482+3 1.000000+0 2.092270-3 4.033990-2-1.816928-1 0.000000+09437 2151     \n"
    " 1.864925+3 1.000000+0 1.781617-3 4.091445-2-1.916259-1 0.000000+09437 2151     \n"
    " 1.868829+3 0.000000+0 2.252456-2 4.667733-2 2.261868+0 1.000000-79437 2151     \n"
    " 1.873767+3 1.000000+0 1.660267-2 4.247915-2-1.106368-1 0.000000+09437 2151     \n"
    " 1.874293+3 1.000000+0 2.347360-2 3.857995-2 1.052983-1 0.000000+09437 2151     \n"
    " 1.876245+3 1.000000+0 2.819939-3 4.184800-2-9.769245-2 0.000000+09437 2151     \n"
    " 1.878438+3 1.000000+0 1.724067-2 3.661624-2-1.309373-2 0.000000+09437 2151     \n"
    " 1.881749+3 0.000000+0 2.379626-2 3.726952-2-1.990200+0 1.000000-79437 2151     \n"
    " 1.886432+3 1.000000+0 5.176145-3 3.970175-2 1.128750-2 0.000000+09437 2151     \n"
    " 1.887712+3 1.000000+0 2.907434-3 4.444454-2 2.014829-2 0.000000+09437 2151     \n"
    " 1.890202+3 1.000000+0 2.089314-3 4.361563-2-1.480986-1 0.000000+09437 2151     \n"
    " 1.892832+3 1.000000+0 1.067590-2 3.502816-2 2.003010-2 0.000000+09437 2151     \n"
    " 1.896809+3 0.000000+0 3.049036-3 3.947625-2 4.857613+0 1.000000-79437 2151     \n"
    " 1.896862+3 1.000000+0 1.833640-2 3.680690-2-1.237599-1 0.000000+09437 2151     \n"
    " 1.898252+3 1.000000+0 4.458617-2 4.062546-2 2.339613-2 0.000000+09437 2151     \n"
    " 1.901814+3 1.000000+0 1.783316-3 3.905936-2-1.070232-1 0.000000+09437 2151     \n"
    " 1.904500+3 1.000000+0 1.435310-2 3.602175-2-7.645693-2 0.000000+09437 2151     \n"
    " 1.905953+3 0.000000+0 6.919762-3 4.064975-2 1.935356+0 1.000000-79437 2151     \n"
    " 1.907851+3 1.000000+0 3.862864-3 3.770130-2 2.974175-2 0.000000+09437 2151     \n"
    " 1.909941+3 1.000000+0 5.007714-4 3.878670-2 6.134440-2 0.000000+09437 2151     \n"
    " 1.910468+3 1.000000+0 1.633536-3 4.180700-2 3.156390-2 0.000000+09437 2151     \n"
    " 1.912310+3 1.000000+0 3.728859-2 3.124288-2-1.325304-2 0.000000+09437 2151     \n"
    " 1.913897+3 0.000000+0 3.686217-3 3.909247-2-2.599948+0 1.000000-79437 2151     \n"
    " 1.919222+3 1.000000+0 1.014184-2 3.622871-2-1.332703-2 0.000000+09437 2151     \n"
    " 1.921456+3 1.000000+0 6.965350-3 3.487346-2 1.433860-1 0.000000+09437 2151     \n"
    " 1.924497+3 1.000000+0 8.697789-3 3.395444-2 1.519711-2 0.000000+09437 2151     \n"
    " 1.929482+3 1.000000+0 1.801005-2 3.954849-2 1.517909-2 0.000000+09437 2151     \n"
    " 1.930696+3 0.000000+0 1.015663-3 3.934661-2 2.095204+0 1.000000-79437 2151     \n"
    " 1.932609+3 1.000000+0 1.275695-2 3.857974-2-2.111439-2 0.000000+09437 2151     \n"
    " 1.933880+3 1.000000+0 7.762908-3 3.814652-2-3.602851-2 0.000000+09437 2151     \n"
    " 1.938739+3 0.000000+0 6.562808-2 4.867894-2 1.112319-1 1.000000-79437 2151     \n"
    " 1.941404+3 1.000000+0 7.600405-3 4.618863-2-1.555852-2 0.000000+09437 2151     \n"
    " 1.942477+3 1.000000+0 1.066440-3 4.563059-2-1.206249-2 0.000000+09437 2151     \n"
    " 1.942649+3 0.000000+0 1.131916-2 3.982420-2 6.164993+0 1.000000-79437 2151     \n"
    " 1.946810+3 1.000000+0 5.016763-4 4.132358-2 8.497657-2 0.000000+09437 2151     \n"
    " 1.949450+3 1.000000+0 3.932724-2 4.041945-2 1.802922-2 0.000000+09437 2151     \n"
    " 1.951738+3 1.000000+0 1.645920-2 3.887614-2-4.367301-2 0.000000+09437 2151     \n"
    " 1.954313+3 1.000000+0 5.706548-2 4.305220-2 5.503020-3 0.000000+09437 2151     \n"
    " 1.957133+3 0.000000+0 8.067423-3 3.958476-2 1.577022+0 1.000000-79437 2151     \n"
    " 1.960508+3 1.000000+0 1.782499-2 4.030867-2-5.518953-2 0.000000+09437 2151     \n"
    " 1.962888+3 1.000000+0 2.043731-3 3.560831-2 4.044925-2 0.000000+09437 2151     \n"
    " 1.965142+3 1.000000+0 1.289593-2 3.961310-2-1.177907-2 0.000000+09437 2151     \n"
    " 1.966565+3 1.000000+0 1.979712-2 3.609119-2 1.550365-2 0.000000+09437 2151     \n"
    " 1.968230+3 1.000000+0 9.685734-4 3.953201-2 1.023339-1 0.000000+09437 2151     \n"
    " 1.970887+3 1.000000+0 5.383918-3 3.064580-2-2.380174-2 0.000000+09437 2151     \n"
    " 1.972393+3 1.000000+0 3.476063-3 4.567003-2 8.708368-3 0.000000+09437 2151     \n"
    " 1.972921+3 1.000000+0 3.171582-2 4.409313-2-5.271527-3 0.000000+09437 2151     \n"
    " 1.976735+3 1.000000+0 4.904093-3 4.216678-2-4.303219-3 0.000000+09437 2151     \n"
    " 1.978102+3 1.000000+0 1.786307-3 4.591341-2-1.227959-2 0.000000+09437 2151     \n"
    " 1.982048+3 1.000000+0 7.325897-4 3.238251-2 3.307333-2 0.000000+09437 2151     \n"
    " 1.982197+3 0.000000+0 8.478486-3 3.612020-2 4.923013+0 1.000000-79437 2151     \n"
    " 1.984908+3 0.000000+0 5.062818-2 2.214625-2 9.552852-1 1.000000-79437 2151     \n"
    " 1.987342+3 1.000000+0 8.927402-3 4.428911-2-3.321284-2 0.000000+09437 2151     \n"
    " 1.991094+3 1.000000+0 3.858645-3 4.248951-2 4.259036-2 0.000000+09437 2151     \n"
    " 1.992806+3 1.000000+0 3.442919-3 3.770163-2-3.974193-2 0.000000+09437 2151     \n"
    " 1.995053+3 1.000000+0 1.195400-2 3.922595-2-7.292491-3 0.000000+09437 2151     \n"
    " 1.996966+3 1.000000+0 2.791911-2 4.081704-2 9.168416-3 0.000000+09437 2151     \n"
    " 1.999382+3 1.000000+0 2.593286-3 3.968152-2 1.060782-2 0.000000+09437 2151     \n"
    " 2.000224+3 1.000000+0 1.208631-2 4.292026-2-6.291131-3 0.000000+09437 2151     \n"
    " 2.002937+3 0.000000+0 6.204478-3 3.756603-2 1.009206-1 1.000000-79437 2151     \n"
    " 2.006801+3 1.000000+0 4.471843-3 4.143149-2 1.587210-2 0.000000+09437 2151     \n"
    " 2.007856+3 1.000000+0 3.901714-3 3.851482-2-1.963365-2 0.000000+09437 2151     \n"
    " 2.008942+3 1.000000+0 1.758549-2 4.126256-2-1.072214-2 0.000000+09437 2151     \n"
    " 2.011326+3 1.000000+0 1.078857-3 4.711429-2 4.995740-2 0.000000+09437 2151     \n"
    " 2.011712+3 0.000000+0 5.494292-3 3.929363-2 3.997577+0 1.000000-79437 2151     \n"
    " 2.013166+3 1.000000+0 3.542766-3 3.600698-2-2.412080-2 0.000000+09437 2151     \n"
    " 2.015261+3 1.000000+0 1.110768-3 1.755669-2 6.170730-3 0.000000+09437 2151     \n"
    " 2.017503+3 0.000000+0 3.176388-2 5.014250-2 1.065723+0 1.000000-79437 2151     \n"
    " 2.019094+3 1.000000+0 3.012926-2 4.332363-2 9.899716-2 0.000000+09437 2151     \n"
    " 2.019524+3 1.000000+0 1.554868-3 4.010445-2 3.260430-2 0.000000+09437 2151     \n"
    " 2.021884+3 1.000000+0 1.733808-2 4.775342-2-2.147677-2 0.000000+09437 2151     \n"
    " 2.025955+3 1.000000+0 1.877880-3 4.426396-2-4.274753-3 0.000000+09437 2151     \n"
    " 2.027722+3 1.000000+0 1.238582-3 3.500693-2 1.805151-2 0.000000+09437 2151     \n"
    " 2.028301+3 1.000000+0 5.221551-3 3.096293-2 3.839535-3 0.000000+09437 2151     \n"
    " 2.032382+3 0.000000+0 2.173619-3 1.783553-2 3.126112-1 1.000000-79437 2151     \n"
    " 2.034429+3 0.000000+0 3.784343-3 4.013090-2-2.794602+0 1.000000-79437 2151     \n"
    " 2.035756+3 1.000000+0 3.146359-3 4.568920-2-2.065112-2 0.000000+09437 2151     \n"
    " 2.035961+3 0.000000+0 1.069828-2 3.882435-2-1.160704+1 1.000000-79437 2151     \n"
    " 2.036997+3 1.000000+0 3.004212-3 4.335961-2-3.821441-2 0.000000+09437 2151     \n"
    " 2.039079+3 1.000000+0 1.550312-2 4.860854-2 1.165717-1 0.000000+09437 2151     \n"
    " 2.040493+3 1.000000+0 8.647200-3 4.140921-2-3.071465-1 0.000000+09437 2151     \n"
    " 2.042985+3 1.000000+0 1.757964-2 4.202762-2 7.632964-3 0.000000+09437 2151     \n"
    " 2.046337+3 1.000000+0 2.488558-3 3.669514-2-2.567506-2 0.000000+09437 2151     \n"
    " 2.051117+3 1.000000+0 1.597228-3 2.499388-2 1.234156-1 0.000000+09437 2151     \n"
    " 2.051284+3 0.000000+0 2.022764-2 3.728365-2 7.469845+0 1.000000-79437 2151     \n"
    " 2.052540+3 0.000000+0 1.402395-3 3.338492-2-7.275207-1 1.000000-79437 2151     \n"
    " 2.054343+3 1.000000+0 3.598474-3 3.462234-2-2.423923-2 0.000000+09437 2151     \n"
    " 2.055301+3 1.000000+0 1.433273-2 4.080422-2 3.101911-3 0.000000+09437 2151     \n"
    " 2.056711+3 1.000000+0 1.669088-3 3.950310-2-1.656836-2 0.000000+09437 2151     \n"
    " 2.059387+3 1.000000+0 6.799629-3 3.967792-2-1.748491-2 0.000000+09437 2151     \n"
    " 2.061483+3 1.000000+0 1.121010-3 4.055615-2 5.100863-2 0.000000+09437 2151     \n"
    " 2.064716+3 1.000000+0 1.934257-3 4.387297-2 1.070339-2 0.000000+09437 2151     \n"
    " 2.066409+3 1.000000+0 4.117391-2 4.248523-2-5.848718-3 0.000000+09437 2151     \n"
    " 2.068316+3 1.000000+0 1.037228-2 4.640571-2 1.141887-2 0.000000+09437 2151     \n"
    " 2.069328+3 1.000000+0 4.430025-3 4.268364-2-1.241931-1 0.000000+09437 2151     \n"
    " 2.073066+3 1.000000+0 2.801642-2 3.873929-2-1.843217-2 0.000000+09437 2151     \n"
    " 2.074483+3 0.000000+0 9.596636-3 3.971789-2-5.979099+0 1.000000-79437 2151     \n"
    " 2.074988+3 1.000000+0 3.671240-3 4.493108-2 1.834708-1 0.000000+09437 2151     \n"
    " 2.078406+3 1.000000+0 2.119101-2 4.259195-2 8.334089-3 0.000000+09437 2151     \n"
    " 2.081419+3 1.000000+0 2.331685-3 4.154518-2-1.979007-1 0.000000+09437 2151     \n"
    " 2.084564+3 1.000000+0 1.122306-2 3.983860-2-4.036004-2 0.000000+09437 2151     \n"
    " 2.085968+3 1.000000+0 2.081318-2 4.110912-2 3.563219-3 0.000000+09437 2151     \n"
    " 2.088471+3 1.000000+0 4.422189-2 4.065505-2 1.014620-2 0.000000+09437 2151     \n"
    " 2.090772+3 0.000000+0 6.156402-3 3.864207-2-2.209422-1 1.000000-79437 2151     \n"
    " 2.093008+3 1.000000+0 2.206815-3 3.557900-2 1.694222-1 0.000000+09437 2151     \n"
    " 2.095339+3 1.000000+0 2.244980-2 3.861175-2-1.033176-2 0.000000+09437 2151     \n"
    " 2.098344+3 1.000000+0 2.646243-2 4.422202-2 6.069811-3 0.000000+09437 2151     \n"
    " 2.101982+3 1.000000+0 1.995764-2 3.641869-2-1.259660-2 0.000000+09437 2151     \n"
    " 2.103303+3 0.000000+0 6.062369-3 3.604469-2 1.505853-1 1.000000-79437 2151     \n"
    " 2.104719+3 1.000000+0 2.413126-3 3.132249-2 8.017105-3 0.000000+09437 2151     \n"
    " 2.108409+3 1.000000+0 1.626485-2 4.248755-2 6.300488-3 0.000000+09437 2151     \n"
    " 2.111106+3 0.000000+0 7.270788-3 3.945213-2-6.429099+0 1.000000-79437 2151     \n"
    " 2.111643+3 1.000000+0 2.690906-2 4.611288-2 3.654087-3 0.000000+09437 2151     \n"
    " 2.113089+3 1.000000+0 9.367311-4 3.608550-2 3.747281-2 0.000000+09437 2151     \n"
    " 2.113174+3 1.000000+0 3.434872-2 3.727510-2-1.059728-2 0.000000+09437 2151     \n"
    " 2.118134+3 0.000000+0 1.217425-2 4.147475-2-8.753085-2 1.000000-79437 2151     \n"
    " 2.121453+3 1.000000+0 6.368312-2 4.304666-2-5.103230-3 0.000000+09437 2151     \n"
    " 2.124165+3 1.000000+0 1.260232-2 4.342024-2 2.740066-2 0.000000+09437 2151     \n"
    " 2.126750+3 1.000000+0 3.359436-2 4.022321-2-1.133049-2 0.000000+09437 2151     \n"
    " 2.130611+3 0.000000+0 4.329111-2 3.944103-2 1.310970+1 1.000000-79437 2151     \n"
    " 2.130841+3 1.000000+0 2.935495-2 4.184943-2 4.880522-3 0.000000+09437 2151     \n"
    " 2.133119+3 0.000000+0 3.223332-3 3.849405-2-2.058017-1 1.000000-79437 2151     \n"
    " 2.137724+3 1.000000+0 1.034220-2 3.821267-2-7.023978-3 0.000000+09437 2151     \n"
    " 2.139222+3 1.000000+0 4.012505-3 3.701797-2 3.055131-3 0.000000+09437 2151     \n"
    " 2.142620+3 1.000000+0 6.948947-4 3.626985-2 3.580697-3 0.000000+09437 2151     \n"
    " 2.147097+3 1.000000+0 5.526107-3 4.229239-2-3.909620-2 0.000000+09437 2151     \n"
    " 2.148809+3 0.000000+0 1.569378-1 4.656564-2-7.770651+0 1.000000-79437 2151     \n"
    " 2.149359+3 1.000000+0 1.340164-2 4.248182-2-3.208952-3 0.000000+09437 2151     \n"
    " 2.151859+3 1.000000+0 1.782580-2 4.048659-2-7.516334-2 0.000000+09437 2151     \n"
    " 2.153253+3 1.000000+0 4.694036-2 4.868454-2 5.664604-2 0.000000+09437 2151     \n"
    " 2.155782+3 0.000000+0 9.328562-4 4.594696-2 1.590105+0 1.000000-79437 2151     \n"
    " 2.157340+3 1.000000+0 2.088267-2 3.788074-2-1.698318-2 0.000000+09437 2151     \n"
    " 2.159320+3 1.000000+0 9.231942-3 4.421580-2 8.996125-2 0.000000+09437 2151     \n"
    " 2.160686+3 1.000000+0 1.197584-2 3.976726-2 2.246279-2 0.000000+09437 2151     \n"
    " 2.164815+3 1.000000+0 2.390328-2 4.174622-2-1.488385-2 0.000000+09437 2151     \n"
    " 2.166778+3 1.000000+0 2.486256-3 3.896671-2-8.682833-2 0.000000+09437 2151     \n"
    " 2.169304+3 1.000000+0 3.315154-2 4.299262-2-1.080015-2 0.000000+09437 2151     \n"
    " 2.171681+3 1.000000+0 1.599791-2 4.191154-2 1.098836-2 0.000000+09437 2151     \n"
    " 2.174323+3 0.000000+0 9.226570-3 3.902123-2-3.713866-2 1.000000-79437 2151     \n"
    " 2.176171+3 1.000000+0 1.597003-2 4.631179-2 1.291736-2 0.000000+09437 2151     \n"
    " 2.178780+3 1.000000+0 2.359097-2 4.294947-2 6.808729-3 0.000000+09437 2151     \n"
    " 2.179792+3 1.000000+0 3.962686-3 3.691358-2-2.664350-2 0.000000+09437 2151     \n"
    " 2.184703+3 1.000000+0 5.681006-3 4.626915-2-2.439285-2 0.000000+09437 2151     \n"
    " 2.186547+3 1.000000+0 5.436049-4 3.409959-2 2.622390-2 0.000000+09437 2151     \n"
    " 2.189803+3 1.000000+0 4.788398-3 4.512204-2-1.606631-1 0.000000+09437 2151     \n"
    " 2.192620+3 1.000000+0 8.875556-3 3.691526-2 1.506114-1 0.000000+09437 2151     \n"
    " 2.194115+3 0.000000+0 7.969122-2 3.742445-2-2.069393-1 1.000000-79437 2151     \n"
    " 2.195611+3 1.000000+0 1.742913-2 4.691567-2-2.763797-1 0.000000+09437 2151     \n"
    " 2.198985+3 1.000000+0 1.316166-2 3.905728-2 6.229653-3 0.000000+09437 2151     \n"
    " 2.201632+3 1.000000+0 3.919791-3 3.891155-2 9.726939-2 0.000000+09437 2151     \n"
    " 2.203380+3 1.000000+0 6.401885-2 4.279449-2 1.071957-2 0.000000+09437 2151     \n"
    " 2.206338+3 0.000000+0 3.798439-2 4.630020-2-3.677324-2 1.000000-79437 2151     \n"
    " 2.206611+3 1.000000+0 3.643981-4 3.922158-2-8.954007-2 0.000000+09437 2151     \n"
    " 2.211037+3 1.000000+0 1.438875-3 3.551663-2 3.193184-3 0.000000+09437 2151     \n"
    " 2.213451+3 1.000000+0 3.131928-2 3.899644-2 4.971989-3 0.000000+09437 2151     \n"
    " 2.215313+3 1.000000+0 9.179737-3 4.319904-2-2.318613-2 0.000000+09437 2151     \n"
    " 2.217343+3 0.000000+0 9.432184-3 4.007261-2 3.144011+0 1.000000-79437 2151     \n"
    " 2.218734+3 1.000000+0 2.502833-3 4.805535-2 6.762811-2 0.000000+09437 2151     \n"
    " 2.220572+3 1.000000+0 8.321700-3 3.705414-2-1.449670-2 0.000000+09437 2151     \n"
    " 2.223537+3 1.000000+0 2.414024-2 3.550126-2-8.936435-3 0.000000+09437 2151     \n"
    " 2.225104+3 1.000000+0 8.560178-5 3.974940-2 4.300977-2 0.000000+09437 2151     \n"
    " 2.226225+3 0.000000+0 1.163234-2 4.005833-2-2.135768-1 1.000000-79437 2151     \n"
    " 2.228507+3 1.000000+0 2.373177-2 4.437303-2 1.093624-2 0.000000+09437 2151     \n"
    " 2.232138+3 1.000000+0 2.158151-3 2.310204-2-1.437166-2 0.000000+09437 2151     \n"
    " 2.234180+3 1.000000+0 7.074490-3 3.540847-2 2.087951-2 0.000000+09437 2151     \n"
    " 2.237278+3 1.000000+0 1.344089-3 3.209526-2 1.421916-2 0.000000+09437 2151     \n"
    " 2.238740+3 1.000000+0 1.428040-2 4.338134-2-9.542316-3 0.000000+09437 2151     \n"
    " 2.241753+3 1.000000+0 2.028986-2 4.443146-2 2.612992-3 0.000000+09437 2151     \n"
    " 2.244602+3 0.000000+0 1.355113-1 4.123770-2 2.292985-1 3.975553+09437 2151     \n"
    " 2.245179+3 1.000000+0 1.015033-3 3.919099-2-3.041409-3 0.000000+09437 2151     \n"
    " 2.249207+3 1.000000+0 3.059499-3 5.099279-2 8.921748-3 0.000000+09437 2151     \n"
    " 2.251437+3 1.000000+0 1.384986-2 4.667644-2 1.084873-2 0.000000+09437 2151     \n"
    " 2.253159+3 0.000000+0 4.989452-2 3.544650-2-2.174407-1 1.000000-79437 2151     \n"
    " 2.254114+3 1.000000+0 3.111839-3 3.984540-2-7.847550-2 0.000000+09437 2151     \n"
    " 2.257791+3 1.000000+0 1.886510-3 4.241857-2-7.435329-2 0.000000+09437 2151     \n"
    " 2.259854+3 1.000000+0 2.051468-2 4.065321-2 1.085016-2 0.000000+09437 2151     \n"
    " 2.262742+3 1.000000+0 4.060248-2 4.364635-2-7.644257-3 0.000000+09437 2151     \n"
    " 2.264511+3 0.000000+0 8.976590-3 3.906788-2 3.030848+0 1.000000-79437 2151     \n"
    " 2.265214+3 1.000000+0 4.313227-3 4.116370-2 6.090438-2 0.000000+09437 2151     \n"
    " 2.268127+3 1.000000+0 1.171167-2 4.165004-2 7.016500-2 0.000000+09437 2151     \n"
    " 2.270079+3 1.000000+0 3.063515-2 3.927836-2-3.934826-2 0.000000+09437 2151     \n"
    " 2.274698+3 0.000000+0 1.614696-2 4.225598-2-1.997802-1 1.000000-79437 2151     \n"
    " 2.276852+3 1.000000+0 1.373407-3 2.999896-2 3.678107-3 0.000000+09437 2151     \n"
    " 2.278821+3 1.000000+0 9.067187-3 4.042639-2 6.261985-2 0.000000+09437 2151     \n"
    " 2.282169+3 1.000000+0 6.799302-3 3.833407-2-7.585087-2 0.000000+09437 2151     \n"
    " 2.284866+3 1.000000+0 1.410530-3 3.095617-2-1.045052-1 0.000000+09437 2151     \n"
    " 2.288072+3 0.000000+0 1.878178-2 3.628915-2 1.174812+0 1.000000-79437 2151     \n"
    " 2.289214+3 1.000000+0 2.866950-2 4.380880-2-1.421272-2 0.000000+09437 2151     \n"
    " 2.292590+3 1.000000+0 5.031047-3 4.366698-2 2.336725-1 0.000000+09437 2151     \n"
    " 2.297360+3 1.000000+0 1.337171-1 4.306273-2 6.342066-3 0.000000+09437 2151     \n"
    " 2.302100+3 0.000000+0 2.426252-3 3.525766-2-1.762211-2 1.000000-79437 2151     \n"
    " 2.304044+3 1.000000+0 7.032930-2 4.739982-2-2.235681-1 0.000000+09437 2151     \n"
    " 2.305977+3 1.000000+0 1.447976-2 4.052549-2 1.004501-1 0.000000+09437 2151     \n"
    " 2.311182+3 1.000000+0 3.247019-2 4.466184-2 1.327924-1 0.000000+09437 2151     \n"
    " 2.316908+3 1.000000+0 1.611031-2 4.025475-2-3.201371-2 0.000000+09437 2151     \n"
    " 2.317339+3 1.000000+0 1.375105-3 3.826571-2-2.648093-2 0.000000+09437 2151     \n"
    " 2.318728+3 0.000000+0 1.524758-1 5.875229-2 9.748594-1 1.000000-79437 2151     \n"
    " 2.321506+3 1.000000+0 4.320247-3 4.852679-2-9.443706-2 0.000000+09437 2151     \n"
    " 2.322106+3 1.000000+0 2.641197-3 4.317863-2 6.560024-2 0.000000+09437 2151     \n"
    " 2.323758+3 1.000000+0 1.431401-2 3.965691-2-2.269174-2 0.000000+09437 2151     \n"
    " 2.325562+3 1.000000+0 9.813499-3 4.008493-2-4.223092-2 0.000000+09437 2151     \n"
    " 2.327372+3 1.000000+0 2.593998-3 3.619644-2-4.311255-2 0.000000+09437 2151     \n"
    " 2.330995+3 1.000000+0 3.534057-2 4.280785-2 4.105250-3 0.000000+09437 2151     \n"
    " 2.333653+3 0.000000+0 1.088940-1 4.034352-2 4.296217+0 1.000000-79437 2151     \n"
    " 2.335568+3 1.000000+0 4.462989-2 4.430467-2-9.663180-3 0.000000+09437 2151     \n"
    " 2.338443+3 1.000000+0 2.684090-2 4.321854-2 4.745834-3 0.000000+09437 2151     \n"
    " 2.340096+3 1.000000+0 1.049295-3 4.036187-2-4.120463-2 0.000000+09437 2151     \n"
    " 2.342587+3 1.000000+0 3.508878-3 3.647737-2-4.908357-2 0.000000+09437 2151     \n"
    " 2.343342+3 0.000000+0 7.549059-3 3.926941-2 4.844125+0 1.000000-79437 2151     \n"
    " 2.344881+3 1.000000+0 5.103697-2 4.262505-2 1.944830-3 0.000000+09437 2151     \n"
    " 2.347394+3 1.000000+0 3.899662-2 4.238115-2 8.750218-3 0.000000+09437 2151     \n"
    " 2.350946+3 1.000000+0 5.563939-3 3.936866-2-4.331138-2 0.000000+09437 2151     \n"
    " 2.351911+3 1.000000+0 8.579570-3 3.269701-2-2.263338-2 0.000000+09437 2151     \n"
    " 2.356348+3 1.000000+0 1.481004-2 3.980822-2 8.908409-3 0.000000+09437 2151     \n"
    " 2.358035+3 1.000000+0 4.600875-2 4.214233-2-8.815548-3 0.000000+09437 2151     \n"
    " 2.363288+3 0.000000+0 3.030832-2 5.141165-2 3.646843-1 1.000000-79437 2151     \n"
    " 2.365913+3 1.000000+0 1.868270-2 4.469454-2-2.579650-2 0.000000+09437 2151     \n"
    " 2.369461+3 1.000000+0 3.873434-3 4.063349-2 1.340909-1 0.000000+09437 2151     \n"
    " 2.372896+3 1.000000+0 1.427538-2 4.307093-2 3.350099-2 0.000000+09437 2151     \n"
    " 2.375592+3 1.000000+0 6.891944-2 4.077528-2 3.312209-2 0.000000+09437 2151     \n"
    " 2.378791+3 1.000000+0 1.486878-2 4.209407-2-2.007397-2 0.000000+09437 2151     \n"
    " 2.378930+3 0.000000+0 9.006650-3 3.943538-2 2.258307+0 1.000000-79437 2151     \n"
    " 2.382096+3 1.000000+0 1.511459-2 4.239330-2-2.144772-2 0.000000+09437 2151     \n"
    " 2.384016+3 1.000000+0 2.182967-2 3.940891-2 1.762001-2 0.000000+09437 2151     \n"
    " 2.385965+3 1.000000+0 5.374500-3 4.431403-2-8.985235-2 0.000000+09437 2151     \n"
    " 2.389117+3 1.000000+0 5.006930-3 4.234616-2 7.605238-2 0.000000+09437 2151     \n"
    " 2.391353+3 1.000000+0 8.046657-3 4.056547-2 1.221953-2 0.000000+09437 2151     \n"
    " 2.395056+3 1.000000+0 5.714356-3 4.143499-2-4.434808-2 0.000000+09437 2151     \n"
    " 2.396711+3 1.000000+0 1.743268-2 4.413565-2-2.742783-2 0.000000+09437 2151     \n"
    " 2.399056+3 0.000000+0 1.469614-2 3.947846-2 2.685494+0 1.000000-79437 2151     \n"
    " 2.399139+3 1.000000+0 7.233525-3 4.791467-2-6.310280-2 0.000000+09437 2151     \n"
    " 2.401663+3 1.000000+0 3.293240-3 4.280458-2 1.445564-1 0.000000+09437 2151     \n"
    " 2.405377+3 1.000000+0 2.866199-2 4.634825-2 8.590478-3 0.000000+09437 2151     \n"
    " 2.407505+3 1.000000+0 1.133941-2 4.188191-2-1.104672-1 0.000000+09437 2151     \n"
    " 2.409798+3 1.000000+0 2.929108-2 4.251974-2-3.311877-2 0.000000+09437 2151     \n"
    " 2.411607+3 0.000000+0 2.656691-2 4.565647-2 1.991763-1 1.000000-79437 2151     \n"
    " 2.414457+3 1.000000+0 7.474390-3 4.357733-2-8.555420-2 0.000000+09437 2151     \n"
    " 2.416545+3 1.000000+0 2.688922-2 4.006468-2 1.214732-2 0.000000+09437 2151     \n"
    " 2.418446+3 1.000000+0 5.775164-3 4.024681-2-1.601892-1 0.000000+09437 2151     \n"
    " 2.420994+3 1.000000+0 5.416815-3 3.862109-2-3.374989-2 0.000000+09437 2151     \n"
    " 2.423309+3 1.000000+0 2.837293-2 4.149984-2 3.721636-2 0.000000+09437 2151     \n"
    " 2.425980+3 1.000000+0 5.135813-3 3.926671-2-6.010302-2 0.000000+09437 2151     \n"
    " 2.429061+3 0.000000+0 6.949047-2 3.775351-2 3.095689+0 1.000000-79437 2151     \n"
    " 2.431032+3 1.000000+0 5.790285-3 2.440887-2 3.572645-2 0.000000+09437 2151     \n"
    " 2.433112+3 0.000000+0 5.787605-3 5.270107-2 9.857250-2 1.000000-79437 2151     \n"
    " 2.436500+3 1.000000+0 4.047958-3 4.539627-2-1.579660-2 0.000000+09437 2151     \n"
    " 2.439757+3 1.000000+0 4.824480-2 4.702040-2-7.361523-3 0.000000+09437 2151     \n"
    " 2.442239+3 1.000000+0 4.350887-2 4.083304-2 2.685357-2 0.000000+09437 2151     \n"
    " 2.446971+3 0.000000+0 1.099686-1 3.999520-2-7.701130+0 1.000000-79437 2151     \n"
    " 2.448527+3 1.000000+0 7.076491-3 4.264229-2 3.730102-2 0.000000+09437 2151     \n"
    " 2.450849+3 1.000000+0 6.815714-2 3.985104-2-7.957897-2 0.000000+09437 2151     \n"
    " 2.456585+3 1.000000+0 2.022294-2 4.290266-2 1.145612-1 0.000000+09437 2151     \n"
    " 2.460223+3 0.000000+0 1.078125-2 4.406361-2 2.653408-1 1.000000-79437 2151     \n"
    " 2.462643+3 1.000000+0 2.507513-2 4.084632-2-5.396590-2 0.000000+09437 2151     \n"
    " 2.466188+3 1.000000+0 4.681939-3 3.025920-2-2.556120-2 0.000000+09437 2151     \n"
    " 2.469624+3 1.000000+0 3.680925-2 3.981501-2 8.971018-3 0.000000+09437 2151     \n"
    " 2.471418+3 1.000000+0 4.835053-3 3.641830-2 8.294363-2 0.000000+09437 2151     \n"
    " 2.473086+3 1.000000+0 8.894406-3 3.534725-2-2.462633-2 0.000000+09437 2151     \n"
    " 2.476036+3 0.000000+0 6.996022-3 4.300423-2 1.453861+0 1.000000-79437 2151     \n"
    " 2.477380+3 1.000000+0 2.277438-2 3.543532-2-2.084932-2 0.000000+09437 2151     \n"
    " 2.480042+3 1.000000+0 5.399367-3 3.204663-2-5.341699-2 0.000000+09437 2151     \n"
    " 2.482107+3 1.000000+0 2.353772-2 4.530128-2 1.693341-1 0.000000+09437 2151     \n"
    " 2.483391+3 0.000000+0 1.505215-3 3.899041-2-1.291420+0 1.000000-79437 2151     \n"
    " 2.484322+3 1.000000+0 4.534185-2 4.024252-2-1.063069-2 0.000000+09437 2151     \n"
    " 2.485523+3 1.000000+0 4.862761-2 6.168987-2 2.027231-1 0.000000+09437 2151     \n"
    " 2.492112+3 1.000000+0 1.524686-2 1.587093-2-1.448784-2 0.000000+09437 2151     \n"
    " 2.494814+3 1.000000+0 5.331384-2 3.178817-2 1.106315-2 0.000000+09437 2151     \n"
    " 2.497162+3 1.000000+0 9.579957-4 3.351754-2-4.020848-3 0.000000+09437 2151     \n"
    " 2.499517+3 1.000000+0 4.093699-3 3.262239-2-7.769965-3 0.000000+09437 2151     \n"
    " 2.504370+3 1.000000+0 3.307642-2 3.927163-2 1.818048-1 0.000000+09437 2151     \n"
    " 2.506800+3 0.000000+0 3.453994-1 3.927646-2 1.326153+0 9.075218-29437 2151     \n"
    " 2.509230+3 1.000000+0 9.130709-2 3.926767-2-1.348610-1 0.000000+09437 2151     \n"
    " 2.511660+3 0.000000+0 3.123558-1 3.926596-2 1.002582+0 1.971943-29437 2151     \n"
    " 2.514090+3 1.000000+0 1.322666-1 3.926598-2-2.281954-1 0.000000+09437 2151     \n"
    " 2.516520+3 1.000000+0 1.436018-1 3.926556-2 3.888109-2 0.000000+09437 2151     \n"
    " 2.518950+3 1.000000+0 1.662716-1 3.926594-2 2.224877-1 0.000000+09437 2151     \n"
    " 2.521380+3 1.000000+0 1.600835-1 3.926564-2 8.565578-2 0.000000+09437 2151     \n"
    " 2.523800+3 1.000000+0 1.829061-1 3.926666-2 4.074383-1 0.000000+09437 2151     \n"
    "                                                                  9437 2  0     \n";
}

std::string Si29() {

  // Si29 ENDF/B-VIII.0 LRF=3 resonance evaluation - with duplicate eliminated
  // channels

  return
    " 1.402900+4 2.872800+1          0          0          1          01428 2151     \n"
    " 1.402900+4 1.000000+0          0          0          1          01428 2151     \n"
    " 1.000000-5 1.300000+6          1          3          0          11428 2151     \n"
    " 5.000000-1 4.400000-1          1          0          3          31428 2151     \n"
    " 2.872800+1 4.400000-1          0          0         54          91428 2151     \n"
    "-2.179600+6 0.000000+0 1.722200+6 4.090800+2 0.000000+0 0.000000+01428 2151     \n"
    "-8.602400+5 0.000000+0 3.417000+4 9.999700-1 0.000000+0 0.000000+01428 2151     \n"
    "-4.312800+5 0.000000+0 2.285100+5 1.005900+0 0.000000+0 0.000000+01428 2151     \n"
    " 3.857642+5 1.000000+0 2.413300+4 4.670000+0 0.000000+0 0.000000+01428 2151     \n"
    " 7.167713+5 0.000000+0 2.193000+5 3.000000+0 0.000000+0 0.000000+01428 2151     \n"
    " 8.620037+5 0.000000+0 4.329300+5 3.000000+0 0.000000+0 0.000000+01428 2151     \n"
    " 1.178602+6 1.000000+0 8.295900+3 3.000000+0 0.000000+0 0.000000+01428 2151     \n"
    " 1.769073+6 0.000000+0 3.213600+1 3.000000+0 0.000000+0 0.000000+01428 2151     \n"
    " 2.248487+6 0.000000+0 1.693200+2 3.000000+0 0.000000+0 0.000000+01428 2151     \n"
    " 2.872800+1 4.400000-1          1          0        102         171428 2151     \n"
    " 1.528200+4-1.000000+0 1.000000+1 1.646000+0 0.000000+0 0.000000+01428 2151     \n"
    " 3.881900+4 2.000000+0 7.592600+1 2.400000+0 0.000000+0 0.000000+01428 2151     \n"
    " 1.596830+5-1.000000+0 1.200300+3 1.900000+0 0.000000+0 0.000000+01428 2151     \n"
    " 1.844565+5-1.000000+0 1.367400+2 1.500000+0 0.000000+0 0.000000+01428 2151     \n"
    " 3.367903+5-1.000000+0 2.512800+3 8.000000-1 0.000000+0 0.000000+01428 2151     \n"
    " 5.522418+5 2.000000+0 1.298900+3 5.700000+0 0.000000+0 0.000000+01428 2151     \n"
    " 5.665584+5-1.000000+0 7.082000+4 3.000000+0 0.000000+0 0.000000+01428 2151     \n"
    " 6.196646+5 2.000000+0 7.259600+2 3.000000+0 0.000000+0 0.000000+01428 2151     \n"
    " 6.497260+5 2.000000+0 1.095900+3 3.000000+0 0.000000+0 0.000000+01428 2151     \n"
    " 6.530646+5 1.000000+0 1.938600+4 6.300000+0 0.000000+0 0.000000+01428 2151     \n"
    " 7.150646+5-1.000000+0 9.785700+2 3.000000-1 0.000000+0 0.000000+01428 2151     \n"
    " 8.724835+5-1.000000+0 1.733500+4 3.000000-1 0.000000+0 0.000000+01428 2151     \n"
    " 9.558912+5 2.000000+0 9.828900+2 3.000000-1 0.000000+0 0.000000+01428 2151     \n"
    " 1.113808+6-1.000000+0 7.653300+4 3.000000-1 0.000000+0 0.000000+01428 2151     \n"
    " 1.122279+6 2.000000+0 4.881600+3 3.000000-1 0.000000+0 0.000000+01428 2151     \n"
    " 1.192268+6 1.000000+0 3.750600+2 3.000000-1 0.000000+0 0.000000+01428 2151     \n"
    " 1.207629+6 2.000000+0 1.979500+4 3.000000-1 0.000000+0 0.000000+01428 2151     \n"
    " 2.872800+1 4.400000-1          2          0         18          31428 2151     \n"
    " 8.022589+5 1.000000+0 9.934900+3 3.000000+0 0.000000+0 0.000000+01428 2151     \n"
    " 1.098426+6 1.000000+0 5.778700+1 3.000000+0 0.000000+0 0.000000+01428 2151     \n"
    " 1.388859+6 1.000000+0 4.271400+3 3.000000+0 0.000000+0 0.000000+01428 2151     \n"
    " 0.000000+0 0.000000+0          0          0          0          01428 2  0     \n";
}
