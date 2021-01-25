void checkCompoundSystem( const CompoundSystem< ReichMoore, ShiftFactor >& );

SCENARIO( "CompoundSystem" ) {

  GIVEN( "valid data for a CompoundSystem" ) {

    // particles
    Particle photon( ParticleID( "g" ), 0.0 * daltons, 0.0 * coulombs, 1., +1);
    Particle neutron( ParticleID( "n" ), neutronMass, 0.0 * coulombs, 0.5, +1);
    Particle fe55( ParticleID( "Fe55" ), 5.446635e+1 * neutronMass,
                   26.0 * coulombs, 0.0, +1);
    Particle fe54( ParticleID( "Fe54" ), 5.347624e+1 * neutronMass,
                   26.0 * coulombs, 0.0, +1);

    // particle pairs
    ParticlePair out( photon, fe55 );
    ParticlePair in( neutron, fe54 );

    // channels
    Channel< Photon > capture1( in, out, 0. * electronVolt, { 0, 0.0, 0.5, -1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic1( in, in, 0. * electronVolt, { 1, 0.5, 0.5, -1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );
    Channel< Photon > capture2( in, out, 0. * electronVolt, { 0, 0.0, 0.5, +1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic2( in, in, 0. * electronVolt, { 0, 0.5, 0.5, +1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );
    Channel< Photon > capture3( in, out, 0. * electronVolt, { 0, 0.0, 1.5, -1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic3( in, in, 0. * electronVolt, { 1, 0.5, 1.5, -1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );
    Channel< Photon > capture4( in, out, 0. * electronVolt, { 0, 0.0, 1.5, +1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic4( in, in, 0. * electronVolt, { 2, 0.5, 1.5, +1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );
    Channel< Photon > capture5( in, out, 0. * electronVolt, { 0, 0.0, 2.5, +1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic5( in, in, 0. * electronVolt, { 2, 0.5, 2.5, +1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );

    // single resonance tables
    ResonanceTable table1(
      { elastic1.channelID() },
      { Resonance( 10. * electronVolt,
                   { 1.0 * rootElectronVolt },
                   2.0 * rootElectronVolt ) } );
    ResonanceTable table2(
      { elastic2.channelID() },
      { Resonance( 100. * electronVolt,
                   { 3.0 * rootElectronVolt },
                   4.0 * rootElectronVolt ) } );
    ResonanceTable table3(
      { elastic3.channelID() },
      { Resonance( 1000. * electronVolt,
                   { 5.0 * rootElectronVolt },
                   6.0 * rootElectronVolt ) } );
    ResonanceTable table4(
      { elastic4.channelID() },
      { Resonance( 10000. * electronVolt,
                   { 7.0 * rootElectronVolt },
                   8.0 * rootElectronVolt ) } );
    ResonanceTable table5(
      { elastic5.channelID() },
      { Resonance( 1e+5 * electronVolt,
                   { 0.0 * rootElectronVolt },
                   10.0 * rootElectronVolt ),
        Resonance( 1.5e+5 * electronVolt,
                   { 9.0 * rootElectronVolt },
                   0. * rootElectronVolt ) } );

    SpinGroup< ReichMoore, ShiftFactor >
        group1( { elastic1 }, std::move( table1 ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group2( { elastic2 }, std::move( table2 ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group3( { elastic3 }, std::move( table3 ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group4( { elastic4 }, std::move( table4 ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group5( { elastic5 }, std::move( table5 ) );

    // particle channel data instances
    ParticleChannelData data11( capture1,
                                { 10. * electronVolt },
                                { 2. * rootElectronVolt }, true );
    ParticleChannelData data12( elastic1,
                                { 10. * electronVolt },
                                { 1. * rootElectronVolt } );
    ParticleChannelData data21( capture2,
                                { 100. * electronVolt },
                                { 4. * rootElectronVolt }, true );
    ParticleChannelData data22( elastic2,
                                { 100. * electronVolt },
                                { 3. * rootElectronVolt } );
    ParticleChannelData data31( capture3,
                                { 1000. * electronVolt },
                                { 6. * rootElectronVolt }, true );
    ParticleChannelData data32( elastic3,
                                { 1000. * electronVolt },
                                { 5. * rootElectronVolt } );
    ParticleChannelData data41( capture4,
                                { 10000. * electronVolt },
                                { 8. * rootElectronVolt }, true );
    ParticleChannelData data42( elastic4,
                                { 10000. * electronVolt },
                                { 7. * rootElectronVolt } );
    ParticleChannelData data51( capture5,
                                { 1e+5 * electronVolt },
                                { 10. * rootElectronVolt }, true );
    ParticleChannelData data52( elastic5,
                                { 1.5e+5 * electronVolt },
                                { 9. * rootElectronVolt } );

    THEN( "a SpinGroup can be constructed using spin groups" ) {

      CompoundSystem< ReichMoore, ShiftFactor >
          system( { group1, group2, group3, group4, group5 } );

      checkCompoundSystem( system );
    } // THEN

    THEN( "a SpinGroup can be constructed using channel data" ) {

      CompoundSystem< ReichMoore, ShiftFactor >
          system( { data11, data12, data21, data22, data31, data32,
                    data41, data42, data51, data52 } );

      checkCompoundSystem( system );
    } // THEN
  } // GIVEN

  GIVEN( "data for a CompoundSystem with errors" ) {

    // particles
    Particle photon( ParticleID( "g" ), 0.0 * daltons, 0.0 * coulombs, 1., +1);
    Particle neutron( ParticleID( "n" ), neutronMass, 0.0 * coulombs, 0.5, +1);
    Particle fe55( ParticleID( "Fe55" ), 5.446635e+1 * neutronMass,
                   26.0 * coulombs, 0.0, +1);
    Particle fe54( ParticleID( "Fe54" ), 5.347624e+1 * neutronMass,
                   26.0 * coulombs, 0.0, +1);

    // particle pairs
    ParticlePair out( photon, fe55 );
    ParticlePair in( neutron, fe54 );

    // channels
    Channel< Photon > capture1( in, out, 0. * electronVolt, { 0, 0.0, 0.5, +1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic1( in, in, 0. * electronVolt, { 0, 0.5, 0.5, +1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );
    Channel< Photon > capture2( in, out, 0. * electronVolt, { 0, 0.0, 0.5, -1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic2( in, in, 0. * electronVolt, { 1, 0.5, 0.5, -1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );
    Channel< Photon > capture3( in, out, 0. * electronVolt, { 0, 0.0, 1.5, -1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic3( in, in, 0. * electronVolt, { 1, 0.5, 1.5, -1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );
    Channel< Photon > capture4( in, out, 0. * electronVolt, { 0, 0.0, 1.5, +1 },
                                { 0.0 * rootBarn } );
    Channel< Neutron > elastic4( in, in, 0. * electronVolt, { 2, 0.5, 1.5, +1 },
                                 { 5.437300e-1 * rootBarn,
                                   5.437300e-1 * rootBarn } );

    // single resonance tables
    ResonanceTable table1(
      { elastic1.channelID() },
      { Resonance( 10. * electronVolt,
                   { 1.0 * rootElectronVolt },
                   2.0 * rootElectronVolt ) } );
    ResonanceTable table2(
      { elastic2.channelID() },
      { Resonance( 100. * electronVolt,
                   { 3.0 * rootElectronVolt },
                   4.0 * rootElectronVolt ) } );
    ResonanceTable table3(
      { elastic3.channelID() },
      { Resonance( 1000. * electronVolt,
                   { 5.0 * rootElectronVolt },
                   6.0 * rootElectronVolt ) } );
    ResonanceTable table4(
      { elastic4.channelID() },
      { Resonance( 10000. * electronVolt,
                   { 7.0 * rootElectronVolt },
                   8.0 * rootElectronVolt ) } );

    SpinGroup< ReichMoore, ShiftFactor >
        group1( { elastic1 }, std::move( table1 ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group2( { elastic2 }, std::move( table2 ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group3( { elastic3 }, std::move( table3 ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group4( { elastic4 }, std::move( table4 ) );

    // no need to test with other template parameter: construction of
    // CompoundSystem is independent of the template parameters
    using ReichMooreShiftFactorCompoundSystem =
        CompoundSystem< ReichMoore, ShiftFactor >;
    using ReichMooreShiftFactorSpinGroup =
        SpinGroup< ReichMoore, ShiftFactor >;

    THEN( "an exception is thrown at construction when there are no spin "
          "groups" ) {

      CHECK_THROWS( ReichMooreShiftFactorCompoundSystem
                          ( std::vector< ReichMooreShiftFactorSpinGroup >{} ) );
    } // THEN

    THEN( "an exception is thrown at construction when there is no channel "
          "data" ) {

      CHECK_THROWS( ReichMooreShiftFactorCompoundSystem
                          ( std::vector< ParticleChannelData >{} ) );
    } // THEN

    THEN( "an exception is thrown at construction when the spin groups are not "
          "unique" ) {

      CHECK_THROWS( ReichMooreShiftFactorCompoundSystem
                          ( { group1, group2, group3, group4, group1 } ) );
    } // THEN
  } // GIVEN
} // SCENARIO

void checkCompoundSystem( const CompoundSystem< ReichMoore, ShiftFactor >& system ) {

  // number of spin groups
  CHECK( 5 == system.spinGroups().size() );

  // reaction identifiers
  auto reactions = system.reactionIDs();
  CHECK( 2 == reactions.size() );
  CHECK( "n,Fe54->capture" == reactions[0].symbol() );
  CHECK( "n,Fe54->n,Fe54" == reactions[1].symbol() );

  // group 1 - Jpi = 0.5-
  auto group = system.spinGroups()[0];
  CHECK( 1 == group.incidentChannels().size() );
  CHECK( 1 == group.channels().size() );
  auto channel = std::get< Channel< Neutron > >( group.channels()[0] );
  CHECK( "n,Fe54{1,1/2,1/2-}" == channel.channelID() );

  auto table = group.resonanceTable();
  CHECK( 1 == table.numberChannels() );
  CHECK( 1 == table.channels().size() );
  CHECK( "n,Fe54{1,1/2,1/2-}" == table.channels()[0] );
  CHECK( 1 == table.numberResonances() );
  CHECK( 1 == table.resonances().size() );
  CHECK( 1 == table.energies().size() );
  CHECK( 10. == Approx( table.energies()[0].value ) );
  auto resonance = table.resonances()[0];
  CHECK( 10. == Approx( resonance.energy().value ) );
  CHECK( 2. == Approx( resonance.eliminatedWidth().value ) );
  CHECK( 1 == resonance.widths().size() );
  CHECK( 1. == Approx( resonance.widths()[0].value ) );

  // group 2 - Jpi = 0.5+
  group = system.spinGroups()[1];
  CHECK( 1 == group.incidentChannels().size() );
  CHECK( 1 == group.channels().size() );
  channel = std::get< Channel< Neutron > >( group.channels()[0] );
  CHECK( "n,Fe54{0,1/2,1/2+}" == channel.channelID() );

  table = group.resonanceTable();
  CHECK( 1 == table.numberChannels() );
  CHECK( 1 == table.channels().size() );
  CHECK( "n,Fe54{0,1/2,1/2+}" == table.channels()[0] );
  CHECK( 1 == table.numberResonances() );
  CHECK( 1 == table.resonances().size() );
  CHECK( 1 == table.energies().size() );
  CHECK( 100. == Approx( table.energies()[0].value ) );
  resonance = table.resonances()[0];
  CHECK( 100. == Approx( resonance.energy().value ) );
  CHECK( 4. == Approx( resonance.eliminatedWidth().value ) );
  CHECK( 1 == resonance.widths().size() );
  CHECK( 3. == Approx( resonance.widths()[0].value ) );

  // group 3 - Jpi = 1.5-
  group = system.spinGroups()[2];
  CHECK( 1 == group.incidentChannels().size() );
  CHECK( 1 == group.channels().size() );
  channel = std::get< Channel< Neutron > >( group.channels()[0] );
  CHECK( "n,Fe54{1,1/2,3/2-}" == channel.channelID() );

  table = group.resonanceTable();
  CHECK( 1 == table.numberChannels() );
  CHECK( 1 == table.channels().size() );
  CHECK( "n,Fe54{1,1/2,3/2-}" == table.channels()[0] );
  CHECK( 1 == table.numberResonances() );
  CHECK( 1 == table.resonances().size() );
  CHECK( 1 == table.energies().size() );
  CHECK( 1000. == Approx( table.energies()[0].value ) );
  resonance = table.resonances()[0];
  CHECK( 1000. == Approx( resonance.energy().value ) );
  CHECK( 6. == Approx( resonance.eliminatedWidth().value ) );
  CHECK( 1 == resonance.widths().size() );
  CHECK( 5. == Approx( resonance.widths()[0].value ) );

  // group 4 - Jpi = 1.5+
  group = system.spinGroups()[3];
  CHECK( 1 == group.incidentChannels().size() );
  CHECK( 1 == group.channels().size() );
  channel = std::get< Channel< Neutron > >( group.channels()[0] );
  CHECK( "n,Fe54{2,1/2,3/2+}" == channel.channelID() );

  table = group.resonanceTable();
  CHECK( 1 == table.numberChannels() );
  CHECK( 1 == table.channels().size() );
  CHECK( "n,Fe54{2,1/2,3/2+}" == table.channels()[0] );
  CHECK( 1 == table.numberResonances() );
  CHECK( 1 == table.resonances().size() );
  CHECK( 1 == table.energies().size() );
  CHECK( 10000. == Approx( table.energies()[0].value ) );
  resonance = table.resonances()[0];
  CHECK( 10000. == Approx( resonance.energy().value ) );
  CHECK( 8. == Approx( resonance.eliminatedWidth().value ) );
  CHECK( 1 == resonance.widths().size() );
  CHECK( 7. == Approx( resonance.widths()[0].value ) );

  // group 5 - Jpi = 2.5+
  group = system.spinGroups()[4];
  CHECK( 1 == group.incidentChannels().size() );
  CHECK( 1 == group.channels().size() );
  channel = std::get< Channel< Neutron > >( group.channels()[0] );
  CHECK( "n,Fe54{2,1/2,5/2+}" == channel.channelID() );

  table = group.resonanceTable();
  CHECK( 1 == table.numberChannels() );
  CHECK( 1 == table.channels().size() );
  CHECK( "n,Fe54{2,1/2,5/2+}" == table.channels()[0] );
  CHECK( 2 == table.numberResonances() );
  CHECK( 2 == table.resonances().size() );
  CHECK( 2 == table.energies().size() );
  CHECK( 1e+5 == Approx( table.energies()[0].value ) );
  resonance = table.resonances()[0];
  CHECK( 1e+5 == Approx( resonance.energy().value ) );
  CHECK( 10. == Approx( resonance.eliminatedWidth().value ) );
  CHECK( 1 == resonance.widths().size() );
  CHECK( 0. == Approx( resonance.widths()[0].value ) );
  resonance = table.resonances()[1];
  CHECK( 1.5e+5 == Approx( resonance.energy().value ) );
  CHECK( 0. == Approx( resonance.eliminatedWidth().value ) );
  CHECK( 1 == resonance.widths().size() );
  CHECK( 9. == Approx( resonance.widths()[0].value ) );

  // grid information
  auto grid = system.grid();
  CHECK( 18 == grid.size() );
  CHECK( 10. == Approx( grid[1].value ) );
  CHECK( 100. == Approx( grid[4].value ) );
  CHECK( 1000. == Approx( grid[7].value ) );
  CHECK( 10000. == Approx( grid[10].value ) );
  CHECK( 100000. == Approx( grid[13].value ) );
  CHECK( 150000. == Approx( grid[16].value ) );
}
