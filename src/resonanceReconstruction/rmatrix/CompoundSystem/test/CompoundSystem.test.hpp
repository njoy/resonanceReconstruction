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

    // conversion from Gamma to gamma
    auto eGamma = [&] ( double width, const Energy& energy,
                        const Channel< Neutron >& elastic ) -> ReducedWidth {
      return std::sqrt( width / 2. / elastic.penetrability( energy ) ) *
             rootElectronVolt;
    };
    auto cGamma = [&] ( double width ) -> ReducedWidth {
      return std::sqrt( width / 2. ) * rootElectronVolt;
    };

    // single resonance tables
    ResonanceTable table1(
      { elastic1.channelID() },
      { Resonance( 5.152000e+4 * electronVolt,
                   { eGamma( 1.600200e+1, 5.152000e+4 * electronVolt,
                             elastic1 ) },
                   cGamma( 3.600000e-1 ) ) } );
    ResonanceTable table2(
      { elastic2.channelID() },
      { Resonance( 7.788000e+3 * electronVolt,
                   { eGamma( 1.187354e+3, 7.788000e+3 * electronVolt,
                             elastic2 ) },
                   cGamma( 1.455000e+0 ) ) } );
    ResonanceTable table3(
      { elastic3.channelID() },
      { Resonance( 3.099000e+3 * electronVolt,
                   { eGamma( 1.400000e-3, 3.099000e+3 * electronVolt,
                             elastic3 ) },
                   cGamma( 5.900000e-1 ) ) } );
    ResonanceTable table4(
      { elastic4.channelID() },
      { Resonance( 9.480000e+3 * electronVolt,
                   { eGamma( 1.200000e+0, 9.480000e+3 * electronVolt,
                             elastic4 ) },
                   cGamma( 2.700000e-1 ) ) } );
    ResonanceTable table5(
      { elastic5.channelID() },
      { Resonance( 1.264000e+5 * electronVolt,
                   { eGamma( 2.900000e+0, 1.264000e+5 * electronVolt,
                             elastic5 ) },
                   cGamma( 1.100000e+0 ) ) } );

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
                                { 7.788000e+3 * electronVolt },
                                { cGamma( 1.455000e+0 ) }, true );
    ParticleChannelData data12( elastic1,
                                { 7.788000e+3 * electronVolt },
                                { eGamma( 1.187354e+3, 7.788000e+3 * electronVolt,
                                          elastic1 ) } );
    ParticleChannelData data21( capture2,
                                { 5.152000e+4 * electronVolt },
                                { cGamma( 3.600000e-1 ) }, true );
    ParticleChannelData data22( elastic2,
                                { 5.152000e+4 * electronVolt },
                                { eGamma( 1.600200e+1, 5.152000e+4 * electronVolt,
                                          elastic2 ) } );
    ParticleChannelData data31( capture3,
                                { 3.099000e+3 * electronVolt },
                                { cGamma( 5.900000e-1 ) }, true );
    ParticleChannelData data32( elastic3,
                                { 3.099000e+3 * electronVolt },
                                { eGamma( 1.400000e-3, 3.099000e+3 * electronVolt,
                                          elastic3 ) } );
    ParticleChannelData data41( capture4,
                                { 9.480000e+3 * electronVolt },
                                { cGamma( 2.700000e-1 ) }, true );
    ParticleChannelData data42( elastic4,
                                { 9.480000e+3 * electronVolt },
                                { eGamma( 1.200000e+0, 9.480000e+3 * electronVolt,
                                          elastic4 ) } );
    ParticleChannelData data51( capture5,
                                { 1.264000e+5 * electronVolt },
                                { cGamma( 1.100000e+0 ) }, true );
    ParticleChannelData data52( elastic5,
                                { 1.264000e+5 * electronVolt },
                                { eGamma( 2.900000e+0, 1.264000e+5 * electronVolt,
                                          elastic5 ) } );

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

    // conversion from Gamma to gamma
    auto eGamma = [&] ( double width, const Energy& energy,
                        const Channel< Neutron >& elastic ) -> ReducedWidth {
      return std::sqrt( width / 2. / elastic.penetrability( energy ) ) *
             rootElectronVolt;
    };
    auto cGamma = [&] ( double width ) -> ReducedWidth {
      return std::sqrt( width / 2. ) * rootElectronVolt;
    };

    // multiple resonance table
    ResonanceTable table1(
      { elastic1.channelID() },
      { Resonance( 7.788000e+3 * electronVolt,
                   { eGamma( 1.187354e+3, 7.788000e+3 * electronVolt,
                             elastic1 ) },
                   cGamma( 1.455000e+0 ) ) } );
    ResonanceTable table2(
      { elastic2.channelID() },
      { Resonance( 5.152000e+4 * electronVolt,
                   { eGamma( 1.600200e+1, 5.152000e+4 * electronVolt,
                             elastic2 ) },
                   cGamma( 3.600000e-1 ) ) } );
    ResonanceTable table3(
      { elastic3.channelID() },
      { Resonance( 3.099000e+3 * electronVolt,
                   { eGamma( 1.400000e-3, 3.099000e+3 * electronVolt,
                             elastic3 ) },
                   cGamma( 5.900000e-1 ) ) } );
    ResonanceTable table4(
      { elastic4.channelID() },
      { Resonance( 9.480000e+3 * electronVolt,
                   { eGamma( 1.200000e+0, 9.480000e+3 * electronVolt,
                             elastic4 ) },
                   cGamma( 2.700000e-1 ) ) } );

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

  CHECK( 5 == system.spinGroups().size() );

  // group 1 - Jpi = 0.5-
  auto group = system.spinGroups()[0];
  CHECK( 1 == group.incidentChannels().size() );
  CHECK( 1 == group.channels().size() );
  auto channel = std::get< Channel< Neutron > >( group.channels()[0] );
  CHECK( "n,Fe54{1,1/2,1/2-}" == channel.channelID() );

  // group 2 - Jpi = 0.5+
  group = system.spinGroups()[1];
  CHECK( 1 == group.incidentChannels().size() );
  CHECK( 1 == group.channels().size() );
  channel = std::get< Channel< Neutron > >( group.channels()[0] );
  CHECK( "n,Fe54{0,1/2,1/2+}" == channel.channelID() );

  // group 3 - Jpi = 1.5-
  group = system.spinGroups()[2];
  CHECK( 1 == group.incidentChannels().size() );
  CHECK( 1 == group.channels().size() );
  channel = std::get< Channel< Neutron > >( group.channels()[0] );
  CHECK( "n,Fe54{1,1/2,3/2-}" == channel.channelID() );

  // group 4 - Jpi = 1.5+
  group = system.spinGroups()[3];
  CHECK( 1 == group.incidentChannels().size() );
  CHECK( 1 == group.channels().size() );
  channel = std::get< Channel< Neutron > >( group.channels()[0] );
  CHECK( "n,Fe54{2,1/2,3/2+}" == channel.channelID() );

  // group 5 - Jpi = 2.5+
  group = system.spinGroups()[4];
  CHECK( 1 == group.incidentChannels().size() );
  CHECK( 1 == group.channels().size() );
  channel = std::get< Channel< Neutron > >( group.channels()[0] );
  CHECK( "n,Fe54{2,1/2,5/2+}" == channel.channelID() );
}
