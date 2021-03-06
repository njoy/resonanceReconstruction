SCENARIO( "evaluate" ) {

  GIVEN( "valid data for a CompoundSystem with only one SpinGroup without "
         "missing J values using the Reich Moore formalism" ) {

    // test based on Fe54 ENDF/B-VIII.0 LRF7 resonance evaluation
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required
    // cross section values extracted from NJOY2016.39

    // because the orbital angular momentum l = 0 for these SpinGroup,
    // SpinGroup< ReichMoore, ShiftFactor > and SpinGroup< ReichMoore, Constant >
    // should give the same results as B = S = 0

    // using SpinGroup< ReichMoore, ShiftFactor > is equivalent to NJOY2016's
    // LRF7 reconstruction

    // particles
    Particle photon( ParticleID( "g" ), 0.0 * daltons, 0.0 * coulombs, 1., +1);
    Particle neutron( ParticleID( "n" ), neutronMass, 0.0 * coulombs, 0.5, +1);
    Particle fe55( ParticleID( "Fe55" ), 5.446635e+1 * neutronMass,
                   26.0 * coulombs, 0.0, +1);
    Particle fe54( ParticleID( "Fe54" ), 5.347624e+1 * neutronMass,
                   26.0 * coulombs, 0.0, +1);

    // particle pairs
    ParticlePair in( neutron, fe54 );
    ParticlePair out( photon, fe55 );

    // channels
    Channel< Photon > capture( in, out, 0. * electronVolt, { 0, 0.0, 0.5, +1 },
                               { 0.0 * rootBarn },
                               0.0 );
    Channel< Neutron > elastic( in, in, 0. * electronVolt, { 0, 0.5, 0.5, +1 },
                                { 5.437300e-1 * rootBarn,
                                  5.437300e-1 * rootBarn },
                                0.0 );

    // conversion from Gamma to gamma
    auto eGamma = [&] ( double width, const Energy& energy ) -> ReducedWidth {
      return std::sqrt( width / 2. / elastic.penetrability( energy ) ) *
             rootElectronVolt;
    };
    auto cGamma = [&] ( double width ) -> ReducedWidth {
      return std::sqrt( width / 2. ) * rootElectronVolt;
    };

    // single resonance table
    ResonanceTable single(
      { elastic.channelID() },
      { Resonance( 7.788000e+3 * electronVolt,
                   { eGamma( 1.187354e+3, 7.788000e+3 * electronVolt ) },
                   cGamma( 1.455000e+0 ) ) } );
    ResonanceTable single2 = single;

    // multiple resonance table
    ResonanceTable multiple(
      { elastic.channelID() },
      { Resonance( 7.788000e+3 * electronVolt,
                   { eGamma( 1.187354e+3, 7.788000e+3 * electronVolt ) },
                   cGamma( 1.455000e+0 ) ),
        Resonance( 5.287200e+4 * electronVolt,
                   { eGamma( 2.000345e+3, 5.287200e+4 * electronVolt ) },
                   cGamma( 2.000000e+0 ) ),
        Resonance( 7.190500e+4 * electronVolt,
                   { eGamma( 1.781791e+3, 7.190500e+4 * electronVolt ) },
                   cGamma( 2.000000e+0 ) ) } );
    ResonanceTable multiple2 = multiple;

    SpinGroup< ReichMoore, ShiftFactor >
        group1( { elastic }, std::move( single ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group2( { elastic }, std::move( multiple ) );
    SpinGroup< ReichMoore, Constant >
        group3( { elastic }, std::move( single2 ) );
    SpinGroup< ReichMoore, Constant >
        group4( { elastic }, std::move( multiple2 ) );

    CompoundSystem< ReichMoore, ShiftFactor > system1( { group1 } );
    CompoundSystem< ReichMoore, ShiftFactor > system2( { group2 } );
    CompoundSystem< ReichMoore, Constant > system3( { group3 } );
    CompoundSystem< ReichMoore, Constant > system4( { group4 } );

    ReactionID elas( "n,Fe54->n,Fe54" );
    ReactionID capt( "n,Fe54->capture" );

    THEN( "cross sections can be calculated for a single resonance using the "
          "ShiftFactor boundary condition" ) {

      std::map< ReactionID, CrossSection > xs;
      system1.evaluate( 1e-5 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.575879e-1 == Approx( xs[ elas ].value ) );
      CHECK( 6.895037e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system1.evaluate( 1e-4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.575879e-1 == Approx( xs[ elas ].value ) );
      CHECK( 2.180403e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system1.evaluate( 1e-3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.575877e-1 == Approx( xs[ elas ].value ) );
      CHECK( 6.895040e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system1.evaluate( 1e-2 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.575860e-1 == Approx( xs[ elas ].value ) );
      CHECK( 2.180408e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system1.evaluate( 1e-1 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.575694e-1 == Approx( xs[ elas ].value ) );
      CHECK( 6.895214e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system1.evaluate( 1e+0 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.574031e-1 == Approx( xs[ elas ].value ) );
      CHECK( 2.180961e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system1.evaluate( 1e+1 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.557415e-1 == Approx( xs[ elas ].value ) );
      CHECK( 6.912724e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system1.evaluate( 1e+2 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.392160e-1 == Approx( xs[ elas ].value ) );
      CHECK( 2.237319e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system1.evaluate( 1e+3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.932778e-2 == Approx( xs[ elas ].value ) );
      CHECK( 9.067251e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      system1.evaluate( 1e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.342056e+1 == Approx( xs[ elas ].value ) );
      CHECK( 2.473549e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system1.evaluate( 1e+5 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.979435e+0 == Approx( xs[ elas ].value ) );
      CHECK( 4.915667e-6 == Approx( xs[ capt ].value ) );
      xs.clear();

      system1.evaluate( 1e+6 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.308817e+0 == Approx( xs[ elas ].value ) );
      CHECK( 1.343259e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      system1.evaluate( 7.788000e+3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.424287e+2 == Approx( xs[ elas ].value ) );
      CHECK( 4.241421e-1 == Approx( xs[ capt ].value ) );
      xs.clear();
    } // THEN

    THEN( "cross sections can be calculated for multiple resonances using the "
          "ShiftFactor boundary condition" ) {

      std::map< ReactionID, CrossSection > xs;
      system2.evaluate( 1e-5 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.781791e-2 == Approx( xs[ elas ].value ) );
      CHECK( 7.082909e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 1e-4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.781790e-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.239813e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 1e-3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.781781e-2 == Approx( xs[ elas ].value ) );
      CHECK( 7.082911e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 1e-2 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.781682e-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.239818e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 1e-1 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.780692e-2 == Approx( xs[ elas ].value ) );
      CHECK( 7.083086e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 1e+0 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.770804e-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.240372e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 1e+1 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.672107e-2 == Approx( xs[ elas ].value ) );
      CHECK( 7.100642e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 1e+2 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 7.704194e-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.296876e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 1e+3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 7.098264e-3 == Approx( xs[ elas ].value ) );
      CHECK( 9.259109e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 1e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.062823e+1 == Approx( xs[ elas ].value ) );
      CHECK( 2.502603e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 1e+5 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 5.330137e+0 == Approx( xs[ elas ].value ) );
      CHECK( 5.716074e-5 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 1e+6 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.324263e+0 == Approx( xs[ elas ].value ) );
      CHECK( 3.695041e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 7.788000e+3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.424287e+2 == Approx( xs[ elas ].value ) );
      CHECK( 4.241421e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 5.287200e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.738763e+1 == Approx( xs[ elas ].value ) );
      CHECK( 5.099758e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system2.evaluate( 7.190500e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.390972e+1 == Approx( xs[ elas ].value ) );
      CHECK( 4.208797e-2 == Approx( xs[ capt ].value ) );
      xs.clear();
    } // THEN

    THEN( "cross sections can be calculated for a single resonance using the "
          "Constant boundary condition" ) {

      std::map< ReactionID, CrossSection > xs;
      system3.evaluate( 1e-5 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.575879e-1 == Approx( xs[ elas ].value ) );
      CHECK( 6.895037e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system3.evaluate( 1e-4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.575879e-1 == Approx( xs[ elas ].value ) );
      CHECK( 2.180403e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system3.evaluate( 1e-3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.575877e-1 == Approx( xs[ elas ].value ) );
      CHECK( 6.895040e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system3.evaluate( 1e-2 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.575860e-1 == Approx( xs[ elas ].value ) );
      CHECK( 2.180408e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system3.evaluate( 1e-1 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.575694e-1 == Approx( xs[ elas ].value ) );
      CHECK( 6.895214e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system3.evaluate( 1e+0 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.574031e-1 == Approx( xs[ elas ].value ) );
      CHECK( 2.180961e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system3.evaluate( 1e+1 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.557415e-1 == Approx( xs[ elas ].value ) );
      CHECK( 6.912724e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system3.evaluate( 1e+2 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.392160e-1 == Approx( xs[ elas ].value ) );
      CHECK( 2.237319e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system3.evaluate( 1e+3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.932778e-2 == Approx( xs[ elas ].value ) );
      CHECK( 9.067251e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      system3.evaluate( 1e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.342056e+1 == Approx( xs[ elas ].value ) );
      CHECK( 2.473549e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system3.evaluate( 1e+5 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.979435e+0 == Approx( xs[ elas ].value ) );
      CHECK( 4.915667e-6 == Approx( xs[ capt ].value ) );
      xs.clear();

      system3.evaluate( 1e+6 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.308817e+0 == Approx( xs[ elas ].value ) );
      CHECK( 1.343259e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      system3.evaluate( 7.788000e+3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.424287e+2 == Approx( xs[ elas ].value ) );
      CHECK( 4.241421e-1 == Approx( xs[ capt ].value ) );
      xs.clear();
    } // THEN

    THEN( "cross sections can be calculated for multiple resonances using the "
          "Constant boundary condition" ) {

      std::map< ReactionID, CrossSection > xs;
      system4.evaluate( 1e-5 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.781791e-2 == Approx( xs[ elas ].value ) );
      CHECK( 7.082909e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 1e-4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.781790e-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.239813e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 1e-3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.781781e-2 == Approx( xs[ elas ].value ) );
      CHECK( 7.082911e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 1e-2 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.781682e-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.239818e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 1e-1 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.780692e-2 == Approx( xs[ elas ].value ) );
      CHECK( 7.083086e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 1e+0 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.770804e-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.240372e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 1e+1 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.672107e-2 == Approx( xs[ elas ].value ) );
      CHECK( 7.100642e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 1e+2 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 7.704194e-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.296876e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 1e+3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 7.098264e-3 == Approx( xs[ elas ].value ) );
      CHECK( 9.259109e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 1e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.062823e+1 == Approx( xs[ elas ].value ) );
      CHECK( 2.502603e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 1e+5 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 5.330137e+0 == Approx( xs[ elas ].value ) );
      CHECK( 5.716074e-5 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 1e+6 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.324263e+0 == Approx( xs[ elas ].value ) );
      CHECK( 3.695041e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 7.788000e+3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.424287e+2 == Approx( xs[ elas ].value ) );
      CHECK( 4.241421e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 5.287200e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.738763e+1 == Approx( xs[ elas ].value ) );
      CHECK( 5.099758e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system4.evaluate( 7.190500e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.390972e+1 == Approx( xs[ elas ].value ) );
      CHECK( 4.208797e-2 == Approx( xs[ capt ].value ) );
      xs.clear();
    } // THEN
  } // GIVEN

  GIVEN( "valid data for a CompoundSystem with five SpinGroup without "
         "missing J values" ) {

    // test based on Fe54 ENDF/B-VIII.0 LRF7 resonance evaluation
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required
    // cross section values extracted from NJOY2016.39

    // because the oribital angular momentum values for these SpinGroup are
    // 0, 1 and 2, CompoundSystem< ShiftFactor > and CompoundSystem< Constant > will
    // give different results

    // using SpinGroup< ShiftFactor > is equivalent to NJOY2016's LRF7 reconstruction

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

    // multiple resonance table
    ResonanceTable table1(
      { elastic1.channelID() },
      { Resonance( 7.788000e+3 * electronVolt,
                   { eGamma( 1.187354e+3, 7.788000e+3 * electronVolt,
                             elastic1 ) },
                   cGamma( 1.455000e+0 ) ),
        Resonance( 5.287200e+4 * electronVolt,
                   { eGamma( 2.000345e+3, 5.287200e+4 * electronVolt,
                             elastic1 ) },
                   cGamma( 2.000000e+0 ) ),
        Resonance( 7.190500e+4 * electronVolt,
                   { eGamma( 1.781791e+3, 7.190500e+4 * electronVolt,
                             elastic1 ) },
                   cGamma( 2.000000e+0 ) ) } );
    ResonanceTable table2(
      { elastic2.channelID() },
      { Resonance( 5.152000e+4 * electronVolt,
                   { eGamma( 1.600200e+1, 5.152000e+4 * electronVolt,
                             elastic2 ) },
                   cGamma( 3.600000e-1 ) ),
        Resonance( 5.359000e+4 * electronVolt,
                   { eGamma( 1.700000e+1, 5.359000e+4 * electronVolt,
                             elastic2 ) },
                   cGamma( 1.500000e+0 ) ),
        Resonance( 5.545900e+4 * electronVolt,
                   { eGamma( 3.200000e+1, 5.545900e+4 * electronVolt,
                             elastic2 ) },
                   cGamma( 5.600000e-1 ) ) } );
    ResonanceTable table3(
      { elastic3.channelID() },
      { Resonance( 3.099000e+3 * electronVolt,
                   { eGamma( 1.400000e-3, 3.099000e+3 * electronVolt,
                             elastic3 ) },
                   cGamma( 5.900000e-1 ) ),
        Resonance( 1.358100e+4 * electronVolt,
                   { eGamma( 1.750000e-2, 1.358100e+4 * electronVolt,
                             elastic3 ) },
                   cGamma( 5.900000e-1 ) ),
        Resonance( 1.927800e+4 * electronVolt,
                   { eGamma( 2.750000e-2, 1.927800e+4 * electronVolt,
                             elastic3 ) },
                   cGamma( 5.900000e-1 ) ) } );
    ResonanceTable table4(
      { elastic4.channelID() },
      { Resonance( 9.480000e+3 * electronVolt,
                   { eGamma( 1.200000e+0, 9.480000e+3 * electronVolt,
                             elastic4 ) },
                   cGamma( 2.700000e-1 ) ),
        Resonance( 1.118000e+4 * electronVolt,
                   { eGamma( 3.850100e+0, 1.118000e+4 * electronVolt,
                             elastic4 ) },
                   cGamma( 3.500000e-1 ) ),
        Resonance( 1.445000e+4 * electronVolt,
                   { eGamma( 7.000200e-1, 1.445000e+4 * electronVolt,
                             elastic4 ) },
                   cGamma( 3.500000e-1 ) ) } );
    ResonanceTable table5(
      { elastic5.channelID() },
      { Resonance( 1.264000e+5 * electronVolt,
                   { eGamma( 2.900000e+0, 1.264000e+5 * electronVolt,
                             elastic5 ) },
                   cGamma( 1.100000e+0 ) ),
        Resonance( 1.504700e+5 * electronVolt,
                   { eGamma( 2.600000e+0, 1.504700e+5 * electronVolt,
                             elastic5 ) },
                   cGamma( 9.600000e-1 ) ),
        Resonance( 1.779400e+5 * electronVolt,
                   { eGamma( 1.400000e+0, 1.779400e+5 * electronVolt,
                             elastic5 ) },
                   cGamma( 9.600000e-1 ) ) } );

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

    CompoundSystem< ReichMoore, ShiftFactor >
        system( { group1, group2, group3, group4, group5 } );

    ReactionID elas( "n,Fe54->n,Fe54" );
    ReactionID capt( "n,Fe54->capture" );

    THEN( "cross sections can be calculated" ) {

      std::map< ReactionID, CrossSection > xs;
      system.evaluate( 1e-5 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.781787e-2 == Approx( xs[ elas ].value ) );
      CHECK( 7.082910e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e-4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.781786e-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.239813e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e-3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.781776e-2 == Approx( xs[ elas ].value ) );
      CHECK( 7.082912e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e-2 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.781677e-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.239818e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e-1 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.780687e-2 == Approx( xs[ elas ].value ) );
      CHECK( 7.083087e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e+0 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.770799e-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.240372e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e+1 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.672103e-2 == Approx( xs[ elas ].value ) );
      CHECK( 7.100643e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e+2 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 7.704191e-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.296878e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e+3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 7.100428e-3 == Approx( xs[ elas ].value ) );
      CHECK( 9.259258e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.062845e+1 == Approx( xs[ elas ].value ) );
      CHECK( 2.531090e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e+5 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 5.355058e+0 == Approx( xs[ elas ].value ) );
      CHECK( 5.836614e-5 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1e+6 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.240021e+0 == Approx( xs[ elas ].value ) );
      CHECK( 1.613911e-7 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 7.788000e+3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.424290e+2 == Approx( xs[ elas ].value ) );
      CHECK( 4.241628e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 5.287200e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.739414e+1 == Approx( xs[ elas ].value ) );
      CHECK( 5.169252e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 7.190500e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.392479e+1 == Approx( xs[ elas ].value ) );
      CHECK( 4.209120e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 5.152000e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 5.688153e+1 == Approx( xs[ elas ].value ) );
      CHECK( 1.147127e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 5.359000e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.699520e+1 == Approx( xs[ elas ].value ) );
      CHECK( 3.790368e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 5.545900e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 6.447108e+1 == Approx( xs[ elas ].value ) );
      CHECK( 8.302520e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 3.099000e+3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.288509e-1 == Approx( xs[ elas ].value ) );
      CHECK( 4.129390e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1.358100e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 1.269370e+1 == Approx( xs[ elas ].value ) );
      CHECK( 1.113654e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1.927800e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 7.072284e+0 == Approx( xs[ elas ].value ) );
      CHECK( 1.192964e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 9.480000e+3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.384506e+2 == Approx( xs[ elas ].value ) );
      CHECK( 8.551894e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1.118000e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.293224e+2 == Approx( xs[ elas ].value ) );
      CHECK( 3.693614e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1.445000e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 1.769195e+2 == Approx( xs[ elas ].value ) );
      CHECK( 8.311454e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1.264000e+5 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.830078e+1 == Approx( xs[ elas ].value ) );
      CHECK( 1.278686e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1.504700e+5 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.301758e+1 == Approx( xs[ elas ].value ) );
      CHECK( 1.061047e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      system.evaluate( 1.779400e+5 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.010266e+1 == Approx( xs[ elas ].value ) );
      CHECK( 1.099364e+1 == Approx( xs[ capt ].value ) );
      xs.clear();
    } // THEN
  } // GIVEN
} // SCENARIO
