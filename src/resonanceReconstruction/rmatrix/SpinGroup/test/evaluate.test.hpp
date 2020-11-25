SCENARIO( "evaluate" ) {

  //! @todo add test with more than one entrance channel
  //! @todo add test with a charged particle channel

  GIVEN( "valid data for a SpinGroup with one eliminated capture channel "
         "and one elastic channel using the Reich Moore formalism" ) {

    // test based on Fe54 ENDF/B-VIII.0 LRF7 resonance evaluation
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required
    // cross section values extracted from NJOY2016.43

    // because the orbital angular momentum l = 0 for these SpinGroup,
    // SpinGroup< ReichMoore, ShiftFactor > and SpinGroup< ReichMoore, Constant >
    // should give the same results as B = S = 0

    // using SpinGroup< ReichMoore, ShiftFactor > is equivalent to NJOY2016's
    // LRF7 reconstruction

    // particles
    Particle photon( ParticleID( "g" ), 0.0 * daltons,
                     0.0 * elementary, 1., +1);
    Particle neutron( ParticleID( "n" ), neutronMass,
                      0.0 * elementary, 0.5, +1);
    Particle fe55( ParticleID( "Fe55" ), 5.446635e+1 * neutronMass,
                   26.0 * elementary, 0.0, +1);
    Particle fe54( ParticleID( "Fe54" ), 5.347624e+1 * neutronMass,
                   26.0 * elementary, 0.0, +1);

    // particle pairs
    ParticlePair in( neutron, fe54 );
    ParticlePair out( photon, fe55 );

    // channels
    Channel< Photon > capture( in, out, 0.0 * electronVolt, { 0, 0.0, 0.5, +1 },
                               { 0.0 * rootBarn },
                               0.0 );
    Channel< Neutron > elastic( in, in, 0.0 * electronVolt, { 0, 0.5, 0.5, +1 },
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

    SpinGroup< ReichMoore, ShiftFactor > group1( { elastic },
                                                 std::move( single ) );
    SpinGroup< ReichMoore, ShiftFactor > group2( { elastic },
                                                 std::move( multiple ) );
    SpinGroup< ReichMoore, Constant > group3( { elastic },
                                              std::move( single2 ) );
    SpinGroup< ReichMoore, Constant > group4( { elastic },
                                              std::move( multiple2 ) );

    ReactionID elas( "n,Fe54->n,Fe54" );
    ReactionID capt( "n,Fe54->capture" );

    THEN( "cross sections can be calculated for a single resonance using the "
          "ShiftFactor boundary condition" ) {

      Map< ReactionID, CrossSection > xs;
      group1.evaluate( 1e-5 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.575879e-1 == Approx( xs[ elas ].value ) );
      CHECK( 6.895037e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e-4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.575879e-1 == Approx( xs[ elas ].value ) );
      CHECK( 2.180403e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e-3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.575877e-1 == Approx( xs[ elas ].value ) );
      CHECK( 6.895040e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e-2 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.575860e-1 == Approx( xs[ elas ].value ) );
      CHECK( 2.180408e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e-1 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.575694e-1 == Approx( xs[ elas ].value ) );
      CHECK( 6.895214e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+0 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.574031e-1 == Approx( xs[ elas ].value ) );
      CHECK( 2.180961e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+1 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.557415e-1 == Approx( xs[ elas ].value ) );
      CHECK( 6.912724e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+2 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.392160e-1 == Approx( xs[ elas ].value ) );
      CHECK( 2.237319e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.932778e-2 == Approx( xs[ elas ].value ) );
      CHECK( 9.067251e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.342056e+1 == Approx( xs[ elas ].value ) );
      CHECK( 2.473549e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+5 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.979435e+0 == Approx( xs[ elas ].value ) );
      CHECK( 4.915667e-6 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+6 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.308817e+0 == Approx( xs[ elas ].value ) );
      CHECK( 1.343259e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 7.788000e+3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.424287e+2 == Approx( xs[ elas ].value ) );
      CHECK( 4.241421e-1 == Approx( xs[ capt ].value ) );
      xs.clear();
    } // THEN

    THEN( "cross sections can be calculated for multiple resonances using the "
          "ShiftFactor boundary condition" ) {

      Map< ReactionID, CrossSection > xs;
      group2.evaluate( 1e-5 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.781791e-2 == Approx( xs[ elas ].value ) );
      CHECK( 7.082909e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e-4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.781790e-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.239813e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e-3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.781781e-2 == Approx( xs[ elas ].value ) );
      CHECK( 7.082911e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e-2 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.781682e-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.239818e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e-1 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.780692e-2 == Approx( xs[ elas ].value ) );
      CHECK( 7.083086e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+0 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.770804e-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.240372e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+1 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.672107e-2 == Approx( xs[ elas ].value ) );
      CHECK( 7.100642e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+2 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 7.704194e-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.296876e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 7.098264e-3 == Approx( xs[ elas ].value ) );
      CHECK( 9.259109e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.062823e+1 == Approx( xs[ elas ].value ) );
      CHECK( 2.502603e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+5 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 5.330137e+0 == Approx( xs[ elas ].value ) );
      CHECK( 5.716074e-5 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+6 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.324263e+0 == Approx( xs[ elas ].value ) );
      CHECK( 3.695041e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 7.788000e+3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.424287e+2 == Approx( xs[ elas ].value ) );
      CHECK( 4.241421e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 5.287200e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.738763e+1 == Approx( xs[ elas ].value ) );
      CHECK( 5.099758e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 7.190500e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.390972e+1 == Approx( xs[ elas ].value ) );
      CHECK( 4.208797e-2 == Approx( xs[ capt ].value ) );
      xs.clear();
    } // THEN

    THEN( "cross sections can be calculated for a single resonance using the "
          "Constant boundary condition" ) {

      Map< ReactionID, CrossSection > xs;
      group3.evaluate( 1e-5 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.575879e-1 == Approx( xs[ elas ].value ) );
      CHECK( 6.895037e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e-4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.575879e-1 == Approx( xs[ elas ].value ) );
      CHECK( 2.180403e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e-3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.575877e-1 == Approx( xs[ elas ].value ) );
      CHECK( 6.895040e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e-2 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.575860e-1 == Approx( xs[ elas ].value ) );
      CHECK( 2.180408e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e-1 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.575694e-1 == Approx( xs[ elas ].value ) );
      CHECK( 6.895214e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+0 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.574031e-1 == Approx( xs[ elas ].value ) );
      CHECK( 2.180961e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+1 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.557415e-1 == Approx( xs[ elas ].value ) );
      CHECK( 6.912724e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+2 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.392160e-1 == Approx( xs[ elas ].value ) );
      CHECK( 2.237319e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.932778e-2 == Approx( xs[ elas ].value ) );
      CHECK( 9.067251e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.342056e+1 == Approx( xs[ elas ].value ) );
      CHECK( 2.473549e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+5 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.979435e+0 == Approx( xs[ elas ].value ) );
      CHECK( 4.915667e-6 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+6 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.308817e+0 == Approx( xs[ elas ].value ) );
      CHECK( 1.343259e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 7.788000e+3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.424287e+2 == Approx( xs[ elas ].value ) );
      CHECK( 4.241421e-1 == Approx( xs[ capt ].value ) );
      xs.clear();
    } // THEN

    THEN( "cross sections can be calculated for multiple resonances using the "
          "Constant boundary condition" ) {

      Map< ReactionID, CrossSection > xs;
      group4.evaluate( 1e-5 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.781791e-2 == Approx( xs[ elas ].value ) );
      CHECK( 7.082909e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e-4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.781790e-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.239813e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e-3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.781781e-2 == Approx( xs[ elas ].value ) );
      CHECK( 7.082911e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e-2 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.781682e-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.239818e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e-1 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.780692e-2 == Approx( xs[ elas ].value ) );
      CHECK( 7.083086e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+0 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.770804e-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.240372e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+1 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.672107e-2 == Approx( xs[ elas ].value ) );
      CHECK( 7.100642e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+2 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 7.704194e-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.296876e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 7.098264e-3 == Approx( xs[ elas ].value ) );
      CHECK( 9.259109e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.062823e+1 == Approx( xs[ elas ].value ) );
      CHECK( 2.502603e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+5 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 5.330137e+0 == Approx( xs[ elas ].value ) );
      CHECK( 5.716074e-5 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+6 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.324263e+0 == Approx( xs[ elas ].value ) );
      CHECK( 3.695041e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 7.788000e+3 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.424287e+2 == Approx( xs[ elas ].value ) );
      CHECK( 4.241421e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 5.287200e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.738763e+1 == Approx( xs[ elas ].value ) );
      CHECK( 5.099758e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 7.190500e+4 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.390972e+1 == Approx( xs[ elas ].value ) );
      CHECK( 4.208797e-2 == Approx( xs[ capt ].value ) );
      xs.clear();
    } // THEN
  } // GIVEN

  GIVEN( "valid data for a SpinGroup with one eliminated capture channel, "
         "one elastic channel and two fission channels" ) {

    // test based on Pu239 ENDF/B-VIII.0 LRF3 resonance evaluation (some of the
    // widths used in this test were negative in the original evaluation)
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required
    // cross section values extracted from NJOY2016.43 (note: LRF7 in NJOY2016
    // doesn't add potential scattering for missing J values)

    // because the oribital angular momentum l = 0 for these SpinGroup,
    // SpinGroup< ReichMoore, ShiftFactor > and SpinGroup< ReichMoore, Constant >
    // should give the same results

    // using SpinGroup< ReichMoore, ShiftFactor > is equivalent to NJOY2016's
    // LRF7 reconstruction

    // particles
    Particle photon( ParticleID( "g" ), 0.0 * daltons,
                     0.0 * elementary, 1., +1);
    Particle neutron( ParticleID( "n" ), neutronMass,
                      0.0 * elementary, 0.5, +1);
    Particle pu240( ParticleID( "Pu240" ), 2.379916e+2 * neutronMass,
                    94.0 * elementary, 0.5, +1);
    Particle pu239( ParticleID( "Pu239" ), 2.369986e+2 * neutronMass,
                    94.0 * elementary, 0.5, +1);

    // particle pairs
    ParticlePair in( neutron, pu239 );
    ParticlePair out1( photon, pu240 );
    ParticlePair out2( neutron, pu239, ParticlePairID( "fission" ) );

    // channels
    Channel< Photon > capture( in, out1, 0. * electronVolt, { 0, 0.0, 0.0, +1 },
                               { 0.0 * rootBarn },
                               0.0 );
    Channel< Neutron > elastic( in, in, 0. * electronVolt, { 0, 0.5, 0.0, +1 },
                                { 9.410000e-1 * rootBarn },
                                0.0 );
    Channel< Fission > fission1( in, out2, "fission1", 0. * electronVolt,
                                 { 0, 0.0, 0.0, +1 },
                                 { 9.410000e-1 * rootBarn },
                                 0.0 );
    Channel< Fission > fission2( in, out2, "fission2", 0. * electronVolt,
                                 { 0, 0.0, 0.0, +1 },
                                 { 9.410000e-1 * rootBarn },
                                 0.0 );

    // conversion from Gamma to gamma
    auto eGamma = [&] ( double width, const Energy& energy ) -> ReducedWidth {
      return std::sqrt( width / 2. / elastic.penetrability( energy ) ) *
             rootElectronVolt;
    };
    auto fGamma = [&] ( double width ) -> ReducedWidth {
      return std::sqrt( width / 2. ) * rootElectronVolt;
    };
    auto cGamma = [&] ( double width ) -> ReducedWidth {
      return std::sqrt( width / 2. ) * rootElectronVolt;
    };

    // single resonance table
    ResonanceTable single(
      { elastic.channelID(), fission1.channelID(), fission2.channelID() },
      { Resonance( 1.541700e+1 * electronVolt,
                   { eGamma( 2.056203e-3, 1.541700e+1 * electronVolt ),
                     fGamma( 1.093928e-6 ),
                     fGamma( 7.550000e-1 ) },
                   cGamma( 4.054259e-2 ) ) } );
    ResonanceTable single2 = single;

    // multiple resonance table
    ResonanceTable multiple(
      { elastic.channelID(), fission1.channelID(), fission2.channelID() },
      { Resonance( 1.541700e+1 * electronVolt,
                   { eGamma( 2.056203e-3, 1.541700e+1 * electronVolt ),
                     fGamma( 1.093928e-6 ),
                     fGamma( 7.550000e-1 ) },
                   cGamma( 4.054259e-2 ) ),
        Resonance( 3.232700e+1 * electronVolt,
                   { eGamma( 8.678823e-4, 3.232700e+1 * electronVolt ),
                     fGamma( 5.235058e-3 ),
                     fGamma( 1.279000e-1 ) },
                   cGamma( 4.182541e-2 ) ),
        Resonance( 4.753400e+1 * electronVolt,
                   { eGamma( 5.171861e-3, 4.753400e+1 * electronVolt ),
                     fGamma( 5.548812e-1 ),
                     fGamma( 1.274000e-7 ) },
                   cGamma( 2.938826e-2 ) ) } );
    ResonanceTable multiple2 = multiple;

    SpinGroup< ReichMoore, ShiftFactor >
        group1( { elastic, fission1, fission2 }, std::move( single ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group2( { elastic, fission1, fission2 }, std::move( multiple ) );
    SpinGroup< ReichMoore, Constant >
        group3( { elastic, fission1, fission2 }, std::move( single2 ) );
    SpinGroup< ReichMoore, Constant >
        group4( { elastic, fission1, fission2 }, std::move( multiple2 ) );

    ReactionID elas( "n,Pu239->n,Pu239" );
    ReactionID fiss( "n,Pu239->fission" );
    ReactionID capt( "n,Pu239->capture" );

    THEN( "cross sections can be calculated for a single resonance using the "
          "ShiftFactor boundary condition" ) {

      Map< ReactionID, CrossSection > xs;
      group1.evaluate( 1e-5 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.736137e+0 == Approx( xs[ elas ].value ) );
      CHECK( 8.627570e+1 == Approx( xs[ fiss ].value ) );
      CHECK( 4.632894e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e-4 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.736137e+0 == Approx( xs[ elas ].value ) );
      CHECK( 2.728309e+1 == Approx( xs[ fiss ].value ) );
      CHECK( 1.465067e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e-3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.736134e+0 == Approx( xs[ elas ].value ) );
      CHECK( 8.628677e+0 == Approx( xs[ fiss ].value ) );
      CHECK( 4.633489e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e-2 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.736107e+0 == Approx( xs[ elas ].value ) );
      CHECK( 2.731814e+0 == Approx( xs[ fiss ].value ) );
      CHECK( 1.466949e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e-1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.735839e+0 == Approx( xs[ elas ].value ) );
      CHECK( 8.740501e-1 == Approx( xs[ fiss ].value ) );
      CHECK( 4.693537e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+0 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.732981e+0 == Approx( xs[ elas ].value ) );
      CHECK( 3.119579e-1 == Approx( xs[ fiss ].value ) );
      CHECK( 1.675174e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.653310e+0 == Approx( xs[ elas ].value ) );
      CHECK( 6.955268e-1 == Approx( xs[ fiss ].value ) );
      CHECK( 3.734890e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+2 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.789797e+0 == Approx( xs[ elas ].value ) );
      CHECK( 9.069864e-4 == Approx( xs[ fiss ].value ) );
      CHECK( 4.870401e-5 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.778610e+0 == Approx( xs[ elas ].value ) );
      CHECK( 2.116762e-6 == Approx( xs[ fiss ].value ) );
      CHECK( 1.136674e-7 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 2e+3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.774325e+0 == Approx( xs[ elas ].value ) );
      CHECK( 3.684030e-7 == Approx( xs[ fiss ].value ) );
      CHECK( 1.978277e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1.541700e+1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 3.050416e+0 == Approx( xs[ elas ].value ) );
      CHECK( 1.039122e+2 == Approx( xs[ fiss ].value ) );
      CHECK( 5.579952e+0 == Approx( xs[ capt ].value ) );
      xs.clear();
    } // THEN

    THEN( "cross sections can be calculated for multiple resonances using the "
          "ShiftFactor boundary condition" ) {

      Map< ReactionID, CrossSection > xs;
      group2.evaluate( 1e-5 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.708738e+0 == Approx( xs[ elas ].value ) );
      CHECK( 1.164550e+2 == Approx( xs[ fiss ].value ) );
      CHECK( 5.455577e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e-4 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.708738e+0 == Approx( xs[ elas ].value ) );
      CHECK( 3.682667e+1 == Approx( xs[ fiss ].value ) );
      CHECK( 1.725223e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e-3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.708734e+0 == Approx( xs[ elas ].value ) );
      CHECK( 1.164683e+1 == Approx( xs[ fiss ].value ) );
      CHECK( 5.456212e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e-2 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.708702e+0 == Approx( xs[ elas ].value ) );
      CHECK( 3.686888e+0 == Approx( xs[ fiss ].value ) );
      CHECK( 1.727233e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e-1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.708378e+0 == Approx( xs[ elas ].value ) );
      CHECK( 1.178142e+0 == Approx( xs[ fiss ].value ) );
      CHECK( 5.520311e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+0 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.704945e+0 == Approx( xs[ elas ].value ) );
      CHECK( 4.151231e-1 == Approx( xs[ fiss ].value ) );
      CHECK( 1.948731e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.618127e+0 == Approx( xs[ elas ].value ) );
      CHECK( 7.902586e-1 == Approx( xs[ fiss ].value ) );
      CHECK( 3.880615e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+2 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.812236e+0 == Approx( xs[ elas ].value ) );
      CHECK( 4.139402e-3 == Approx( xs[ fiss ].value ) );
      CHECK( 2.030171e-4 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.779884e+0 == Approx( xs[ elas ].value ) );
      CHECK( 5.781516e-6 == Approx( xs[ fiss ].value ) );
      CHECK( 2.751724e-7 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 2e+3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.774946e+0 == Approx( xs[ elas ].value ) );
      CHECK( 9.897643e-7 == Approx( xs[ fiss ].value ) );
      CHECK( 4.705688e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1.541700e+1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 3.018352e+0 == Approx( xs[ elas ].value ) );
      CHECK( 1.039298e+2 == Approx( xs[ fiss ].value ) );
      CHECK( 5.580787e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 3.232700e+1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 3.200788e+0 == Approx( xs[ elas ].value ) );
      CHECK( 7.594133e+1 == Approx( xs[ fiss ].value ) );
      CHECK( 2.384032e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 4.753400e+1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 3.826411e+0 == Approx( xs[ elas ].value ) );
      CHECK( 1.140901e+2 == Approx( xs[ fiss ].value ) );
      CHECK( 6.042426e+0 == Approx( xs[ capt ].value ) );
      xs.clear();
    } // THEN

    THEN( "cross sections can be calculated for a single resonance using the "
          "Constant boundary condition" ) {

      Map< ReactionID, CrossSection > xs;
      group3.evaluate( 1e-5 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.736137e+0 == Approx( xs[ elas ].value ) );
      CHECK( 8.627570e+1 == Approx( xs[ fiss ].value ) );
      CHECK( 4.632894e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e-4 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.736137e+0 == Approx( xs[ elas ].value ) );
      CHECK( 2.728309e+1 == Approx( xs[ fiss ].value ) );
      CHECK( 1.465067e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e-3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.736134e+0 == Approx( xs[ elas ].value ) );
      CHECK( 8.628677e+0 == Approx( xs[ fiss ].value ) );
      CHECK( 4.633489e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e-2 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.736107e+0 == Approx( xs[ elas ].value ) );
      CHECK( 2.731814e+0 == Approx( xs[ fiss ].value ) );
      CHECK( 1.466949e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e-1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.735839e+0 == Approx( xs[ elas ].value ) );
      CHECK( 8.740501e-1 == Approx( xs[ fiss ].value ) );
      CHECK( 4.693537e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+0 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.732981e+0 == Approx( xs[ elas ].value ) );
      CHECK( 3.119579e-1 == Approx( xs[ fiss ].value ) );
      CHECK( 1.675174e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.653310e+0 == Approx( xs[ elas ].value ) );
      CHECK( 6.955268e-1 == Approx( xs[ fiss ].value ) );
      CHECK( 3.734890e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+2 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.789797e+0 == Approx( xs[ elas ].value ) );
      CHECK( 9.069864e-4 == Approx( xs[ fiss ].value ) );
      CHECK( 4.870401e-5 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.778610e+0 == Approx( xs[ elas ].value ) );
      CHECK( 2.116762e-6 == Approx( xs[ fiss ].value ) );
      CHECK( 1.136674e-7 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 2e+3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.774325e+0 == Approx( xs[ elas ].value ) );
      CHECK( 3.684030e-7 == Approx( xs[ fiss ].value ) );
      CHECK( 1.978277e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1.541700e+1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 3.050416e+0 == Approx( xs[ elas ].value ) );
      CHECK( 1.039122e+2 == Approx( xs[ fiss ].value ) );
      CHECK( 5.579952e+0 == Approx( xs[ capt ].value ) );
      xs.clear();
    } // THEN

    THEN( "cross sections can be calculated for multiple resonances using the "
          "Constant boundary condition" ) {

      Map< ReactionID, CrossSection > xs;
      group4.evaluate( 1e-5 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.708738e+0 == Approx( xs[ elas ].value ) );
      CHECK( 1.164550e+2 == Approx( xs[ fiss ].value ) );
      CHECK( 5.455577e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e-4 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.708738e+0 == Approx( xs[ elas ].value ) );
      CHECK( 3.682667e+1 == Approx( xs[ fiss ].value ) );
      CHECK( 1.725223e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e-3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.708734e+0 == Approx( xs[ elas ].value ) );
      CHECK( 1.164683e+1 == Approx( xs[ fiss ].value ) );
      CHECK( 5.456212e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e-2 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.708702e+0 == Approx( xs[ elas ].value ) );
      CHECK( 3.686888e+0 == Approx( xs[ fiss ].value ) );
      CHECK( 1.727233e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e-1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.708378e+0 == Approx( xs[ elas ].value ) );
      CHECK( 1.178142e+0 == Approx( xs[ fiss ].value ) );
      CHECK( 5.520311e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+0 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.704945e+0 == Approx( xs[ elas ].value ) );
      CHECK( 4.151231e-1 == Approx( xs[ fiss ].value ) );
      CHECK( 1.948731e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.618127e+0 == Approx( xs[ elas ].value ) );
      CHECK( 7.902586e-1 == Approx( xs[ fiss ].value ) );
      CHECK( 3.880615e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+2 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.812236e+0 == Approx( xs[ elas ].value ) );
      CHECK( 4.139402e-3 == Approx( xs[ fiss ].value ) );
      CHECK( 2.030171e-4 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.779884e+0 == Approx( xs[ elas ].value ) );
      CHECK( 5.781516e-6 == Approx( xs[ fiss ].value ) );
      CHECK( 2.751724e-7 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 2e+3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.774946e+0 == Approx( xs[ elas ].value ) );
      CHECK( 9.897643e-7 == Approx( xs[ fiss ].value ) );
      CHECK( 4.705688e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1.541700e+1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 3.018352e+0 == Approx( xs[ elas ].value ) );
      CHECK( 1.039298e+2 == Approx( xs[ fiss ].value ) );
      CHECK( 5.580787e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 3.232700e+1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 3.200788e+0 == Approx( xs[ elas ].value ) );
      CHECK( 7.594133e+1 == Approx( xs[ fiss ].value ) );
      CHECK( 2.384032e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 4.753400e+1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 3.826411e+0 == Approx( xs[ elas ].value ) );
      CHECK( 1.140901e+2 == Approx( xs[ fiss ].value ) );
      CHECK( 6.042426e+0 == Approx( xs[ capt ].value ) );
      xs.clear();
    } // THEN
  } // GIVEN

  GIVEN( "valid data for a SpinGroup with a resonance at a negative energy" ) {

    // test based on Pu239 ENDF/B-VIII.0 LRF3 resonance evaluation (some of the
    // widths used in this test were negative in the original evaluation)
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required
    // cross section values extracted from NJOY2016.43 (note: LRF7 in NJOY2016
    // doesn't add potential scattering for missing J values)

    // particles
    Particle photon( ParticleID( "g" ), 0.0 * daltons, 0.0 * elementary, 1., +1);
    Particle neutron( ParticleID( "n" ), neutronMass, 0.0 * elementary, 0.5, +1);
    Particle pu240( ParticleID( "Pu240" ), 2.379916e+2 * neutronMass,
                    94.0 * elementary, 0.5, +1);
    Particle pu239( ParticleID( "Pu239" ), 2.369986e+2 * neutronMass,
                    94.0 * elementary, 0.5, +1);

    // particle pairs
    ParticlePair in( neutron, pu239 );
    ParticlePair out1( photon, pu240 );
    ParticlePair out2( neutron, pu239, ParticlePairID( "fission" ) );

    // channels
    Channel< Photon > capture( in, out1, 0. * electronVolt, { 0, 0.0, 0.0, +1 },
                               { 0.0 * rootBarn },
                               0.0 );
    Channel< Neutron > elastic( in, in, 0. * electronVolt, { 0, 0.5, 0.0, +1 },
                                { 9.410000e-1 * rootBarn },
                                0.0 );
    Channel< Fission > fission1( in, out2, "fission1", 0. * electronVolt,
                                 { 0, 0.0, 0.0, +1 },
                                 { 9.410000e-1 * rootBarn },
                                 0.0 );
    Channel< Fission > fission2( in, out2, "fission2", 0. * electronVolt,
                                 { 0, 0.0, 0.0, +1 },
                                 { 9.410000e-1 * rootBarn },
                                 0.0 );

    // conversion from Gamma to gamma
    auto eGamma = [&] ( double width, const Energy& energy ) -> ReducedWidth {
      return std::sqrt( width / 2. / elastic.penetrability( energy ) ) *
             rootElectronVolt;
    };
    auto fGamma = [&] ( double width ) -> ReducedWidth {
      return std::sqrt( width / 2. ) * rootElectronVolt;
    };
    auto cGamma = [&] ( double width ) -> ReducedWidth {
      return std::sqrt( width / 2. ) * rootElectronVolt;
    };

    // multiple resonance table
    ResonanceTable table(
      { elastic.channelID(), fission1.channelID(), fission2.channelID() },
      { Resonance( -1.541700e+1 * electronVolt,
                   { eGamma( 2.056203e-3, 1.541700e+1 * electronVolt ),
                     fGamma( 1.093928e-6 ),
                     fGamma( 7.550000e-1 ) },
                   cGamma( 4.054259e-2 ) ),
        Resonance( 3.232700e+1 * electronVolt,
                   { eGamma( 8.678823e-4, 3.232700e+1 * electronVolt ),
                     fGamma( 5.235058e-3 ),
                     fGamma( 1.279000e-1 ) },
                   cGamma( 4.182541e-2 ) ),
        Resonance( 4.753400e+1 * electronVolt,
                   { eGamma( 5.171861e-3, 4.753400e+1 * electronVolt ),
                     fGamma( 5.548812e-1 ),
                     fGamma( 1.274000e-7 ) },
                   cGamma( 2.938826e-2 ) ) } );

    SpinGroup< ReichMoore, ShiftFactor >
        group( { elastic, fission1, fission2 }, std::move( table ) );

    ReactionID elas( "n,Pu239->n,Pu239" );
    ReactionID fiss( "n,Pu239->fission" );
    ReactionID capt( "n,Pu239->capture" );

    THEN( "cross sections can be calculated" ) {

      Map< ReactionID, CrossSection > xs;
      group.evaluate( 1e-5 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.800024e+0 == Approx( xs[ elas ].value ) );
      CHECK( 7.969823e+1 == Approx( xs[ fiss ].value ) );
      CHECK( 5.456509e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e-4 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.800023e+0 == Approx( xs[ elas ].value ) );
      CHECK( 2.520251e+1 == Approx( xs[ fiss ].value ) );
      CHECK( 1.725484e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e-3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.800020e+0 == Approx( xs[ elas ].value ) );
      CHECK( 7.968829e+0 == Approx( xs[ fiss ].value ) );
      CHECK( 5.455955e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e-2 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.799987e+0 == Approx( xs[ elas ].value ) );
      CHECK( 2.517113e+0 == Approx( xs[ fiss ].value ) );
      CHECK( 1.723733e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e-1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.799661e+0 == Approx( xs[ elas ].value ) );
      CHECK( 7.870531e-1 == Approx( xs[ fiss ].value ) );
      CHECK( 5.401134e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e+0 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.796554e+0 == Approx( xs[ elas ].value ) );
      CHECK( 2.233478e-1 == Approx( xs[ fiss ].value ) );
      CHECK( 1.566002e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e+1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.773379e+0 == Approx( xs[ elas ].value ) );
      CHECK( 3.532512e-2 == Approx( xs[ fiss ].value ) );
      CHECK( 3.181981e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e+2 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.809990e+0 == Approx( xs[ elas ].value ) );
      CHECK( 3.583571e-3 == Approx( xs[ fiss ].value ) );
      CHECK( 1.804712e-4 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e+3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.779863e+0 == Approx( xs[ elas ].value ) );
      CHECK( 5.625611e-6 == Approx( xs[ fiss ].value ) );
      CHECK( 2.683740e-7 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 2e+3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.774941e+0 == Approx( xs[ elas ].value ) );
      CHECK( 9.760313e-7 == Approx( xs[ fiss ].value ) );
      CHECK( 4.645620e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 3.232700e+1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 3.200741e+0 == Approx( xs[ elas ].value ) );
      CHECK( 7.593905e+1 == Approx( xs[ fiss ].value ) );
      CHECK( 2.384228e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 4.753400e+1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 3.815638e+0 == Approx( xs[ elas ].value ) );
      CHECK( 1.140803e+2 == Approx( xs[ fiss ].value ) );
      CHECK( 6.042125e+0 == Approx( xs[ capt ].value ) );
      xs.clear();
    } // THEN
  } // GIVEN

  GIVEN( "valid data for a SpinGroup with resonances using negative widths" ) {

    // test based on Pu239 ENDF/B-VIII.0 LRF3 resonance evaluation (all
    // widths used in this test are from the original evaluation)
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required
    // cross section values extracted from NJOY2016.43 (note: LRF7 in NJOY2016
    // doesn't add potential scattering for missing J values)

    // particles
    Particle photon( ParticleID( "g" ), 0.0 * daltons,
                     0.0 * elementary, 1., +1);
    Particle neutron( ParticleID( "n" ), neutronMass,
                      0.0 * elementary, 0.5, +1);
    Particle pu240( ParticleID( "Pu240" ), 2.379916e+2 * neutronMass,
                    94.0 * elementary, 0.5, +1);
    Particle pu239( ParticleID( "Pu239" ), 2.369986e+2 * neutronMass,
                    94.0 * elementary, 0.5, +1);

    // particle pairs
    ParticlePair in( neutron, pu239 );
    ParticlePair out1( photon, pu240 );
    ParticlePair out2( neutron, pu239, ParticlePairID( "fission" ) );

    // channels
    Channel< Photon > capture( in, out1, 0. * electronVolt, { 0, 0.0, 0.0, +1 },
                               { 0.0 * rootBarn },
                               0.0 );
    Channel< Neutron > elastic( in, in, 0. * electronVolt, { 0, 0.5, 0.0, +1 },
                                { 9.410000e-1 * rootBarn },
                                0.0 );
    Channel< Fission > fission1( in, out2, "fission1", 0. * electronVolt,
                                 { 0, 0.0, 0.0, +1 },
                                 { 9.410000e-1 * rootBarn },
                                 0.0 );
    Channel< Fission > fission2( in, out2, "fission2", 0. * electronVolt,
                                 { 0, 0.0, 0.0, +1 },
                                 { 9.410000e-1 * rootBarn },
                                 0.0 );

    // conversion from Gamma to gamma
    auto eGamma = [&] ( double width, const Energy& energy ) -> ReducedWidth {
      double sign = 1.0;
      if ( width < 0 ) sign = -1.0;
      return sign * std::sqrt( std::abs( width ) / 2. /
                               elastic.penetrability( energy ) ) *
             rootElectronVolt;
    };
    auto fGamma = [&] ( double width ) -> ReducedWidth {
      double sign = 1.0;
      if ( width < 0 ) sign = -1.0;
      return sign *std::sqrt( std::abs( width ) / 2. ) * rootElectronVolt;
    };
    auto cGamma = [&] ( double width ) -> ReducedWidth {
      return std::sqrt( width / 2. ) * rootElectronVolt;
    };

    // multiple resonance table
    ResonanceTable table(
      { elastic.channelID(), fission1.channelID(), fission2.channelID() },
      { Resonance( 1.541700e+1 * electronVolt,
                   { eGamma( 2.056203e-3, 1.541700e+1 * electronVolt ),
                     fGamma( -1.093928e-6 ),
                     fGamma( 7.550000e-1 ) },
                   cGamma( 4.054259e-2 ) ),
        Resonance( 3.232700e+1 * electronVolt,
                   { eGamma( 8.678823e-4, 3.232700e+1 * electronVolt ),
                     fGamma( 5.235058e-3 ),
                     fGamma( -1.279000e-1 ) },
                   cGamma( 4.182541e-2 ) ),
        Resonance( 4.753400e+1 * electronVolt,
                   { eGamma( 5.171861e-3, 4.753400e+1 * electronVolt ),
                     fGamma( 5.548812e-1 ),
                     fGamma( -1.274000e-7 ) },
                   cGamma( 2.938826e-2 ) ) } );

    SpinGroup< ReichMoore, ShiftFactor >
        group( { elastic, fission1, fission2 }, std::move( table ) );

    ReactionID elas( "n,Pu239->n,Pu239" );
    ReactionID fiss( "n,Pu239->fission" );
    ReactionID capt( "n,Pu239->capture" );

    THEN( "cross sections can be calculated" ) {

      Map< ReactionID, CrossSection > xs;
      group.evaluate( 1e-5 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.708724e+0 == Approx( xs[ elas ].value ) );
      CHECK( 7.968614e+1 == Approx( xs[ fiss ].value ) );
      CHECK( 5.456366e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e-4 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.708723e+0 == Approx( xs[ elas ].value ) );
      CHECK( 2.519925e+1 == Approx( xs[ fiss ].value ) );
      CHECK( 1.725473e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e-3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.708720e+0 == Approx( xs[ elas ].value ) );
      CHECK( 7.969600e+0 == Approx( xs[ fiss ].value ) );
      CHECK( 5.457001e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e-2 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.708688e+0 == Approx( xs[ elas ].value ) );
      CHECK( 2.523044e+0 == Approx( xs[ fiss ].value ) );
      CHECK( 1.727483e-1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e-1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.708364e+0 == Approx( xs[ elas ].value ) );
      CHECK( 8.069116e-1 == Approx( xs[ fiss ].value ) );
      CHECK( 5.521116e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e+0 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.704930e+0 == Approx( xs[ elas ].value ) );
      CHECK( 2.868371e-1 == Approx( xs[ fiss ].value ) );
      CHECK( 1.949035e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e+1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.617986e+0 == Approx( xs[ elas ].value ) );
      CHECK( 6.396225e-1 == Approx( xs[ fiss ].value ) );
      CHECK( 3.881964e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e+2 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.812237e+0 == Approx( xs[ elas ].value ) );
      CHECK( 3.121466e-3 == Approx( xs[ fiss ].value ) );
      CHECK( 2.030192e-4 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1e+3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.779884e+0 == Approx( xs[ elas ].value ) );
      CHECK( 3.851580e-6 == Approx( xs[ fiss ].value ) );
      CHECK( 2.751724e-7 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 2e+3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.774946e+0 == Approx( xs[ elas ].value ) );
      CHECK( 6.568294e-7 == Approx( xs[ fiss ].value ) );
      CHECK( 4.705688e-8 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 1.541700e+1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.983437e+0 == Approx( xs[ elas ].value ) );
      CHECK( 1.039338e+2 == Approx( xs[ fiss ].value ) );
      CHECK( 5.584260e+0 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 3.232700e+1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 3.360875e+0 == Approx( xs[ elas ].value ) );
      CHECK( 7.591453e+1 == Approx( xs[ fiss ].value ) );
      CHECK( 2.389296e+1 == Approx( xs[ capt ].value ) );
      xs.clear();

      group.evaluate( 4.753400e+1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 3.826384e+0 == Approx( xs[ elas ].value ) );
      CHECK( 1.140768e+2 == Approx( xs[ fiss ].value ) );
      CHECK( 6.042620e+0 == Approx( xs[ capt ].value ) );
      xs.clear();
    } // THEN
  } // GIVEN

  GIVEN( "valid data for a SpinGroup with one eliminated capture channel, "
         "one elastic channel and a proton channel" ) {

    // test based on Cl35 ENDF/B-VIII.0 LRF7 resonance evaluation
    // data given in Gamma = 2 gamma^2 P(Er) so conversion is required
    // cross section values extracted from NJOY2016.43

    // because the oribital angular momentum l = 0 for these SpinGroup,
    // SpinGroup< ReichMoore, ShiftFactor > and SpinGroup< ReichMoore, Constant >
    // should give the same results

    // using SpinGroup< ReichMoore, ShiftFactor > is equivalent to NJOY2016's
    // LRF7 reconstruction

    // particles
    Particle photon( ParticleID( "g" ), 0.0 * daltons,
                     0.0 * elementary, 1., +1);
    Particle neutron( ParticleID( "n" ), neutronMass,
                      0.0 * elementary, 0.5, +1);
    Particle proton( ParticleID( "p" ), 9.986235e-1 * neutronMass,
                     elementary, 0.5, +1);
    Particle cl36( ParticleID( "Cl36" ), 3.565932e+1 * neutronMass,
                   17.0 * elementary, 0., +1);
    Particle cl35( ParticleID( "Cl35" ), 3.466845e+1 * neutronMass,
                   17.0 * elementary, 1.5, +1);
    Particle s36( ParticleID( "S36" ), 3.466863e+1 * neutronMass,
                  16.0 * elementary, 1.5, +1);

    // particle pairs
    ParticlePair in( neutron, cl35 );
    ParticlePair out1( photon, cl36 );
    ParticlePair out2( proton, s36 );

    // channels
    Channel< Photon > capture( in, out1, 0. * electronVolt, { 0, 0.0, 1.0, +1 },
                               { 0.0 * rootBarn },
                               0.0 );
    Channel< Neutron > elastic( in, in, 0. * electronVolt, { 0, 1.0, 1.0, +1 },
                                { 4.822220e-1 * rootBarn,
                                  3.667980e-1 * rootBarn },
                                0.0 );
    Channel< ChargedParticle > protonemission( in, out2,
                                               6.152200e+5 * electronVolt,
                                               { 0, 1.0, 1.0, +1 },
                                               { 4.822220e-1 * rootBarn,
                                                 3.667980e-1 * rootBarn },
                                               0.0 );

    // conversion from Gamma to gamma
    auto eGamma = [&] ( double width, const Energy& energy ) -> ReducedWidth {
      return std::sqrt( width / 2. / elastic.penetrability( energy ) ) *
             rootElectronVolt;
    };
    auto pGamma = [&] ( double width, const Energy& energy ) -> ReducedWidth {
      return std::sqrt( width / 2. / protonemission.penetrability( energy ) ) *
             rootElectronVolt;
    };
    auto cGamma = [&] ( double width ) -> ReducedWidth {
      return std::sqrt( width / 2. ) * rootElectronVolt;
    };

    // single resonance table
    ResonanceTable single(
      { elastic.channelID(), protonemission.channelID() },
      { Resonance( 6.823616e+4 * electronVolt,
                   { eGamma( 2.179040e+2, 6.823616e+4 * electronVolt ),
                     pGamma( 1.000000e-5, 6.823616e+4 * electronVolt ) },
                   cGamma( 3.933600e-1 ) ) } );
    ResonanceTable single2 = single;

    // multiple resonance table
    ResonanceTable multiple(
      { elastic.channelID(), protonemission.channelID() },
      { Resonance( 6.823616e+4 * electronVolt,
                   { eGamma( 2.179040e+2, 6.823616e+4 * electronVolt ),
                     pGamma( 1.000000e-5, 6.823616e+4 * electronVolt ) },
                   cGamma( 3.933600e-1 ) ),
        Resonance( 1.825230e+5 * electronVolt,
                   { eGamma( 1.759740e+3, 1.825230e+5 * electronVolt ),
                     pGamma( 4.000000e-1, 1.825230e+5 * electronVolt ) },
                   cGamma( 7.451500e-1 ) ),
        Resonance( 2.397427e+5 * electronVolt,
                   { eGamma( 2.685470e+2, 2.397427e+5 * electronVolt ),
                     pGamma( 0.0, 2.397427e+5 * electronVolt ) },
                   cGamma( 6.871600e-1 ) ) } );
    ResonanceTable multiple2 = multiple;

    SpinGroup< ReichMoore, ShiftFactor >
        group1( { elastic, protonemission }, std::move( single ) );
    SpinGroup< ReichMoore, ShiftFactor >
        group2( { elastic, protonemission }, std::move( multiple ) );
    SpinGroup< ReichMoore, Constant >
        group3( { elastic, protonemission }, std::move( single2 ) );
    SpinGroup< ReichMoore, Constant >
        group4( { elastic, protonemission }, std::move( multiple2 ) );

    ReactionID elas( "n,Cl35->n,Cl35" );
    ReactionID pemi( "n,Cl35->p,S36" );
    ReactionID capt( "n,Cl35->capture" );

    THEN( "cross sections can be calculated for a single resonance using the "
          "ShiftFactor boundary condition" ) {

      Map< ReactionID, CrossSection > xs;
      group1.evaluate( 1e-5 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.241509e-1 == Approx( xs[ elas ].value ) );
      CHECK( 5.576194e-8 == Approx( xs[ pemi ].value ) );
      CHECK( 5.758583e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e-4 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.241509e-1 == Approx( xs[ elas ].value ) );
      CHECK( 1.763347e-8 == Approx( xs[ pemi ].value ) );
      CHECK( 1.821024e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e-3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.241509e-1 == Approx( xs[ elas ].value ) );
      CHECK( 5.576194e-9 == Approx( xs[ pemi ].value ) );
      CHECK( 5.758583e-4 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e-2 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.241509e-1 == Approx( xs[ elas ].value ) );
      CHECK( 1.763348e-9 == Approx( xs[ pemi ].value ) );
      CHECK( 1.821024e-4 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e-1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.241508e-1 == Approx( xs[ elas ].value ) );
      CHECK( 5.57622e-10 == Approx( xs[ pemi ].value ) );
      CHECK( 5.758600e-5 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+0 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.241506e-1 == Approx( xs[ elas ].value ) );
      CHECK( 1.76343e-10 == Approx( xs[ pemi ].value ) );
      CHECK( 1.821077e-5 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.241482e-1 == Approx( xs[ elas ].value ) );
      CHECK( 5.57868e-11 == Approx( xs[ pemi ].value ) );
      CHECK( 5.760271e-6 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+2 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.241239e-1 == Approx( xs[ elas ].value ) );
      CHECK( 1.77124e-11 == Approx( xs[ pemi ].value ) );
      CHECK( 1.826373e-6 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.238793e-1 == Approx( xs[ elas ].value ) );
      CHECK( 5.83185e-12 == Approx( xs[ pemi ].value ) );
      CHECK( 5.931151e-7 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+4 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.212174e-1 == Approx( xs[ elas ].value ) );
      CHECK( 2.81649e-12 == Approx( xs[ pemi ].value ) );
      CHECK( 2.500110e-7 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+5 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.417078e-1 == Approx( xs[ elas ].value ) );
      CHECK( 1.01231e-12 == Approx( xs[ pemi ].value ) );
      CHECK( 2.657488e-7 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 1e+6 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 5.149891e-1 == Approx( xs[ elas ].value ) );
      CHECK( 1.06136e-11 == Approx( xs[ pemi ].value ) );
      CHECK( 9.76637e-11 == Approx( xs[ capt ].value ) );
      xs.clear();

      group1.evaluate( 6.823616e+4 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 1.447223e+1 == Approx( xs[ elas ].value ) );
      CHECK( 6.926440e-7 == Approx( xs[ pemi ].value ) );
      CHECK( 2.724583e-2 == Approx( xs[ capt ].value ) );
      xs.clear();
    } // THEN

    THEN( "cross sections can be calculated for multiple resonances using the "
          "ShiftFactor boundary condition" ) {

      Map< ReactionID, CrossSection > xs;
      group2.evaluate( 1e-5 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.043473e-1 == Approx( xs[ elas ].value ) );
      CHECK( 4.169719e-4 == Approx( xs[ pemi ].value ) );
      CHECK( 1.382266e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e-4 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.043473e-1 == Approx( xs[ elas ].value ) );
      CHECK( 1.318581e-4 == Approx( xs[ pemi ].value ) );
      CHECK( 4.371108e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e-3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.043473e-1 == Approx( xs[ elas ].value ) );
      CHECK( 4.169719e-5 == Approx( xs[ pemi ].value ) );
      CHECK( 1.382266e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e-2 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.043473e-1 == Approx( xs[ elas ].value ) );
      CHECK( 1.318581e-5 == Approx( xs[ pemi ].value ) );
      CHECK( 4.371109e-4 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e-1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.043472e-1 == Approx( xs[ elas ].value ) );
      CHECK( 4.169730e-6 == Approx( xs[ pemi ].value ) );
      CHECK( 1.382268e-4 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+0 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.043469e-1 == Approx( xs[ elas ].value ) );
      CHECK( 1.318616e-6 == Approx( xs[ pemi ].value ) );
      CHECK( 4.371189e-5 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.043436e-1 == Approx( xs[ elas ].value ) );
      CHECK( 4.170824e-7 == Approx( xs[ pemi ].value ) );
      CHECK( 1.382522e-5 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+2 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.043108e-1 == Approx( xs[ elas ].value ) );
      CHECK( 1.322079e-7 == Approx( xs[ pemi ].value ) );
      CHECK( 4.379210e-6 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.039801e-1 == Approx( xs[ elas ].value ) );
      CHECK( 4.281706e-8 == Approx( xs[ pemi ].value ) );
      CHECK( 1.408289e-6 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+4 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.004127e-1 == Approx( xs[ elas ].value ) );
      CHECK( 1.721283e-8 == Approx( xs[ pemi ].value ) );
      CHECK( 5.349224e-7 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+5 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 5.999024e-1 == Approx( xs[ elas ].value ) );
      CHECK( 7.662583e-8 == Approx( xs[ pemi ].value ) );
      CHECK( 6.497930e-7 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 1e+6 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 5.179658e-1 == Approx( xs[ elas ].value ) );
      CHECK( 7.263042e-8 == Approx( xs[ pemi ].value ) );
      CHECK( 1.452944e-9 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 6.823616e+4 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 1.447223e+1 == Approx( xs[ elas ].value ) );
      CHECK( 6.925359e-7 == Approx( xs[ pemi ].value ) );
      CHECK( 2.724593e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group2.evaluate( 2.397427e+5 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 3.680009e+0 == Approx( xs[ elas ].value ) );
      CHECK( 3.03460e-12 == Approx( xs[ pemi ].value ) );
      CHECK( 1.097562e-2 == Approx( xs[ capt ].value ) );
      xs.clear();
    } // THEN

    THEN( "cross sections can be calculated for a single resonance using the "
          "Constant boundary condition" ) {

      Map< ReactionID, CrossSection > xs;
      group3.evaluate( 1e-5 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.241509e-1 == Approx( xs[ elas ].value ) );
      CHECK( 5.576194e-8 == Approx( xs[ pemi ].value ) );
      CHECK( 5.758583e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e-4 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.241509e-1 == Approx( xs[ elas ].value ) );
      CHECK( 1.763347e-8 == Approx( xs[ pemi ].value ) );
      CHECK( 1.821024e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e-3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.241509e-1 == Approx( xs[ elas ].value ) );
      CHECK( 5.576194e-9 == Approx( xs[ pemi ].value ) );
      CHECK( 5.758583e-4 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e-2 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.241509e-1 == Approx( xs[ elas ].value ) );
      CHECK( 1.763348e-9 == Approx( xs[ pemi ].value ) );
      CHECK( 1.821024e-4 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e-1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.241508e-1 == Approx( xs[ elas ].value ) );
      CHECK( 5.57622e-10 == Approx( xs[ pemi ].value ) );
      CHECK( 5.758600e-5 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+0 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.241506e-1 == Approx( xs[ elas ].value ) );
      CHECK( 1.76343e-10 == Approx( xs[ pemi ].value ) );
      CHECK( 1.821077e-5 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.241482e-1 == Approx( xs[ elas ].value ) );
      CHECK( 5.57868e-11 == Approx( xs[ pemi ].value ) );
      CHECK( 5.760271e-6 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+2 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.241239e-1 == Approx( xs[ elas ].value ) );
      CHECK( 1.77124e-11 == Approx( xs[ pemi ].value ) );
      CHECK( 1.826373e-6 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.238793e-1 == Approx( xs[ elas ].value ) );
      CHECK( 5.83185e-12 == Approx( xs[ pemi ].value ) );
      CHECK( 5.931151e-7 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+4 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.212174e-1 == Approx( xs[ elas ].value ) );
      CHECK( 2.81649e-12 == Approx( xs[ pemi ].value ) );
      CHECK( 2.500110e-7 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+5 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.417078e-1 == Approx( xs[ elas ].value ) );
      CHECK( 1.01231e-12 == Approx( xs[ pemi ].value ) );
      CHECK( 2.657488e-7 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 1e+6 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 5.149891e-1 == Approx( xs[ elas ].value ) );
      CHECK( 1.06136e-11 == Approx( xs[ pemi ].value ) );
      CHECK( 9.76637e-11 == Approx( xs[ capt ].value ) );
      xs.clear();

      group3.evaluate( 6.823616e+4 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 1.447223e+1 == Approx( xs[ elas ].value ) );
      CHECK( 6.926440e-7 == Approx( xs[ pemi ].value ) );
      CHECK( 2.724583e-2 == Approx( xs[ capt ].value ) );
      xs.clear();
    } // THEN

    THEN( "cross sections can be calculated for multiple resonances using the "
          "Constant boundary condition" ) {

      Map< ReactionID, CrossSection > xs;
      group4.evaluate( 1e-5 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.043473e-1 == Approx( xs[ elas ].value ) );
      CHECK( 4.169719e-4 == Approx( xs[ pemi ].value ) );
      CHECK( 1.382266e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e-4 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.043473e-1 == Approx( xs[ elas ].value ) );
      CHECK( 1.318581e-4 == Approx( xs[ pemi ].value ) );
      CHECK( 4.371108e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e-3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.043473e-1 == Approx( xs[ elas ].value ) );
      CHECK( 4.169719e-5 == Approx( xs[ pemi ].value ) );
      CHECK( 1.382266e-3 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e-2 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.043473e-1 == Approx( xs[ elas ].value ) );
      CHECK( 1.318581e-5 == Approx( xs[ pemi ].value ) );
      CHECK( 4.371109e-4 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e-1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.043472e-1 == Approx( xs[ elas ].value ) );
      CHECK( 4.169730e-6 == Approx( xs[ pemi ].value ) );
      CHECK( 1.382268e-4 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+0 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.043469e-1 == Approx( xs[ elas ].value ) );
      CHECK( 1.318616e-6 == Approx( xs[ pemi ].value ) );
      CHECK( 4.371189e-5 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+1 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.043436e-1 == Approx( xs[ elas ].value ) );
      CHECK( 4.170824e-7 == Approx( xs[ pemi ].value ) );
      CHECK( 1.382522e-5 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+2 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.043108e-1 == Approx( xs[ elas ].value ) );
      CHECK( 1.322079e-7 == Approx( xs[ pemi ].value ) );
      CHECK( 4.379210e-6 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+3 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.039801e-1 == Approx( xs[ elas ].value ) );
      CHECK( 4.281706e-8 == Approx( xs[ pemi ].value ) );
      CHECK( 1.408289e-6 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+4 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.004127e-1 == Approx( xs[ elas ].value ) );
      CHECK( 1.721283e-8 == Approx( xs[ pemi ].value ) );
      CHECK( 5.349224e-7 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+5 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 5.999024e-1 == Approx( xs[ elas ].value ) );
      CHECK( 7.662583e-8 == Approx( xs[ pemi ].value ) );
      CHECK( 6.497930e-7 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 1e+6 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 5.179658e-1 == Approx( xs[ elas ].value ) );
      CHECK( 7.263042e-8 == Approx( xs[ pemi ].value ) );
      CHECK( 1.452944e-9 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 6.823616e+4 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 1.447223e+1 == Approx( xs[ elas ].value ) );
      CHECK( 6.925359e-7 == Approx( xs[ pemi ].value ) );
      CHECK( 2.724593e-2 == Approx( xs[ capt ].value ) );
      xs.clear();

      group4.evaluate( 2.397427e+5 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 3.680009e+0 == Approx( xs[ elas ].value ) );
      CHECK( 3.03460e-12 == Approx( xs[ pemi ].value ) );
      CHECK( 1.097562e-2 == Approx( xs[ capt ].value ) );
      xs.clear();
    } // THEN
  } // GIVEN*/
} // SCENARIO
