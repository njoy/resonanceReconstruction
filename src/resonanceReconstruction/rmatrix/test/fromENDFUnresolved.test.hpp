std::string Na22();

SCENARIO( "fromENDF - legacy unresolved resonances" ) {

  GIVEN( "valid ENDF data for Na22" ) {

    std::string string = Na22();
    auto begin = string.begin();
    auto end = string.end();
    long lineNumber = 1;

    njoy::ENDFtk::HeadRecord head( begin, end, lineNumber );
    njoy::ENDFtk::section::Type< 2, 151 > endf( head, begin, end, lineNumber, 1122 );
    ResonanceRange endfResonanceRange = endf.isotopes().front().resonanceRanges().front();

    auto resonances = fromENDF( endfResonanceRange, neutronMass, elementaryCharge, ParticleID( "n" ), ParticleID( "Na22" ) );

    THEN( "the appropriate CompoundSystem is returned" ) {

      auto compoundsystem = std::get< legacy::unresolved::CompoundSystem >( resonances.compoundSystem() );

      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // content verification
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      // spin groups
      auto spingroups = compoundsystem.spinGroups();
      CHECK( 12 == spingroups.size() );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup0 = spingroups[0];
    } // THEN

    THEN( "cross sections can be reconstructed" ) {

      ReactionID elas( "n,Na22->n,Na22" );
      ReactionID capt( "n,Na22->capture" );
      std::map< ReactionID, CrossSection > xs;

      xs = resonances( 1.5e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 9.3027405025178052 == Approx( xs[ elas ].value ) );
      CHECK( 4.0639184850639387E-002 == Approx( xs[ capt ].value ) );

      xs = resonances( 2e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 8.7051878520691872 == Approx( xs[ elas ].value ) );
      CHECK( 3.2650798234536135E-002 == Approx( xs[ capt ].value ) );

      xs = resonances( 3e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 8.0164891832868452 == Approx( xs[ elas ].value ) );
      CHECK( 2.3797442539692449E-002 == Approx( xs[ capt ].value ) );

      xs = resonances( 4e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.6215788391323818 == Approx( xs[ elas ].value ) );
      CHECK( 1.8870640640168301E-002 == Approx( xs[ capt ].value ) );

      xs = resonances( 5e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.3617143134432972 == Approx( xs[ elas ].value ) );
      CHECK( 1.5670838905457892E-002 == Approx( xs[ capt ].value ) );

      xs = resonances( 6e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.1763623147011355 == Approx( xs[ elas ].value ) );
      CHECK( 1.3406527410576555E-002 == Approx( xs[ capt ].value ) );

      xs = resonances( 7e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.0369090178438851 == Approx( xs[ elas ].value ) );
      CHECK( 1.1714340063536310E-002 == Approx( xs[ capt ].value ) );

      xs = resonances( 8e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 6.9279078516440951 == Approx( xs[ elas ].value ) );
      CHECK( 1.0400317466079910E-002 == Approx( xs[ capt ].value ) );

      xs = resonances( 9e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 6.8402239876247961 == Approx( xs[ elas ].value ) );
      CHECK( 9.3501533213065698E-003 == Approx( xs[ capt ].value ) );

      xs = resonances( 1e+5 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 6.7680817547271817 == Approx( xs[ elas ].value ) );
      CHECK( 8.4916589760332239E-003 == Approx( xs[ capt ].value ) );
    } // THEN
  } // GIVEN
} // SCENARIO

std::string Na22() {

  // Fe54 ENDF/B-VIII.0 LRU=2 resonance evaluation

  return
    " 1.102200+4 2.180550+1          0          0          1          01122 2151     \n"
    " 1.102200+4 1.000000+0          0          0          1          01122 2151     \n"
    " 1.500000+4 1.000000+5          2          2          0          01122 2151     \n"
    " 3.000000+0 5.700000-1          0          0          3          01122 2151     \n"
    " 2.180550+1 0.000000+0          0          0          2          01122 2151     \n"
    " 2.500000+0 0.000000+0          2          0         66         101122 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 1.000000+0 0.000000+0 0.000000+01122 2151     \n"
    " 1.500000+4 5.579370+4 0.000000+0 7.476360+0 1.081650+0 0.000000+01122 2151     \n"
    " 2.000000+4 5.579370+4 0.000000+0 7.476360+0 1.081650+0 0.000000+01122 2151     \n"
    " 3.000000+4 5.579370+4 0.000000+0 7.476360+0 1.081650+0 0.000000+01122 2151     \n"
    " 4.000000+4 5.579370+4 0.000000+0 7.476360+0 1.081650+0 0.000000+01122 2151     \n"
    " 5.000000+4 5.579370+4 0.000000+0 7.476360+0 1.081650+0 0.000000+01122 2151     \n"
    " 6.000000+4 5.579370+4 0.000000+0 7.476360+0 1.081650+0 0.000000+01122 2151     \n"
    " 7.000000+4 5.579370+4 0.000000+0 7.476360+0 1.081650+0 0.000000+01122 2151     \n"
    " 8.000000+4 5.579370+4 0.000000+0 7.476360+0 1.081650+0 0.000000+01122 2151     \n"
    " 9.000000+4 5.579370+4 0.000000+0 7.476360+0 1.081650+0 0.000000+01122 2151     \n"
    " 1.000000+5 5.579370+4 0.000000+0 7.476360+0 1.081650+0 0.000000+01122 2151     \n"
    " 3.500000+0 0.000000+0          2          0         66         101122 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 1.000000+0 0.000000+0 0.000000+01122 2151     \n"
    " 1.500000+4 6.840720+4 0.000000+0 9.166560+0 1.081650+0 0.000000+01122 2151     \n"
    " 2.000000+4 6.840720+4 0.000000+0 9.166560+0 1.081650+0 0.000000+01122 2151     \n"
    " 3.000000+4 6.840720+4 0.000000+0 9.166560+0 1.081650+0 0.000000+01122 2151     \n"
    " 4.000000+4 6.840720+4 0.000000+0 9.166560+0 1.081650+0 0.000000+01122 2151     \n"
    " 5.000000+4 6.840720+4 0.000000+0 9.166560+0 1.081650+0 0.000000+01122 2151     \n"
    " 6.000000+4 6.840720+4 0.000000+0 9.166560+0 1.081650+0 0.000000+01122 2151     \n"
    " 7.000000+4 6.840720+4 0.000000+0 9.166560+0 1.081650+0 0.000000+01122 2151     \n"
    " 8.000000+4 6.840720+4 0.000000+0 9.166560+0 1.081650+0 0.000000+01122 2151     \n"
    " 9.000000+4 6.840720+4 0.000000+0 9.166560+0 1.081650+0 0.000000+01122 2151     \n"
    " 1.000000+5 6.840720+4 0.000000+0 9.166560+0 1.081650+0 0.000000+01122 2151     \n"
    " 2.180550+1 0.000000+0          1          0          4          01122 2151     \n"
    " 1.500000+0 0.000000+0          2          0         66         101122 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 1.000000+0 0.000000+0 0.000000+01122 2151     \n"
    " 1.500000+4 5.891260+4 0.000000+0 1.643660+1 5.444310+0 0.000000+01122 2151     \n"
    " 2.000000+4 5.891260+4 0.000000+0 1.643660+1 5.444310+0 0.000000+01122 2151     \n"
    " 3.000000+4 5.891260+4 0.000000+0 1.643660+1 5.444310+0 0.000000+01122 2151     \n"
    " 4.000000+4 5.891260+4 0.000000+0 1.643660+1 5.444310+0 0.000000+01122 2151     \n"
    " 5.000000+4 5.891260+4 0.000000+0 1.643660+1 5.444310+0 0.000000+01122 2151     \n"
    " 6.000000+4 5.891260+4 0.000000+0 1.643660+1 5.444310+0 0.000000+01122 2151     \n"
    " 7.000000+4 5.891260+4 0.000000+0 1.643660+1 5.444310+0 0.000000+01122 2151     \n"
    " 8.000000+4 5.891260+4 0.000000+0 1.643660+1 5.444310+0 0.000000+01122 2151     \n"
    " 9.000000+4 5.891260+4 0.000000+0 1.643660+1 5.444310+0 0.000000+01122 2151     \n"
    " 1.000000+5 5.891260+4 0.000000+0 1.643660+1 5.444310+0 0.000000+01122 2151     \n"
    " 2.500000+0 0.000000+0          2          0         66         101122 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 2.000000+0 0.000000+0 0.000000+01122 2151     \n"
    " 1.500000+4 5.579370+4 0.000000+0 3.113290+1 5.444310+0 0.000000+01122 2151     \n"
    " 2.000000+4 5.579370+4 0.000000+0 3.113290+1 5.444310+0 0.000000+01122 2151     \n"
    " 3.000000+4 5.579370+4 0.000000+0 3.113290+1 5.444310+0 0.000000+01122 2151     \n"
    " 4.000000+4 5.579370+4 0.000000+0 3.113290+1 5.444310+0 0.000000+01122 2151     \n"
    " 5.000000+4 5.579370+4 0.000000+0 3.113290+1 5.444310+0 0.000000+01122 2151     \n"
    " 6.000000+4 5.579370+4 0.000000+0 3.113290+1 5.444310+0 0.000000+01122 2151     \n"
    " 7.000000+4 5.579370+4 0.000000+0 3.113290+1 5.444310+0 0.000000+01122 2151     \n"
    " 8.000000+4 5.579370+4 0.000000+0 3.113290+1 5.444310+0 0.000000+01122 2151     \n"
    " 9.000000+4 5.579370+4 0.000000+0 3.113290+1 5.444310+0 0.000000+01122 2151     \n"
    " 1.000000+5 5.579370+4 0.000000+0 3.113290+1 5.444310+0 0.000000+01122 2151     \n"
    " 3.500000+0 0.000000+0          2          0         66         101122 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 2.000000+0 0.000000+0 0.000000+01122 2151     \n"
    " 1.500000+4 6.840720+4 0.000000+0 3.817120+1 5.444310+0 0.000000+01122 2151     \n"
    " 2.000000+4 6.840720+4 0.000000+0 3.817120+1 5.444310+0 0.000000+01122 2151     \n"
    " 3.000000+4 6.840720+4 0.000000+0 3.817120+1 5.444310+0 0.000000+01122 2151     \n"
    " 4.000000+4 6.840720+4 0.000000+0 3.817120+1 5.444310+0 0.000000+01122 2151     \n"
    " 5.000000+4 6.840720+4 0.000000+0 3.817120+1 5.444310+0 0.000000+01122 2151     \n"
    " 6.000000+4 6.840720+4 0.000000+0 3.817120+1 5.444310+0 0.000000+01122 2151     \n"
    " 7.000000+4 6.840720+4 0.000000+0 3.817120+1 5.444310+0 0.000000+01122 2151     \n"
    " 8.000000+4 6.840720+4 0.000000+0 3.817120+1 5.444310+0 0.000000+01122 2151     \n"
    " 9.000000+4 6.840720+4 0.000000+0 3.817120+1 5.444310+0 0.000000+01122 2151     \n"
    " 1.000000+5 6.840720+4 0.000000+0 3.817120+1 5.444310+0 0.000000+01122 2151     \n"
    " 4.500000+0 0.000000+0          2          0         66         101122 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 1.000000+0 0.000000+0 0.000000+01122 2151     \n"
    " 1.500000+4 1.029520+5 0.000000+0 2.872350+1 5.444310+0 0.000000+01122 2151     \n"
    " 2.000000+4 1.029520+5 0.000000+0 2.872350+1 5.444310+0 0.000000+01122 2151     \n"
    " 3.000000+4 1.029520+5 0.000000+0 2.872350+1 5.444310+0 0.000000+01122 2151     \n"
    " 4.000000+4 1.029520+5 0.000000+0 2.872350+1 5.444310+0 0.000000+01122 2151     \n"
    " 5.000000+4 1.029520+5 0.000000+0 2.872350+1 5.444310+0 0.000000+01122 2151     \n"
    " 6.000000+4 1.029520+5 0.000000+0 2.872350+1 5.444310+0 0.000000+01122 2151     \n"
    " 7.000000+4 1.029520+5 0.000000+0 2.872350+1 5.444310+0 0.000000+01122 2151     \n"
    " 8.000000+4 1.029520+5 0.000000+0 2.872350+1 5.444310+0 0.000000+01122 2151     \n"
    " 9.000000+4 1.029520+5 0.000000+0 2.872350+1 5.444310+0 0.000000+01122 2151     \n"
    " 1.000000+5 1.029520+5 0.000000+0 2.872350+1 5.444310+0 0.000000+01122 2151     \n"
    " 2.180550+1 0.000000+0          2          0          6          01122 2151     \n"
    " 5.000000-1 0.000000+0          2          0         66         101122 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 1.000000+0 0.000000+0 0.000000+01122 2151     \n"
    " 1.500000+4 9.544600+4 0.000000+0 3.149720+1 1.081650+0 0.000000+01122 2151     \n"
    " 2.000000+4 9.544600+4 0.000000+0 3.149720+1 1.081650+0 0.000000+01122 2151     \n"
    " 3.000000+4 9.544600+4 0.000000+0 3.149720+1 1.081650+0 0.000000+01122 2151     \n"
    " 4.000000+4 9.544600+4 0.000000+0 3.149720+1 1.081650+0 0.000000+01122 2151     \n"
    " 5.000000+4 9.544600+4 0.000000+0 3.149720+1 1.081650+0 0.000000+01122 2151     \n"
    " 6.000000+4 9.544600+4 0.000000+0 3.149720+1 1.081650+0 0.000000+01122 2151     \n"
    " 7.000000+4 9.544600+4 0.000000+0 3.149720+1 1.081650+0 0.000000+01122 2151     \n"
    " 8.000000+4 9.544600+4 0.000000+0 3.149720+1 1.081650+0 0.000000+01122 2151     \n"
    " 9.000000+4 9.544600+4 0.000000+0 3.149720+1 1.081650+0 0.000000+01122 2151     \n"
    " 1.000000+5 9.544600+4 0.000000+0 3.149720+1 1.081650+0 0.000000+01122 2151     \n"
    " 1.500000+0 0.000000+0          2          0         66         101122 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 2.000000+0 0.000000+0 0.000000+01122 2151     \n"
    " 1.500000+4 5.891260+4 0.000000+0 3.888230+1 1.081650+0 0.000000+01122 2151     \n"
    " 2.000000+4 5.891260+4 0.000000+0 3.888230+1 1.081650+0 0.000000+01122 2151     \n"
    " 3.000000+4 5.891260+4 0.000000+0 3.888230+1 1.081650+0 0.000000+01122 2151     \n"
    " 4.000000+4 5.891260+4 0.000000+0 3.888230+1 1.081650+0 0.000000+01122 2151     \n"
    " 5.000000+4 5.891260+4 0.000000+0 3.888230+1 1.081650+0 0.000000+01122 2151     \n"
    " 6.000000+4 5.891260+4 0.000000+0 3.888230+1 1.081650+0 0.000000+01122 2151     \n"
    " 7.000000+4 5.891260+4 0.000000+0 3.888230+1 1.081650+0 0.000000+01122 2151     \n"
    " 8.000000+4 5.891260+4 0.000000+0 3.888230+1 1.081650+0 0.000000+01122 2151     \n"
    " 9.000000+4 5.891260+4 0.000000+0 3.888230+1 1.081650+0 0.000000+01122 2151     \n"
    " 1.000000+5 5.891260+4 0.000000+0 3.888230+1 1.081650+0 0.000000+01122 2151     \n"
    " 2.500000+0 0.000000+0          2          0         66         101122 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 2.000000+0 0.000000+0 0.000000+01122 2151     \n"
    " 1.500000+4 5.579370+4 0.000000+0 3.682380+1 1.081650+0 0.000000+01122 2151     \n"
    " 2.000000+4 5.579370+4 0.000000+0 3.682380+1 1.081650+0 0.000000+01122 2151     \n"
    " 3.000000+4 5.579370+4 0.000000+0 3.682380+1 1.081650+0 0.000000+01122 2151     \n"
    " 4.000000+4 5.579370+4 0.000000+0 3.682380+1 1.081650+0 0.000000+01122 2151     \n"
    " 5.000000+4 5.579370+4 0.000000+0 3.682380+1 1.081650+0 0.000000+01122 2151     \n"
    " 6.000000+4 5.579370+4 0.000000+0 3.682380+1 1.081650+0 0.000000+01122 2151     \n"
    " 7.000000+4 5.579370+4 0.000000+0 3.682380+1 1.081650+0 0.000000+01122 2151     \n"
    " 8.000000+4 5.579370+4 0.000000+0 3.682380+1 1.081650+0 0.000000+01122 2151     \n"
    " 9.000000+4 5.579370+4 0.000000+0 3.682380+1 1.081650+0 0.000000+01122 2151     \n"
    " 1.000000+5 5.579370+4 0.000000+0 3.682380+1 1.081650+0 0.000000+01122 2151     \n"
    " 3.500000+0 0.000000+0          2          0         66         101122 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 2.000000+0 0.000000+0 0.000000+01122 2151     \n"
    " 1.500000+4 6.840720+4 0.000000+0 4.514870+1 1.081650+0 0.000000+01122 2151     \n"
    " 2.000000+4 6.840720+4 0.000000+0 4.514870+1 1.081650+0 0.000000+01122 2151     \n"
    " 3.000000+4 6.840720+4 0.000000+0 4.514870+1 1.081650+0 0.000000+01122 2151     \n"
    " 4.000000+4 6.840720+4 0.000000+0 4.514870+1 1.081650+0 0.000000+01122 2151     \n"
    " 5.000000+4 6.840720+4 0.000000+0 4.514870+1 1.081650+0 0.000000+01122 2151     \n"
    " 6.000000+4 6.840720+4 0.000000+0 4.514870+1 1.081650+0 0.000000+01122 2151     \n"
    " 7.000000+4 6.840720+4 0.000000+0 4.514870+1 1.081650+0 0.000000+01122 2151     \n"
    " 8.000000+4 6.840720+4 0.000000+0 4.514870+1 1.081650+0 0.000000+01122 2151     \n"
    " 9.000000+4 6.840720+4 0.000000+0 4.514870+1 1.081650+0 0.000000+01122 2151     \n"
    " 1.000000+5 6.840720+4 0.000000+0 4.514870+1 1.081650+0 0.000000+01122 2151     \n"
    " 4.500000+0 0.000000+0          2          0         66         101122 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 2.000000+0 0.000000+0 0.000000+01122 2151     \n"
    " 1.500000+4 1.029520+5 0.000000+0 6.794820+1 1.081650+0 0.000000+01122 2151     \n"
    " 2.000000+4 1.029520+5 0.000000+0 6.794820+1 1.081650+0 0.000000+01122 2151     \n"
    " 3.000000+4 1.029520+5 0.000000+0 6.794820+1 1.081650+0 0.000000+01122 2151     \n"
    " 4.000000+4 1.029520+5 0.000000+0 6.794820+1 1.081650+0 0.000000+01122 2151     \n"
    " 5.000000+4 1.029520+5 0.000000+0 6.794820+1 1.081650+0 0.000000+01122 2151     \n"
    " 6.000000+4 1.029520+5 0.000000+0 6.794820+1 1.081650+0 0.000000+01122 2151     \n"
    " 7.000000+4 1.029520+5 0.000000+0 6.794820+1 1.081650+0 0.000000+01122 2151     \n"
    " 8.000000+4 1.029520+5 0.000000+0 6.794820+1 1.081650+0 0.000000+01122 2151     \n"
    " 9.000000+4 1.029520+5 0.000000+0 6.794820+1 1.081650+0 0.000000+01122 2151     \n"
    " 1.000000+5 1.029520+5 0.000000+0 6.794820+1 1.081650+0 0.000000+01122 2151     \n"
    " 5.500000+0 0.000000+0          2          0         66         101122 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 1.000000+0 0.000000+0 0.000000+01122 2151     \n"
    " 1.500000+4 1.857300+5 0.000000+0 6.129090+1 1.081650+0 0.000000+01122 2151     \n"
    " 2.000000+4 1.857300+5 0.000000+0 6.129090+1 1.081650+0 0.000000+01122 2151     \n"
    " 3.000000+4 1.857300+5 0.000000+0 6.129090+1 1.081650+0 0.000000+01122 2151     \n"
    " 4.000000+4 1.857300+5 0.000000+0 6.129090+1 1.081650+0 0.000000+01122 2151     \n"
    " 5.000000+4 1.857300+5 0.000000+0 6.129090+1 1.081650+0 0.000000+01122 2151     \n"
    " 6.000000+4 1.857300+5 0.000000+0 6.129090+1 1.081650+0 0.000000+01122 2151     \n"
    " 7.000000+4 1.857300+5 0.000000+0 6.129090+1 1.081650+0 0.000000+01122 2151     \n"
    " 8.000000+4 1.857300+5 0.000000+0 6.129090+1 1.081650+0 0.000000+01122 2151     \n"
    " 9.000000+4 1.857300+5 0.000000+0 6.129090+1 1.081650+0 0.000000+01122 2151     \n"
    " 1.000000+5 1.857300+5 0.000000+0 6.129090+1 1.081650+0 0.000000+01122 2151     \n"
    " 0.000000+0 0.000000+0          0          0          0          01122 2  0     \n";
}
