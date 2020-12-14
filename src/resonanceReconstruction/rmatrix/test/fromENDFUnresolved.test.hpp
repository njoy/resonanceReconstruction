std::string Na22();
std::string Pu239();
std::string Er167();
std::string Au197();
std::string Ag107();

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

      CHECK( false == resonances.isResolved() );
      CHECK( true == resonances.isUnresolved() );
      CHECK( 15000. == Approx( resonances.lowerEnergy().value ) );
      CHECK( 100000. == Approx( resonances.upperEnergy().value ) );

      auto compoundsystem = std::get< legacy::unresolved::CompoundSystem >( resonances.compoundSystem() );

      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // content verification
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      // spin groups
      auto spingroups = compoundsystem.spinGroups();
      CHECK( 12 == spingroups.size() );

      auto grid = compoundsystem.grid();
      CHECK( 13 == grid.size() );
      CHECK( 15000. == Approx( grid[0].value ) );
      CHECK( 17000. == Approx( grid[1].value ) );
      CHECK( 20000. == Approx( grid[2].value ) );
      CHECK( 25000. == Approx( grid[3].value ) );
      CHECK( 30000. == Approx( grid[4].value ) );
      CHECK( 35000. == Approx( grid[5].value ) );
      CHECK( 40000. == Approx( grid[6].value ) );
      CHECK( 50000. == Approx( grid[7].value ) );
      CHECK( 60000. == Approx( grid[8].value ) );
      CHECK( 70000. == Approx( grid[9].value ) );
      CHECK( 80000. == Approx( grid[10].value ) );
      CHECK( 90000. == Approx( grid[11].value ) );
      CHECK( 100000. == Approx( grid[12].value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup0 = spingroups[0];

      auto channel0 = spingroup0.incidentChannel();
      CHECK( "n,Na22{0,1/2,5/2+}" == channel0.channelID() );
      CHECK( "n,Na22->n,Na22" == channel0.reactionID().symbol() );

      // resonance table
      auto table0 = spingroup0.resonanceTable();

      CHECK( 10 == table0.numberResonances() );

      auto energies0 = table0.energies();
      CHECK( 15000. == Approx( energies0.front().value ) );
      CHECK( 100000. == Approx( energies0.back().value ) );

      auto resonances0 = table0.resonances();
      CHECK( 15000. == Approx( resonances0.front().energy().value ) );
      CHECK( 100000. == Approx( resonances0.back().energy().value ) );
      CHECK( 55793.70 == Approx( resonances0.front().levelSpacing().value ) );
      CHECK( 55793.70 == Approx( resonances0.back().levelSpacing().value ) );
      CHECK( 7.476360 == Approx( resonances0.front().elastic().value ) );
      CHECK( 7.476360 == Approx( resonances0.back().elastic().value ) );
      CHECK( 1.081650 == Approx( resonances0.front().capture().value ) );
      CHECK( 1.081650 == Approx( resonances0.back().capture().value ) );
      CHECK( 0. == Approx( resonances0.front().fission().value ) );
      CHECK( 0. == Approx( resonances0.back().fission().value ) );
      CHECK( 0. == Approx( resonances0.front().competition().value ) );
      CHECK( 0. == Approx( resonances0.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup1 = spingroups[1];

      auto channel1 = spingroup1.incidentChannel();
      CHECK( "n,Na22{0,1/2,7/2+}" == channel1.channelID() );
      CHECK( "n,Na22->n,Na22" == channel1.reactionID().symbol() );

      // resonance table
      auto table1 = spingroup1.resonanceTable();

      CHECK( 10 == table1.numberResonances() );

      auto energies1 = table1.energies();
      CHECK( 15000. == Approx( energies1.front().value ) );
      CHECK( 100000. == Approx( energies1.back().value ) );

      auto resonances1 = table1.resonances();
      CHECK( 15000. == Approx( resonances1.front().energy().value ) );
      CHECK( 100000. == Approx( resonances1.back().energy().value ) );
      CHECK( 68407.20 == Approx( resonances1.front().levelSpacing().value ) );
      CHECK( 68407.20 == Approx( resonances1.back().levelSpacing().value ) );
      CHECK( 9.166560 == Approx( resonances1.front().elastic().value ) );
      CHECK( 9.166560 == Approx( resonances1.back().elastic().value ) );
      CHECK( 1.081650 == Approx( resonances1.front().capture().value ) );
      CHECK( 1.081650 == Approx( resonances1.back().capture().value ) );
      CHECK( 0. == Approx( resonances1.front().fission().value ) );
      CHECK( 0. == Approx( resonances1.back().fission().value ) );
      CHECK( 0. == Approx( resonances1.front().competition().value ) );
      CHECK( 0. == Approx( resonances1.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 2
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup2 = spingroups[2];

      auto channel2 = spingroup2.incidentChannel();
      CHECK( "n,Na22{1,1/2,3/2-}" == channel2.channelID() );
      CHECK( "n,Na22->n,Na22" == channel2.reactionID().symbol() );

      // resonance table
      auto table2 = spingroup2.resonanceTable();

      CHECK( 10 == table2.numberResonances() );

      auto energies2 = table2.energies();
      CHECK( 15000. == Approx( energies2.front().value ) );
      CHECK( 100000. == Approx( energies2.back().value ) );

      auto resonances2 = table2.resonances();
      CHECK( 15000. == Approx( resonances2.front().energy().value ) );
      CHECK( 100000. == Approx( resonances2.back().energy().value ) );
      CHECK( 58912.60 == Approx( resonances2.front().levelSpacing().value ) );
      CHECK( 58912.60 == Approx( resonances2.back().levelSpacing().value ) );
      CHECK( 16.43660 == Approx( resonances2.front().elastic().value ) );
      CHECK( 16.43660 == Approx( resonances2.back().elastic().value ) );
      CHECK( 5.444310 == Approx( resonances2.front().capture().value ) );
      CHECK( 5.444310 == Approx( resonances2.back().capture().value ) );
      CHECK( 0. == Approx( resonances2.front().fission().value ) );
      CHECK( 0. == Approx( resonances2.back().fission().value ) );
      CHECK( 0. == Approx( resonances2.front().competition().value ) );
      CHECK( 0. == Approx( resonances2.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 3
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup3 = spingroups[3];

      auto channel3 = spingroup3.incidentChannel();
      CHECK( "n,Na22{1,1/2,5/2-}" == channel3.channelID() );
      CHECK( "n,Na22->n,Na22" == channel3.reactionID().symbol() );

      // resonance table
      auto table3 = spingroup3.resonanceTable();

      CHECK( 10 == table3.numberResonances() );

      auto energies3 = table3.energies();
      CHECK( 15000. == Approx( energies3.front().value ) );
      CHECK( 100000. == Approx( energies3.back().value ) );

      auto resonances3 = table3.resonances();
      CHECK( 15000. == Approx( resonances3.front().energy().value ) );
      CHECK( 100000. == Approx( resonances3.back().energy().value ) );
      CHECK( 55793.70 == Approx( resonances3.front().levelSpacing().value ) );
      CHECK( 55793.70 == Approx( resonances3.back().levelSpacing().value ) );
      CHECK( 31.13290 == Approx( resonances3.front().elastic().value ) );
      CHECK( 31.13290 == Approx( resonances3.back().elastic().value ) );
      CHECK( 5.444310 == Approx( resonances3.front().capture().value ) );
      CHECK( 5.444310 == Approx( resonances3.back().capture().value ) );
      CHECK( 0. == Approx( resonances3.front().fission().value ) );
      CHECK( 0. == Approx( resonances3.back().fission().value ) );
      CHECK( 0. == Approx( resonances3.front().competition().value ) );
      CHECK( 0. == Approx( resonances3.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 4
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup4 = spingroups[4];

      auto channel4 = spingroup4.incidentChannel();
      CHECK( "n,Na22{1,1/2,7/2-}" == channel4.channelID() );
      CHECK( "n,Na22->n,Na22" == channel4.reactionID().symbol() );

      // resonance table
      auto table4 = spingroup4.resonanceTable();

      CHECK( 10 == table4.numberResonances() );

      auto energies4 = table4.energies();
      CHECK( 15000. == Approx( energies4.front().value ) );
      CHECK( 100000. == Approx( energies4.back().value ) );

      auto resonances4 = table4.resonances();
      CHECK( 15000. == Approx( resonances4.front().energy().value ) );
      CHECK( 100000. == Approx( resonances4.back().energy().value ) );
      CHECK( 68407.20 == Approx( resonances4.front().levelSpacing().value ) );
      CHECK( 68407.20 == Approx( resonances4.back().levelSpacing().value ) );
      CHECK( 38.17120 == Approx( resonances4.front().elastic().value ) );
      CHECK( 38.17120 == Approx( resonances4.back().elastic().value ) );
      CHECK( 5.444310 == Approx( resonances4.front().capture().value ) );
      CHECK( 5.444310 == Approx( resonances4.back().capture().value ) );
      CHECK( 0. == Approx( resonances4.front().fission().value ) );
      CHECK( 0. == Approx( resonances4.back().fission().value ) );
      CHECK( 0. == Approx( resonances4.front().competition().value ) );
      CHECK( 0. == Approx( resonances4.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 5
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup5 = spingroups[5];

      auto channel5 = spingroup5.incidentChannel();
      CHECK( "n,Na22{1,1/2,9/2-}" == channel5.channelID() );
      CHECK( "n,Na22->n,Na22" == channel5.reactionID().symbol() );

      // resonance table
      auto table5 = spingroup5.resonanceTable();

      CHECK( 10 == table5.numberResonances() );

      auto energies5 = table5.energies();
      CHECK( 15000. == Approx( energies5.front().value ) );
      CHECK( 100000. == Approx( energies5.back().value ) );

      auto resonances5 = table5.resonances();
      CHECK( 15000. == Approx( resonances5.front().energy().value ) );
      CHECK( 100000. == Approx( resonances5.back().energy().value ) );
      CHECK( 102952.0 == Approx( resonances5.front().levelSpacing().value ) );
      CHECK( 102952.0 == Approx( resonances5.back().levelSpacing().value ) );
      CHECK( 28.72350 == Approx( resonances5.front().elastic().value ) );
      CHECK( 28.72350 == Approx( resonances5.back().elastic().value ) );
      CHECK( 5.444310 == Approx( resonances5.front().capture().value ) );
      CHECK( 5.444310 == Approx( resonances5.back().capture().value ) );
      CHECK( 0. == Approx( resonances5.front().fission().value ) );
      CHECK( 0. == Approx( resonances5.back().fission().value ) );
      CHECK( 0. == Approx( resonances5.front().competition().value ) );
      CHECK( 0. == Approx( resonances5.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 6
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup6 = spingroups[6];

      auto channel6 = spingroup6.incidentChannel();
      CHECK( "n,Na22{2,1/2,1/2+}" == channel6.channelID() );
      CHECK( "n,Na22->n,Na22" == channel6.reactionID().symbol() );

      // resonance table
      auto table6 = spingroup6.resonanceTable();

      CHECK( 10 == table6.numberResonances() );

      auto energies6 = table6.energies();
      CHECK( 15000. == Approx( energies6.front().value ) );
      CHECK( 100000. == Approx( energies6.back().value ) );

      auto resonances6 = table6.resonances();
      CHECK( 15000. == Approx( resonances6.front().energy().value ) );
      CHECK( 100000. == Approx( resonances6.back().energy().value ) );
      CHECK( 95446.00 == Approx( resonances6.front().levelSpacing().value ) );
      CHECK( 95446.00 == Approx( resonances6.back().levelSpacing().value ) );
      CHECK( 31.49720 == Approx( resonances6.front().elastic().value ) );
      CHECK( 31.49720 == Approx( resonances6.back().elastic().value ) );
      CHECK( 1.081650 == Approx( resonances6.front().capture().value ) );
      CHECK( 1.081650 == Approx( resonances6.back().capture().value ) );
      CHECK( 0. == Approx( resonances6.front().fission().value ) );
      CHECK( 0. == Approx( resonances6.back().fission().value ) );
      CHECK( 0. == Approx( resonances6.front().competition().value ) );
      CHECK( 0. == Approx( resonances6.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 7
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup7 = spingroups[7];

      auto channel7 = spingroup7.incidentChannel();
      CHECK( "n,Na22{2,1/2,3/2+}" == channel7.channelID() );
      CHECK( "n,Na22->n,Na22" == channel7.reactionID().symbol() );

      // resonance table
      auto table7 = spingroup7.resonanceTable();

      CHECK( 10 == table7.numberResonances() );

      auto energies7 = table0.energies();
      CHECK( 15000. == Approx( energies7.front().value ) );
      CHECK( 100000. == Approx( energies7.back().value ) );

      auto resonances7 = table7.resonances();
      CHECK( 15000. == Approx( resonances7.front().energy().value ) );
      CHECK( 100000. == Approx( resonances7.back().energy().value ) );
      CHECK( 58912.60 == Approx( resonances7.front().levelSpacing().value ) );
      CHECK( 58912.60 == Approx( resonances7.back().levelSpacing().value ) );
      CHECK( 38.88230 == Approx( resonances7.front().elastic().value ) );
      CHECK( 38.88230 == Approx( resonances7.back().elastic().value ) );
      CHECK( 1.081650 == Approx( resonances7.front().capture().value ) );
      CHECK( 1.081650 == Approx( resonances7.back().capture().value ) );
      CHECK( 0. == Approx( resonances7.front().fission().value ) );
      CHECK( 0. == Approx( resonances7.back().fission().value ) );
      CHECK( 0. == Approx( resonances7.front().competition().value ) );
      CHECK( 0. == Approx( resonances7.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 8
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup8 = spingroups[8];

      auto channel8 = spingroup8.incidentChannel();
      CHECK( "n,Na22{2,1/2,5/2+}" == channel8.channelID() );
      CHECK( "n,Na22->n,Na22" == channel8.reactionID().symbol() );

      // resonance table
      auto table8 = spingroup8.resonanceTable();

      CHECK( 10 == table8.numberResonances() );

      auto energies8 = table8.energies();
      CHECK( 15000. == Approx( energies8.front().value ) );
      CHECK( 100000. == Approx( energies8.back().value ) );

      auto resonances8 = table8.resonances();
      CHECK( 15000. == Approx( resonances8.front().energy().value ) );
      CHECK( 100000. == Approx( resonances8.back().energy().value ) );
      CHECK( 55793.70 == Approx( resonances8.front().levelSpacing().value ) );
      CHECK( 55793.70 == Approx( resonances8.back().levelSpacing().value ) );
      CHECK( 36.82380 == Approx( resonances8.front().elastic().value ) );
      CHECK( 36.82380 == Approx( resonances8.back().elastic().value ) );
      CHECK( 1.081650 == Approx( resonances8.front().capture().value ) );
      CHECK( 1.081650 == Approx( resonances8.back().capture().value ) );
      CHECK( 0. == Approx( resonances8.front().fission().value ) );
      CHECK( 0. == Approx( resonances8.back().fission().value ) );
      CHECK( 0. == Approx( resonances8.front().competition().value ) );
      CHECK( 0. == Approx( resonances8.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 9
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup9 = spingroups[9];

      auto channel9 = spingroup9.incidentChannel();
      CHECK( "n,Na22{2,1/2,7/2+}" == channel9.channelID() );
      CHECK( "n,Na22->n,Na22" == channel9.reactionID().symbol() );

      // resonance table
      auto table9 = spingroup9.resonanceTable();

      CHECK( 10 == table9.numberResonances() );

      auto energies9 = table9.energies();
      CHECK( 15000. == Approx( energies9.front().value ) );
      CHECK( 100000. == Approx( energies9.back().value ) );

      auto resonances9 = table9.resonances();
      CHECK( 15000. == Approx( resonances9.front().energy().value ) );
      CHECK( 100000. == Approx( resonances9.back().energy().value ) );
      CHECK( 68407.20 == Approx( resonances9.front().levelSpacing().value ) );
      CHECK( 68407.20 == Approx( resonances9.back().levelSpacing().value ) );
      CHECK( 45.14870 == Approx( resonances9.front().elastic().value ) );
      CHECK( 45.14870 == Approx( resonances9.back().elastic().value ) );
      CHECK( 1.081650 == Approx( resonances9.front().capture().value ) );
      CHECK( 1.081650 == Approx( resonances9.back().capture().value ) );
      CHECK( 0. == Approx( resonances9.front().fission().value ) );
      CHECK( 0. == Approx( resonances9.back().fission().value ) );
      CHECK( 0. == Approx( resonances9.front().competition().value ) );
      CHECK( 0. == Approx( resonances9.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 10
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup10 = spingroups[10];

      auto channel10 = spingroup10.incidentChannel();
      CHECK( "n,Na22{2,1/2,9/2+}" == channel10.channelID() );
      CHECK( "n,Na22->n,Na22" == channel10.reactionID().symbol() );

      // resonance table
      auto table10 = spingroup10.resonanceTable();

      CHECK( 10 == table10.numberResonances() );

      auto energies10 = table10.energies();
      CHECK( 15000. == Approx( energies10.front().value ) );
      CHECK( 100000. == Approx( energies10.back().value ) );

      auto resonances10 = table10.resonances();
      CHECK( 15000. == Approx( resonances10.front().energy().value ) );
      CHECK( 100000. == Approx( resonances10.back().energy().value ) );
      CHECK( 102952.0 == Approx( resonances10.front().levelSpacing().value ) );
      CHECK( 102952.0 == Approx( resonances10.back().levelSpacing().value ) );
      CHECK( 67.94820 == Approx( resonances10.front().elastic().value ) );
      CHECK( 67.94820 == Approx( resonances10.back().elastic().value ) );
      CHECK( 1.081650 == Approx( resonances10.front().capture().value ) );
      CHECK( 1.081650 == Approx( resonances10.back().capture().value ) );
      CHECK( 0. == Approx( resonances10.front().fission().value ) );
      CHECK( 0. == Approx( resonances10.back().fission().value ) );
      CHECK( 0. == Approx( resonances10.front().competition().value ) );
      CHECK( 0. == Approx( resonances10.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 11
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup11 = spingroups[11];

      auto channel11 = spingroup11.incidentChannel();
      CHECK( "n,Na22{2,1/2,11/2+}" == channel11.channelID() );
      CHECK( "n,Na22->n,Na22" == channel11.reactionID().symbol() );

      // resonance table
      auto table11 = spingroup11.resonanceTable();

      CHECK( 10 == table11.numberResonances() );

      auto energies11 = table11.energies();
      CHECK( 15000. == Approx( energies11.front().value ) );
      CHECK( 100000. == Approx( energies11.back().value ) );

      auto resonances11 = table11.resonances();
      CHECK( 15000. == Approx( resonances11.front().energy().value ) );
      CHECK( 100000. == Approx( resonances11.back().energy().value ) );
      CHECK( 185730.0 == Approx( resonances11.front().levelSpacing().value ) );
      CHECK( 185730.0 == Approx( resonances11.back().levelSpacing().value ) );
      CHECK( 61.29090 == Approx( resonances11.front().elastic().value ) );
      CHECK( 61.29090 == Approx( resonances11.back().elastic().value ) );
      CHECK( 1.081650 == Approx( resonances11.front().capture().value ) );
      CHECK( 1.081650 == Approx( resonances11.back().capture().value ) );
      CHECK( 0. == Approx( resonances11.front().fission().value ) );
      CHECK( 0. == Approx( resonances11.back().fission().value ) );
      CHECK( 0. == Approx( resonances11.front().competition().value ) );
      CHECK( 0. == Approx( resonances11.back().competition().value ) );
    } // THEN

    THEN( "cross sections can be reconstructed" ) {

      ReactionID elas( "n,Na22->n,Na22" );
      ReactionID capt( "n,Na22->capture" );
      Map< ReactionID, CrossSection > xs;

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

  GIVEN( "valid ENDF data for Pu239" ) {

    std::string string = Pu239();
    auto begin = string.begin();
    auto end = string.end();
    long lineNumber = 1;

    njoy::ENDFtk::HeadRecord head( begin, end, lineNumber );
    njoy::ENDFtk::section::Type< 2, 151 > endf( head, begin, end, lineNumber, 9437 );
    ResonanceRange endfResonanceRange = endf.isotopes().front().resonanceRanges().front();

    auto resonances = fromENDF( endfResonanceRange, neutronMass, elementaryCharge, ParticleID( "n" ), ParticleID( "Pu239" ) );

    THEN( "the appropriate CompoundSystem is returned" ) {

      CHECK( false == resonances.isResolved() );
      CHECK( true == resonances.isUnresolved() );
      CHECK( 2500. == Approx( resonances.lowerEnergy().value ) );
      CHECK( 30000. == Approx( resonances.upperEnergy().value ) );

      auto compoundsystem = std::get< legacy::unresolved::CompoundSystem >( resonances.compoundSystem() );

      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // content verification
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      // spin groups
      auto spingroups = compoundsystem.spinGroups();
      CHECK( 5 == spingroups.size() );

      auto grid = compoundsystem.grid();
      CHECK( 71 == grid.size() );
      CHECK(  2500. == Approx( grid[0].value ) );
      CHECK(  2550. == Approx( grid[1].value ) );
      CHECK(  2650. == Approx( grid[2].value ) );
      CHECK(  2750. == Approx( grid[3].value ) );
      CHECK(  2850. == Approx( grid[4].value ) );
      CHECK(  2950. == Approx( grid[5].value ) );
      CHECK(  3050. == Approx( grid[6].value ) );
      CHECK(  3150. == Approx( grid[7].value ) );
      CHECK(  3250. == Approx( grid[8].value ) );
      CHECK(  3350. == Approx( grid[9].value ) );
      CHECK(  3450. == Approx( grid[10].value ) );
      CHECK(  3550. == Approx( grid[11].value ) );
      CHECK(  3650. == Approx( grid[12].value ) );
      CHECK(  3750. == Approx( grid[13].value ) );
      CHECK(  3850. == Approx( grid[14].value ) );
      CHECK(  3950. == Approx( grid[15].value ) );
      CHECK(  4125. == Approx( grid[16].value ) );
      CHECK(  4375. == Approx( grid[17].value ) );
      CHECK(  4625. == Approx( grid[18].value ) );
      CHECK(  4875. == Approx( grid[19].value ) );
      CHECK(  5125. == Approx( grid[20].value ) );
      CHECK(  5375. == Approx( grid[21].value ) );
      CHECK(  5625. == Approx( grid[22].value ) );
      CHECK(  5875. == Approx( grid[23].value ) );
      CHECK(  6125. == Approx( grid[24].value ) );
      CHECK(  6375. == Approx( grid[25].value ) );
      CHECK(  6625. == Approx( grid[26].value ) );
      CHECK(  6875. == Approx( grid[27].value ) );
      CHECK(  7125. == Approx( grid[28].value ) );
      CHECK(  7375. == Approx( grid[29].value ) );
      CHECK(  7625. == Approx( grid[30].value ) );
      CHECK(  7875. == Approx( grid[31].value ) );
      CHECK(  8125. == Approx( grid[32].value ) );
      CHECK(  8375. == Approx( grid[33].value ) );
      CHECK(  8625. == Approx( grid[34].value ) );
      CHECK(  8875. == Approx( grid[35].value ) );
      CHECK(  9125. == Approx( grid[36].value ) );
      CHECK(  9375. == Approx( grid[37].value ) );
      CHECK(  9625. == Approx( grid[38].value ) );
      CHECK(  9875. == Approx( grid[39].value ) );
      CHECK( 10250. == Approx( grid[40].value ) );
      CHECK( 10750. == Approx( grid[41].value ) );
      CHECK( 11250. == Approx( grid[42].value ) );
      CHECK( 11750. == Approx( grid[43].value ) );
      CHECK( 12250. == Approx( grid[44].value ) );
      CHECK( 12750. == Approx( grid[45].value ) );
      CHECK( 13250. == Approx( grid[46].value ) );
      CHECK( 13750. == Approx( grid[47].value ) );
      CHECK( 14250. == Approx( grid[48].value ) );
      CHECK( 14750. == Approx( grid[49].value ) );
      CHECK( 15250. == Approx( grid[50].value ) );
      CHECK( 15750. == Approx( grid[51].value ) );
      CHECK( 16250. == Approx( grid[52].value ) );
      CHECK( 16750. == Approx( grid[53].value ) );
      CHECK( 17250. == Approx( grid[54].value ) );
      CHECK( 17750. == Approx( grid[55].value ) );
      CHECK( 18250. == Approx( grid[56].value ) );
      CHECK( 18750. == Approx( grid[57].value ) );
      CHECK( 19250. == Approx( grid[58].value ) );
      CHECK( 19750. == Approx( grid[59].value ) );
      CHECK( 20500. == Approx( grid[60].value ) );
      CHECK( 21500. == Approx( grid[61].value ) );
      CHECK( 22500. == Approx( grid[62].value ) );
      CHECK( 23500. == Approx( grid[63].value ) );
      CHECK( 24500. == Approx( grid[64].value ) );
      CHECK( 25500. == Approx( grid[65].value ) );
      CHECK( 26500. == Approx( grid[66].value ) );
      CHECK( 27500. == Approx( grid[67].value ) );
      CHECK( 28500. == Approx( grid[68].value ) );
      CHECK( 29500. == Approx( grid[69].value ) );
      CHECK( 30000. == Approx( grid[70].value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup0 = spingroups[0];

      auto channel0 = spingroup0.incidentChannel();
      CHECK( "n,Pu239{0,1/2,0+}" == channel0.channelID() );
      CHECK( "n,Pu239->n,Pu239" == channel0.reactionID().symbol() );

      // resonance table
      auto table0 = spingroup0.resonanceTable();

      CHECK( 71 == table0.numberResonances() );

      auto energies0 = table0.energies();
      CHECK( 2500. == Approx( energies0.front().value ) );
      CHECK( 30000. == Approx( energies0.back().value ) );

      auto resonances0 = table0.resonances();
      CHECK( 2500. == Approx( resonances0.front().energy().value ) );
      CHECK( 30000. == Approx( resonances0.back().energy().value ) );
      CHECK( 8.917200 == Approx( resonances0.front().levelSpacing().value ) );
      CHECK( 8.465900 == Approx( resonances0.back().levelSpacing().value ) );
      CHECK( 9.508100e-4 == Approx( resonances0.front().elastic().value ) );
      CHECK( 9.465000e-4 == Approx( resonances0.back().elastic().value ) );
      CHECK( 4.070000e-2 == Approx( resonances0.front().capture().value ) );
      CHECK( 4.070000e-2 == Approx( resonances0.back().capture().value ) );
      CHECK( 2.842000 == Approx( resonances0.front().fission().value ) );
      CHECK( 2.693000 == Approx( resonances0.back().fission().value ) );
      CHECK( 0. == Approx( resonances0.front().competition().value ) );
      CHECK( 0. == Approx( resonances0.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup1 = spingroups[1];

      auto channel1 = spingroup1.incidentChannel();
      CHECK( "n,Pu239{0,1/2,1+}" == channel1.channelID() );
      CHECK( "n,Pu239->n,Pu239" == channel1.reactionID().symbol() );

      // resonance table
      auto table1 = spingroup1.resonanceTable();

      CHECK( 71 == table1.numberResonances() );

      auto energies1 = table1.energies();
      CHECK( 2500. == Approx( energies1.front().value ) );
      CHECK( 30000. == Approx( energies1.back().value ) );

      auto resonances1 = table1.resonances();
      CHECK( 2500. == Approx( resonances1.front().energy().value ) );
      CHECK( 30000. == Approx( resonances1.back().energy().value ) );
      CHECK( 3.044300 == Approx( resonances1.front().levelSpacing().value ) );
      CHECK( 2.890100 == Approx( resonances1.back().levelSpacing().value ) );
      CHECK( 3.246100e-4 == Approx( resonances1.front().elastic().value ) );
      CHECK( 3.231000e-4 == Approx( resonances1.back().elastic().value ) );
      CHECK( 4.030000e-2 == Approx( resonances1.front().capture().value ) );
      CHECK( 4.030000e-2 == Approx( resonances1.back().capture().value ) );
      CHECK( 7.100000e-2 == Approx( resonances1.front().fission().value ) );
      CHECK( 7.647000e-2 == Approx( resonances1.back().fission().value ) );
      CHECK( 0. == Approx( resonances1.front().competition().value ) );
      CHECK( 0. == Approx( resonances1.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 2
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup2 = spingroups[2];

      auto channel2 = spingroup2.incidentChannel();
      CHECK( "n,Pu239{1,1/2,0-}" == channel2.channelID() );
      CHECK( "n,Pu239->n,Pu239" == channel2.reactionID().symbol() );

      // resonance table
      auto table2 = spingroup2.resonanceTable();

      CHECK( 71 == table2.numberResonances() );

      auto energies2 = table2.energies();
      CHECK( 2500. == Approx( energies2.front().value ) );
      CHECK( 30000. == Approx( energies2.back().value ) );

      auto resonances2 = table2.resonances();
      CHECK( 2500. == Approx( resonances2.front().energy().value ) );
      CHECK( 30000. == Approx( resonances2.back().energy().value ) );
      CHECK( 8.917200 == Approx( resonances2.front().levelSpacing().value ) );
      CHECK( 8.465900 == Approx( resonances2.back().levelSpacing().value ) );
      CHECK( 1.573800e-3 == Approx( resonances2.front().elastic().value ) );
      CHECK( 1.411300e-3 == Approx( resonances2.back().elastic().value ) );
      CHECK( 1.150000e-2 == Approx( resonances2.front().capture().value ) );
      CHECK( 1.150000e-2 == Approx( resonances2.back().capture().value ) );
      CHECK( 0. == Approx( resonances2.front().fission().value ) );
      CHECK( 0. == Approx( resonances2.back().fission().value ) );
      CHECK( 0. == Approx( resonances2.front().competition().value ) );
      CHECK( 0. == Approx( resonances2.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 3
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup3 = spingroups[3];

      auto channel3 = spingroup3.incidentChannel();
      CHECK( "n,Pu239{1,1/2,1-}" == channel3.channelID() );
      CHECK( "n,Pu239->n,Pu239" == channel3.reactionID().symbol() );

      // resonance table
      auto table3 = spingroup3.resonanceTable();

      CHECK( 71 == table3.numberResonances() );

      auto energies3 = table3.energies();
      CHECK( 2500. == Approx( energies3.front().value ) );
      CHECK( 30000. == Approx( energies3.back().value ) );

      auto resonances3 = table3.resonances();
      CHECK( 2500. == Approx( resonances3.front().energy().value ) );
      CHECK( 30000. == Approx( resonances3.back().energy().value ) );
      CHECK( 3.044300 == Approx( resonances3.front().levelSpacing().value ) );
      CHECK( 2.890100 == Approx( resonances3.back().levelSpacing().value ) );
      CHECK( 2.686450e-4 == Approx( resonances3.front().elastic().value ) );
      CHECK( 2.408900e-4 == Approx( resonances3.back().elastic().value ) );
      CHECK( 3.035000e-2 == Approx( resonances3.front().capture().value ) );
      CHECK( 3.03000e-2 == Approx( resonances3.back().capture().value ) );
      CHECK( 9.720000e-1 == Approx( resonances3.front().fission().value ) );
      CHECK( 9.180000e-1 == Approx( resonances3.back().fission().value ) );
      CHECK( 0. == Approx( resonances3.front().competition().value ) );
      CHECK( 0. == Approx( resonances3.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 4
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup4 = spingroups[4];

      auto channel4 = spingroup4.incidentChannel();
      CHECK( "n,Pu239{1,1/2,2-}" == channel4.channelID() );
      CHECK( "n,Pu239->n,Pu239" == channel4.reactionID().symbol() );

      // resonance table
      auto table4 = spingroup4.resonanceTable();

      CHECK( 71 == table4.numberResonances() );

      auto energies4 = table4.energies();
      CHECK( 2500. == Approx( energies4.front().value ) );
      CHECK( 30000. == Approx( energies4.back().value ) );

      auto resonances4 = table4.resonances();
      CHECK( 2500. == Approx( resonances4.front().energy().value ) );
      CHECK( 30000. == Approx( resonances4.back().energy().value ) );
      CHECK( 1.915900 == Approx( resonances4.front().levelSpacing().value ) );
      CHECK( 1.818700 == Approx( resonances4.back().levelSpacing().value ) );
      CHECK( 3.381400e-4 == Approx( resonances4.front().elastic().value ) );
      CHECK( 3.032000e-4 == Approx( resonances4.back().elastic().value ) );
      CHECK( 3.335000e-2 == Approx( resonances4.front().capture().value ) );
      CHECK( 3.335000e-2 == Approx( resonances4.back().capture().value ) );
      CHECK( 6.130000e-1 == Approx( resonances4.front().fission().value ) );
      CHECK( 5.770000e-1 == Approx( resonances4.back().fission().value ) );
      CHECK( 0. == Approx( resonances4.front().competition().value ) );
      CHECK( 0. == Approx( resonances4.back().competition().value ) );
    } // THEN

    THEN( "cross sections can be reconstructed" ) {

      ReactionID elas( "n,Pu239->n,Pu239" );
      ReactionID capt( "n,Pu239->capture" );
      ReactionID fiss( "n,Pu239->fission" );
      Map< ReactionID, CrossSection > xs;

      xs = resonances( 2500. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 13.431867799494993 == Approx( xs[ elas ].value ) );
      CHECK( 2.4246294118305469 == Approx( xs[ capt ].value ) );
      CHECK( 4.2347077829718600 == Approx( xs[ fiss ].value ) );

      xs = resonances( 2550. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 13.757049331751272 == Approx( xs[ elas ].value ) );
      CHECK( 2.7535588417225401 == Approx( xs[ capt ].value ) );
      CHECK( 2.7250738402813770 == Approx( xs[ fiss ].value ) );

      xs = resonances( 2650. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 14.785619814411167 == Approx( xs[ elas ].value ) );
      CHECK( 3.4250992815937775 == Approx( xs[ capt ].value ) );
      CHECK( 3.1034768421167902 == Approx( xs[ fiss ].value ) );

      xs = resonances( 2750. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 13.015052442868628 == Approx( xs[ elas ].value ) );
      CHECK( 2.0094674966235324 == Approx( xs[ capt ].value ) );
      CHECK( 4.1691468530984688 == Approx( xs[ fiss ].value ) );

      xs = resonances( 2850. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 13.166067715187607 == Approx( xs[ elas ].value ) );
      CHECK( 2.0766921192491146 == Approx( xs[ capt ].value ) );
      CHECK( 4.1257369135080824 == Approx( xs[ fiss ].value ) );

      xs = resonances( 2950. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 16.630485325284845 == Approx( xs[ elas ].value ) );
      CHECK( 3.7100196889326873 == Approx( xs[ capt ].value ) );
      CHECK( 3.3622150057314930 == Approx( xs[ fiss ].value ) );

      xs = resonances( 3050. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.828566425344780 == Approx( xs[ elas ].value ) );
      CHECK( 1.9979093154182013 == Approx( xs[ capt ].value ) );
      CHECK( 3.0165662797466308 == Approx( xs[ fiss ].value ) );

      xs = resonances( 3150. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 13.383826625287547 == Approx( xs[ elas ].value ) );
      CHECK( 1.9338653188386783 == Approx( xs[ capt ].value ) );
      CHECK( 4.8963477437143519 == Approx( xs[ fiss ].value ) );

      xs = resonances( 3250. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 13.735444943400172 == Approx( xs[ elas ].value ) );
      CHECK( 2.2765782247111428 == Approx( xs[ capt ].value ) );
      CHECK( 3.9544111307709917 == Approx( xs[ fiss ].value ) );

      xs = resonances( 3350. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.828163098228208 == Approx( xs[ elas ].value ) );
      CHECK( 2.1659848701277475 == Approx( xs[ capt ].value ) );
      CHECK( 1.7100938479382908 == Approx( xs[ fiss ].value ) );

      xs = resonances( 3450. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 13.817894408227316 == Approx( xs[ elas ].value ) );
      CHECK( 2.5718330146392030 == Approx( xs[ capt ].value ) );
      CHECK( 2.1984392121425937 == Approx( xs[ fiss ].value ) );

      xs = resonances( 3550. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.646757244603497 == Approx( xs[ elas ].value ) );
      CHECK( 1.8854034302751068 == Approx( xs[ capt ].value ) );
      CHECK( 2.2139612408447187 == Approx( xs[ fiss ].value ) );

      xs = resonances( 3650. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 15.310909041429509 == Approx( xs[ elas ].value ) );
      CHECK( 2.9480222219153105 == Approx( xs[ capt ].value ) );
      CHECK( 2.3944569670663327 == Approx( xs[ fiss ].value ) );

      xs = resonances( 3750. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.640575959272745 == Approx( xs[ elas ].value ) );
      CHECK( 1.6237580373380918 == Approx( xs[ capt ].value ) );
      CHECK( 3.0669295520029860 == Approx( xs[ fiss ].value ) );

      xs = resonances( 3850. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 13.722029715032093 == Approx( xs[ elas ].value ) );
      CHECK( 2.1218992499625551 == Approx( xs[ capt ].value ) );
      CHECK( 3.5565563633157784 == Approx( xs[ fiss ].value ) );

      xs = resonances( 3950. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 14.059124787355131 == Approx( xs[ elas ].value ) );
      CHECK( 2.3965504607649137 == Approx( xs[ capt ].value ) );
      CHECK( 2.9310589935473237 == Approx( xs[ fiss ].value ) );

      xs = resonances( 4125. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 13.528961154554985 == Approx( xs[ elas ].value ) );
      CHECK( 2.2704057957425268 == Approx( xs[ capt ].value ) );
      CHECK( 2.1142575574972580 == Approx( xs[ fiss ].value ) );

      xs = resonances( 4375. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 13.609780185300750 == Approx( xs[ elas ].value ) );
      CHECK( 2.1293843563176686 == Approx( xs[ capt ].value ) );
      CHECK( 2.5096480505990488 == Approx( xs[ fiss ].value ) );

      xs = resonances( 4625. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 13.052069928804340 == Approx( xs[ elas ].value ) );
      CHECK( 1.7147979400601743 == Approx( xs[ capt ].value ) );
      CHECK( 2.7719514755312127 == Approx( xs[ fiss ].value ) );

      xs = resonances( 4875. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 13.719677116547597 == Approx( xs[ elas ].value ) );
      CHECK( 2.1863558896560784 == Approx( xs[ capt ].value ) );
      CHECK( 1.9797770726152923 == Approx( xs[ fiss ].value ) );

      xs = resonances( 5125. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 13.502623926973984 == Approx( xs[ elas ].value ) );
      CHECK( 1.9162545111510034 == Approx( xs[ capt ].value ) );
      CHECK( 2.4065505126237392 == Approx( xs[ fiss ].value ) );

      xs = resonances( 5375. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 13.562053958643903 == Approx( xs[ elas ].value ) );
      CHECK( 1.9532105636617532 == Approx( xs[ capt ].value ) );
      CHECK( 2.1533363795826519 == Approx( xs[ fiss ].value ) );

      xs = resonances( 5625. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 13.422695511184223 == Approx( xs[ elas ].value ) );
      CHECK( 1.8068452088555180 == Approx( xs[ capt ].value ) );
      CHECK( 2.2939498019909874 == Approx( xs[ fiss ].value ) );

      xs = resonances( 5875. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 13.390192639539050 == Approx( xs[ elas ].value ) );
      CHECK( 1.7627114148361509 == Approx( xs[ capt ].value ) );
      CHECK( 2.2332467864787167 == Approx( xs[ fiss ].value ) );

      xs = resonances( 6125. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 13.198804329211262 == Approx( xs[ elas ].value ) );
      CHECK( 1.5843835161262598 == Approx( xs[ capt ].value ) );
      CHECK( 2.4741616941079800 == Approx( xs[ fiss ].value ) );

      xs = resonances( 6375. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 13.424706176346824 == Approx( xs[ elas ].value ) );
      CHECK( 1.7678398184607547 == Approx( xs[ capt ].value ) );
      CHECK( 1.9446042404062127 == Approx( xs[ fiss ].value ) );

      xs = resonances( 6625. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 13.455973570069764 == Approx( xs[ elas ].value ) );
      CHECK( 1.6667307768985120 == Approx( xs[ capt ].value ) );
      CHECK( 1.8811526705591437 == Approx( xs[ fiss ].value ) );

      xs = resonances( 6875. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 13.467089082469917 == Approx( xs[ elas ].value ) );
      CHECK( 1.6814347745297269 == Approx( xs[ capt ].value ) );
      CHECK( 1.7470805116987502 == Approx( xs[ fiss ].value ) );

      xs = resonances( 7125. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 13.488548279720980 == Approx( xs[ elas ].value ) );
      CHECK( 1.5784312329223060 == Approx( xs[ capt ].value ) );
      CHECK( 1.8044139743918477 == Approx( xs[ fiss ].value ) );

      xs = resonances( 7375. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 13.020788627011836 == Approx( xs[ elas ].value ) );
      CHECK( 1.3053812558314954 == Approx( xs[ capt ].value ) );
      CHECK( 2.4205352196525434 == Approx( xs[ fiss ].value ) );

      xs = resonances( 7625. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.934090129198637 == Approx( xs[ elas ].value ) );
      CHECK( 1.2483261737077900 == Approx( xs[ capt ].value ) );
      CHECK( 2.4706296545479458 == Approx( xs[ fiss ].value ) );

      xs = resonances( 7875. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 13.246671145851739 == Approx( xs[ elas ].value ) );
      CHECK( 1.5201992860824949 == Approx( xs[ capt ].value ) );
      CHECK( 1.7997565204799648 == Approx( xs[ fiss ].value ) );

      xs = resonances( 8125. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.958211772802807 == Approx( xs[ elas ].value ) );
      CHECK( 1.3036570594495487 == Approx( xs[ capt ].value ) );
      CHECK( 2.1867221998170163 == Approx( xs[ fiss ].value ) );

      xs = resonances( 8375. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.770586335936564 == Approx( xs[ elas ].value ) );
      CHECK( 1.1804309950377845 == Approx( xs[ capt ].value ) );
      CHECK( 2.4176316769531470 == Approx( xs[ fiss ].value ) );

      xs = resonances( 8625. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.947806523947804 == Approx( xs[ elas ].value ) );
      CHECK( 1.3050046097936769 == Approx( xs[ capt ].value ) );
      CHECK( 2.0381641692884105 == Approx( xs[ fiss ].value ) );

      xs = resonances(  8875. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.760432686695797 == Approx( xs[ elas ].value ) );
      CHECK( 1.1640784418020065 == Approx( xs[ capt ].value ) );
      CHECK( 2.2929613239197519 == Approx( xs[ fiss ].value ) );

      xs = resonances(  9125. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.559170184948195 == Approx( xs[ elas ].value ) );
      CHECK( 1.1655220356423068 == Approx( xs[ capt ].value ) );
      CHECK( 1.8347520770291275 == Approx( xs[ fiss ].value ) );

      xs = resonances(  9375. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.592262067702467 == Approx( xs[ elas ].value ) );
      CHECK( 1.1595571193848930 == Approx( xs[ capt ].value ) );
      CHECK( 1.6663681889910804 == Approx( xs[ fiss ].value ) );

      xs = resonances(  9625. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.635026246096960 == Approx( xs[ elas ].value ) );
      CHECK( 1.0544019009565162 == Approx( xs[ capt ].value ) );
      CHECK( 2.1334064846186003 == Approx( xs[ fiss ].value ) );

      xs = resonances(  9875. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.782067401286186 == Approx( xs[ elas ].value ) );
      CHECK( 1.1648992178296613 == Approx( xs[ capt ].value ) );
      CHECK( 1.7784561774099026 == Approx( xs[ fiss ].value ) );

      xs = resonances( 10250. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.819335366868941 == Approx( xs[ elas ].value ) );
      CHECK( 1.1507794220832457 == Approx( xs[ capt ].value ) );
      CHECK( 1.7951686104929920 == Approx( xs[ fiss ].value ) );

      xs = resonances( 10750. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.611896798348319 == Approx( xs[ elas ].value ) );
      CHECK( 1.0596920794676248 == Approx( xs[ capt ].value ) );
      CHECK( 1.9818319032543621 == Approx( xs[ fiss ].value ) );

      xs = resonances( 11250. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.633103421665119 == Approx( xs[ elas ].value ) );
      CHECK( 1.0711604035852116 == Approx( xs[ capt ].value ) );
      CHECK( 1.8434964727682279 == Approx( xs[ fiss ].value ) );

      xs = resonances( 11750. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.594176815248932 == Approx( xs[ elas ].value ) );
      CHECK( 1.0419874775673770 == Approx( xs[ capt ].value ) );
      CHECK( 1.8133049719717498 == Approx( xs[ fiss ].value ) );

      xs = resonances( 12250. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.486667608800810 == Approx( xs[ elas ].value ) );
      CHECK( 0.96866065027378934 == Approx( xs[ capt ].value ) );
      CHECK(  1.901355820300749 == Approx( xs[ fiss ].value ) );

      xs = resonances( 12750. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.455471187853776 == Approx( xs[ elas ].value ) );
      CHECK( 0.94742843284310452 == Approx( xs[ capt ].value ) );
      CHECK(  1.865477512557371 == Approx( xs[ fiss ].value ) );

      xs = resonances( 13250. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.410143794353386 == Approx( xs[ elas ].value ) );
      CHECK( 0.91755631630504819 == Approx( xs[ capt ].value ) );
      CHECK(  1.858572995399656 == Approx( xs[ fiss ].value ) );

      xs = resonances( 13750. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.449135602633906 == Approx( xs[ elas ].value ) );
      CHECK( 0.94201061234645433 == Approx( xs[ capt ].value ) );
      CHECK(  1.716542591252183 == Approx( xs[ fiss ].value ) );

      xs = resonances( 14250. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.591634237352824 == Approx( xs[ elas ].value ) );
      CHECK( 0.94850656944138445 == Approx( xs[ capt ].value ) );
      CHECK(  1.492806266937279 == Approx( xs[ fiss ].value ) );

      xs = resonances( 14750. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.308022371357769 == Approx( xs[ elas ].value ) );
      CHECK( 0.85450571851375257 == Approx( xs[ capt ].value ) );
      CHECK(  1.798497061940447 == Approx( xs[ fiss ].value ) );

      xs = resonances( 15250. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.210845244801529 == Approx( xs[ elas ].value ) );
      CHECK( 0.79752517089583075 == Approx( xs[ capt ].value ) );
      CHECK(  1.884810485963839 == Approx( xs[ fiss ].value ) );

      xs = resonances( 15750. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.286799583842889 == Approx( xs[ elas ].value ) );
      CHECK( 0.84320542472853222 == Approx( xs[ capt ].value ) );
      CHECK(  1.698221819747725 == Approx( xs[ fiss ].value ) );

      xs = resonances( 16250. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.180623078524269 == Approx( xs[ elas ].value ) );
      CHECK( 0.78202638203132013 == Approx( xs[ capt ].value ) );
      CHECK(  1.802674092760821 == Approx( xs[ fiss ].value ) );

      xs = resonances( 16750. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.251610953569406 == Approx( xs[ elas ].value ) );
      CHECK( 0.82483223820061735 == Approx( xs[ capt ].value ) );
      CHECK(  1.629080574217000 == Approx( xs[ fiss ].value ) );

      xs = resonances( 17250. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.329807047967776 == Approx( xs[ elas ].value ) );
      CHECK( 0.81924952671935924 == Approx( xs[ capt ].value ) );
      CHECK(  1.500099191463520 == Approx( xs[ fiss ].value ) );

      xs = resonances( 17750. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.028767354762378 == Approx( xs[ elas ].value ) );
      CHECK( 0.70155388127063223 == Approx( xs[ capt ].value ) );
      CHECK(  1.863657774693272 == Approx( xs[ fiss ].value ) );

      xs = resonances( 18250. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.090362225314085 == Approx( xs[ elas ].value ) );
      CHECK( 0.73696423297333669 == Approx( xs[ capt ].value ) );
      CHECK(  1.712835190216758 == Approx( xs[ fiss ].value ) );

      xs = resonances( 18750. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.107455154053609 == Approx( xs[ elas ].value ) );
      CHECK( 0.74833832603990780 == Approx( xs[ capt ].value ) );
      CHECK(  1.633891650953657 == Approx( xs[ fiss ].value ) );

      xs = resonances( 19250. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 12.005872458454432 == Approx( xs[ elas ].value ) );
      CHECK( 0.69448668036096795 == Approx( xs[ capt ].value ) );
      CHECK(  1.739938668328656 == Approx( xs[ fiss ].value ) );

      xs = resonances( 19750. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 11.969798878422685 == Approx( xs[ elas ].value ) );
      CHECK( 0.67732933688504371 == Approx( xs[ capt ].value ) );
      CHECK(  1.745465190754702 == Approx( xs[ fiss ].value ) );

      xs = resonances( 20500. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 11.968800519801665 == Approx( xs[ elas ].value ) );
      CHECK( 0.68008741025293007 == Approx( xs[ capt ].value ) );
      CHECK(  1.674538629296025 == Approx( xs[ fiss ].value ) );

      xs = resonances( 21500. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 11.925981702969445 == Approx( xs[ elas ].value ) );
      CHECK( 0.66213720180738678 == Approx( xs[ capt ].value ) );
      CHECK(  1.648644743940811 == Approx( xs[ fiss ].value ) );

      xs = resonances( 22500. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 11.983238513168740 == Approx( xs[ elas ].value ) );
      CHECK( 0.69807106458315982 == Approx( xs[ capt ].value ) );
      CHECK(  1.473934886435610 == Approx( xs[ fiss ].value ) );

      xs = resonances( 23500. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 11.823209511276264 == Approx( xs[ elas ].value ) );
      CHECK( 0.61943892610228057 == Approx( xs[ capt ].value ) );
      CHECK(  1.635054817328119 == Approx( xs[ fiss ].value ) );

      xs = resonances( 24500. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 11.768399134539791 == Approx( xs[ elas ].value ) );
      CHECK( 0.59731025362879375 == Approx( xs[ capt ].value ) );
      CHECK(  1.639098826246428 == Approx( xs[ fiss ].value ) );

      xs = resonances( 25500. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 11.778121013407796 == Approx( xs[ elas ].value ) );
      CHECK( 0.60740269130167657 == Approx( xs[ capt ].value ) );
      CHECK(  1.550000939633398 == Approx( xs[ fiss ].value ) );

      xs = resonances( 26500. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 11.676168727767354 == Approx( xs[ elas ].value ) );
      CHECK( 0.56319077036259035 == Approx( xs[ capt ].value ) );
      CHECK(  1.630580402041155 == Approx( xs[ fiss ].value ) );

      xs = resonances( 27500. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 11.685793766139076 == Approx( xs[ elas ].value ) );
      CHECK( 0.57319551308004868 == Approx( xs[ capt ].value ) );
      CHECK(  1.547042997823886 == Approx( xs[ fiss ].value ) );

      xs = resonances( 28500. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 11.624645091802503 == Approx( xs[ elas ].value ) );
      CHECK( 0.54965485457342067 == Approx( xs[ capt ].value ) );
      CHECK(  1.571659583301897 == Approx( xs[ fiss ].value ) );

      xs = resonances( 29500. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 11.553222436781994 == Approx( xs[ elas ].value ) );
      CHECK( 0.52209263866606681 == Approx( xs[ capt ].value ) );
      CHECK(  1.612989558297598 == Approx( xs[ fiss ].value ) );

      xs = resonances( 30000. * electronVolt );
      CHECK( 3 == xs.size() );
      CHECK( 11.560334256440548 == Approx( xs[ elas ].value ) );
      CHECK( 0.52805263879204856 == Approx( xs[ capt ].value ) );
      CHECK(  1.572007382214403 == Approx( xs[ fiss ].value ) );
    } // THEN
  } // GIVEN

  GIVEN( "valid ENDF data for Er167" ) {

    std::string string = Er167();
    auto begin = string.begin();
    auto end = string.end();
    long lineNumber = 1;

    njoy::ENDFtk::HeadRecord head( begin, end, lineNumber );
    njoy::ENDFtk::section::Type< 2, 151 > endf( head, begin, end, lineNumber, 6840 );
    ResonanceRange endfResonanceRange = endf.isotopes().front().resonanceRanges().front();

    auto resonances = fromENDF( endfResonanceRange, neutronMass, elementaryCharge, ParticleID( "n" ), ParticleID( "Er167" ) );

    THEN( "the appropriate CompoundSystem is returned" ) {

      auto compoundsystem = std::get< legacy::unresolved::CompoundSystem >( resonances.compoundSystem() );

      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // content verification
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      // spin groups
      auto spingroups = compoundsystem.spinGroups();
      CHECK( 6 == spingroups.size() );

      auto grid = compoundsystem.grid();
      CHECK( 11 == grid.size() );
      CHECK( 1750. == Approx( grid[0].value ) );
      CHECK( 2000. == Approx( grid[1].value ) );
      CHECK( 2500. == Approx( grid[2].value ) );
      CHECK( 3000. == Approx( grid[3].value ) );
      CHECK( 3500. == Approx( grid[4].value ) );
      CHECK( 4000. == Approx( grid[5].value ) );
      CHECK( 5000. == Approx( grid[6].value ) );
      CHECK( 6000. == Approx( grid[7].value ) );
      CHECK( 7200. == Approx( grid[8].value ) );
      CHECK( 8500. == Approx( grid[9].value ) );
      CHECK( 10000. == Approx( grid[10].value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup0 = spingroups[0];

      auto channel0 = spingroup0.incidentChannel();
      CHECK( "n,Er167{0,1/2,3+}" == channel0.channelID() );
      CHECK( "n,Er167->n,Er167" == channel0.reactionID().symbol() );

      // resonance table
      auto table0 = spingroup0.resonanceTable();

      CHECK( 2 == table0.numberResonances() );

      auto energies0 = table0.energies();
      CHECK( 1750. == Approx( energies0.front().value ) );
      CHECK( 10000. == Approx( energies0.back().value ) );

      auto resonances0 = table0.resonances();
      CHECK( 1750. == Approx( resonances0.front().energy().value ) );
      CHECK( 10000. == Approx( resonances0.back().energy().value ) );
      CHECK( 9.142900 == Approx( resonances0.front().levelSpacing().value ) );
      CHECK( 9.142900 == Approx( resonances0.back().levelSpacing().value ) );
      CHECK( 1.748200e-3 == Approx( resonances0.front().elastic().value ) );
      CHECK( 1.748200e-3 == Approx( resonances0.back().elastic().value ) );
      CHECK( 1.126380e-1 == Approx( resonances0.front().capture().value ) );
      CHECK( 1.126380e-1 == Approx( resonances0.back().capture().value ) );
      CHECK( 0. == Approx( resonances0.front().fission().value ) );
      CHECK( 0. == Approx( resonances0.back().fission().value ) );
      CHECK( 0. == Approx( resonances0.front().competition().value ) );
      CHECK( 0. == Approx( resonances0.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup1 = spingroups[1];

      auto channel1 = spingroup1.incidentChannel();
      CHECK( "n,Er167{0,1/2,4+}" == channel1.channelID() );
      CHECK( "n,Er167->n,Er167" == channel1.reactionID().symbol() );

      // resonance table
      auto table1 = spingroup1.resonanceTable();

      CHECK( 2 == table1.numberResonances() );

      auto energies1 = table1.energies();
      CHECK( 1750. == Approx( energies1.front().value ) );
      CHECK( 10000. == Approx( energies1.back().value ) );

      auto resonances1 = table1.resonances();
      CHECK( 1750. == Approx( resonances1.front().energy().value ) );
      CHECK( 10000. == Approx( resonances1.back().energy().value ) );
      CHECK( 7.111100 == Approx( resonances1.front().levelSpacing().value ) );
      CHECK( 7.111100 == Approx( resonances1.back().levelSpacing().value ) );
      CHECK( 1.321300e-3 == Approx( resonances1.front().elastic().value ) );
      CHECK( 1.321300e-3 == Approx( resonances1.back().elastic().value ) );
      CHECK( 1.126380e-1 == Approx( resonances1.front().capture().value ) );
      CHECK( 1.126380e-1 == Approx( resonances1.back().capture().value ) );
      CHECK( 0. == Approx( resonances1.front().fission().value ) );
      CHECK( 0. == Approx( resonances1.back().fission().value ) );
      CHECK( 0. == Approx( resonances1.front().competition().value ) );
      CHECK( 0. == Approx( resonances1.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 2
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup2 = spingroups[2];

      auto channel2 = spingroup2.incidentChannel();
      CHECK( "n,Er167{1,1/2,2-}" == channel2.channelID() );
      CHECK( "n,Er167->n,Er167" == channel2.reactionID().symbol() );

      // resonance table
      auto table2 = spingroup2.resonanceTable();

      CHECK( 2 == table2.numberResonances() );

      auto energies2 = table2.energies();
      CHECK( 1750. == Approx( energies2.front().value ) );
      CHECK( 10000. == Approx( energies2.back().value ) );

      auto resonances2 = table2.resonances();
      CHECK( 1750. == Approx( resonances2.front().energy().value ) );
      CHECK( 10000. == Approx( resonances2.back().energy().value ) );
      CHECK( 12.8 == Approx( resonances2.front().levelSpacing().value ) );
      CHECK( 12.8 == Approx( resonances2.back().levelSpacing().value ) );
      CHECK( 1.236000e-3 == Approx( resonances2.front().elastic().value ) );
      CHECK( 1.236000e-3 == Approx( resonances2.back().elastic().value ) );
      CHECK( 1.126380e-1 == Approx( resonances2.front().capture().value ) );
      CHECK( 1.126380e-1 == Approx( resonances2.back().capture().value ) );
      CHECK( 0. == Approx( resonances2.front().fission().value ) );
      CHECK( 0. == Approx( resonances2.back().fission().value ) );
      CHECK( 0. == Approx( resonances2.front().competition().value ) );
      CHECK( 0. == Approx( resonances2.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 3
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup3 = spingroups[3];

      auto channel3 = spingroup3.incidentChannel();
      CHECK( "n,Er167{1,1/2,3-}" == channel3.channelID() );
      CHECK( "n,Er167->n,Er167" == channel3.reactionID().symbol() );

      // resonance table
      auto table3 = spingroup3.resonanceTable();

      CHECK( 2 == table3.numberResonances() );

      auto energies3 = table3.energies();
      CHECK( 1750. == Approx( energies3.front().value ) );
      CHECK( 10000. == Approx( energies3.back().value ) );

      auto resonances3 = table3.resonances();
      CHECK( 1750. == Approx( resonances3.front().energy().value ) );
      CHECK( 10000. == Approx( resonances3.back().energy().value ) );
      CHECK( 9.142900 == Approx( resonances3.front().levelSpacing().value ) );
      CHECK( 9.142900 == Approx( resonances3.back().levelSpacing().value ) );
      CHECK( 1.894300e-3 == Approx( resonances3.front().elastic().value ) );
      CHECK( 1.894300e-3 == Approx( resonances3.back().elastic().value ) );
      CHECK( 1.126380e-1 == Approx( resonances3.front().capture().value ) );
      CHECK( 1.126380e-1 == Approx( resonances3.back().capture().value ) );
      CHECK( 0. == Approx( resonances3.front().fission().value ) );
      CHECK( 0. == Approx( resonances3.back().fission().value ) );
      CHECK( 0. == Approx( resonances3.front().competition().value ) );
      CHECK( 0. == Approx( resonances3.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 4
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup4 = spingroups[4];

      auto channel4 = spingroup4.incidentChannel();
      CHECK( "n,Er167{1,1/2,4-}" == channel4.channelID() );
      CHECK( "n,Er167->n,Er167" == channel4.reactionID().symbol() );

      // resonance table
      auto table4 = spingroup4.resonanceTable();

      CHECK( 2 == table4.numberResonances() );

      auto energies4 = table4.energies();
      CHECK( 1750. == Approx( energies4.front().value ) );
      CHECK( 10000. == Approx( energies4.back().value ) );

      auto resonances4 = table4.resonances();
      CHECK( 1750. == Approx( resonances4.front().energy().value ) );
      CHECK( 10000. == Approx( resonances4.back().energy().value ) );
      CHECK( 7.1111 == Approx( resonances4.front().levelSpacing().value ) );
      CHECK( 7.1111 == Approx( resonances4.back().levelSpacing().value ) );
      CHECK( 1.406600e-3 == Approx( resonances4.front().elastic().value ) );
      CHECK( 1.406600e-3 == Approx( resonances4.back().elastic().value ) );
      CHECK( 1.126380e-1 == Approx( resonances4.front().capture().value ) );
      CHECK( 1.126380e-1 == Approx( resonances4.back().capture().value ) );
      CHECK( 0. == Approx( resonances4.front().fission().value ) );
      CHECK( 0. == Approx( resonances4.back().fission().value ) );
      CHECK( 0. == Approx( resonances4.front().competition().value ) );
      CHECK( 0. == Approx( resonances4.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 5
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup5 = spingroups[5];

      auto channel5 = spingroup5.incidentChannel();
      CHECK( "n,Er167{1,1/2,5-}" == channel5.channelID() );
      CHECK( "n,Er167->n,Er167" == channel5.reactionID().symbol() );

      // resonance table
      auto table5 = spingroup5.resonanceTable();

      CHECK( 2 == table5.numberResonances() );

      auto energies5 = table5.energies();
      CHECK( 1750. == Approx( energies5.front().value ) );
      CHECK( 10000. == Approx( energies5.back().value ) );

      auto resonances5 = table5.resonances();
      CHECK( 1750. == Approx( resonances5.front().energy().value ) );
      CHECK( 10000. == Approx( resonances5.back().energy().value ) );
      CHECK( 5.8182 == Approx( resonances5.front().levelSpacing().value ) );
      CHECK( 5.8182 == Approx( resonances5.back().levelSpacing().value ) );
      CHECK( 6.681800e-4 == Approx( resonances5.front().elastic().value ) );
      CHECK( 6.681800e-4 == Approx( resonances5.back().elastic().value ) );
      CHECK( 1.126380e-1 == Approx( resonances5.front().capture().value ) );
      CHECK( 1.126380e-1 == Approx( resonances5.back().capture().value ) );
      CHECK( 0. == Approx( resonances5.front().fission().value ) );
      CHECK( 0. == Approx( resonances5.back().fission().value ) );
      CHECK( 0. == Approx( resonances5.front().competition().value ) );
      CHECK( 0. == Approx( resonances5.back().competition().value ) );
    } // THEN

    THEN( "cross sections can be reconstructed" ) {

      ReactionID elas( "n,Er167->n,Er167" );
      ReactionID capt( "n,Er167->capture" );
      Map< ReactionID, CrossSection > xs;

      xs = resonances( 1750. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 18.241628628010506 == Approx( xs[ elas ].value ) );
      CHECK( 8.8470221152128197 == Approx( xs[ capt ].value ) );

      xs = resonances(  2000. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 17.820760274730802 == Approx( xs[ elas ].value ) );
      CHECK( 8.0651842000962901 == Approx( xs[ capt ].value ) );

      xs = resonances(  2500. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 17.135195181202420 == Approx( xs[ elas ].value ) );
      CHECK( 6.9113059738704008 == Approx( xs[ capt ].value ) );

      xs = resonances(  3000. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 16.592348035233275 == Approx( xs[ elas ].value ) );
      CHECK( 6.0957654963405616 == Approx( xs[ capt ].value ) );

      xs = resonances(  3500. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 16.146013683922963 == Approx( xs[ elas ].value ) );
      CHECK( 5.4858210475023412 == Approx( xs[ capt ].value ) );

      xs = resonances(  4000. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 15.768997902900036 == Approx( xs[ elas ].value ) );
      CHECK( 5.0108853038756882 == Approx( xs[ capt ].value ) );

      xs = resonances(  5000. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 15.159379589050337 == Approx( xs[ elas ].value ) );
      CHECK( 4.3164971574593238 == Approx( xs[ capt ].value ) );

      xs = resonances(  6000. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 14.680652351726371 == Approx( xs[ elas ].value ) );
      CHECK( 3.8310594535618745 == Approx( xs[ capt ].value ) );

      xs = resonances(  7200. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 14.219700422962447 == Approx( xs[ elas ].value ) );
      CHECK( 3.4100036907692624 == Approx( xs[ capt ].value ) );

      xs = resonances(  8500. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 13.815771022819700 == Approx( xs[ elas ].value ) );
      CHECK( 3.0759470502728590 == Approx( xs[ capt ].value ) );

      xs = resonances( 10000. * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 13.435021117381460 == Approx( xs[ elas ].value ) );
      CHECK( 2.7888116700218495 == Approx( xs[ capt ].value ) );
    } // THEN
  } // GIVEN

  GIVEN( "valid ENDF data for Au197" ) {

    std::string string = Au197();
    auto begin = string.begin();
    auto end = string.end();
    long lineNumber = 1;

    njoy::ENDFtk::HeadRecord head( begin, end, lineNumber );
    njoy::ENDFtk::section::Type< 2, 151 > endf( head, begin, end, lineNumber, 7925 );
    ResonanceRange endfResonanceRange = endf.isotopes().front().resonanceRanges().front();

    auto resonances = fromENDF( endfResonanceRange, neutronMass, elementaryCharge, ParticleID( "n" ), ParticleID( "Au197" ) );

    THEN( "the appropriate CompoundSystem is returned" ) {

      auto compoundsystem = std::get< legacy::unresolved::CompoundSystem >( resonances.compoundSystem() );

      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // content verification
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      // spin groups
      auto spingroups = compoundsystem.spinGroups();
      CHECK( 11 == spingroups.size() );

      auto grid = compoundsystem.grid();
      CHECK( 100 == grid.size() );
      CHECK( 2.000000e+3 == Approx( grid[0].value ) );
      CHECK( 2.250000e+3 == Approx( grid[1].value ) );
      CHECK( 2.500000e+3 == Approx( grid[2].value ) );
      CHECK( 2.750000e+3 == Approx( grid[3].value ) );
      CHECK( 3.000000e+3 == Approx( grid[4].value ) );
      CHECK( 3.250000e+3 == Approx( grid[5].value ) );
      CHECK( 3.500000e+3 == Approx( grid[6].value ) );
      CHECK( 3.750000e+3 == Approx( grid[7].value ) );
      CHECK( 4.000000e+3 == Approx( grid[8].value ) );
      CHECK( 4.250000e+3 == Approx( grid[9].value ) );
      CHECK( 4.500000e+3 == Approx( grid[10].value ) );
      CHECK( 4.750000e+3 == Approx( grid[11].value ) );
      CHECK( 5.000000e+3 == Approx( grid[12].value ) );
      CHECK( 5.250000e+3 == Approx( grid[13].value ) );
      CHECK( 5.500000e+3 == Approx( grid[14].value ) );
      CHECK( 6.000000e+3 == Approx( grid[15].value ) );
      CHECK( 6.500000e+3 == Approx( grid[16].value ) );
      CHECK( 7.000000e+3 == Approx( grid[17].value ) );
      CHECK( 7.500000e+3 == Approx( grid[18].value ) );
      CHECK( 8.000000e+3 == Approx( grid[19].value ) );
      CHECK( 8.500000e+3 == Approx( grid[20].value ) );
      CHECK( 9.000000e+3 == Approx( grid[21].value ) );
      CHECK( 9.500000e+3 == Approx( grid[22].value ) );
      CHECK( 1.000000e+4 == Approx( grid[23].value ) );
      CHECK( 1.100000e+4 == Approx( grid[24].value ) );
      CHECK( 1.200000e+4 == Approx( grid[25].value ) );
      CHECK( 1.300000e+4 == Approx( grid[26].value ) );
      CHECK( 1.400000e+4 == Approx( grid[27].value ) );
      CHECK( 1.500000e+4 == Approx( grid[28].value ) );
      CHECK( 1.600000e+4 == Approx( grid[29].value ) );
      CHECK( 1.700000e+4 == Approx( grid[30].value ) );
      CHECK( 1.800000e+4 == Approx( grid[31].value ) );
      CHECK( 1.900000e+4 == Approx( grid[32].value ) );
      CHECK( 2.000000e+4 == Approx( grid[33].value ) );
      CHECK( 2.100000e+4 == Approx( grid[34].value ) );
      CHECK( 2.200000e+4 == Approx( grid[35].value ) );
      CHECK( 2.300000e+4 == Approx( grid[36].value ) );
      CHECK( 2.400000e+4 == Approx( grid[37].value ) );
      CHECK( 2.500000e+4 == Approx( grid[38].value ) );
      CHECK( 2.600000e+4 == Approx( grid[39].value ) );
      CHECK( 2.700000e+4 == Approx( grid[40].value ) );
      CHECK( 2.800000e+4 == Approx( grid[41].value ) );
      CHECK( 2.900000e+4 == Approx( grid[42].value ) );
      CHECK( 3.000000e+4 == Approx( grid[43].value ) );
      CHECK( 3.100000e+4 == Approx( grid[44].value ) );
      CHECK( 3.200000e+4 == Approx( grid[45].value ) );
      CHECK( 3.300000e+4 == Approx( grid[46].value ) );
      CHECK( 3.400000e+4 == Approx( grid[47].value ) );
      CHECK( 3.500000e+4 == Approx( grid[48].value ) );
      CHECK( 3.600000e+4 == Approx( grid[49].value ) );
      CHECK( 3.700000e+4 == Approx( grid[50].value ) );
      CHECK( 3.800000e+4 == Approx( grid[51].value ) );
      CHECK( 3.900000e+4 == Approx( grid[52].value ) );
      CHECK( 4.000000e+4 == Approx( grid[53].value ) );
      CHECK( 4.100000e+4 == Approx( grid[54].value ) );
      CHECK( 4.200000e+4 == Approx( grid[55].value ) );
      CHECK( 4.300000e+4 == Approx( grid[56].value ) );
      CHECK( 4.400000e+4 == Approx( grid[57].value ) );
      CHECK( 4.500000e+4 == Approx( grid[58].value ) );
      CHECK( 4.600000e+4 == Approx( grid[59].value ) );
      CHECK( 4.700000e+4 == Approx( grid[60].value ) );
      CHECK( 4.800000e+4 == Approx( grid[61].value ) );
      CHECK( 5.000000e+4 == Approx( grid[62].value ) );
      CHECK( 5.250000e+4 == Approx( grid[63].value ) );
      CHECK( 5.500000e+4 == Approx( grid[64].value ) );
      CHECK( 5.750000e+4 == Approx( grid[65].value ) );
      CHECK( 6.000000e+4 == Approx( grid[66].value ) );
      CHECK( 6.250000e+4 == Approx( grid[67].value ) );
      CHECK( 6.500000e+4 == Approx( grid[68].value ) );
      CHECK( 6.750000e+4 == Approx( grid[69].value ) );
      CHECK( 7.000000e+4 == Approx( grid[70].value ) );
      CHECK( 7.200000e+4 == Approx( grid[71].value ) );
      CHECK( 7.400000e+4 == Approx( grid[72].value ) );
      CHECK( 7.500000e+4 == Approx( grid[73].value ) );
      CHECK( 7.775000e+4 == Approx( grid[74].value ) );
      CHECK( 7.779600e+4 == Approx( grid[75].value ) );
      CHECK( 7.781000e+4 == Approx( grid[76].value ) );
      CHECK( 7.783000e+4 == Approx( grid[77].value ) );
      CHECK( 7.786000e+4 == Approx( grid[78].value ) );
      CHECK( 7.789000e+4 == Approx( grid[79].value ) );
      CHECK( 7.793000e+4 == Approx( grid[80].value ) );
      CHECK( 7.798000e+4 == Approx( grid[81].value ) );
      CHECK( 7.807000e+4 == Approx( grid[82].value ) );
      CHECK( 7.820000e+4 == Approx( grid[83].value ) );
      CHECK( 7.835000e+4 == Approx( grid[84].value ) );
      CHECK( 7.862000e+4 == Approx( grid[85].value ) );
      CHECK( 7.890000e+4 == Approx( grid[86].value ) );
      CHECK( 7.945000e+4 == Approx( grid[87].value ) );
      CHECK( 8.000000e+4 == Approx( grid[88].value ) );
      CHECK( 8.100000e+4 == Approx( grid[89].value ) );
      CHECK( 8.200000e+4 == Approx( grid[90].value ) );
      CHECK( 8.300000e+4 == Approx( grid[91].value ) );
      CHECK( 8.500000e+4 == Approx( grid[92].value ) );
      CHECK( 8.750000e+4 == Approx( grid[93].value ) );
      CHECK( 9.000000e+4 == Approx( grid[94].value ) );
      CHECK( 9.200000e+4 == Approx( grid[95].value ) );
      CHECK( 9.350000e+4 == Approx( grid[96].value ) );
      CHECK( 9.500000e+4 == Approx( grid[97].value ) );
      CHECK( 9.750000e+4 == Approx( grid[98].value ) );
      CHECK( 1.000000e+5 == Approx( grid[99].value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup0 = spingroups[0];

      auto channel0 = spingroup0.incidentChannel();
      CHECK( "n,Au197{0,1/2,1+}" == channel0.channelID() );
      CHECK( "n,Au197->n,Au197" == channel0.reactionID().symbol() );

      // resonance table
      auto table0 = spingroup0.resonanceTable();

      CHECK( 100 == table0.numberResonances() );

      auto energies0 = table0.energies();
      CHECK( 2000. == Approx( energies0.front().value ) );
      CHECK( 100000. == Approx( energies0.back().value ) );

      auto resonances0 = table0.resonances();
      CHECK( 2000. == Approx( resonances0.front().energy().value ) );
      CHECK( 100000. == Approx( resonances0.back().energy().value ) );
      CHECK( 38.48710 == Approx( resonances0.front().levelSpacing().value ) );
      CHECK( 33.54027 == Approx( resonances0.back().levelSpacing().value ) );
      CHECK( 6.450977e-3 == Approx( resonances0.front().elastic().value ) );
      CHECK( 5.129039e-3  == Approx( resonances0.back().elastic().value ) );
      CHECK( 1.292688e-1 == Approx( resonances0.front().capture().value ) );
      CHECK( 1.337047e-1 == Approx( resonances0.back().capture().value ) );
      CHECK( 0. == Approx( resonances0.front().fission().value ) );
      CHECK( 0. == Approx( resonances0.back().fission().value ) );
      CHECK( 0. == Approx( resonances0.front().competition().value ) );
      CHECK( .8007053 == Approx( resonances0.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 1
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup1 = spingroups[1];

      auto channel1 = spingroup1.incidentChannel();
      CHECK( "n,Au197{0,1/2,2+}" == channel1.channelID() );
      CHECK( "n,Au197->n,Au197" == channel1.reactionID().symbol() );

      // resonance table
      auto table1 = spingroup1.resonanceTable();

      CHECK( 100 == table1.numberResonances() );

      auto energies1 = table1.energies();
      CHECK( 2000. == Approx( energies1.front().value ) );
      CHECK( 100000. == Approx( energies1.back().value ) );

      auto resonances1 = table1.resonances();
      CHECK( 2000. == Approx( resonances1.front().energy().value ) );
      CHECK( 100000. == Approx( resonances1.back().energy().value ) );
      CHECK( 23.93741 == Approx( resonances1.front().levelSpacing().value ) );
      CHECK( 20.85515 == Approx( resonances1.back().levelSpacing().value ) );
      CHECK( 4.012245e-3 == Approx( resonances1.front().elastic().value ) );
      CHECK( 3.189207e-3 == Approx( resonances1.back().elastic().value ) );
      CHECK( 1.292688e-1 == Approx( resonances1.front().capture().value ) );
      CHECK( 1.337047e-1 == Approx( resonances1.back().capture().value ) );
      CHECK( 0. == Approx( resonances1.front().fission().value ) );
      CHECK( 0. == Approx( resonances1.back().fission().value ) );
      CHECK( 0. == Approx( resonances1.front().competition().value ) );
      CHECK( 1.065243e-3 == Approx( resonances1.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 2
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup2 = spingroups[2];

      auto channel2 = spingroup2.incidentChannel();
      CHECK( "n,Au197{1,1/2,0-}" == channel2.channelID() );
      CHECK( "n,Au197->n,Au197" == channel2.reactionID().symbol() );

      // resonance table
      auto table2 = spingroup2.resonanceTable();

      CHECK( 100 == table2.numberResonances() );

      auto energies2 = table2.energies();
      CHECK( 2000. == Approx( energies2.front().value ) );
      CHECK( 100000. == Approx( energies2.back().value ) );

      auto resonances2 = table2.resonances();
      CHECK( 2000. == Approx( resonances2.front().energy().value ) );
      CHECK( 100000. == Approx( resonances2.back().energy().value ) );
      CHECK( 113.4047 == Approx( resonances2.front().levelSpacing().value ) );
      CHECK( 98.84168 == Approx( resonances2.back().levelSpacing().value ) );
      CHECK( 6.028556e-3 == Approx( resonances2.front().elastic().value ) );
      CHECK( 5.721812e-3 == Approx( resonances2.back().elastic().value ) );
      CHECK( 4.205961e-2 == Approx( resonances2.front().capture().value ) );
      CHECK( 4.350289e-2 == Approx( resonances2.back().capture().value ) );
      CHECK( 0. == Approx( resonances2.front().fission().value ) );
      CHECK( 0. == Approx( resonances2.back().fission().value ) );
      CHECK( 0. == Approx( resonances2.front().competition().value ) );
      CHECK( 5.033883e-2 == Approx( resonances2.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 3
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup3 = spingroups[3];

      auto channel3 = spingroup3.incidentChannel();
      CHECK( "n,Au197{1,1/2,1-}" == channel3.channelID() );
      CHECK( "n,Au197->n,Au197" == channel3.reactionID().symbol() );

      // resonance table
      auto table3 = spingroup3.resonanceTable();

      CHECK( 100 == table3.numberResonances() );

      auto energies3 = table3.energies();
      CHECK( 2000. == Approx( energies3.front().value ) );
      CHECK( 100000. == Approx( energies3.back().value ) );

      auto resonances3 = table3.resonances();
      CHECK( 2000. == Approx( resonances3.front().energy().value ) );
      CHECK( 100000. == Approx( resonances3.back().energy().value ) );
      CHECK( 38.48710 == Approx( resonances3.front().levelSpacing().value ) );
      CHECK( 33.54027 == Approx( resonances3.back().levelSpacing().value ) );
      CHECK( 2.045961e-3 == Approx( resonances3.front().elastic().value ) );
      CHECK( 1.941601e-3 == Approx( resonances3.back().elastic().value ) );
      CHECK( 4.205961e-2 == Approx( resonances3.front().capture().value ) );
      CHECK( 4.350289e-2 == Approx( resonances3.back().capture().value ) );
      CHECK( 0. == Approx( resonances3.front().fission().value ) );
      CHECK( 0. == Approx( resonances3.back().fission().value ) );
      CHECK( 0. == Approx( resonances3.front().competition().value ) );
      CHECK( 3.416328e-2 == Approx( resonances3.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 4
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup4 = spingroups[4];

      auto channel4 = spingroup4.incidentChannel();
      CHECK( "n,Au197{1,1/2,2-}" == channel4.channelID() );
      CHECK( "n,Au197->n,Au197" == channel4.reactionID().symbol() );

      // resonance table
      auto table4 = spingroup4.resonanceTable();

      CHECK( 100 == table4.numberResonances() );

      auto energies4 = table4.energies();
      CHECK( 2000. == Approx( energies4.front().value ) );
      CHECK( 100000. == Approx( energies4.back().value ) );

      auto resonances4 = table4.resonances();
      CHECK( 2000. == Approx( resonances4.front().energy().value ) );
      CHECK( 100000. == Approx( resonances4.back().energy().value ) );
      CHECK( 23.93741  == Approx( resonances4.front().levelSpacing().value ) );
      CHECK( 20.85515 == Approx( resonances4.back().levelSpacing().value ) );
      CHECK( 1.272505e-3 == Approx( resonances4.front().elastic().value ) );
      CHECK( 1.207276e-3 == Approx( resonances4.back().elastic().value ) );
      CHECK( 4.205961e-2 == Approx( resonances4.front().capture().value ) );
      CHECK( 4.350289e-2 == Approx( resonances4.back().capture().value ) );
      CHECK( 0. == Approx( resonances4.front().fission().value ) );
      CHECK( 0. == Approx( resonances4.back().fission().value ) );
      CHECK( 0. == Approx( resonances4.front().competition().value ) );
      CHECK( 1.062126e-2 == Approx( resonances4.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 5
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup5 = spingroups[5];

      auto channel5 = spingroup5.incidentChannel();
      CHECK( "n,Au197{1,1/2,3-}" == channel5.channelID() );
      CHECK( "n,Au197->n,Au197" == channel5.reactionID().symbol() );

      // resonance table
      auto table5 = spingroup5.resonanceTable();

      CHECK( 100 == table5.numberResonances() );

      auto energies5 = table5.energies();
      CHECK( 2000. == Approx( energies5.front().value ) );
      CHECK( 100000. == Approx( energies5.back().value ) );

      auto resonances5 = table5.resonances();
      CHECK( 2000. == Approx( resonances5.front().energy().value ) );
      CHECK( 100000. == Approx( resonances5.back().energy().value ) );
      CHECK( 18.04535 == Approx( resonances5.front().levelSpacing().value ) );
      CHECK( 15.71550 == Approx( resonances5.back().levelSpacing().value ) );
      CHECK( 9.592847e-4 == Approx( resonances5.front().elastic().value ) );
      CHECK( 9.097494e-4 == Approx( resonances5.back().elastic().value ) );
      CHECK( 4.205961e-2 == Approx( resonances5.front().capture().value ) );
      CHECK( 4.350289e-2 == Approx( resonances5.back().capture().value ) );
      CHECK( 0. == Approx( resonances5.front().fission().value ) );
      CHECK( 0. == Approx( resonances5.back().fission().value ) );
      CHECK( 0. == Approx( resonances5.front().competition().value ) );
      CHECK( 0. == Approx( resonances5.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 6
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup6 = spingroups[6];

      auto channel6 = spingroup6.incidentChannel();
      CHECK( "n,Au197{2,1/2,0+}" == channel6.channelID() );
      CHECK( "n,Au197->n,Au197" == channel6.reactionID().symbol() );

      // resonance table
      auto table6 = spingroup6.resonanceTable();

      CHECK( 100 == table6.numberResonances() );

      auto energies6 = table6.energies();
      CHECK( 2000. == Approx( energies6.front().value ) );
      CHECK( 100000. == Approx( energies6.back().value ) );

      auto resonances6 = table6.resonances();
      CHECK( 2000. == Approx( resonances6.front().energy().value ) );
      CHECK( 100000. == Approx( resonances6.back().energy().value ) );
      CHECK( 113.4047 == Approx( resonances6.front().levelSpacing().value ) );
      CHECK( 98.84168 == Approx( resonances6.back().levelSpacing().value ) );
      CHECK( 3.981167e-2 == Approx( resonances6.front().elastic().value ) );
      CHECK( 3.380303e-2 == Approx( resonances6.back().elastic().value ) );
      CHECK( 1.292688e-1 == Approx( resonances6.front().capture().value ) );
      CHECK( 1.337047e-1 == Approx( resonances6.back().capture().value ) );
      CHECK( 0. == Approx( resonances6.front().fission().value ) );
      CHECK( 0. == Approx( resonances6.back().fission().value ) );
      CHECK( 0. == Approx( resonances6.front().competition().value ) );
      CHECK( 2.359643 == Approx( resonances6.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 7
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup7 = spingroups[7];

      auto channel7 = spingroup7.incidentChannel();
      CHECK( "n,Au197{2,1/2,1+}" == channel7.channelID() );
      CHECK( "n,Au197->n,Au197" == channel7.reactionID().symbol() );

      // resonance table
      auto table7 = spingroup7.resonanceTable();

      CHECK( 100 == table7.numberResonances() );

      auto energies7 = table7.energies();
      CHECK( 2000. == Approx( energies7.front().value ) );
      CHECK( 100000. == Approx( energies7.back().value ) );

      auto resonances7 = table7.resonances();
      CHECK( 2000. == Approx( resonances7.front().energy().value ) );
      CHECK( 100000. == Approx( resonances7.back().energy().value ) );
      CHECK( 38.48710 == Approx( resonances7.front().levelSpacing().value ) );
      CHECK( 33.54027 == Approx( resonances7.back().levelSpacing().value ) );
      CHECK( 1.351122e-2 == Approx( resonances7.front().elastic().value ) );
      CHECK( 1.147049e-2 == Approx( resonances7.back().elastic().value ) );
      CHECK( 1e-9 == Approx( resonances7.front().capture().value ) );
      CHECK( 1e-9 == Approx( resonances7.back().capture().value ) );
      CHECK( 0. == Approx( resonances7.front().fission().value ) );
      CHECK( 0. == Approx( resonances7.back().fission().value ) );
      CHECK( 0. == Approx( resonances7.front().competition().value ) );
      CHECK( 1e-9 == Approx( resonances7.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 8
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup8 = spingroups[8];

      auto channel8 = spingroup8.incidentChannel();
      CHECK( "n,Au197{2,1/2,2+}" == channel8.channelID() );
      CHECK( "n,Au197->n,Au197" == channel8.reactionID().symbol() );

      // resonance table
      auto table8 = spingroup8.resonanceTable();

      CHECK( 100 == table8.numberResonances() );

      auto energies8 = table8.energies();
      CHECK( 2000. == Approx( energies8.front().value ) );
      CHECK( 100000. == Approx( energies8.back().value ) );

      auto resonances8 = table8.resonances();
      CHECK( 2000. == Approx( resonances8.front().energy().value ) );
      CHECK( 100000. == Approx( resonances8.back().energy().value ) );
      CHECK( 23.93741 == Approx( resonances8.front().levelSpacing().value ) );
      CHECK( 20.85515 == Approx( resonances8.back().levelSpacing().value ) );
      CHECK( 8.403428e-3 == Approx( resonances8.front().elastic().value ) );
      CHECK( 7.132285e-3 == Approx( resonances8.back().elastic().value ) );
      CHECK( 1e-9 == Approx( resonances8.front().capture().value ) );
      CHECK( 1e-9 == Approx( resonances8.back().capture().value ) );
      CHECK( 0. == Approx( resonances8.front().fission().value ) );
      CHECK( 0. == Approx( resonances8.back().fission().value ) );
      CHECK( 0. == Approx( resonances8.front().competition().value ) );
      CHECK( 1e-9 == Approx( resonances8.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 9
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup9 = spingroups[9];

      auto channel9 = spingroup9.incidentChannel();
      CHECK( "n,Au197{2,1/2,3+}" == channel9.channelID() );
      CHECK( "n,Au197->n,Au197" == channel9.reactionID().symbol() );

      // resonance table
      auto table9 = spingroup9.resonanceTable();

      CHECK( 100 == table9.numberResonances() );

      auto energies9 = table9.energies();
      CHECK( 2000. == Approx( energies9.front().value ) );
      CHECK( 100000. == Approx( energies9.back().value ) );

      auto resonances9 = table9.resonances();
      CHECK( 2000. == Approx( resonances9.front().energy().value ) );
      CHECK( 100000. == Approx( resonances9.back().energy().value ) );
      CHECK( 18.04535 == Approx( resonances9.front().levelSpacing().value ) );
      CHECK( 15.71550 == Approx( resonances9.back().levelSpacing().value ) );
      CHECK( 6.334971e-3 == Approx( resonances9.front().elastic().value ) );
      CHECK( 5.374571e-3 == Approx( resonances9.back().elastic().value ) );
      CHECK( 1.292688e-1 == Approx( resonances9.front().capture().value ) );
      CHECK( 1.337047e-1 == Approx( resonances9.back().capture().value ) );
      CHECK( 0. == Approx( resonances9.front().fission().value ) );
      CHECK( 0. == Approx( resonances9.back().fission().value ) );
      CHECK( 0. == Approx( resonances9.front().competition().value ) );
      CHECK( 4.013598e-4 == Approx( resonances9.back().competition().value ) );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 10
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup10 = spingroups[10];

      auto channel10 = spingroup10.incidentChannel();
      CHECK( "n,Au197{2,1/2,4+}" == channel10.channelID() );
      CHECK( "n,Au197->n,Au197" == channel10.reactionID().symbol() );

      // resonance table
      auto table10 = spingroup10.resonanceTable();

      CHECK( 100 == table10.numberResonances() );

      auto energies10 = table10.energies();
      CHECK( 2000. == Approx( energies10.front().value ) );
      CHECK( 100000. == Approx( energies10.back().value ) );

      auto resonances10 = table10.resonances();
      CHECK( 2000. == Approx( resonances10.front().energy().value ) );
      CHECK( 100000. == Approx( resonances10.back().energy().value ) );
      CHECK( 15.08143 == Approx( resonances10.front().levelSpacing().value ) );
      CHECK( 13.12728 == Approx( resonances10.back().levelSpacing().value ) );
      CHECK( 5.294460e-3  == Approx( resonances10.front().elastic().value ) );
      CHECK( 4.489419e-3 == Approx( resonances10.back().elastic().value ) );
      CHECK( 1.292688e-1 == Approx( resonances10.front().capture().value ) );
      CHECK( 1.337047e-1 == Approx( resonances10.back().capture().value ) );
      CHECK( 0. == Approx( resonances10.front().fission().value ) );
      CHECK( 0. == Approx( resonances10.back().fission().value ) );
      CHECK( 0. == Approx( resonances10.front().competition().value ) );
      CHECK( 0. == Approx( resonances10.back().competition().value ) );
    } // THEN

    THEN( "cross sections can be reconstructed" ) {

      ReactionID elas( "n,Au197->n,Au197" );
      ReactionID capt( "n,Au197->capture" );
      Map< ReactionID, CrossSection > xs;

      xs = resonances( 2.000000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 21.717463944441640 == Approx( xs[ elas ].value ) );
      CHECK( 4.0755783688968394 == Approx( xs[ capt ].value ) );

      xs = resonances( 2.250000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 21.158864319324596 == Approx( xs[ elas ].value ) );
      CHECK( 3.7158831891462314 == Approx( xs[ capt ].value ) );

      xs = resonances( 2.500000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 20.675250662181227 == Approx( xs[ elas ].value ) );
      CHECK( 3.4212736235276160 == Approx( xs[ capt ].value ) );

      xs = resonances( 2.750000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 20.250584093091774 == Approx( xs[ elas ].value ) );
      CHECK( 3.1752190383786236 == Approx( xs[ capt ].value ) );

      xs = resonances( 3.000000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 19.873308261469475 == Approx( xs[ elas ].value ) );
      CHECK( 2.9664032270683394 == Approx( xs[ capt ].value ) );

      xs = resonances( 3.250000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 19.534856996177762 == Approx( xs[ elas ].value ) );
      CHECK( 2.7868057954688967 == Approx( xs[ capt ].value ) );

      xs = resonances( 3.500000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 19.228721582261858 == Approx( xs[ elas ].value ) );
      CHECK( 2.6305798446432438 == Approx( xs[ capt ].value ) );

      xs = resonances( 3.750000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 18.949847738094306 == Approx( xs[ elas ].value ) );
      CHECK( 2.4933566437618562 == Approx( xs[ capt ].value ) );

      xs = resonances( 4.000000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 18.694231175777620 == Approx( xs[ elas ].value ) );
      CHECK( 2.3717993044001462 == Approx( xs[ capt ].value ) );

      xs = resonances( 4.250000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 18.458664493051280 == Approx( xs[ elas ].value ) );
      CHECK( 2.2633198567515889 == Approx( xs[ capt ].value ) );

      xs = resonances( 4.500000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 18.240535580356699 == Approx( xs[ elas ].value ) );
      CHECK( 2.1658743202365174 == Approx( xs[ capt ].value ) );

      xs = resonances( 4.750000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 18.037695288729775 == Approx( xs[ elas ].value ) );
      CHECK( 2.0778259237275489 == Approx( xs[ capt ].value ) );

      xs = resonances( 5.000000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 17.848346256202223 == Approx( xs[ elas ].value ) );
      CHECK( 1.9978510296165168 == Approx( xs[ capt ].value ) );

      xs = resonances( 5.250000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 17.670983314728886 == Approx( xs[ elas ].value ) );
      CHECK( 1.9248644823191390 == Approx( xs[ capt ].value ) );

      xs = resonances( 5.500000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 17.504330554149128 == Approx( xs[ elas ].value ) );
      CHECK( 1.8579699828272891 == Approx( xs[ capt ].value ) );

      xs = resonances( 6.000000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 17.198934575078447 == Approx( xs[ elas ].value ) );
      CHECK( 1.7395805199602514 == Approx( xs[ capt ].value ) );

      xs = resonances( 6.500000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 16.925094020175777 == Approx( xs[ elas ].value ) );
      CHECK( 1.6379890049983796 == Approx( xs[ capt ].value ) );

      xs = resonances( 7.000000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 16.677457045942749 == Approx( xs[ elas ].value ) );
      CHECK( 1.5497831902779751 == Approx( xs[ capt ].value ) );

      xs = resonances( 7.500000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 16.451873438905494 == Approx( xs[ elas ].value ) );
      CHECK( 1.4724238432815815 == Approx( xs[ capt ].value ) );

      xs = resonances( 8.000000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 16.245079381249582 == Approx( xs[ elas ].value ) );
      CHECK( 1.4039818616357036 == Approx( xs[ capt ].value ) );

      xs = resonances( 8.500000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 16.054448067400635 == Approx( xs[ elas ].value ) );
      CHECK( 1.3429599210586312 == Approx( xs[ capt ].value ) );

      xs = resonances( 9.000000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 15.877848567333569 == Approx( xs[ elas ].value ) );
      CHECK( 1.2881828494631309 == Approx( xs[ capt ].value ) );

      xs = resonances( 9.500000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 15.713536471730269 == Approx( xs[ elas ].value ) );
      CHECK( 1.2387143805384573 == Approx( xs[ capt ].value ) );

      xs = resonances( 1.000000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 15.560052842156665 == Approx( xs[ elas ].value ) );
      CHECK( 1.1937949933438716 == Approx( xs[ capt ].value ) );

      xs = resonances( 1.100000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 15.280863352795697 == Approx( xs[ elas ].value ) );
      CHECK( 1.1152412333681043 == Approx( xs[ capt ].value ) );

      xs = resonances( 1.200000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 15.032547678790767 == Approx( xs[ elas ].value ) );
      CHECK( 1.0487421277131179 == Approx( xs[ capt ].value ) );

      xs = resonances( 1.300000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 14.809401184463059 == Approx( xs[ elas ].value ) );
      CHECK( .99164473585799573 == Approx( xs[ capt ].value ) );

      xs = resonances( 1.400000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 14.607109423464392 == Approx( xs[ elas ].value ) );
      CHECK( .94202874095227407 == Approx( xs[ capt ].value ) );

      xs = resonances( 1.500000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 14.422339701040150 == Approx( xs[ elas ].value ) );
      CHECK( .89847028489672787 == Approx( xs[ capt ].value ) );

      xs = resonances( 1.600000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 14.252463192348928 == Approx( xs[ elas ].value ) );
      CHECK( .85989073009023231 == Approx( xs[ capt ].value ) );

      xs = resonances( 1.700000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 14.095385829483586 == Approx( xs[ elas ].value ) );
      CHECK( .82545743403548077 == Approx( xs[ capt ].value ) );

      xs = resonances( 1.800000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 13.949403660587343 == Approx( xs[ elas ].value ) );
      CHECK( .79451597508763327 == Approx( xs[ capt ].value ) );

      xs = resonances( 1.900000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 13.813118722052089 == Approx( xs[ elas ].value ) );
      CHECK( .76654722097535133 == Approx( xs[ capt ].value ) );

      xs = resonances( 2.000000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 13.685372007352729 == Approx( xs[ elas ].value ) );
      CHECK( .74113073929425466 == Approx( xs[ capt ].value ) );

      xs = resonances( 2.100000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 13.565192872332991 == Approx( xs[ elas ].value ) );
      CHECK( .71792524679036973 == Approx( xs[ capt ].value ) );

      xs = resonances( 2.200000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 13.451756379977457 == Approx( xs[ elas ].value ) );
      CHECK( .69664766898018293 == Approx( xs[ capt ].value ) );

      xs = resonances( 2.300000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 13.344365904626565 == Approx( xs[ elas ].value ) );
      CHECK( .67706331584012869 == Approx( xs[ capt ].value ) );

      xs = resonances( 2.400000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 13.242420223864354 == Approx( xs[ elas ].value ) );
      CHECK( .65897534932829682 == Approx( xs[ capt ].value ) );

      xs = resonances( 2.500000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 13.145400474084278 == Approx( xs[ elas ].value ) );
      CHECK( .64221701974202139 == Approx( xs[ capt ].value ) );

      xs = resonances( 2.600000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 13.052853533016076 == Approx( xs[ elas ].value ) );
      CHECK( .62664580590970842 == Approx( xs[ capt ].value ) );

      xs = resonances( 2.700000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 12.964388970908846 == Approx( xs[ elas ].value ) );
      CHECK( .61214007962122519 == Approx( xs[ capt ].value ) );

      xs = resonances( 2.800000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 12.879655664965817 == Approx( xs[ elas ].value ) );
      CHECK( .59859451415292686 == Approx( xs[ capt ].value ) );

      xs = resonances( 2.900000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 12.798353821023653 == Approx( xs[ elas ].value ) );
      CHECK( .58591716944452343 == Approx( xs[ capt ].value ) );

      xs = resonances( 3.000000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 12.720209993718839 == Approx( xs[ elas ].value ) );
      CHECK( .57402899567880261 == Approx( xs[ capt ].value ) );

      xs = resonances( 3.100000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 12.644983348389196 == Approx( xs[ elas ].value ) );
      CHECK( .56285946715693913 == Approx( xs[ capt ].value ) );

      xs = resonances( 3.200000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 12.572456575897453 == Approx( xs[ elas ].value ) );
      CHECK( .55234704257736711 == Approx( xs[ capt ].value ) );

      xs = resonances( 3.300000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 12.502441049928386 == Approx( xs[ elas ].value ) );
      CHECK( .54243708160879345 == Approx( xs[ capt ].value ) );

      xs = resonances( 3.400000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 12.434759865569720 == Approx( xs[ elas ].value ) );
      CHECK( .53308068757792715 == Approx( xs[ capt ].value ) );

      xs = resonances( 3.500000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 12.369256185398221 == Approx( xs[ elas ].value ) );
      CHECK( .52423503833685547 == Approx( xs[ capt ].value ) );

      xs = resonances( 3.600000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 12.305788123703831 == Approx( xs[ elas ].value ) );
      CHECK( .51586080871686835 == Approx( xs[ capt ].value ) );

      xs = resonances( 3.700000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 12.244230125691011 == Approx( xs[ elas ].value ) );
      CHECK( .50792323174994325 == Approx( xs[ capt ].value ) );

      xs = resonances( 3.800000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 12.184460418436371 == Approx( xs[ elas ].value ) );
      CHECK( .50039098362396517 == Approx( xs[ capt ].value ) );

      xs = resonances( 3.900000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 12.126376686163534 == Approx( xs[ elas ].value ) );
      CHECK( .49323562972159007 == Approx( xs[ capt ].value ) );

      xs = resonances( 4.000000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 12.069876312417096 == Approx( xs[ elas ].value ) );
      CHECK( .48643104620781485 == Approx( xs[ capt ].value ) );

      xs = resonances( 4.100000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 12.014871495720216 == Approx( xs[ elas ].value ) );
      CHECK( .47995407192215767 == Approx( xs[ capt ].value ) );

      xs = resonances( 4.200000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 11.961280418275395 == Approx( xs[ elas ].value ) );
      CHECK( .47378322564598246 == Approx( xs[ capt ].value ) );

      xs = resonances( 4.300000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 11.909025218200741 == Approx( xs[ elas ].value ) );
      CHECK( .46789880230941822 == Approx( xs[ capt ].value ) );

      xs = resonances( 4.400000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 11.858039516347056 == Approx( xs[ elas ].value ) );
      CHECK( .46228312444519093 == Approx( xs[ capt ].value ) );

      xs = resonances( 4.500000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 11.808254819676240 == Approx( xs[ elas ].value ) );
      CHECK( .45691952085470089 == Approx( xs[ capt ].value ) );

      xs = resonances( 4.600000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 11.759612070177585 == Approx( xs[ elas ].value ) );
      CHECK( .45179302231672125 == Approx( xs[ capt ].value ) );

      xs = resonances( 4.700000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 11.712058346276915 == Approx( xs[ elas ].value ) );
      CHECK( .44688976930800184 == Approx( xs[ capt ].value ) );

      xs = resonances( 4.800000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 11.665539240319038 == Approx( xs[ elas ].value ) );
      CHECK( .44219634199449981 == Approx( xs[ capt ].value ) );

      xs = resonances( 5.000000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 11.575418026679081 == Approx( xs[ elas ].value ) );
      CHECK( .43339378347413404 == Approx( xs[ capt ].value ) );

      xs = resonances( 5.250000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 11.467800542038042 == Approx( xs[ elas ].value ) );
      CHECK( .42337771879260866 == Approx( xs[ capt ].value ) );

      xs = resonances( 5.500000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 11.365236558312915 == Approx( xs[ elas ].value ) );
      CHECK( .41432761858167666 == Approx( xs[ capt ].value ) );

      xs = resonances( 5.750000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 11.267231306400037 == Approx( xs[ elas ].value ) );
      CHECK( .40612384667209511 == Approx( xs[ capt ].value ) );

      xs = resonances( 6.000000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 11.173358376643804 == Approx( xs[ elas ].value ) );
      CHECK( .39866533103724966 == Approx( xs[ capt ].value ) );

      xs = resonances( 6.250000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 11.083250748982454 == Approx( xs[ elas ].value ) );
      CHECK( .39186548194296128 == Approx( xs[ capt ].value ) );

      xs = resonances( 6.500000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.996591797856304 == Approx( xs[ elas ].value ) );
      CHECK( .38564973620515730 == Approx( xs[ capt ].value ) );

      xs = resonances( 6.750000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.913098569184749 == Approx( xs[ elas ].value ) );
      CHECK( .37995428701372730 == Approx( xs[ capt ].value ) );

      xs = resonances( 7.000000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.832528762563291 == Approx( xs[ elas ].value ) );
      CHECK( .37472274211778367 == Approx( xs[ capt ].value ) );

      xs = resonances( 7.200000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.770030184550793 == Approx( xs[ elas ].value ) );
      CHECK( .37083893705970072 == Approx( xs[ capt ].value ) );

      xs = resonances( 7.400000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.709162951003291 == Approx( xs[ elas ].value ) );
      CHECK( .36719851443808604 == Approx( xs[ capt ].value ) );

      xs = resonances( 7.500000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.679311982327039 == Approx( xs[ elas ].value ) );
      CHECK( .36546344665758973 == Approx( xs[ capt ].value ) );

      xs = resonances( 7.775000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.599124738439293 == Approx( xs[ elas ].value ) );
      CHECK( .36096132033780043 == Approx( xs[ capt ].value ) );

      xs = resonances( 7.779600e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.597806013479772 == Approx( xs[ elas ].value ) );
      CHECK( .36088924245341347 == Approx( xs[ capt ].value ) );

      xs = resonances( 7.781000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.590144144170473 == Approx( xs[ elas ].value ) );
      CHECK( .35922309551236625 == Approx( xs[ capt ].value ) );

      xs = resonances( 7.783000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.585874533146217 == Approx( xs[ elas ].value ) );
      CHECK( .35847539590795097 == Approx( xs[ capt ].value ) );

      xs = resonances( 7.786000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.581374281798254 == Approx( xs[ elas ].value ) );
      CHECK( .35776995184313243 == Approx( xs[ capt ].value ) );

      xs = resonances( 7.789000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.577788971531614 == Approx( xs[ elas ].value ) );
      CHECK( .35725415241543368 == Approx( xs[ capt ].value ) );

      xs = resonances( 7.793000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.573728544577582 == Approx( xs[ elas ].value ) );
      CHECK( .35670935548455218 == Approx( xs[ capt ].value ) );

      xs = resonances( 7.798000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.569328123841821 == Approx( xs[ elas ].value ) );
      CHECK( .35615730221482944 == Approx( xs[ capt ].value ) );

      xs = resonances( 7.807000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.562367638551017 == Approx( xs[ elas ].value ) );
      CHECK( .35535029540629282 == Approx( xs[ capt ].value ) );

      xs = resonances( 7.820000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.553855525071533 == Approx( xs[ elas ].value ) );
      CHECK( .35444371200797808 == Approx( xs[ capt ].value ) );

      xs = resonances( 7.835000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.545168766148270 == Approx( xs[ elas ].value ) );
      CHECK( .35358937715373634 == Approx( xs[ capt ].value ) );

      xs = resonances( 7.862000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.531256719810232 == Approx( xs[ elas ].value ) );
      CHECK( .35232812691864424 == Approx( xs[ capt ].value ) );

      xs = resonances( 7.890000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.518207086136112 == Approx( xs[ elas ].value ) );
      CHECK( .35123198331241617 == Approx( xs[ capt ].value ) );

      xs = resonances( 7.945000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.494873327949790 == Approx( xs[ elas ].value ) );
      CHECK( .34941374387606017 == Approx( xs[ capt ].value ) );

      xs = resonances( 8.000000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.473367161114492 == Approx( xs[ elas ].value ) );
      CHECK( .34784907453939978 == Approx( xs[ capt ].value ) );

      xs = resonances( 8.100000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.436985846496702 == Approx( xs[ elas ].value ) );
      CHECK( .34536205456504243 == Approx( xs[ capt ].value ) );

      xs = resonances( 8.200000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.402818918876486 == Approx( xs[ elas ].value ) );
      CHECK( .34315393355474422 == Approx( xs[ capt ].value ) );

      xs = resonances( 8.300000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.370142011134579 == Approx( xs[ elas ].value ) );
      CHECK( .34112734352676644 == Approx( xs[ capt ].value ) );

      xs = resonances( 8.500000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.307913234827785 == Approx( xs[ elas ].value ) );
      CHECK( .33744533944640664 == Approx( xs[ capt ].value ) );

      xs = resonances( 8.750000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.234301122794545 == Approx( xs[ elas ].value ) );
      CHECK( .33333221045367423 == Approx( xs[ capt ].value ) );

      xs = resonances( 9.000000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.164125045943718 == Approx( xs[ elas ].value ) );
      CHECK( .32961757753543425 == Approx( xs[ capt ].value ) );

      xs = resonances( 9.200000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.109983604516954 == Approx( xs[ elas ].value ) );
      CHECK( .32687573457210550 == Approx( xs[ capt ].value ) );

      xs = resonances( 9.350000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.070394126277289 == Approx( xs[ elas ].value ) );
      CHECK( .32493459235800237 == Approx( xs[ capt ].value ) );

      xs = resonances( 9.500000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 10.031596813680748 == Approx( xs[ elas ].value ) );
      CHECK( .32308194131781010 == Approx( xs[ capt ].value ) );

      xs = resonances( 9.750000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 9.9685514411795690 == Approx( xs[ elas ].value ) );
      CHECK( .32017160113668208 == Approx( xs[ capt ].value ) );

      xs = resonances( 1.000000e+5 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 9.9073506528349853 == Approx( xs[ elas ].value ) );
      CHECK( .31745831628731552 == Approx( xs[ capt ].value ) );
    } // THEN
  } // GIVEN

  GIVEN( "valid ENDF data for Ag107" ) {

    std::string string = Ag107();
    auto begin = string.begin();
    auto end = string.end();
    long lineNumber = 1;

    njoy::ENDFtk::HeadRecord head( begin, end, lineNumber );
    njoy::ENDFtk::section::Type< 2, 151 > endf( head, begin, end, lineNumber, 4725 );
    ResonanceRange endfResonanceRange = endf.isotopes().front().resonanceRanges().front();

    auto resonances = fromENDF( endfResonanceRange, neutronMass, elementaryCharge, ParticleID( "n" ), ParticleID( "Ag107" ) );

    THEN( "cross sections can be reconstructed" ) {

      ReactionID elas( "n,Ag107->n,Ag107" );
      ReactionID capt( "n,Ag107->capture" );
      Map< ReactionID, CrossSection > xs;

      xs = resonances( 6.500000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.399512 == Approx( xs[ elas ].value ) );
      CHECK( 1.570063 == Approx( xs[ capt ].value ) );

      xs = resonances( 7.400000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.363600 == Approx( xs[ elas ].value ) );
      CHECK( 1.503679 == Approx( xs[ capt ].value ) );

      xs = resonances( 8.500000e+3 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.311894 == Approx( xs[ elas ].value ) );
      CHECK( 1.433198 == Approx( xs[ capt ].value ) );

      xs = resonances( 1.000000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.260546 == Approx( xs[ elas ].value ) );
      CHECK( 1.354022 == Approx( xs[ capt ].value ) );

      xs = resonances( 1.067000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.242543 == Approx( xs[ elas ].value ) );
      CHECK( 1.323163 == Approx( xs[ capt ].value ) );

      xs = resonances( 1.25000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.147463 == Approx( xs[ elas ].value ) );
      CHECK( 1.277123 == Approx( xs[ capt ].value ) );

      xs = resonances( 1.394000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.070861 == Approx( xs[ elas ].value ) );
      CHECK( 1.247018 == Approx( xs[ capt ].value ) );

      xs = resonances( 1.721000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 6.938058 == Approx( xs[ elas ].value ) );
      CHECK( 1.137130 == Approx( xs[ capt ].value ) );

      xs = resonances( 2.000000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 6.946458 == Approx( xs[ elas ].value ) );
      CHECK( 1.072939 == Approx( xs[ capt ].value ) );

      xs = resonances( 2.500000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 6.994424 == Approx( xs[ elas ].value ) );
      CHECK( 9.752918e-1 == Approx( xs[ capt ].value ) );

      xs = resonances( 3.000000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.066323 == Approx( xs[ elas ].value ) );
      CHECK( 8.927133e-1 == Approx( xs[ capt ].value ) );

      xs = resonances( 3.181000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.095747 == Approx( xs[ elas ].value ) );
      CHECK( 8.656563e-1 == Approx( xs[ capt ].value ) );

      xs = resonances( 3.500000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.151891 == Approx( xs[ elas ].value ) );
      CHECK( 8.192528e-1 == Approx( xs[ capt ].value ) );

      xs = resonances( 4.000000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.228131 == Approx( xs[ elas ].value ) );
      CHECK( 7.541488e-1 == Approx( xs[ capt ].value ) );

      xs = resonances( 4.641000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.307337 == Approx( xs[ elas ].value ) );
      CHECK( 6.821784e-1 == Approx( xs[ capt ].value ) );

      xs = resonances( 5.371000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.350422 == Approx( xs[ elas ].value ) );
      CHECK( 6.355932e-1 == Approx( xs[ capt ].value ) );

      xs = resonances( 6.100000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.384957 == Approx( xs[ elas ].value ) );
      CHECK( 5.764472e-1 == Approx( xs[ capt ].value ) );

      xs = resonances( 6.831000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.440526 == Approx( xs[ elas ].value ) );
      CHECK( 5.324228e-1 == Approx( xs[ capt ].value ) );

      xs = resonances( 7.561000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.435362 == Approx( xs[ elas ].value ) );
      CHECK( 5.064053e-1 == Approx( xs[ capt ].value ) );

      xs = resonances( 9.174000e+4 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.370678 == Approx( xs[ elas ].value ) );
      CHECK( 4.829130e-1 == Approx( xs[ capt ].value ) );

      xs = resonances( 1.000000e+5 * electronVolt );
      CHECK( 2 == xs.size() );
      CHECK( 7.377410 == Approx( xs[ elas ].value ) );
      CHECK( 4.705200e-1 == Approx( xs[ capt ].value ) );
    } // THEN
  } // GIVEN
} // SCENARIO

std::string Na22() {

  // Na22 ENDF/B-VIII.0 LRU=2 resonance evaluation

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

std::string Pu239() {

  // Pu239 ENDF/B-VIII.0 LRU=2 resonance evaluation

  return
    " 9.423900+4 2.369986+2          0          0          1          09437 2151     \n"
    " 9.423900+4 1.000000+0          0          1          1          09437 2151     \n"
    " 2.500000+3 3.000000+4          2          2          0          09437 2151     \n"
    " 5.000000-1 9.460000-1          1          0          2          09437 2151     \n"
    " 2.369986+2 0.000000+0          0          0          2          09437 2151     \n"
    " 0.000000+0 0.000000+0          2          0        432         719437 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 1.000000+0 0.000000+0 2.000000+09437 2151     \n"
    " 2.500000+3 8.917200+0 0.000000+0 9.508100-4 4.070000-2 2.842000+09437 2151     \n"
    " 2.550000+3 8.915500+0 0.000000+0 8.673000-4 4.070000-2 4.020000-19437 2151     \n"
    " 2.650000+3 8.913900+0 0.000000+0 1.113200-3 4.070000-2 2.841000+09437 2151     \n"
    " 2.750000+3 8.912200+0 0.000000+0 8.952900-4 4.070000-2 2.840000+09437 2151     \n"
    " 2.850000+3 8.910500+0 0.000000+0 9.309800-4 4.070000-2 2.840000+09437 2151     \n"
    " 2.950000+3 8.908800+0 0.000000+0 1.450600-3 4.070000-2 7.840000-19437 2151     \n"
    " 3.050000+3 8.907200+0 0.000000+0 7.845300-4 4.070000-2 2.839000+09437 2151     \n"
    " 3.150000+3 8.905500+0 0.000000+0 1.079100-3 4.070000-2 2.838000+09437 2151     \n"
    " 3.250000+3 8.903800+0 0.000000+0 1.065500-3 4.070000-2 2.838000+09437 2151     \n"
    " 3.350000+3 8.902200+0 0.000000+0 6.778100-4 4.070000-2 1.150000+09437 2151     \n"
    " 3.450000+3 8.900500+0 0.000000+0 9.239600-4 4.070000-2 9.960000-19437 2151     \n"
    " 3.550000+3 8.898800+0 0.000000+0 7.025700-4 4.070000-2 2.836000+09437 2151     \n"
    " 3.650000+3 8.897100+0 0.000000+0 1.215900-3 4.070000-2 3.930000-19437 2151     \n"
    " 3.750000+3 8.895500+0 0.000000+0 7.980200-4 4.070000-2 2.835000+09437 2151     \n"
    " 3.850000+3 8.893800+0 0.000000+0 1.082000-3 4.070000-2 2.835000+09437 2151     \n"
    " 3.950000+3 8.892100+0 0.000000+0 1.093700-3 4.070000-2 2.834000+09437 2151     \n"
    " 4.125000+3 8.889100+0 0.000000+0 9.226000-4 4.070000-2 2.829000+09437 2151     \n"
    " 4.375000+3 8.884900+0 0.000000+0 9.985000-4 4.070000-2 2.828000+09437 2151     \n"
    " 4.625000+3 8.880800+0 0.000000+0 9.217000-4 4.070000-2 2.827000+09437 2151     \n"
    " 4.875000+3 8.876600+0 0.000000+0 9.974000-4 4.070000-2 2.826000+09437 2151     \n"
    " 5.125000+3 8.872400+0 0.000000+0 1.011400-3 4.070000-2 2.820000+09437 2151     \n"
    " 5.375000+3 8.868200+0 0.000000+0 1.011000-3 4.070000-2 2.820000+09437 2151     \n"
    " 5.625000+3 8.864100+0 0.000000+0 1.010500-3 4.070000-2 2.820000+09437 2151     \n"
    " 5.875000+3 8.859900+0 0.000000+0 1.010000-3 4.070000-2 2.820000+09437 2151     \n"
    " 6.125000+3 8.855700+0 0.000000+0 1.009500-3 4.070000-2 2.820000+09437 2151     \n"
    " 6.375000+3 8.851500+0 0.000000+0 1.009100-3 4.070000-2 2.820000+09437 2151     \n"
    " 6.625000+3 8.847300+0 0.000000+0 1.008600-3 4.070000-2 9.320000-19437 2151     \n"
    " 6.875000+3 8.843100+0 0.000000+0 1.008100-3 4.070000-2 1.056500+09437 2151     \n"
    " 7.125000+3 8.839000+0 0.000000+0 1.019800-3 4.070000-2 4.814000-19437 2151     \n"
    " 7.375000+3 8.834800+0 0.000000+0 1.014200-3 4.070000-2 1.130000+09437 2151     \n"
    " 7.625000+3 8.830700+0 0.000000+0 1.013100-3 4.070000-2 1.309500+09437 2151     \n"
    " 7.875000+3 8.826500+0 0.000000+0 1.013400-3 4.070000-2 2.078000+09437 2151     \n"
    " 8.125000+3 8.822300+0 0.000000+0 1.005700-3 4.070000-2 2.140000+09437 2151     \n"
    " 8.375000+3 8.818200+0 0.000000+0 1.005300-3 4.070000-2 2.807000+09437 2151     \n"
    " 8.625000+3 8.814000+0 0.000000+0 1.004800-3 4.070000-2 2.806000+09437 2151     \n"
    " 8.875000+3 8.809900+0 0.000000+0 1.004300-3 4.070000-2 2.804000+09437 2151     \n"
    " 9.125000+3 8.805700+0 0.000000+0 8.810000-4 4.070000-2 2.803000+09437 2151     \n"
    " 9.375000+3 8.801400+0 0.000000+0 8.622000-4 4.070000-2 1.249000+09437 2151     \n"
    " 9.625000+3 8.797400+0 0.000000+0 9.643000-4 4.070000-2 1.250000+09437 2151     \n"
    " 9.875000+3 8.793200+0 0.000000+0 9.609000-4 4.070000-2 1.250000+09437 2151     \n"
    " 1.025000+4 8.787100+0 0.000000+0 9.823000-4 4.070000-2 1.182000+09437 2151     \n"
    " 1.075000+4 8.778900+0 0.000000+0 9.815000-4 4.070000-2 2.794000+09437 2151     \n"
    " 1.125000+4 8.770500+0 0.000000+0 9.805000-4 4.070000-2 2.792000+09437 2151     \n"
    " 1.175000+4 8.762200+0 0.000000+0 9.797000-4 4.070000-2 2.789000+09437 2151     \n"
    " 1.225000+4 8.754100+0 0.000000+0 9.787000-4 4.070000-2 2.786000+09437 2151     \n"
    " 1.275000+4 8.745800+0 0.000000+0 9.778000-4 4.070000-2 2.784000+09437 2151     \n"
    " 1.325000+4 8.737600+0 0.000000+0 9.768000-4 4.070000-2 2.781000+09437 2151     \n"
    " 1.375000+4 8.729400+0 0.000000+0 9.759000-4 4.070000-2 2.779000+09437 2151     \n"
    " 1.425000+4 8.721200+0 0.000000+0 9.750000-4 4.070000-2 6.620000-19437 2151     \n"
    " 1.475000+4 8.712900+0 0.000000+0 9.740000-4 4.070000-2 2.773000+09437 2151     \n"
    " 1.525000+4 8.704700+0 0.000000+0 9.732000-4 4.070000-2 2.771000+09437 2151     \n"
    " 1.575000+4 8.696500+0 0.000000+0 9.723000-4 4.070000-2 2.768000+09437 2151     \n"
    " 1.625000+4 8.688300+0 0.000000+0 9.714000-4 4.070000-2 2.766000+09437 2151     \n"
    " 1.675000+4 8.680100+0 0.000000+0 9.694000-4 4.070000-2 2.763000+09437 2151     \n"
    " 1.725000+4 8.672000+0 0.000000+0 9.697000-4 4.070000-2 9.465000-19437 2151     \n"
    " 1.775000+4 8.663800+0 0.000000+0 9.687000-4 4.070000-2 2.758000+09437 2151     \n"
    " 1.825000+4 8.655700+0 0.000000+0 9.678000-4 4.070000-2 2.755000+09437 2151     \n"
    " 1.875000+4 8.647500+0 0.000000+0 9.669000-4 4.070000-2 2.753000+09437 2151     \n"
    " 1.925000+4 8.639300+0 0.000000+0 9.660000-4 4.070000-2 2.750000+09437 2151     \n"
    " 1.975000+4 8.631200+0 0.000000+0 9.651000-4 4.070000-2 2.747000+09437 2151     \n"
    " 2.050000+4 8.619000+0 0.000000+0 9.636000-4 4.070000-2 2.743000+09437 2151     \n"
    " 2.150000+4 8.602800+0 0.000000+0 9.618000-4 4.070000-2 2.738000+09437 2151     \n"
    " 2.250000+4 8.586700+0 0.000000+0 9.600000-4 4.070000-2 2.733000+09437 2151     \n"
    " 2.350000+4 8.570500+0 0.000000+0 9.582000-4 4.070000-2 2.728000+09437 2151     \n"
    " 2.450000+4 8.554400+0 0.000000+0 9.564000-4 4.070000-2 2.723000+09437 2151     \n"
    " 2.550000+4 8.538200+0 0.000000+0 9.546000-4 4.070000-2 2.718000+09437 2151     \n"
    " 2.650000+4 8.522200+0 0.000000+0 9.528000-4 4.070000-2 2.713000+09437 2151     \n"
    " 2.750000+4 8.506200+0 0.000000+0 9.510000-4 4.070000-2 2.708000+09437 2151     \n"
    " 2.850000+4 8.490100+0 0.000000+0 9.492000-4 4.070000-2 2.703000+09437 2151     \n"
    " 2.950000+4 8.474100+0 0.000000+0 9.474000-4 4.070000-2 2.697000+09437 2151     \n"
    " 3.000000+4 8.465900+0 0.000000+0 9.465000-4 4.070000-2 2.693000+09437 2151     \n"
    " 1.000000+0 0.000000+0          2          0        432         719437 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 1.000000+0 0.000000+0 1.000000+09437 2151     \n"
    " 2.500000+3 3.044300+0 0.000000+0 3.246100-4 4.030000-2 7.100000-29437 2151     \n"
    " 2.550000+3 3.043700+0 0.000000+0 2.960900-4 4.030000-2 3.400000-29437 2151     \n"
    " 2.650000+3 3.043100+0 0.000000+0 3.800500-4 4.030000-2 1.120000-29437 2151     \n"
    " 2.750000+3 3.042500+0 0.000000+0 3.056100-4 4.030000-2 1.028000-19437 2151     \n"
    " 2.850000+3 3.042000+0 0.000000+0 3.178300-4 4.030000-2 9.160000-29437 2151     \n"
    " 2.950000+3 3.041400+0 0.000000+0 4.954000-4 4.030000-2 1.580000-29437 2151     \n"
    " 3.050000+3 3.040800+0 0.000000+0 2.666500-4 4.030000-2 5.200000-29437 2151     \n"
    " 3.150000+3 3.040200+0 0.000000+0 3.684000-4 4.030000-2 1.420000-19437 2151     \n"
    " 3.250000+3 3.039700+0 0.000000+0 3.637500-4 4.030000-2 6.230000-29437 2151     \n"
    " 3.350000+3 3.039100+0 0.000000+0 2.314000-4 4.030000-2 1.203000-29437 2151     \n"
    " 3.450000+3 3.038500+0 0.000000+0 3.154300-4 4.030000-2 1.373000-29437 2151     \n"
    " 3.550000+3 3.038000+0 0.000000+0 2.398600-4 4.030000-2 2.810000-29437 2151     \n"
    " 3.650000+3 3.037400+0 0.000000+0 4.151200-4 4.030000-2 1.660000-29437 2151     \n"
    " 3.750000+3 3.036800+0 0.000000+0 2.724300-4 4.030000-2 8.070000-29437 2151     \n"
    " 3.850000+3 3.036200+0 0.000000+0 3.693900-4 4.030000-2 5.380000-29437 2151     \n"
    " 3.950000+3 3.035700+0 0.000000+0 3.733900-4 4.030000-2 2.340000-29437 2151     \n"
    " 4.125000+3 3.034600+0 0.000000+0 3.148000-4 4.030000-2 1.098000-29437 2151     \n"
    " 4.375000+3 3.033200+0 0.000000+0 3.400000-4 4.030000-2 2.170000-29437 2151     \n"
    " 4.625000+3 3.031800+0 0.000000+0 3.145000-4 4.030000-2 5.070000-29437 2151     \n"
    " 4.875000+3 3.030400+0 0.000000+0 3.403000-4 4.030000-2 7.500000-39437 2151     \n"
    " 5.125000+3 3.029000+0 0.000000+0 3.453000-4 4.030000-2 2.360000-29437 2151     \n"
    " 5.375000+3 3.027600+0 0.000000+0 3.451000-4 4.030000-2 1.500000-29437 2151     \n"
    " 5.625000+3 3.026100+0 0.000000+0 3.450000-4 4.030000-2 2.300000-29437 2151     \n"
    " 5.875000+3 3.024700+0 0.000000+0 3.448000-4 4.030000-2 2.220000-29437 2151     \n"
    " 6.125000+3 3.023200+0 0.000000+0 3.446000-4 4.030000-2 3.910000-29437 2151     \n"
    " 6.375000+3 3.021800+0 0.000000+0 3.445000-4 4.030000-2 1.270000-29437 2151     \n"
    " 6.625000+3 3.020400+0 0.000000+0 3.443400-4 4.030000-2 2.210000-29437 2151     \n"
    " 6.875000+3 3.018900+0 0.000000+0 3.441700-4 4.030000-2 1.555000-29437 2151     \n"
    " 7.125000+3 3.017500+0 0.000000+0 3.483000-4 4.030000-2 2.920000-29437 2151     \n"
    " 7.375000+3 3.016100+0 0.000000+0 3.464000-4 4.030000-2 6.930000-29437 2151     \n"
    " 7.625000+3 3.014700+0 0.000000+0 3.463000-4 4.030000-2 7.648000-29437 2151     \n"
    " 7.875000+3 3.013300+0 0.000000+0 3.463000-4 4.030000-2 1.544000-29437 2151     \n"
    " 8.125000+3 3.011900+0 0.000000+0 3.433000-4 4.030000-2 4.330000-29437 2151     \n"
    " 8.375000+3 3.010500+0 0.000000+0 3.432000-4 4.030000-2 6.650000-29437 2151     \n"
    " 8.625000+3 3.009000+0 0.000000+0 3.430000-4 4.030000-2 3.180000-29437 2151     \n"
    " 8.875000+3 3.007600+0 0.000000+0 3.429000-4 4.030000-2 5.860000-29437 2151     \n"
    " 9.125000+3 3.006200+0 0.000000+0 2.986000-4 4.030000-2 3.264000-29437 2151     \n"
    " 9.375000+3 3.004800+0 0.000000+0 2.923000-4 4.030000-2 3.136000-29437 2151     \n"
    " 9.625000+3 3.003300+0 0.000000+0 3.292000-4 4.030000-2 7.470000-29437 2151     \n"
    " 9.875000+3 3.001900+0 0.000000+0 3.260600-4 4.030000-2 3.570000-29437 2151     \n"
    " 1.025000+4 2.999800+0 0.000000+0 3.354000-4 4.030000-2 3.720000-29437 2151     \n"
    " 1.075000+4 2.997000+0 0.000000+0 3.351000-4 4.030000-2 4.610000-29437 2151     \n"
    " 1.125000+4 2.994100+0 0.000000+0 3.347000-4 4.030000-2 3.510000-29437 2151     \n"
    " 1.175000+4 2.991300+0 0.000000+0 3.344000-4 4.030000-2 3.501100-29437 2151     \n"
    " 1.225000+4 2.988500+0 0.000000+0 3.341000-4 4.030000-2 4.810000-29437 2151     \n"
    " 1.275000+4 2.985700+0 0.000000+0 3.337000-4 4.030000-2 4.710000-29437 2151     \n"
    " 1.325000+4 2.982900+0 0.000000+0 3.335000-4 4.030000-2 4.960000-29437 2151     \n"
    " 1.375000+4 2.980100+0 0.000000+0 3.332000-4 4.030000-2 3.510000-29437 2151     \n"
    " 1.425000+4 2.977300+0 0.000000+0 3.329000-4 4.030000-2 3.660000-29437 2151     \n"
    " 1.475000+4 2.974400+0 0.000000+0 3.325000-4 4.030000-2 5.090000-29437 2151     \n"
    " 1.525000+4 2.971600+0 0.000000+0 3.322000-4 4.030000-2 6.860000-29437 2151     \n"
    " 1.575000+4 2.968800+0 0.000000+0 3.319000-4 4.030000-2 4.220000-29437 2151     \n"
    " 1.625000+4 2.966000+0 0.000000+0 3.315000-4 4.030000-2 6.110000-29437 2151     \n"
    " 1.675000+4 2.963200+0 0.000000+0 3.313000-4 4.030000-2 3.740000-29437 2151     \n"
    " 1.725000+4 2.960400+0 0.000000+0 3.310000-4 4.030000-2 4.130000-29437 2151     \n"
    " 1.775000+4 2.957600+0 0.000000+0 3.307000-4 4.030000-2 8.470000-29437 2151     \n"
    " 1.825000+4 2.954900+0 0.000000+0 3.303000-4 4.030000-2 5.690000-29437 2151     \n"
    " 1.875000+4 2.952100+0 0.000000+0 3.301000-4 4.030000-2 4.590000-29437 2151     \n"
    " 1.925000+4 2.949300+0 0.000000+0 3.298000-4 4.030000-2 6.770000-29437 2151     \n"
    " 1.975000+4 2.946500+0 0.000000+0 3.295000-4 4.030000-2 7.190000-29437 2151     \n"
    " 2.050000+4 2.942400+0 0.000000+0 3.290000-4 4.030000-2 6.100000-29437 2151     \n"
    " 2.150000+4 2.936900+0 0.000000+0 3.284000-4 4.030000-2 6.040000-29437 2151     \n"
    " 2.250000+4 2.931300+0 0.000000+0 3.278000-4 4.030000-2 3.320000-29437 2151     \n"
    " 2.350000+4 2.925800+0 0.000000+0 3.271000-4 4.030000-2 6.640000-29437 2151     \n"
    " 2.450000+4 2.920200+0 0.000000+0 3.265000-4 4.030000-2 7.190000-29437 2151     \n"
    " 2.550000+4 2.914700+0 0.000000+0 3.259000-4 4.030000-2 5.500000-29437 2151     \n"
    " 2.650000+4 2.909300+0 0.000000+0 3.253000-4 4.030000-2 7.870000-29437 2151     \n"
    " 2.750000+4 2.903800+0 0.000000+0 3.246000-4 4.030000-2 6.090000-29437 2151     \n"
    " 2.850000+4 2.898400+0 0.000000+0 3.240000-4 4.030000-2 7.070000-29437 2151     \n"
    " 2.950000+4 2.892900+0 0.000000+0 3.234000-4 4.030000-2 8.710000-29437 2151     \n"
    " 3.000000+4 2.890100+0 0.000000+0 3.231000-4 4.030000-2 7.647000-29437 2151     \n"
    " 2.369986+2 0.000000+0          1          0          3          09437 2151     \n"
    " 0.000000+0 0.000000+0          2          0        432         719437 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 1.000000+0 0.000000+0 0.000000+09437 2151     \n"
    " 2.500000+3 8.917200+0 0.000000+0 1.573800-3 1.150000-2 0.000000+09437 2151     \n"
    " 2.550000+3 8.915500+0 0.000000+0 1.449800-3 1.150000-2 0.000000+09437 2151     \n"
    " 2.650000+3 8.913900+0 0.000000+0 1.824400-3 1.150000-2 0.000000+09437 2151     \n"
    " 2.750000+3 8.912200+0 0.000000+0 1.497800-3 1.150000-2 0.000000+09437 2151     \n"
    " 2.850000+3 8.910500+0 0.000000+0 1.555000-3 1.150000-2 0.000000+09437 2151     \n"
    " 2.950000+3 8.908800+0 0.000000+0 2.338500-3 1.150000-2 0.000000+09437 2151     \n"
    " 3.050000+3 8.907200+0 0.000000+0 1.337900-3 1.150000-2 0.000000+09437 2151     \n"
    " 3.150000+3 8.905500+0 0.000000+0 1.785600-3 1.150000-2 0.000000+09437 2151     \n"
    " 3.250000+3 8.903800+0 0.000000+0 1.767600-3 1.150000-2 0.000000+09437 2151     \n"
    " 3.350000+3 8.902200+0 0.000000+0 1.182400-3 1.150000-2 0.000000+09437 2151     \n"
    " 3.450000+3 8.900500+0 0.000000+0 1.558300-3 1.150000-2 0.000000+09437 2151     \n"
    " 3.550000+3 8.898800+0 0.000000+0 1.224400-3 1.150000-2 0.000000+09437 2151     \n"
    " 3.650000+3 8.897100+0 0.000000+0 2.003600-3 1.150000-2 0.000000+09437 2151     \n"
    " 3.750000+3 8.895500+0 0.000000+0 1.373800-3 1.150000-2 0.000000+09437 2151     \n"
    " 3.850000+3 8.893800+0 0.000000+0 1.806300-3 1.150000-2 0.000000+09437 2151     \n"
    " 3.950000+3 8.892100+0 0.000000+0 1.826000-3 1.150000-2 0.000000+09437 2151     \n"
    " 4.125000+3 8.889100+0 0.000000+0 1.349000-3 1.150000-2 0.000000+09437 2151     \n"
    " 4.375000+3 8.884900+0 0.000000+0 1.458600-3 1.150000-2 0.000000+09437 2151     \n"
    " 4.625000+3 8.880800+0 0.000000+0 1.348000-3 1.150000-2 0.000000+09437 2151     \n"
    " 4.875000+3 8.876600+0 0.000000+0 1.459600-3 1.150000-2 0.000000+09437 2151     \n"
    " 5.125000+3 8.872400+0 0.000000+0 1.508300-3 1.150000-2 0.000000+09437 2151     \n"
    " 5.375000+3 8.868200+0 0.000000+0 1.507600-3 1.150000-2 0.000000+09437 2151     \n"
    " 5.625000+3 8.864100+0 0.000000+0 1.506900-3 1.150000-2 0.000000+09437 2151     \n"
    " 5.875000+3 8.859900+0 0.000000+0 1.506200-3 1.150000-2 0.000000+09437 2151     \n"
    " 6.125000+3 8.855700+0 0.000000+0 1.505500-3 1.150000-2 0.000000+09437 2151     \n"
    " 6.375000+3 8.851500+0 0.000000+0 1.504800-3 1.150000-2 0.000000+09437 2151     \n"
    " 6.625000+3 8.847300+0 0.000000+0 1.437500-3 1.150000-2 0.000000+09437 2151     \n"
    " 6.875000+3 8.843100+0 0.000000+0 1.437500-3 1.150000-2 0.000000+09437 2151     \n"
    " 7.125000+3 8.839000+0 0.000000+0 1.484100-3 1.170000-2 0.000000+09437 2151     \n"
    " 7.375000+3 8.834800+0 0.000000+0 1.483400-3 1.150000-2 0.000000+09437 2151     \n"
    " 7.625000+3 8.830700+0 0.000000+0 1.482800-3 1.150000-2 0.000000+09437 2151     \n"
    " 7.875000+3 8.826500+0 0.000000+0 1.481100-3 1.150000-2 0.000000+09437 2151     \n"
    " 8.125000+3 8.822300+0 0.000000+0 1.499800-3 1.150000-2 0.000000+09437 2151     \n"
    " 8.375000+3 8.818200+0 0.000000+0 1.499100-3 1.150000-2 0.000000+09437 2151     \n"
    " 8.625000+3 8.814000+0 0.000000+0 1.498400-3 1.150000-2 0.000000+09437 2151     \n"
    " 8.875000+3 8.809900+0 0.000000+0 1.497700-3 1.150000-2 0.000000+09437 2151     \n"
    " 9.125000+3 8.805700+0 0.000000+0 1.479500-3 1.150000-2 0.000000+09437 2151     \n"
    " 9.375000+3 8.801400+0 0.000000+0 1.466200-3 1.150000-2 0.000000+09437 2151     \n"
    " 9.625000+3 8.797400+0 0.000000+0 1.438700-3 1.150000-2 0.000000+09437 2151     \n"
    " 9.875000+3 8.793200+0 0.000000+0 1.424500-3 1.150000-2 0.000000+09437 2151     \n"
    " 1.025000+4 8.787100+0 0.000000+0 1.463500-3 1.150000-2 0.000000+09437 2151     \n"
    " 1.075000+4 8.778900+0 0.000000+0 1.462300-3 1.150000-2 0.000000+09437 2151     \n"
    " 1.125000+4 8.770500+0 0.000000+0 1.461100-3 1.150000-2 0.000000+09437 2151     \n"
    " 1.175000+4 8.762200+0 0.000000+0 1.459800-3 1.150000-2 0.000000+09437 2151     \n"
    " 1.225000+4 8.754100+0 0.000000+0 1.458600-3 1.150000-2 0.000000+09437 2151     \n"
    " 1.275000+4 8.745800+0 0.000000+0 1.457400-3 1.150000-2 0.000000+09437 2151     \n"
    " 1.325000+4 8.737600+0 0.000000+0 1.456200-3 1.150000-2 0.000000+09437 2151     \n"
    " 1.375000+4 8.729400+0 0.000000+0 1.455000-3 1.150000-2 0.000000+09437 2151     \n"
    " 1.425000+4 8.721200+0 0.000000+0 1.453700-3 1.150000-2 0.000000+09437 2151     \n"
    " 1.475000+4 8.712900+0 0.000000+0 1.452400-3 1.150000-2 0.000000+09437 2151     \n"
    " 1.525000+4 8.704700+0 0.000000+0 1.451200-3 1.150000-2 0.000000+09437 2151     \n"
    " 1.575000+4 8.696500+0 0.000000+0 1.449900-3 1.150000-2 0.000000+09437 2151     \n"
    " 1.625000+4 8.688300+0 0.000000+0 1.448500-3 1.150000-2 0.000000+09437 2151     \n"
    " 1.675000+4 8.680100+0 0.000000+0 1.447400-3 1.150000-2 0.000000+09437 2151     \n"
    " 1.725000+4 8.672000+0 0.000000+0 1.445700-3 1.150000-2 0.000000+09437 2151     \n"
    " 1.775000+4 8.663800+0 0.000000+0 1.444400-3 1.150000-2 0.000000+09437 2151     \n"
    " 1.825000+4 8.655700+0 0.000000+0 1.443000-3 1.150000-2 0.000000+09437 2151     \n"
    " 1.875000+4 8.647500+0 0.000000+0 1.441700-3 1.150000-2 0.000000+09437 2151     \n"
    " 1.925000+4 8.639300+0 0.000000+0 1.440200-3 1.150000-2 0.000000+09437 2151     \n"
    " 1.975000+4 8.631200+0 0.000000+0 1.438900-3 1.150000-2 0.000000+09437 2151     \n"
    " 2.050000+4 8.619000+0 0.000000+0 1.436800-3 1.150000-2 0.000000+09437 2151     \n"
    " 2.150000+4 8.602800+0 0.000000+0 1.434100-3 1.150000-2 0.000000+09437 2151     \n"
    " 2.250000+4 8.586700+0 0.000000+0 1.431400-3 1.150000-2 0.000000+09437 2151     \n"
    " 2.350000+4 8.570500+0 0.000000+0 1.428700-3 1.150000-2 0.000000+09437 2151     \n"
    " 2.450000+4 8.554400+0 0.000000+0 1.426000-3 1.150000-2 0.000000+09437 2151     \n"
    " 2.550000+4 8.538200+0 0.000000+0 1.423300-3 1.150000-2 0.000000+09437 2151     \n"
    " 2.650000+4 8.522200+0 0.000000+0 1.420600-3 1.150000-2 0.000000+09437 2151     \n"
    " 2.750000+4 8.506200+0 0.000000+0 1.418000-3 1.150000-2 0.000000+09437 2151     \n"
    " 2.850000+4 8.490100+0 0.000000+0 1.415300-3 1.150000-2 0.000000+09437 2151     \n"
    " 2.950000+4 8.474100+0 0.000000+0 1.412600-3 1.150000-2 0.000000+09437 2151     \n"
    " 3.000000+4 8.465900+0 0.000000+0 1.411300-3 1.150000-2 0.000000+09437 2151     \n"
    " 1.000000+0 0.000000+0          2          0        432         719437 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 2.000000+0 0.000000+0 2.000000+09437 2151     \n"
    " 2.500000+3 3.044300+0 0.000000+0 2.686450-4 3.035000-2 9.720000-19437 2151     \n"
    " 2.550000+3 3.043700+0 0.000000+0 2.474850-4 3.035000-2 9.720000-19437 2151     \n"
    " 2.650000+3 3.043100+0 0.000000+0 3.114150-4 3.035000-2 9.720000-19437 2151     \n"
    " 2.750000+3 3.042500+0 0.000000+0 2.556600-4 3.035000-2 9.710000-19437 2151     \n"
    " 2.850000+3 3.042000+0 0.000000+0 2.654300-4 3.035000-2 9.710000-19437 2151     \n"
    " 2.950000+3 3.041400+0 0.000000+0 3.991700-4 3.035000-2 9.710000-19437 2151     \n"
    " 3.050000+3 3.040800+0 0.000000+0 2.283650-4 3.035000-2 9.710000-19437 2151     \n"
    " 3.150000+3 3.040200+0 0.000000+0 3.047950-4 3.035000-2 9.710000-19437 2151     \n"
    " 3.250000+3 3.039700+0 0.000000+0 3.017200-4 3.035000-2 9.700000-19437 2151     \n"
    " 3.350000+3 3.039100+0 0.000000+0 2.018250-4 3.035000-2 9.700000-19437 2151     \n"
    " 3.450000+3 3.038500+0 0.000000+0 2.659900-4 3.035000-2 9.700000-19437 2151     \n"
    " 3.550000+3 3.038000+0 0.000000+0 2.090000-4 3.035000-2 9.700000-19437 2151     \n"
    " 3.650000+3 3.037400+0 0.000000+0 3.420100-4 3.035000-2 9.700000-19437 2151     \n"
    " 3.750000+3 3.036800+0 0.000000+0 2.345000-4 3.035000-2 9.700000-19437 2151     \n"
    " 3.850000+3 3.036200+0 0.000000+0 3.083250-4 3.035000-2 9.700000-19437 2151     \n"
    " 3.950000+3 3.035700+0 0.000000+0 3.116900-4 3.035000-2 9.700000-19437 2151     \n"
    " 4.125000+3 3.034600+0 0.000000+0 2.302000-4 3.030000-2 9.660000-19437 2151     \n"
    " 4.375000+3 3.033200+0 0.000000+0 2.488800-4 3.030000-2 9.650000-19437 2151     \n"
    " 4.625000+3 3.031800+0 0.000000+0 2.300000-4 3.030000-2 9.650000-19437 2151     \n"
    " 4.875000+3 3.030400+0 0.000000+0 2.490000-4 3.030000-2 9.640000-19437 2151     \n"
    " 5.125000+3 3.029000+0 0.000000+0 2.574500-4 3.030000-2 9.640000-19437 2151     \n"
    " 5.375000+3 3.027600+0 0.000000+0 2.573500-4 3.030000-2 9.640000-19437 2151     \n"
    " 5.625000+3 3.026100+0 0.000000+0 2.572000-4 3.030000-2 9.630000-19437 2151     \n"
    " 5.875000+3 3.024700+0 0.000000+0 2.571000-4 3.030000-2 9.630000-19437 2151     \n"
    " 6.125000+3 3.023200+0 0.000000+0 2.570000-4 3.030000-2 9.620000-19437 2151     \n"
    " 6.375000+3 3.021800+0 0.000000+0 2.568500-4 3.030000-2 9.620000-19437 2151     \n"
    " 6.625000+3 3.020400+0 0.000000+0 2.452800-4 3.030000-2 9.610000-19437 2151     \n"
    " 6.875000+3 3.018900+0 0.000000+0 2.452800-4 3.030000-2 9.610000-19437 2151     \n"
    " 7.125000+3 3.017500+0 0.000000+0 2.532500-4 3.030000-2 9.600000-19437 2151     \n"
    " 7.375000+3 3.016100+0 0.000000+0 2.531000-4 3.030000-2 9.600000-19437 2151     \n"
    " 7.625000+3 3.014700+0 0.000000+0 2.529500-4 3.030000-2 9.600000-19437 2151     \n"
    " 7.875000+3 3.013300+0 0.000000+0 2.528500-4 3.030000-2 9.600000-19437 2151     \n"
    " 8.125000+3 3.011900+0 0.000000+0 2.560000-4 3.030000-2 9.590000-19437 2151     \n"
    " 8.375000+3 3.010500+0 0.000000+0 2.559000-4 3.030000-2 9.580000-19437 2151     \n"
    " 8.625000+3 3.009000+0 0.000000+0 2.557500-4 3.030000-2 9.580000-19437 2151     \n"
    " 8.875000+3 3.007600+0 0.000000+0 2.556500-4 3.030000-2 9.570000-19437 2151     \n"
    " 9.125000+3 3.006200+0 0.000000+0 2.472250-4 3.030000-2 9.570000-19437 2151     \n"
    " 9.375000+3 3.004800+0 0.000000+0 2.502650-4 3.030000-2 9.560000-19437 2151     \n"
    " 9.625000+3 3.003300+0 0.000000+0 2.455650-4 3.030000-2 9.560000-19437 2151     \n"
    " 9.875000+3 3.001900+0 0.000000+0 2.431600-4 3.030000-2 9.550000-19437 2151     \n"
    " 1.025000+4 2.999800+0 0.000000+0 2.500600-4 3.030000-2 9.550000-19437 2151     \n"
    " 1.075000+4 2.997000+0 0.000000+0 2.498200-4 3.030000-2 9.540000-19437 2151     \n"
    " 1.125000+4 2.994100+0 0.000000+0 2.495900-4 3.030000-2 9.530000-19437 2151     \n"
    " 1.175000+4 2.991300+0 0.000000+0 2.493550-4 3.030000-2 9.520000-19437 2151     \n"
    " 1.225000+4 2.988500+0 0.000000+0 2.491150-4 3.030000-2 9.510000-19437 2151     \n"
    " 1.275000+4 2.985700+0 0.000000+0 2.488850-4 3.030000-2 9.500000-19437 2151     \n"
    " 1.325000+4 2.982900+0 0.000000+0 2.486500-4 3.030000-2 9.490000-19437 2151     \n"
    " 1.375000+4 2.980100+0 0.000000+0 2.484200-4 3.030000-2 9.480000-19437 2151     \n"
    " 1.425000+4 2.977300+0 0.000000+0 2.481750-4 3.030000-2 9.480000-19437 2151     \n"
    " 1.475000+4 2.974400+0 0.000000+0 2.479450-4 3.030000-2 9.470000-19437 2151     \n"
    " 1.525000+4 2.971600+0 0.000000+0 2.477100-4 3.030000-2 9.460000-19437 2151     \n"
    " 1.575000+4 2.968800+0 0.000000+0 2.474800-4 3.030000-2 9.450000-19437 2151     \n"
    " 1.625000+4 2.966000+0 0.000000+0 2.472400-4 3.030000-2 9.440000-19437 2151     \n"
    " 1.675000+4 2.963200+0 0.000000+0 2.470100-4 3.030000-2 9.430000-19437 2151     \n"
    " 1.725000+4 2.960400+0 0.000000+0 2.467800-4 3.030000-2 9.420000-19437 2151     \n"
    " 1.775000+4 2.957600+0 0.000000+0 2.465450-4 3.030000-2 9.410000-19437 2151     \n"
    " 1.825000+4 2.954900+0 0.000000+0 2.463100-4 3.030000-2 9.410000-19437 2151     \n"
    " 1.875000+4 2.952100+0 0.000000+0 2.460800-4 3.030000-2 9.400000-19437 2151     \n"
    " 1.925000+4 2.949300+0 0.000000+0 2.458450-4 3.030000-2 9.390000-19437 2151     \n"
    " 1.975000+4 2.946500+0 0.000000+0 2.456150-4 3.030000-2 9.380000-19437 2151     \n"
    " 2.050000+4 2.942400+0 0.000000+0 2.452500-4 3.030000-2 9.360000-19437 2151     \n"
    " 2.150000+4 2.936900+0 0.000000+0 2.448000-4 3.030000-2 9.350000-19437 2151     \n"
    " 2.250000+4 2.931300+0 0.000000+0 2.443500-4 3.030000-2 9.330000-19437 2151     \n"
    " 2.350000+4 2.925800+0 0.000000+0 2.438500-4 3.030000-2 9.310000-19437 2151     \n"
    " 2.450000+4 2.920200+0 0.000000+0 2.434000-4 3.030000-2 9.290000-19437 2151     \n"
    " 2.550000+4 2.914700+0 0.000000+0 2.429500-4 3.030000-2 9.280000-19437 2151     \n"
    " 2.650000+4 2.904300+0 0.000000+0 2.425000-4 3.030000-2 9.260000-19437 2151     \n"
    " 2.750000+4 2.903800+0 0.000000+0 2.420500-4 3.030000-2 9.240000-19437 2151     \n"
    " 2.850000+4 2.898400+0 0.000000+0 2.416000-4 3.030000-2 9.220000-19437 2151     \n"
    " 2.950000+4 2.892900+0 0.000000+0 2.411500-4 3.030000-2 9.210000-19437 2151     \n"
    " 3.000000+4 2.890100+0 0.000000+0 2.408900-4 3.030000-2 9.180000-19437 2151     \n"
    " 2.000000+0 0.000000+0          2          0        432         719437 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 1.000000+0 0.000000+0 2.000000+09437 2151     \n"
    " 2.500000+3 1.915900+0 0.000000+0 3.381400-4 3.335000-2 6.130000-19437 2151     \n"
    " 2.550000+3 1.915600+0 0.000000+0 3.115100-4 3.335000-2 6.130000-19437 2151     \n"
    " 2.650000+3 1.915200+0 0.000000+0 3.919800-4 3.335000-2 6.120000-19437 2151     \n"
    " 2.750000+3 1.914800+0 0.000000+0 3.218000-4 3.335000-2 6.120000-19437 2151     \n"
    " 2.850000+3 1.914500+0 0.000000+0 3.341000-4 3.335000-2 6.120000-19437 2151     \n"
    " 2.950000+3 1.914100+0 0.000000+0 5.024300-4 3.335000-2 6.120000-19437 2151     \n"
    " 3.050000+3 1.913800+0 0.000000+0 2.874500-4 3.335000-2 6.120000-19437 2151     \n"
    " 3.150000+3 1.913400+0 0.000000+0 3.836600-4 3.335000-2 6.120000-19437 2151     \n"
    " 3.250000+3 1.913000+0 0.000000+0 3.797700-4 3.335000-2 6.120000-19437 2151     \n"
    " 3.350000+3 1.912700+0 0.000000+0 2.540400-4 3.335000-2 6.120000-19437 2151     \n"
    " 3.450000+3 1.912200+0 0.000000+0 3.347900-4 3.335000-2 6.120000-19437 2151     \n"
    " 3.550000+3 1.911900+0 0.000000+0 2.630600-4 3.335000-2 6.110000-19437 2151     \n"
    " 3.650000+3 1.911600+0 0.000000+0 4.304900-4 3.335000-2 6.110000-19437 2151     \n"
    " 3.750000+3 1.911200+0 0.000000+0 2.951700-4 3.335000-2 6.110000-19437 2151     \n"
    " 3.850000+3 1.910900+0 0.000000+0 3.881000-4 3.335000-2 6.110000-19437 2151     \n"
    " 3.950000+3 1.910500+0 0.000000+0 3.923200-4 3.335000-2 6.110000-19437 2151     \n"
    " 4.125000+3 1.909900+0 0.000000+0 2.899000-4 3.335000-2 6.080000-19437 2151     \n"
    " 4.375000+3 1.909000+0 0.000000+0 3.127200-4 3.335000-2 6.080000-19437 2151     \n"
    " 4.625000+3 1.908100+0 0.000000+0 2.893000-4 3.335000-2 6.070000-19437 2151     \n"
    " 4.875000+3 1.907200+0 0.000000+0 3.134000-4 3.335000-2 6.070000-19437 2151     \n"
    " 5.125000+3 1.906300+0 0.000000+0 3.241000-4 3.335000-2 6.070000-19437 2151     \n"
    " 5.375000+3 1.905400+0 0.000000+0 3.239000-4 3.335000-2 6.070000-19437 2151     \n"
    " 5.625000+3 1.904500+0 0.000000+0 3.238000-4 3.335000-2 6.070000-19437 2151     \n"
    " 5.875000+3 1.903600+0 0.000000+0 3.236000-4 3.335000-2 6.060000-19437 2151     \n"
    " 6.125000+3 1.902700+0 0.000000+0 3.235000-4 3.335000-2 6.060000-19437 2151     \n"
    " 6.375000+3 1.901800+0 0.000000+0 3.233000-4 3.335000-2 6.060000-19437 2151     \n"
    " 6.625000+3 1.900900+0 0.000000+0 3.087800-4 3.335000-2 6.050000-19437 2151     \n"
    " 6.875000+3 1.900000+0 0.000000+0 3.087800-4 3.335000-2 6.050000-19437 2151     \n"
    " 7.125000+3 1.899100+0 0.000000+0 3.187000-4 3.335000-2 6.040000-19437 2151     \n"
    " 7.375000+3 1.898200+0 0.000000+0 3.186000-4 3.335000-2 6.040000-19437 2151     \n"
    " 7.625000+3 1.897300+0 0.000000+0 3.184000-4 3.335000-2 6.040000-19437 2151     \n"
    " 7.875000+3 1.896400+0 0.000000+0 3.183000-4 3.335000-2 6.040000-19437 2151     \n"
    " 8.125000+3 1.895500+0 0.000000+0 3.222000-4 3.335000-2 6.030000-19437 2151     \n"
    " 8.375000+3 1.894600+0 0.000000+0 3.221000-4 3.335000-2 6.030000-19437 2151     \n"
    " 8.625000+3 1.893700+0 0.000000+0 3.219000-4 3.335000-2 6.020000-19437 2151     \n"
    " 8.875000+3 1.892800+0 0.000000+0 3.218000-4 3.335000-2 6.020000-19437 2151     \n"
    " 9.125000+3 1.891900+0 0.000000+0 3.181200-4 3.335000-2 6.020000-19437 2151     \n"
    " 9.375000+3 1.891000+0 0.000000+0 3.150100-4 3.335000-2 6.020000-19437 2151     \n"
    " 9.625000+3 1.890100+0 0.000000+0 3.091200-4 3.335000-2 6.010000-19437 2151     \n"
    " 9.875000+3 1.889200+0 0.000000+0 3.060200-4 3.335000-2 6.010000-19437 2151     \n"
    " 1.025000+4 1.887900+0 0.000000+0 3.147400-4 3.335000-2 6.010000-19437 2151     \n"
    " 1.075000+4 1.886100+0 0.000000+0 3.144500-4 3.335000-2 6.000000-19437 2151     \n"
    " 1.125000+4 1.884300+0 0.000000+0 3.141500-4 3.335000-2 6.000000-19437 2151     \n"
    " 1.175000+4 1.882500+0 0.000000+0 3.138500-4 3.335000-2 5.990000-19437 2151     \n"
    " 1.225000+4 1.880800+0 0.000000+0 3.135600-4 3.335000-2 5.990000-19437 2151     \n"
    " 1.275000+4 1.879000+0 0.000000+0 3.132600-4 3.335000-2 5.980000-19437 2151     \n"
    " 1.325000+4 1.877200+0 0.000000+0 3.129600-4 3.335000-2 5.970000-19437 2151     \n"
    " 1.375000+4 1.875400+0 0.000000+0 3.126600-4 3.335000-2 5.970000-19437 2151     \n"
    " 1.425000+4 1.873700+0 0.000000+0 3.123700-4 3.335000-2 5.960000-19437 2151     \n"
    " 1.475000+4 1.871900+0 0.000000+0 3.120700-4 3.335000-2 5.960000-19437 2151     \n"
    " 1.525000+4 1.870100+0 0.000000+0 3.117700-4 3.335000-2 5.950000-19437 2151     \n"
    " 1.575000+4 1.868400+0 0.000000+0 3.114800-4 3.335000-2 5.950000-19437 2151     \n"
    " 1.625000+4 1.866600+0 0.000000+0 3.111900-4 3.335000-2 5.940000-19437 2151     \n"
    " 1.675000+4 1.864800+0 0.000000+0 3.108900-4 3.335000-2 5.940000-19437 2151     \n"
    " 1.725000+4 1.863100+0 0.000000+0 3.106100-4 3.335000-2 5.930000-19437 2151     \n"
    " 1.775000+4 1.861300+0 0.000000+0 3.103100-4 3.335000-2 5.920000-19437 2151     \n"
    " 1.825000+4 1.859600+0 0.000000+0 3.100200-4 3.335000-2 5.920000-19437 2151     \n"
    " 1.875000+4 1.857800+0 0.000000+0 3.097300-4 3.335000-2 5.910000-19437 2151     \n"
    " 1.925000+4 1.856100+0 0.000000+0 3.094400-4 3.335000-2 5.910000-19437 2151     \n"
    " 1.975000+4 1.854300+0 0.000000+0 3.091400-4 3.335000-2 5.900000-19437 2151     \n"
    " 2.050000+4 1.851700+0 0.000000+0 3.087000-4 3.335000-2 5.890000-19437 2151     \n"
    " 2.150000+4 1.848200+0 0.000000+0 3.081000-4 3.335000-2 5.880000-19437 2151     \n"
    " 2.250000+4 1.844700+0 0.000000+0 3.075000-4 3.335000-2 5.870000-19437 2151     \n"
    " 2.350000+4 1.841300+0 0.000000+0 3.070000-4 3.335000-2 5.860000-19437 2151     \n"
    " 2.450000+4 1.837800+0 0.000000+0 3.064000-4 3.335000-2 5.850000-19437 2151     \n"
    " 2.550000+4 1.834300+0 0.000000+0 3.058000-4 3.335000-2 5.840000-19437 2151     \n"
    " 2.650000+4 1.830800+0 0.000000+0 3.052000-4 3.335000-2 5.830000-19437 2151     \n"
    " 2.750000+4 1.827400+0 0.000000+0 3.046000-4 3.335000-2 5.820000-19437 2151     \n"
    " 2.850000+4 1.823900+0 0.000000+0 3.041000-4 3.335000-2 5.800000-19437 2151     \n"
    " 2.950000+4 1.820500+0 0.000000+0 3.035000-4 3.335000-2 5.790000-19437 2151     \n"
    " 3.000000+4 1.818700+0 0.000000+0 3.032000-4 3.335000-2 5.770000-19437 2151     \n"
    "                                                                  9437 2  0     \n";
}

std::string Er167() {

  // Er167 ENDF/B-VIII.0 LRU=2 resonance evaluation

  return
    " 6.816700+4 1.654980+2          0          0          1          06840 2151     \n"
    " 6.816700+4 1.000000+0          0          0          1          06840 2151     \n"
    " 1.750000+3 1.000000+4          2          1          0          06840 2151     \n"
    " 3.500000+0 8.200000-1          0          0          2          06840 2151     \n"
    " 1.654980+2 0.000000+0          0          0         12          26840 2151     \n"
    " 9.142900+0 3.000000+0 1.000000+0 1.748200-3 1.126380-1 0.000000+06840 2151     \n"
    " 7.111100+0 4.000000+0 1.000000+0 1.321300-3 1.126380-1 0.000000+06840 2151     \n"
    " 1.654980+2 0.000000+0          1          0         24          46840 2151     \n"
    " 1.280000+1 2.000000+0 1.000000+0 1.236000-3 1.126380-1 0.000000+06840 2151     \n"
    " 9.142900+0 3.000000+0 2.000000+0 1.894300-3 1.126380-1 0.000000+06840 2151     \n"
    " 7.111100+0 4.000000+0 2.000000+0 1.406600-3 1.126380-1 0.000000+06840 2151     \n"
    " 5.818200+0 5.000000+0 1.000000+0 6.681800-4 1.126380-1 0.000000+06840 2151     \n"
    " 0.000000+0 0.000000+0          0          0          0          06840 2  0     \n";
}

std::string Au197() {

  // Au197 ENDF/B-VIII.0 LRU=2 resonance evaluation - energy dependent radius

  return
  std::string(
    " 7.919700+4 1.952740+2          0          0          1          07925 2151     \n"
    " 7.919700+4 1.000000+0          0          0          1          07925 2151     \n"
    " 2.000000+3 1.000000+5          2          2          1          07925 2151     \n"
    " 0.000000+0 0.000000+0          0          0          1        1007925 2151     \n"
    "        100          2                                            7925 2151     \n"
    " 2.000000+3 9.136729-1 2.250000+3 9.136377-1 2.500000+3 9.136022-17925 2151     \n"
    " 2.750000+3 9.135664-1 3.000000+3 9.135304-1 3.250000+3 9.134940-17925 2151     \n"
    " 3.500000+3 9.134574-1 3.750000+3 9.134205-1 4.000000+3 9.133834-17925 2151     \n"
    " 4.250000+3 9.133459-1 4.500000+3 9.133082-1 4.750000+3 9.132703-17925 2151     \n"
    " 5.000000+3 9.132320-1 5.250000+3 9.131935-1 5.500000+3 9.131547-17925 2151     \n"
    " 6.000000+3 9.130763-1 6.500000+3 9.129968-1 7.000000+3 9.129162-17925 2151     \n"
    " 7.500000+3 9.128345-1 8.000000+3 9.127518-1 8.500000+3 9.126680-17925 2151     \n"
    " 9.000000+3 9.125831-1 9.500000+3 9.124972-1 1.000000+4 9.124103-17925 2151     \n"
    " 1.100000+4 9.122332-1 1.200000+4 9.120519-1 1.300000+4 9.118665-17925 2151     \n"
    " 1.400000+4 9.116770-1 1.500000+4 9.114834-1 1.600000+4 9.112857-17925 2151     \n"
    " 1.700000+4 9.110840-1 1.800000+4 9.108784-1 1.900000+4 9.106688-17925 2151     \n"
    " 2.000000+4 9.104553-1 2.100000+4 9.102379-1 2.200000+4 9.100166-17925 2151     \n"
    " 2.300000+4 9.097915-1 2.400000+4 9.095627-1 2.500000+4 9.093301-17925 2151     \n"
    " 2.600000+4 9.090937-1 2.700000+4 9.088537-1 2.800000+4 9.086099-17925 2151     \n"
    " 2.900000+4 9.083626-1 3.000000+4 9.081117-1 3.100000+4 9.078572-17925 2151     \n"
    " 3.200000+4 9.075991-1 3.300000+4 9.073376-1 3.400000+4 9.070726-17925 2151     \n"
    " 3.500000+4 9.068041-1 3.600000+4 9.065322-1 3.700000+4 9.062570-17925 2151     \n"
    " 3.800000+4 9.059784-1 3.900000+4 9.056966-1 4.000000+4 9.054114-17925 2151     \n"
    " 4.100000+4 9.051230-1 4.200000+4 9.048314-1 4.300000+4 9.045365-17925 2151     \n"
    " 4.400000+4 9.042386-1 4.500000+4 9.039375-1 4.600000+4 9.036333-17925 2151     \n"
    " 4.700000+4 9.033261-1 4.800000+4 9.030159-1 5.000000+4 9.023865-17925 2151     \n"
    " 5.250000+4 9.015833-1 5.500000+4 9.007623-1 5.750000+4 8.999240-17925 2151     \n"
    " 6.000000+4 8.990689-1 6.250000+4 8.981975-1 6.500000+4 8.973104-17925 2151     \n"
    " 6.750000+4 8.964079-1 7.000000+4 8.954907-1 7.200000+4 8.947466-17925 2151     \n"
    " 7.400000+4 8.939936-1 7.500000+4 8.936138-1 7.775000+4 8.925586-17925 2151     \n"
    " 7.779600+4 8.925408-1 7.781000+4 8.925354-1 7.783000+4 8.925277-17925 2151     \n"
    " 7.786000+4 8.925161-1 7.789000+4 8.925045-1 7.793000+4 8.924890-17925 2151     \n"
    " 7.798000+4 8.924697-1 7.807000+4 8.924348-1 7.820000+4 8.923845-17925 2151     \n"
    " 7.835000+4 8.923263-1 7.862000+4 8.922216-1 7.890000+4 8.921128-17925 2151     \n"
    " 7.945000+4 8.918986-1 8.000000+4 8.916838-1 8.100000+4 8.912918-17925 2151     \n"
    " 8.200000+4 8.908978-1 8.300000+4 8.905019-1 8.500000+4 8.897046-17925 2151     \n"
    " 8.750000+4 8.886978-1 9.000000+4 8.876803-1 9.200000+4 8.868588-17925 2151     \n"
    " 9.350000+4 8.862385-1 9.500000+4 8.856147-1 9.750000+4 8.845678-17925 2151     \n"
    " 1.000000+5 8.835121-1                                            7925 2151     \n"
    " 1.500000+0 0.000000+0          1          0          3          07925 2151     \n"
    " 1.952740+2 0.000000+0          0          0          2          07925 2151     \n"
    " 1.000000+0 0.000000+0          2          0        606        1007925 2151     \n"
    " 0.000000+0 0.000000+0 1.000000+0 1.000000+0 0.000000+0 0.000000+07925 2151     \n"
    " 2.000000+3 3.848710+1 0.000000+0 6.450977-3 1.292688-1 0.000000+07925 2151     \n"
    " 2.250000+3 3.847356+1 0.000000+0 6.440182-3 1.292801-1 0.000000+07925 2151     \n"
    " 2.500000+3 3.846002+1 0.000000+0 6.429913-3 1.292914-1 0.000000+07925 2151     \n"
    " 2.750000+3 3.844648+1 0.000000+0 6.420093-3 1.293027-1 0.000000+07925 2151     \n"
    " 3.000000+3 3.843295+1 0.000000+0 6.410663-3 1.293141-1 0.000000+07925 2151     \n"
    " 3.250000+3 3.841943+1 0.000000+0 6.401574-3 1.293254-1 0.000000+07925 2151     \n"
    " 3.500000+3 3.840591+1 0.000000+0 6.392788-3 1.293367-1 0.000000+07925 2151     \n"
    " 3.750000+3 3.839239+1 0.000000+0 6.384274-3 1.293480-1 0.000000+07925 2151     \n"
    " 4.000000+3 3.837888+1 0.000000+0 6.376004-3 1.293593-1 0.000000+07925 2151     \n"
    " 4.250000+3 3.836538+1 0.000000+0 6.367956-3 1.293706-1 0.000000+07925 2151     \n"
    " 4.500000+3 3.835188+1 0.000000+0 6.360112-3 1.293820-1 0.000000+07925 2151     \n"
    " 4.750000+3 3.833838+1 0.000000+0 6.352453-3 1.293933-1 0.000000+07925 2151     \n"
    " 5.000000+3 3.832489+1 0.000000+0 6.344967-3 1.294046-1 0.000000+07925 2151     \n"
    " 5.250000+3 3.831141+1 0.000000+0 6.337640-3 1.294159-1 0.000000+07925 2151     \n"
    " 5.500000+3 3.829793+1 0.000000+0 6.330461-3 1.294272-1 0.000000+07925 2151     \n"
    " 6.000000+3 3.827098+1 0.000000+0 6.316509-3 1.294498-1 0.000000+07925 2151     \n"
    " 6.500000+3 3.824406+1 0.000000+0 6.303043-3 1.294725-1 0.000000+07925 2151     \n"
    " 7.000000+3 3.821715+1 0.000000+0 6.290010-3 1.294951-1 0.000000+07925 2151     \n"
    " 7.500000+3 3.819027+1 0.000000+0 6.277363-3 1.295177-1 0.000000+07925 2151     \n"
    " 8.000000+3 3.816340+1 0.000000+0 6.265066-3 1.295404-1 0.000000+07925 2151     \n"
    " 8.500000+3 3.813655+1 0.000000+0 6.253087-3 1.295630-1 0.000000+07925 2151     \n"
    " 9.000000+3 3.810973+1 0.000000+0 6.241398-3 1.295856-1 0.000000+07925 2151     \n"
    " 9.500000+3 3.808292+1 0.000000+0 6.229975-3 1.296083-1 0.000000+07925 2151     \n"
    " 1.000000+4 3.805613+1 0.000000+0 6.218799-3 1.296309-1 0.000000+07925 2151     \n"
    " 1.100000+4 3.800262+1 0.000000+0 6.197116-3 1.296762-1 0.000000+07925 2151     \n"
    " 1.200000+4 3.794918+1 0.000000+0 6.176230-3 1.297214-1 0.000000+07925 2151     \n"
    " 1.300000+4 3.789582+1 0.000000+0 6.156043-3 1.297667-1 0.000000+07925 2151     \n"
    " 1.400000+4 3.784254+1 0.000000+0 6.136480-3 1.298120-1 0.000000+07925 2151     \n"
    " 1.500000+4 3.778933+1 0.000000+0 6.117476-3 1.298572-1 0.000000+07925 2151     \n"
    " 1.600000+4 3.773621+1 0.000000+0 6.098978-3 1.299025-1 0.000000+07925 2151     \n"
    " 1.700000+4 3.768316+1 0.000000+0 6.080941-3 1.299478-1 0.000000+07925 2151     \n"
    " 1.800000+4 3.763019+1 0.000000+0 6.063325-3 1.299930-1 0.000000+07925 2151     \n"
    " 1.900000+4 3.757730+1 0.000000+0 6.046098-3 1.300383-1 0.000000+07925 2151     \n"
    " 2.000000+4 3.752449+1 0.000000+0 6.029231-3 1.300835-1 0.000000+07925 2151     \n"
    " 2.100000+4 3.747175+1 0.000000+0 6.012699-3 1.301288-1 0.000000+07925 2151     \n"
    " 2.200000+4 3.741909+1 0.000000+0 5.996478-3 1.301741-1 0.000000+07925 2151     \n"
    " 2.300000+4 3.736651+1 0.000000+0 5.980549-3 1.302193-1 0.000000+07925 2151     \n"
    " 2.400000+4 3.731400+1 0.000000+0 5.964894-3 1.302646-1 0.000000+07925 2151     \n"
    " 2.500000+4 3.726157+1 0.000000+0 5.949497-3 1.303099-1 0.000000+07925 2151     \n"
    " 2.600000+4 3.720922+1 0.000000+0 5.934343-3 1.303551-1 0.000000+07925 2151     \n"
    " 2.700000+4 3.715694+1 0.000000+0 5.919421-3 1.304004-1 0.000000+07925 2151     \n"
    " 2.800000+4 3.710475+1 0.000000+0 5.904717-3 1.304457-1 0.000000+07925 2151     \n"
    " 2.900000+4 3.705262+1 0.000000+0 5.890221-3 1.304909-1 0.000000+07925 2151     \n"
    " 3.000000+4 3.700058+1 0.000000+0 5.875923-3 1.305362-1 0.000000+07925 2151     \n"
    " 3.100000+4 3.694861+1 0.000000+0 5.861814-3 1.305814-1 0.000000+07925 2151     \n"
    " 3.200000+4 3.689671+1 0.000000+0 5.847885-3 1.306267-1 0.000000+07925 2151     \n"
    " 3.300000+4 3.684490+1 0.000000+0 5.834130-3 1.306720-1 0.000000+07925 2151     \n"
    " 3.400000+4 3.679315+1 0.000000+0 5.820540-3 1.307172-1 0.000000+07925 2151     \n"
    " 3.500000+4 3.674149+1 0.000000+0 5.807109-3 1.307625-1 0.000000+07925 2151     \n"
    " 3.600000+4 3.668990+1 0.000000+0 5.793832-3 1.308078-1 0.000000+07925 2151     \n"
    " 3.700000+4 3.663838+1 0.000000+0 5.780701-3 1.308530-1 0.000000+07925 2151     \n"
    " 3.800000+4 3.658694+1 0.000000+0 5.767712-3 1.308983-1 0.000000+07925 2151     \n"
    " 3.900000+4 3.653558+1 0.000000+0 5.754860-3 1.309436-1 0.000000+07925 2151     \n"
    " 4.000000+4 3.648429+1 0.000000+0 5.742141-3 1.309888-1 0.000000+07925 2151     \n"
    " 4.100000+4 3.643308+1 0.000000+0 5.729548-3 1.310341-1 0.000000+07925 2151     \n"
    " 4.200000+4 3.638194+1 0.000000+0 5.717080-3 1.310794-1 0.000000+07925 2151     \n"
    " 4.300000+4 3.633087+1 0.000000+0 5.704731-3 1.311246-1 0.000000+07925 2151     \n"
    " 4.400000+4 3.627988+1 0.000000+0 5.692498-3 1.311699-1 0.000000+07925 2151     \n"
    " 4.500000+4 3.622897+1 0.000000+0 5.680377-3 1.312151-1 0.000000+07925 2151     \n"
    " 4.600000+4 3.617813+1 0.000000+0 5.668365-3 1.312604-1 0.000000+07925 2151     \n"
    " 4.700000+4 3.612736+1 0.000000+0 5.656460-3 1.313057-1 0.000000+07925 2151     \n"
    " 4.800000+4 3.607667+1 0.000000+0 5.644657-3 1.313509-1 0.000000+07925 2151     \n"
    " 5.000000+4 3.597551+1 0.000000+0 5.621350-3 1.314415-1 0.000000+07925 2151     \n"
    " 5.250000+4 3.584947+1 0.000000+0 5.592748-3 1.315546-1 0.000000+07925 2151     \n"
    " 5.500000+4 3.572390+1 0.000000+0 5.564704-3 1.316678-1 0.000000+07925 2151     \n"
    " 5.750000+4 3.559878+1 0.000000+0 5.537184-3 1.317809-1 0.000000+07925 2151     \n"
    " 6.000000+4 3.547412+1 0.000000+0 5.510159-3 1.318941-1 0.000000+07925 2151     \n"
    " 6.250000+4 3.534992+1 0.000000+0 5.483602-3 1.320073-1 0.000000+07925 2151     \n"
    " 6.500000+4 3.522616+1 0.000000+0 5.457491-3 1.321204-1 0.000000+07925 2151     \n"
    " 6.750000+4 3.510286+1 0.000000+0 5.431802-3 1.322336-1 0.000000+07925 2151     \n"
    " 7.000000+4 3.498001+1 0.000000+0 5.406517-3 1.323467-1 0.000000+07925 2151     \n"
    " 7.200000+4 3.488205+1 0.000000+0 5.386566-3 1.324373-1 0.000000+07925 2151     \n"
    " 7.400000+4 3.478438+1 0.000000+0 5.366853-3 1.325278-1 0.000000+07925 2151     \n"
    " 7.500000+4 3.473565+1 0.000000+0 5.357083-3 1.325731-1 0.000000+07925 2151     \n"
    " 7.775000+4 3.460201+1 0.000000+0 5.330503-3 1.326975-1 0.000000+07925 2151     \n"
    " 7.779600+4 3.459978+1 0.000000+0 5.330062-3 1.326996-1 0.000000+07925 2151     \n"
    " 7.781000+4 3.459910+1 2.215873-2 5.329928-3 1.327003-1 0.000000+07925 2151     \n"
    " 7.783000+4 3.459813+1 3.449246-2 5.329736-3 1.327012-1 0.000000+07925 2151     \n"
    " 7.786000+4 3.459668+1 4.726681-2 5.329448-3 1.327025-1 0.000000+07925 2151     \n"
    " 7.789000+4 3.459523+1 5.722986-2 5.329161-3 1.327039-1 0.000000+07925 2151     \n"
    " 7.793000+4 3.459329+1 6.825848-2 5.328778-3 1.327057-1 0.000000+07925 2151     \n"
    " 7.798000+4 3.459086+1 7.989684-2 5.328299-3 1.327080-1 0.000000+07925 2151     \n"
    " 7.807000+4 3.458650+1 9.733388-2 5.327437-3 1.327120-1 0.000000+07925 2151     \n"
    " 7.820000+4 3.458020+1 1.179512-1 5.326193-3 1.327179-1 0.000000+07925 2151     \n"
    " 7.835000+4 3.457293+1 1.378502-1 5.324758-3 1.327247-1 0.000000+07925 2151     \n"
    " 7.862000+4 3.455985+1 1.676133-1 5.322180-3 1.327369-1 0.000000+07925 2151     \n"
    " 7.890000+4 3.454629+1 1.934947-1 5.319509-3 1.327496-1 0.000000+07925 2151     \n"
    " 7.945000+4 3.451967+1 2.357812-1 5.314276-3 1.327745-1 0.000000+07925 2151     \n"
    " 8.000000+4 3.449307+1 2.711293-1 5.309059-3 1.327994-1 0.000000+07925 2151     \n"
    " 8.100000+4 3.444477+1 3.249295-1 5.299614-3 1.328447-1 0.000000+07925 2151     \n"
    " 8.200000+4 3.439654+1 3.702358-1 5.290219-3 1.328899-1 0.000000+07925 2151     \n"
    " 8.300000+4 3.434837+1 4.099564-1 5.280875-3 1.329352-1 0.000000+07925 2151     \n"
    " 8.500000+4 3.425226+1 4.782272-1 5.262334-3 1.330257-1 0.000000+07925 2151     \n"
    " 8.750000+4 3.413251+1 5.497995-1 5.239428-3 1.331389-1 0.000000+07925 2151     \n"
    " 9.000000+4 3.401320+1 6.112949-1 5.216811-3 1.332520-1 0.000000+07925 2151     \n"
    " 9.200000+4 3.391806+1 6.552733-1 5.198919-3 1.333426-1 0.000000+07925 2151     \n"
    " 9.350000+4 3.384689+1 6.858482-1 5.185614-3 1.334105-1 0.000000+07925 2151     \n"
    " 9.500000+4 3.377587+1 7.146777-1 5.172404-3 1.334783-1 0.000000+07925 2151     \n"
    " 9.750000+4 3.365786+1 7.594286-1 5.150596-3 1.335915-1 0.000000+07925 2151     \n"
    " 1.000000+5 3.354027+1 8.007053-1 5.129039-3 1.337047-1 0.000000+07925 2151     \n"
    " 2.000000+0 0.000000+0          2          0        606        1007925 2151     \n"
    " 0.000000+0 0.000000+0 2.000000+0 1.000000+0 0.000000+0 0.000000+07925 2151     \n"
    " 2.000000+3 2.393741+1 0.000000+0 4.012245-3 1.292688-1 0.000000+07925 2151     \n"
    " 2.250000+3 2.392897+1 0.000000+0 4.005529-3 1.292801-1 0.000000+07925 2151     \n"
    " 2.500000+3 2.392054+1 0.000000+0 3.999139-3 1.292914-1 0.000000+07925 2151     \n"
    " 2.750000+3 2.391210+1 0.000000+0 3.993029-3 1.293027-1 0.000000+07925 2151     \n"
    " 3.000000+3 2.390367+1 0.000000+0 3.987161-3 1.293141-1 0.000000+07925 2151     \n"
    " 3.250000+3 2.389524+1 0.000000+0 3.981505-3 1.293254-1 0.000000+07925 2151     \n"
    " 3.500000+3 2.388682+1 0.000000+0 3.976038-3 1.293367-1 0.000000+07925 2151     \n"
    " 3.750000+3 2.387839+1 0.000000+0 3.970740-3 1.293480-1 0.000000+07925 2151     \n"
    " 4.000000+3 2.386998+1 0.000000+0 3.965594-3 1.293593-1 0.000000+07925 2151     \n"
    " 4.250000+3 2.386156+1 0.000000+0 3.960586-3 1.293706-1 0.000000+07925 2151     \n"
    " 4.500000+3 2.385315+1 0.000000+0 3.955704-3 1.293820-1 0.000000+07925 2151     \n"
    " 4.750000+3 2.384474+1 0.000000+0 3.950938-3 1.293933-1 0.000000+07925 2151     \n"
    " 5.000000+3 2.383633+1 0.000000+0 3.946279-3 1.294046-1 0.000000+07925 2151     \n"
    " 5.250000+3 2.382793+1 0.000000+0 3.941719-3 1.294159-1 0.000000+07925 2151     \n"
    " 5.500000+3 2.381953+1 0.000000+0 3.937252-3 1.294272-1 0.000000+07925 2151     \n"
    " 6.000000+3 2.380274+1 0.000000+0 3.928569-3 1.294498-1 0.000000+07925 2151     \n"
    " 6.500000+3 2.378596+1 0.000000+0 3.920189-3 1.294725-1 0.000000+07925 2151     \n"
    " 7.000000+3 2.376919+1 0.000000+0 3.912077-3 1.294951-1 0.000000+07925 2151     \n"
    " 7.500000+3 2.375244+1 0.000000+0 3.904206-3 1.295177-1 0.000000+07925 2151     \n"
    " 8.000000+3 2.373569+1 0.000000+0 3.896553-3 1.295404-1 0.000000+07925 2151     \n"
    " 8.500000+3 2.371896+1 0.000000+0 3.889097-3 1.295630-1 0.000000+07925 2151     \n"
    " 9.000000+3 2.370225+1 0.000000+0 3.881821-3 1.295856-1 0.000000+07925 2151     \n"
    " 9.500000+3 2.368554+1 0.000000+0 3.874712-3 1.296083-1 0.000000+07925 2151     \n"
    " 1.000000+4 2.366885+1 0.000000+0 3.867756-3 1.296309-1 0.000000+07925 2151     \n"
    " 1.100000+4 2.363550+1 0.000000+0 3.854260-3 1.296762-1 0.000000+07925 2151     \n"
    " 1.200000+4 2.360220+1 0.000000+0 3.841259-3 1.297214-1 0.000000+07925 2151     \n"
    " 1.300000+4 2.356895+1 0.000000+0 3.828694-3 1.297667-1 0.000000+07925 2151     \n"
    " 1.400000+4 2.353575+1 0.000000+0 3.816516-3 1.298120-1 0.000000+07925 2151     \n"
    " 1.500000+4 2.350259+1 0.000000+0 3.804686-3 1.298572-1 0.000000+07925 2151     \n"
    " 1.600000+4 2.346949+1 0.000000+0 3.793171-3 1.299025-1 0.000000+07925 2151     \n"
    " 1.700000+4 2.343643+1 0.000000+0 3.781943-3 1.299478-1 0.000000+07925 2151     \n"
    " 1.800000+4 2.340343+1 0.000000+0 3.770977-3 1.299930-1 0.000000+07925 2151     \n"
    " 1.900000+4 2.337047+1 0.000000+0 3.760253-3 1.300383-1 0.000000+07925 2151     \n"
    " 2.000000+4 2.333756+1 0.000000+0 3.749752-3 1.300835-1 0.000000+07925 2151     \n"
    " 2.100000+4 2.330469+1 0.000000+0 3.739460-3 1.301288-1 0.000000+07925 2151     \n"
    " 2.200000+4 2.327188+1 0.000000+0 3.729362-3 1.301741-1 0.000000+07925 2151     \n"
    " 2.300000+4 2.323911+1 0.000000+0 3.719445-3 1.302193-1 0.000000+07925 2151     \n"
    " 2.400000+4 2.320640+1 0.000000+0 3.709698-3 1.302646-1 0.000000+07925 2151     \n"
    " 2.500000+4 2.317373+1 0.000000+0 3.700113-3 1.303099-1 0.000000+07925 2151     \n"
    " 2.600000+4 2.314110+1 0.000000+0 3.690678-3 1.303551-1 0.000000+07925 2151     \n"
    " 2.700000+4 2.310853+1 0.000000+0 3.681388-3 1.304004-1 0.000000+07925 2151     \n"
    " 2.800000+4 2.307600+1 0.000000+0 3.672233-3 1.304457-1 0.000000+07925 2151     \n"
    " 2.900000+4 2.304353+1 0.000000+0 3.663208-3 1.304909-1 0.000000+07925 2151     \n"
    " 3.000000+4 2.301110+1 0.000000+0 3.654306-3 1.305362-1 0.000000+07925 2151     \n"
    " 3.100000+4 2.297871+1 0.000000+0 3.645521-3 1.305814-1 0.000000+07925 2151     \n"
    " 3.200000+4 2.294638+1 0.000000+0 3.636849-3 1.306267-1 0.000000+07925 2151     \n"
    " 3.300000+4 2.291409+1 0.000000+0 3.628285-3 1.306720-1 0.000000+07925 2151     \n"
    " 3.400000+4 2.288185+1 0.000000+0 3.619823-3 1.307172-1 0.000000+07925 2151     \n"
    " 3.500000+4 2.284965+1 0.000000+0 3.611461-3 1.307625-1 0.000000+07925 2151     \n"
    " 3.600000+4 2.281751+1 0.000000+0 3.603193-3 1.308078-1 0.000000+07925 2151     \n"
    " 3.700000+4 2.278541+1 0.000000+0 3.595018-3 1.308530-1 0.000000+07925 2151     \n"
    " 3.800000+4 2.275336+1 0.000000+0 3.586930-3 1.308983-1 0.000000+07925 2151     \n"
    " 3.900000+4 2.272135+1 0.000000+0 3.578928-3 1.309436-1 0.000000+07925 2151     \n"
    " 4.000000+4 2.268939+1 0.000000+0 3.571008-3 1.309888-1 0.000000+07925 2151     \n"
    " 4.100000+4 2.265748+1 0.000000+0 3.563167-3 1.310341-1 0.000000+07925 2151     \n"
    " 4.200000+4 2.262562+1 0.000000+0 3.555403-3 1.310794-1 0.000000+07925 2151     \n"
    " 4.300000+4 2.259380+1 0.000000+0 3.547714-3 1.311246-1 0.000000+07925 2151     \n"
    " 4.400000+4 2.256203+1 0.000000+0 3.540097-3 1.311699-1 0.000000+07925 2151     \n"
    " 4.500000+4 2.253030+1 0.000000+0 3.532549-3 1.312151-1 0.000000+07925 2151     \n"
    " 4.600000+4 2.249863+1 0.000000+0 3.525070-3 1.312604-1 0.000000+07925 2151     \n"
    " 4.700000+4 2.246699+1 0.000000+0 3.517657-3 1.313057-1 0.000000+07925 2151     \n"
    " 4.800000+4 2.243541+1 0.000000+0 3.510307-3 1.313509-1 0.000000+07925 2151     \n"
    " 5.000000+4 2.237238+1 0.000000+0 3.495794-3 1.314415-1 0.000000+07925 2151     \n"
    " 5.250000+4 2.229385+1 0.000000+0 3.477984-3 1.315546-1 0.000000+07925 2151     \n"
    " 5.500000+4 2.221561+1 0.000000+0 3.460520-3 1.316678-1 0.000000+07925 2151     \n"
    " 5.750000+4 2.213765+1 0.000000+0 3.443383-3 1.317809-1 0.000000+07925 2151     \n"
    " 6.000000+4 2.205998+1 0.000000+0 3.426554-3 1.318941-1 0.000000+07925 2151     \n"
    " 6.250000+4 2.198259+1 0.000000+0 3.410016-3 1.320073-1 0.000000+07925 2151     \n"
    " 6.500000+4 2.190549+1 0.000000+0 3.393756-3 1.321204-1 0.000000+07925 2151     \n"
    " 6.750000+4 2.182866+1 0.000000+0 3.377759-3 1.322336-1 0.000000+07925 2151     \n"
    " 7.000000+4 2.175212+1 0.000000+0 3.362012-3 1.323467-1 0.000000+07925 2151     \n"
    " 7.200000+4 2.169109+1 0.000000+0 3.349588-3 1.324373-1 0.000000+07925 2151     \n"
    " 7.400000+4 2.163024+1 0.000000+0 3.337311-3 1.325278-1 0.000000+07925 2151     \n"
    " 7.500000+4 2.159988+1 0.000000+0 3.331227-3 1.325731-1 0.000000+07925 2151     \n"
    " 7.775000+4 2.151662+1 0.000000+0 3.314674-3 1.326975-1 0.000000+07925 2151     \n"
    " 7.779600+4 2.151523+1 0.000000+0 3.314399-3 1.326996-1 0.000000+07925 2151     \n"
    " 7.781000+4 2.151481+1 0.000000+0 3.314316-3 1.327003-1 0.000000+07925 2151     \n"
    " 7.783000+4 2.151420+1 0.000000+0 3.314196-3 1.327012-1 0.000000+07925 2151     \n"
    " 7.786000+4 2.151329+1 0.000000+0 3.314017-3 1.327025-1 0.000000+07925 2151     \n"
    " 7.789000+4 2.151239+1 1.315010-9 3.313838-3 1.327039-1 0.000000+07925 2151     \n"
    " 7.793000+4 2.151118+1 3.190269-9 3.313599-3 1.327057-1 0.000000+07925 2151     \n"
    " 7.798000+4 2.150967+1 7.047819-9 3.313301-3 1.327080-1 0.000000+07925 2151     \n"
    " 7.807000+4 2.150695+1 1.906722-8 3.312764-3 1.327120-1 0.000000+07925 2151     \n"
    " 7.820000+4 2.150303+1 5.031771-8 3.311990-3 1.327179-1 0.000000+07925 2151     \n"
    " 7.835000+4 2.149850+1 1.107583-7 3.311096-3 1.327247-1 0.000000+07925 2151     \n"
    " 7.862000+4 2.149035+1 2.986221-7 3.309490-3 1.327369-1 0.000000+07925 2151     \n"
    " 7.890000+4 2.148190+1 6.200368-7 3.307828-3 1.327496-1 0.000000+07925 2151     \n"
    " 7.945000+4 2.146532+1 1.701085-6 3.304569-3 1.327745-1 0.000000+07925 2151     \n"
    " 8.000000+4 2.144874+1 3.481835-6 3.301320-3 1.327994-1 0.000000+07925 2151     \n"
    " 8.100000+4 2.141865+1 8.849198-6 3.295437-3 1.328447-1 0.000000+07925 2151     \n"
    " 8.200000+4 2.138860+1 1.740690-5 3.289587-3 1.328899-1 0.000000+07925 2151     \n"
    " 8.300000+4 2.135859+1 2.960047-5 3.283767-3 1.329352-1 0.000000+07925 2151     \n"
    " 8.500000+4 2.129871+1 6.639977-5 3.272221-3 1.330257-1 0.000000+07925 2151     \n"
    " 8.750000+4 2.122411+1 1.389384-4 3.257955-3 1.331389-1 0.000000+07925 2151     \n"
    " 9.000000+4 2.114977+1 2.448526-4 3.243870-3 1.332520-1 0.000000+07925 2151     \n"
    " 9.200000+4 2.109050+1 3.559854-4 3.232727-3 1.333426-1 0.000000+07925 2151     \n"
    " 9.350000+4 2.104616+1 4.557688-4 3.224441-3 1.334105-1 0.000000+07925 2151     \n"
    " 9.500000+4 2.100192+1 5.702988-4 3.216214-3 1.334783-1 0.000000+07925 2151     \n"
    " 9.750000+4 2.092840+1 7.954091-4 3.202632-3 1.335915-1 0.000000+07925 2151     \n"
    " 1.000000+5 2.085515+1 1.065243-3 3.189207-3 1.337047-1 0.000000+07925 2151     \n"
    " 1.952740+2 0.000000+0          1          0          4          07925 2151     \n"
    " 0.000000+0 0.000000+0          2          0        606        1007925 2151     \n"
    " 0.000000+0 0.000000+0 1.000000+0 1.000000+0 0.000000+0 0.000000+07925 2151     \n"
    " 2.000000+3 1.134047+2 0.000000+0 6.028556-3 4.205961-2 0.000000+07925 2151     \n"
    " 2.250000+3 1.133648+2 0.000000+0 6.028166-3 4.206329-2 0.000000+07925 2151     \n"
    " 2.500000+3 1.133250+2 0.000000+0 6.027772-3 4.206698-2 0.000000+07925 2151     \n"
    " 2.750000+3 1.132851+2 0.000000+0 6.027374-3 4.207066-2 0.000000+07925 2151     \n"
    " 3.000000+3 1.132453+2 0.000000+0 6.026972-3 4.207434-2 0.000000+07925 2151     \n"
    " 3.250000+3 1.132055+2 0.000000+0 6.026566-3 4.207802-2 0.000000+07925 2151     \n"
    " 3.500000+3 1.131657+2 0.000000+0 6.026156-3 4.208170-2 0.000000+07925 2151     \n"
    " 3.750000+3 1.131259+2 0.000000+0 6.025742-3 4.208538-2 0.000000+07925 2151     \n"
    " 4.000000+3 1.130861+2 0.000000+0 6.025325-3 4.208907-2 0.000000+07925 2151     \n"
    " 4.250000+3 1.130464+2 0.000000+0 6.024904-3 4.209275-2 0.000000+07925 2151     \n"
    " 4.500000+3 1.130067+2 0.000000+0 6.024479-3 4.209643-2 0.000000+07925 2151     \n"
    " 4.750000+3 1.129669+2 0.000000+0 6.024051-3 4.210011-2 0.000000+07925 2151     \n"
    " 5.000000+3 1.129272+2 0.000000+0 6.023620-3 4.210379-2 0.000000+07925 2151     \n"
    " 5.250000+3 1.128875+2 0.000000+0 6.023186-3 4.210748-2 0.000000+07925 2151     \n"
    " 5.500000+3 1.128478+2 0.000000+0 6.022748-3 4.211116-2 0.000000+07925 2151     \n"
    " 6.000000+3 1.127685+2 0.000000+0 6.021862-3 4.211852-2 0.000000+07925 2151     \n"
    " 6.500000+3 1.126893+2 0.000000+0 6.020965-3 4.212589-2 0.000000+07925 2151     \n"
    " 7.000000+3 1.126101+2 0.000000+0 6.020055-3 4.213325-2 0.000000+07925 2151     \n"
    " 7.500000+3 1.125309+2 0.000000+0 6.019133-3 4.214061-2 0.000000+07925 2151     \n"
    " 8.000000+3 1.124518+2 0.000000+0 6.018199-3 4.214798-2 0.000000+07925 2151     \n"
    " 8.500000+3 1.123728+2 0.000000+0 6.017254-3 4.215534-2 0.000000+07925 2151     \n"
    " 9.000000+3 1.122938+2 0.000000+0 6.016298-3 4.216270-2 0.000000+07925 2151     \n"
    " 9.500000+3 1.122149+2 0.000000+0 6.015330-3 4.217007-2 0.000000+07925 2151     \n"
    " 1.000000+4 1.121361+2 0.000000+0 6.014352-3 4.217743-2 0.000000+07925 2151     \n"
    " 1.100000+4 1.119785+2 0.000000+0 6.012364-3 4.219216-2 0.000000+07925 2151     \n"
    " 1.200000+4 1.118212+2 0.000000+0 6.010334-3 4.220689-2 0.000000+07925 2151     \n"
    " 1.300000+4 1.116641+2 0.000000+0 6.008264-3 4.222161-2 0.000000+07925 2151     \n"
    " 1.400000+4 1.115073+2 0.000000+0 6.006154-3 4.223634-2 0.000000+07925 2151     \n"
    " 1.500000+4 1.113507+2 0.000000+0 6.004007-3 4.225107-2 0.000000+07925 2151     \n"
    " 1.600000+4 1.111943+2 0.000000+0 6.001821-3 4.226580-2 0.000000+07925 2151     \n"
    " 1.700000+4 1.110381+2 0.000000+0 5.999598-3 4.228052-2 0.000000+07925 2151     \n"
    " 1.800000+4 1.108822+2 0.000000+0 5.997340-3 4.229525-2 0.000000+07925 2151     \n"
    " 1.900000+4 1.107265+2 0.000000+0 5.995045-3 4.230998-2 0.000000+07925 2151     \n"
    " 2.000000+4 1.105710+2 0.000000+0 5.992716-3 4.232470-2 0.000000+07925 2151     \n"
    " 2.100000+4 1.104158+2 0.000000+0 5.990352-3 4.233943-2 0.000000+07925 2151     \n"
    " 2.200000+4 1.102608+2 0.000000+0 5.987955-3 4.235416-2 0.000000+07925 2151     \n"
    " 2.300000+4 1.101060+2 0.000000+0 5.985524-3 4.236889-2 0.000000+07925 2151     \n"
    " 2.400000+4 1.099514+2 0.000000+0 5.983060-3 4.238361-2 0.000000+07925 2151     \n"
    " 2.500000+4 1.097971+2 0.000000+0 5.980564-3 4.239834-2 0.000000+07925 2151     \n"
    " 2.600000+4 1.096429+2 0.000000+0 5.978036-3 4.241307-2 0.000000+07925 2151     \n"
    " 2.700000+4 1.094891+2 0.000000+0 5.975476-3 4.242780-2 0.000000+07925 2151     \n"
    " 2.800000+4 1.093354+2 0.000000+0 5.972885-3 4.244252-2 0.000000+07925 2151     \n"
    " 2.900000+4 1.091820+2 0.000000+0 5.970264-3 4.245725-2 0.000000+07925 2151     \n"
    " 3.000000+4 1.090287+2 0.000000+0 5.967612-3 4.247198-2 0.000000+07925 2151     \n"
    " 3.100000+4 1.088757+2 0.000000+0 5.964930-3 4.248671-2 0.000000+07925 2151     \n"
    " 3.200000+4 1.087230+2 0.000000+0 5.962218-3 4.250143-2 0.000000+07925 2151     \n"
    " 3.300000+4 1.085704+2 0.000000+0 5.959477-3 4.251616-2 0.000000+07925 2151     \n"
    " 3.400000+4 1.084181+2 0.000000+0 5.956707-3 4.253089-2 0.000000+07925 2151     \n"
    " 3.500000+4 1.082660+2 0.000000+0 5.953909-3 4.254562-2 0.000000+07925 2151     \n"
    " 3.600000+4 1.081141+2 0.000000+0 5.951082-3 4.256034-2 0.000000+07925 2151     \n"
    " 3.700000+4 1.079625+2 0.000000+0 5.948227-3 4.257507-2 0.000000+07925 2151     \n"
    " 3.800000+4 1.078111+2 0.000000+0 5.945344-3 4.258980-2 0.000000+07925 2151     \n"
    " 3.900000+4 1.076599+2 0.000000+0 5.942434-3 4.260452-2 0.000000+07925 2151     \n"
    " 4.000000+4 1.075089+2 0.000000+0 5.939496-3 4.261925-2 0.000000+07925 2151     \n"
    " 4.100000+4 1.073581+2 0.000000+0 5.936532-3 4.263398-2 0.000000+07925 2151     \n"
    " 4.200000+4 1.072076+2 0.000000+0 5.933541-3 4.264871-2 0.000000+07925 2151     \n"
    " 4.300000+4 1.070572+2 0.000000+0 5.930524-3 4.266343-2 0.000000+07925 2151     \n"
    " 4.400000+4 1.069071+2 0.000000+0 5.927480-3 4.267816-2 0.000000+07925 2151     \n"
    " 4.500000+4 1.067572+2 0.000000+0 5.924411-3 4.269289-2 0.000000+07925 2151     \n"
    " 4.600000+4 1.066076+2 0.000000+0 5.921316-3 4.270762-2 0.000000+07925 2151     \n"
    " 4.700000+4 1.064581+2 0.000000+0 5.918196-3 4.272234-2 0.000000+07925 2151     \n"
    " 4.800000+4 1.063089+2 0.000000+0 5.915050-3 4.273707-2 0.000000+07925 2151     \n"
    " 5.000000+4 1.060111+2 0.000000+0 5.908685-3 4.276653-2 0.000000+07925 2151     \n"
    " 5.250000+4 1.056400+2 0.000000+0 5.900591-3 4.280334-2 0.000000+07925 2151     \n"
    " 5.500000+4 1.052704+2 0.000000+0 5.892348-3 4.284016-2 0.000000+07925 2151     \n"
    " 5.750000+4 1.049020+2 0.000000+0 5.883959-3 4.287698-2 0.000000+07925 2151     \n"
    " 6.000000+4 1.045350+2 0.000000+0 5.875428-3 4.291380-2 0.000000+07925 2151     \n"
    " 6.250000+4 1.041694+2 0.000000+0 5.866756-3 4.295062-2 0.000000+07925 2151     \n"
    " 6.500000+4 1.038050+2 0.000000+0 5.857949-3 4.298744-2 0.000000+07925 2151     \n"
    " 6.750000+4 1.034420+2 0.000000+0 5.849007-3 4.302425-2 0.000000+07925 2151     \n"
    " 7.000000+4 1.030804+2 0.000000+0 5.839935-3 4.306107-2 0.000000+07925 2151     \n"
    " 7.200000+4 1.027920+2 0.000000+0 5.832585-3 4.309053-2 0.000000+07925 2151     \n"
    " 7.400000+4 1.025044+2 0.000000+0 5.825155-3 4.311998-2 0.000000+07925 2151     \n"
    " 7.500000+4 1.023610+2 0.000000+0 5.821410-3 4.313471-2 0.000000+07925 2151     \n"
    " 7.775000+4 1.019675+2 0.000000+0 5.811011-3 4.317521-2 0.000000+07925 2151     \n"
    " 7.779600+4 1.019610+2 0.000000+0 5.810835-3 4.317589-2 0.000000+07925 2151     \n"
    " 7.781000+4 1.019590+2 8.565862-7 5.810782-3 4.317609-2 0.000000+07925 2151     \n"
    " 7.783000+4 1.019561+2 3.241673-6 5.810706-3 4.317639-2 0.000000+07925 2151     \n"
    " 7.786000+4 1.019518+2 8.371025-6 5.810592-3 4.317683-2 0.000000+07925 2151     \n"
    " 7.789000+4 1.019476+2 1.489904-5 5.810477-3 4.317727-2 0.000000+07925 2151     \n"
    " 7.793000+4 1.019418+2 2.535522-5 5.810325-3 4.317786-2 0.000000+07925 2151     \n"
    " 7.798000+4 1.019347+2 4.079120-5 5.810134-3 4.317860-2 0.000000+07925 2151     \n"
    " 7.807000+4 1.019219+2 7.410352-5 5.809791-3 4.317992-2 0.000000+07925 2151     \n"
    " 7.820000+4 1.019033+2 1.326178-4 5.809295-3 4.318184-2 0.000000+07925 2151     \n"
    " 7.835000+4 1.018819+2 2.128545-4 5.808722-3 4.318405-2 0.000000+07925 2151     \n"
    " 7.862000+4 1.018434+2 3.857696-4 5.807690-3 4.318802-2 0.000000+07925 2151     \n"
    " 7.890000+4 1.018035+2 5.977167-4 5.806618-3 4.319215-2 0.000000+07925 2151     \n"
    " 7.945000+4 1.017251+2 1.094127-3 5.804509-3 4.320025-2 0.000000+07925 2151     \n"
    " 8.000000+4 1.016468+2 1.679980-3 5.802394-3 4.320835-2 0.000000+07925 2151     \n"
    " 8.100000+4 1.015046+2 2.935017-3 5.798534-3 4.322307-2 0.000000+07925 2151     \n"
    " 8.200000+4 1.013626+2 4.396938-3 5.794655-3 4.323780-2 0.000000+07925 2151     \n"
    " 8.300000+4 1.012208+2 6.035961-3 5.790757-3 4.325253-2 0.000000+07925 2151     \n"
    " 8.500000+4 1.009378+2 9.767191-3 5.782907-3 4.328198-2 0.000000+07925 2151     \n"
    " 8.750000+4 1.005853+2 1.514601-2 5.772994-3 4.331880-2 0.000000+07925 2151     \n"
    " 9.000000+4 1.002340+2 2.118809-2 5.762969-3 4.335562-2 0.000000+07925 2151     \n"
    " 9.200000+4 9.995393+1 2.643206-2 5.754872-3 4.338508-2 0.000000+07925 2151     \n"
    " 9.350000+4 9.974439+1 3.057842-2 5.748754-3 4.340717-2 0.000000+07925 2151     \n"
    " 9.500000+4 9.953531+1 3.489235-2 5.742598-3 4.342926-2 0.000000+07925 2151     \n"
    " 9.750000+4 9.918787+1 4.242284-2 5.732256-3 4.346608-2 0.000000+07925 2151     \n"
    " 1.000000+5 9.884168+1 5.033883-2 5.721812-3 4.350289-2 0.000000+07925 2151     \n"
    " 1.000000+0 0.000000+0          2          0        606        1007925 2151     \n"
    " 0.000000+0 0.000000+0 2.000000+0 2.000000+0 0.000000+0 0.000000+07925 2151     \n"
    " 2.000000+3 3.848710+1 0.000000+0 2.045961-3 4.205961-2 0.000000+07925 2151     \n"
    " 2.250000+3 3.847356+1 0.000000+0 2.045828-3 4.206329-2 0.000000+07925 2151     \n"
    " 2.500000+3 3.846002+1 0.000000+0 2.045694-3 4.206698-2 0.000000+07925 2151     \n"
    " 2.750000+3 3.844648+1 0.000000+0 2.045558-3 4.207066-2 0.000000+07925 2151     \n"
    " 3.000000+3 3.843295+1 0.000000+0 2.045421-3 4.207434-2 0.000000+07925 2151     \n"
    " 3.250000+3 3.841943+1 0.000000+0 2.045282-3 4.207802-2 0.000000+07925 2151     \n"
    " 3.500000+3 3.840591+1 0.000000+0 2.045143-3 4.208170-2 0.000000+07925 2151     \n"
    " 3.750000+3 3.839239+1 0.000000+0 2.045001-3 4.208538-2 0.000000+07925 2151     \n"
    " 4.000000+3 3.837888+1 0.000000+0 2.044859-3 4.208907-2 0.000000+07925 2151     \n"
    " 4.250000+3 3.836538+1 0.000000+0 2.044716-3 4.209275-2 0.000000+07925 2151     \n"
    " 4.500000+3 3.835188+1 0.000000+0 2.044571-3 4.209643-2 0.000000+07925 2151     \n"
    " 4.750000+3 3.833838+1 0.000000+0 2.044425-3 4.210011-2 0.000000+07925 2151     \n"
    " 5.000000+3 3.832489+1 0.000000+0 2.044278-3 4.210379-2 0.000000+07925 2151     \n"
    " 5.250000+3 3.831141+1 0.000000+0 2.044130-3 4.210748-2 0.000000+07925 2151     \n"
    " 5.500000+3 3.829793+1 0.000000+0 2.043980-3 4.211116-2 0.000000+07925 2151     \n"
    " 6.000000+3 3.827098+1 0.000000+0 2.043679-3 4.211852-2 0.000000+07925 2151     \n"
    " 6.500000+3 3.824406+1 0.000000+0 2.043372-3 4.212589-2 0.000000+07925 2151     \n"
    " 7.000000+3 3.821715+1 0.000000+0 2.043062-3 4.213325-2 0.000000+07925 2151     \n"
    " 7.500000+3 3.819027+1 0.000000+0 2.042748-3 4.214061-2 0.000000+07925 2151     \n"
    " 8.000000+3 3.816340+1 0.000000+0 2.042430-3 4.214798-2 0.000000+07925 2151     \n"
    " 8.500000+3 3.813655+1 0.000000+0 2.042107-3 4.215534-2 0.000000+07925 2151     \n"
    " 9.000000+3 3.810973+1 0.000000+0 2.041782-3 4.216270-2 0.000000+07925 2151     \n"
    " 9.500000+3 3.808292+1 0.000000+0 2.041452-3 4.217007-2 0.000000+07925 2151     \n"
    " 1.000000+4 3.805613+1 0.000000+0 2.041118-3 4.217743-2 0.000000+07925 2151     \n"
    " 1.100000+4 3.800262+1 0.000000+0 2.040441-3 4.219216-2 0.000000+07925 2151     \n"
    " 1.200000+4 3.794918+1 0.000000+0 2.039749-3 4.220689-2 0.000000+07925 2151     \n"
    " 1.300000+4 3.789582+1 0.000000+0 2.039044-3 4.222161-2 0.000000+07925 2151     \n"
    " 1.400000+4 3.784254+1 0.000000+0 2.038325-3 4.223634-2 0.000000+07925 2151     \n"
    " 1.500000+4 3.778933+1 0.000000+0 2.037594-3 4.225107-2 0.000000+07925 2151     \n"
    " 1.600000+4 3.773621+1 0.000000+0 2.036849-3 4.226580-2 0.000000+07925 2151     \n"
    " 1.700000+4 3.768316+1 0.000000+0 2.036092-3 4.228052-2 0.000000+07925 2151     \n"
    " 1.800000+4 3.763019+1 0.000000+0 2.035323-3 4.229525-2 0.000000+07925 2151     \n"
    " 1.900000+4 3.757730+1 0.000000+0 2.034541-3 4.230998-2 0.000000+07925 2151     \n"
    " 2.000000+4 3.752449+1 0.000000+0 2.033748-3 4.232470-2 0.000000+07925 2151     \n"
    " 2.100000+4 3.747175+1 0.000000+0 2.032943-3 4.233943-2 0.000000+07925 2151     \n"
    " 2.200000+4 3.741909+1 0.000000+0 2.032127-3 4.235416-2 0.000000+07925 2151     \n"
    " 2.300000+4 3.736651+1 0.000000+0 2.031299-3 4.236889-2 0.000000+07925 2151     \n"
    " 2.400000+4 3.731400+1 0.000000+0 2.030460-3 4.238361-2 0.000000+07925 2151     \n"
    " 2.500000+4 3.726157+1 0.000000+0 2.029610-3 4.239834-2 0.000000+07925 2151     \n"
    " 2.600000+4 3.720922+1 0.000000+0 2.028749-3 4.241307-2 0.000000+07925 2151     \n"
    " 2.700000+4 3.715694+1 0.000000+0 2.027878-3 4.242780-2 0.000000+07925 2151     \n"
    " 2.800000+4 3.710475+1 0.000000+0 2.026996-3 4.244252-2 0.000000+07925 2151     \n"
    " 2.900000+4 3.705262+1 0.000000+0 2.026103-3 4.245725-2 0.000000+07925 2151     \n"
    " 3.000000+4 3.700058+1 0.000000+0 2.025201-3 4.247198-2 0.000000+07925 2151     \n"
    " 3.100000+4 3.694861+1 0.000000+0 2.024288-3 4.248671-2 0.000000+07925 2151     \n"
    " 3.200000+4 3.689671+1 0.000000+0 2.023365-3 4.250143-2 0.000000+07925 2151     \n"
    " 3.300000+4 3.684490+1 0.000000+0 2.022432-3 4.251616-2 0.000000+07925 2151     \n"
    " 3.400000+4 3.679315+1 0.000000+0 2.021489-3 4.253089-2 0.000000+07925 2151     \n"
    " 3.500000+4 3.674149+1 0.000000+0 2.020537-3 4.254562-2 0.000000+07925 2151     \n"
    " 3.600000+4 3.668990+1 0.000000+0 2.019575-3 4.256034-2 0.000000+07925 2151     \n"
    " 3.700000+4 3.663838+1 0.000000+0 2.018603-3 4.257507-2 0.000000+07925 2151     \n"
    " 3.800000+4 3.658694+1 0.000000+0 2.017622-3 4.258980-2 0.000000+07925 2151     \n"
    " 3.900000+4 3.653558+1 0.000000+0 2.016631-3 4.260452-2 0.000000+07925 2151     \n"
    " 4.000000+4 3.648429+1 0.000000+0 2.015632-3 4.261925-2 0.000000+07925 2151     \n"
    " 4.100000+4 3.643308+1 0.000000+0 2.014623-3 4.263398-2 0.000000+07925 2151     \n"
    " 4.200000+4 3.638194+1 0.000000+0 2.013605-3 4.264871-2 0.000000+07925 2151     \n"
    " 4.300000+4 3.633087+1 0.000000+0 2.012579-3 4.266343-2 0.000000+07925 2151     \n"
    " 4.400000+4 3.627988+1 0.000000+0 2.011543-3 4.267816-2 0.000000+07925 2151     \n"
    " 4.500000+4 3.622897+1 0.000000+0 2.010499-3 4.269289-2 0.000000+07925 2151     \n"
    " 4.600000+4 3.617813+1 0.000000+0 2.009446-3 4.270762-2 0.000000+07925 2151     \n"
    " 4.700000+4 3.612736+1 0.000000+0 2.008384-3 4.272234-2 0.000000+07925 2151     \n"
    " 4.800000+4 3.607667+1 0.000000+0 2.007314-3 4.273707-2 0.000000+07925 2151     \n"
    " 5.000000+4 3.597551+1 0.000000+0 2.005149-3 4.276653-2 0.000000+07925 2151     \n"
    " 5.250000+4 3.584947+1 0.000000+0 2.002395-3 4.280334-2 0.000000+07925 2151     \n"
    " 5.500000+4 3.572390+1 0.000000+0 1.999591-3 4.284016-2 0.000000+07925 2151     \n"
    " 5.750000+4 3.559878+1 0.000000+0 1.996737-3 4.287698-2 0.000000+07925 2151     \n"
    " 6.000000+4 3.547412+1 0.000000+0 1.993835-3 4.291380-2 0.000000+07925 2151     \n"
    " 6.250000+4 3.534992+1 0.000000+0 1.990886-3 4.295062-2 0.000000+07925 2151     \n"
    " 6.500000+4 3.522616+1 0.000000+0 1.987891-3 4.298744-2 0.000000+07925 2151     \n"
    " 6.750000+4 3.510286+1 0.000000+0 1.984850-3 4.302425-2 0.000000+07925 2151     \n"
    " 7.000000+4 3.498001+1 0.000000+0 1.981764-3 4.306107-2 0.000000+07925 2151     \n"
    " 7.200000+4 3.488205+1 0.000000+0 1.979265-3 4.309053-2 0.000000+07925 2151     \n"
    " 7.400000+4 3.478438+1 0.000000+0 1.976738-3 4.311998-2 0.000000+07925 2151     \n"
    " 7.500000+4 3.473565+1 0.000000+0 1.975464-3 4.313471-2 0.000000+07925 2151     \n"
    " 7.775000+4 3.460201+1 0.000000+0 1.971928-3 4.317521-2 0.000000+07925 2151     \n"
    " 7.779600+4 3.459978+1 0.000000+0 1.971869-3 4.317589-2 0.000000+07925 2151     \n"
    " 7.781000+4 3.459910+1 5.813537-7 1.971850-3 4.317609-2 0.000000+07925 2151     \n"
    " 7.783000+4 3.459813+1 2.200081-6 1.971825-3 4.317639-2 0.000000+07925 2151     \n"
    " 7.786000+4 3.459668+1 5.681304-6 1.971786-3 4.317683-2 0.000000+07925 2151     \n"
    " 7.789000+4 3.459523+1 1.011178-5 1.971747-3 4.317727-2 0.000000+07925 2151     \n"
    " 7.793000+4 3.459329+1 1.720825-5 1.971695-3 4.317786-2 0.000000+07925 2151     \n"
    " 7.798000+4 3.459086+1 2.768444-5 1.971630-3 4.317860-2 0.000000+07925 2151     \n"
    " 7.807000+4 3.458650+1 5.029306-5 1.971513-3 4.317992-2 0.000000+07925 2151     \n"
    " 7.820000+4 3.458020+1 9.000590-5 1.971345-3 4.318184-2 0.000000+07925 2151     \n"
    " 7.835000+4 3.457293+1 1.444614-4 1.971150-3 4.318405-2 0.000000+07925 2151     \n"
    " 7.862000+4 3.455985+1 2.618165-4 1.970799-3 4.318802-2 0.000000+07925 2151     \n"
    " 7.890000+4 3.454629+1 4.056618-4 1.970435-3 4.319215-2 0.000000+07925 2151     \n"
    " 7.945000+4 3.451967+1 7.425680-4 1.969717-3 4.320025-2 0.000000+07925 2151     \n"
    " 8.000000+4 3.449307+1 1.140177-3 1.968998-3 4.320835-2 0.000000+07925 2151     \n"
    " 8.100000+4 3.444477+1 1.991948-3 1.967686-3 4.322307-2 0.000000+07925 2151     \n"
    " 8.200000+4 3.439654+1 2.984127-3 1.966367-3 4.323780-2 0.000000+07925 2151     \n"
    " 8.300000+4 3.434837+1 4.096498-3 1.965041-3 4.325253-2 0.000000+07925 2151     \n"
    " 8.500000+4 3.425226+1 6.628800-3 1.962372-3 4.328198-2 0.000000+07925 2151     \n"
    " 8.750000+4 3.413251+1 1.027926-2 1.959002-3 4.331880-2 0.000000+07925 2151     \n"
    " 9.000000+4 3.401320+1 1.437984-2 1.955593-3 4.335562-2 0.000000+07925 2151     \n"
    " 9.200000+4 3.391806+1 1.793875-2 1.952840-3 4.338508-2 0.000000+07925 2151     \n"
    " 9.350000+4 3.384689+1 2.075273-2 1.950761-3 4.340717-2 0.000000+07925 2151     \n"
    " 9.500000+4 3.377587+1 2.368043-2 1.948668-3 4.342926-2 0.000000+07925 2151     \n"
    " 9.750000+4 3.365786+1 2.879106-2 1.945152-3 4.346608-2 0.000000+07925 2151     \n"
    " 1.000000+5 3.354027+1 3.416328-2 1.941601-3 4.350289-2 0.000000+07925 2151     \n"
    " 2.000000+0 0.000000+0          2          0        606        1007925 2151     \n"
    " 0.000000+0 0.000000+0 1.000000+0 2.000000+0 0.000000+0 0.000000+07925 2151     \n"
    " 2.000000+3 2.393741+1 0.000000+0 1.272505-3 4.205961-2 0.000000+07925 2151     \n"
    " 2.250000+3 2.392897+1 0.000000+0 1.272421-3 4.206329-2 0.000000+07925 2151     \n"
    " 2.500000+3 2.392054+1 0.000000+0 1.272337-3 4.206698-2 0.000000+07925 2151     \n"
    " 2.750000+3 2.391210+1 0.000000+0 1.272251-3 4.207066-2 0.000000+07925 2151     \n"
    " 3.000000+3 2.390367+1 0.000000+0 1.272165-3 4.207434-2 0.000000+07925 2151     \n"
    " 3.250000+3 2.389524+1 0.000000+0 1.272078-3 4.207802-2 0.000000+07925 2151     \n"
    " 3.500000+3 2.388682+1 0.000000+0 1.271990-3 4.208170-2 0.000000+07925 2151     \n"
    " 3.750000+3 2.387839+1 0.000000+0 1.271902-3 4.208538-2 0.000000+07925 2151     \n"
    " 4.000000+3 2.386998+1 0.000000+0 1.271812-3 4.208907-2 0.000000+07925 2151     \n"
    " 4.250000+3 2.386156+1 0.000000+0 1.271722-3 4.209275-2 0.000000+07925 2151     \n"
    " 4.500000+3 2.385315+1 0.000000+0 1.271631-3 4.209643-2 0.000000+07925 2151     \n"
    " 4.750000+3 2.384474+1 0.000000+0 1.271540-3 4.210011-2 0.000000+07925 2151     \n"
    " 5.000000+3 2.383633+1 0.000000+0 1.271447-3 4.210379-2 0.000000+07925 2151     \n"
    " 5.250000+3 2.382793+1 0.000000+0 1.271354-3 4.210748-2 0.000000+07925 2151     \n"
    " 5.500000+3 2.381953+1 0.000000+0 1.271261-3 4.211116-2 0.000000+07925 2151     \n"
    " 6.000000+3 2.380274+1 0.000000+0 1.271071-3 4.211852-2 0.000000+07925 2151     \n"
    " 6.500000+3 2.378596+1 0.000000+0 1.270879-3 4.212589-2 0.000000+07925 2151     \n"
    " 7.000000+3 2.376919+1 0.000000+0 1.270684-3 4.213325-2 0.000000+07925 2151     \n"
    " 7.500000+3 2.375244+1 0.000000+0 1.270487-3 4.214061-2 0.000000+07925 2151     \n"
    " 8.000000+3 2.373569+1 0.000000+0 1.270287-3 4.214798-2 0.000000+07925 2151     \n"
    " 8.500000+3 2.371896+1 0.000000+0 1.270085-3 4.215534-2 0.000000+07925 2151     \n"
    " 9.000000+3 2.370225+1 0.000000+0 1.269881-3 4.216270-2 0.000000+07925 2151     \n"
    " 9.500000+3 2.368554+1 0.000000+0 1.269674-3 4.217007-2 0.000000+07925 2151     \n"
    " 1.000000+4 2.366885+1 0.000000+0 1.269465-3 4.217743-2 0.000000+07925 2151     \n"
    " 1.100000+4 2.363550+1 0.000000+0 1.269040-3 4.219216-2 0.000000+07925 2151     \n"
    " 1.200000+4 2.360220+1 0.000000+0 1.268606-3 4.220689-2 0.000000+07925 2151     \n"
    " 1.300000+4 2.356895+1 0.000000+0 1.268164-3 4.222161-2 0.000000+07925 2151     \n"
    " 1.400000+4 2.353575+1 0.000000+0 1.267714-3 4.223634-2 0.000000+07925 2151     \n"
    " 1.500000+4 2.350259+1 0.000000+0 1.267255-3 4.225107-2 0.000000+07925 2151     \n"
    " 1.600000+4 2.346949+1 0.000000+0 1.266789-3 4.226580-2 0.000000+07925 2151     \n"
    " 1.700000+4 2.343643+1 0.000000+0 1.266315-3 4.228052-2 0.000000+07925 2151     \n"
    " 1.800000+4 2.340343+1 0.000000+0 1.265833-3 4.229525-2 0.000000+07925 2151     \n"
    " 1.900000+4 2.337047+1 0.000000+0 1.265343-3 4.230998-2 0.000000+07925 2151     \n"
    " 2.000000+4 2.333756+1 0.000000+0 1.264846-3 4.232470-2 0.000000+07925 2151     \n"
    " 2.100000+4 2.330469+1 0.000000+0 1.264342-3 4.233943-2 0.000000+07925 2151     \n"
    " 2.200000+4 2.327188+1 0.000000+0 1.263831-3 4.235416-2 0.000000+07925 2151     \n"
    " 2.300000+4 2.323911+1 0.000000+0 1.263313-3 4.236889-2 0.000000+07925 2151     \n"
    " 2.400000+4 2.320640+1 0.000000+0 1.262788-3 4.238361-2 0.000000+07925 2151     \n"
    " 2.500000+4 2.317373+1 0.000000+0 1.262256-3 4.239834-2 0.000000+07925 2151     \n"
    " 2.600000+4 2.314110+1 0.000000+0 1.261717-3 4.241307-2 0.000000+07925 2151     \n"
    " 2.700000+4 2.310853+1 0.000000+0 1.261171-3 4.242780-2 0.000000+07925 2151     \n"
    " 2.800000+4 2.307600+1 0.000000+0 1.260620-3 4.244252-2 0.000000+07925 2151     \n"
    " 2.900000+4 2.304353+1 0.000000+0 1.260061-3 4.245725-2 0.000000+07925 2151     \n"
    " 3.000000+4 2.301110+1 0.000000+0 1.259496-3 4.247198-2 0.000000+07925 2151     \n"
    " 3.100000+4 2.297871+1 0.000000+0 1.258925-3 4.248671-2 0.000000+07925 2151     \n"
    " 3.200000+4 2.294638+1 0.000000+0 1.258348-3 4.250143-2 0.000000+07925 2151     \n"
    " 3.300000+4 2.291409+1 0.000000+0 1.257764-3 4.251616-2 0.000000+07925 2151     \n"
    " 3.400000+4 2.288185+1 0.000000+0 1.257174-3 4.253089-2 0.000000+07925 2151     \n"
    " 3.500000+4 2.284965+1 0.000000+0 1.256578-3 4.254562-2 0.000000+07925 2151     \n"
    " 3.600000+4 2.281751+1 0.000000+0 1.255977-3 4.256034-2 0.000000+07925 2151     \n"
    " 3.700000+4 2.278541+1 0.000000+0 1.255369-3 4.257507-2 0.000000+07925 2151     \n"
    " 3.800000+4 2.275336+1 0.000000+0 1.254756-3 4.258980-2 0.000000+07925 2151     \n"
    " 3.900000+4 2.272135+1 0.000000+0 1.254136-3 4.260452-2 0.000000+07925 2151     \n"
    " 4.000000+4 2.268939+1 0.000000+0 1.253511-3 4.261925-2 0.000000+07925 2151     \n"
    " 4.100000+4 2.265748+1 0.000000+0 1.252880-3 4.263398-2 0.000000+07925 2151     \n"
    " 4.200000+4 2.262562+1 0.000000+0 1.252244-3 4.264871-2 0.000000+07925 2151     \n"
    " 4.300000+4 2.259380+1 0.000000+0 1.251602-3 4.266343-2 0.000000+07925 2151     \n"
    " 4.400000+4 2.256203+1 0.000000+0 1.250955-3 4.267816-2 0.000000+07925 2151     \n"
    " 4.500000+4 2.253030+1 0.000000+0 1.250302-3 4.269289-2 0.000000+07925 2151     \n"
    " 4.600000+4 2.249863+1 0.000000+0 1.249644-3 4.270762-2 0.000000+07925 2151     \n"
    " 4.700000+4 2.246699+1 0.000000+0 1.248980-3 4.272234-2 0.000000+07925 2151     \n"
    " 4.800000+4 2.243541+1 0.000000+0 1.248311-3 4.273707-2 0.000000+07925 2151     \n"
    " 5.000000+4 2.237238+1 0.000000+0 1.246958-3 4.276653-2 0.000000+07925 2151     \n"
    " 5.250000+4 2.229385+1 0.000000+0 1.245237-3 4.280334-2 0.000000+07925 2151     \n"
    " 5.500000+4 2.221561+1 0.000000+0 1.243485-3 4.284016-2 0.000000+07925 2151     \n"
    " 5.750000+4 2.213765+1 0.000000+0 1.241702-3 4.287698-2 0.000000+07925 2151     \n"
    " 6.000000+4 2.205998+1 0.000000+0 1.239889-3 4.291380-2 0.000000+07925 2151     \n"
    " 6.250000+4 2.198259+1 0.000000+0 1.238046-3 4.295062-2 0.000000+07925 2151     \n"
    " 6.500000+4 2.190549+1 0.000000+0 1.236175-3 4.298744-2 0.000000+07925 2151     \n"
    " 6.750000+4 2.182866+1 0.000000+0 1.234276-3 4.302425-2 0.000000+07925 2151     \n"
    " 7.000000+4 2.175212+1 0.000000+0 1.232349-3 4.306107-2 0.000000+07925 2151     \n"
    " 7.200000+4 2.169109+1 0.000000+0 1.230788-3 4.309053-2 0.000000+07925 2151     \n"
    " 7.400000+4 2.163024+1 0.000000+0 1.229210-3 4.311998-2 0.000000+07925 2151     \n"
    " 7.500000+4 2.159988+1 0.000000+0 1.228415-3 4.313471-2 0.000000+07925 2151     \n"
    " 7.775000+4 2.151662+1 0.000000+0 1.226207-3 4.317521-2 0.000000+07925 2151     \n" )
  + std::string(
    " 7.779600+4 2.151523+1 0.000000+0 1.226170-3 4.317589-2 0.000000+07925 2151     \n"
    " 7.781000+4 2.151481+1 1.807520-7 1.226158-3 4.317609-2 0.000000+07925 2151     \n"
    " 7.783000+4 2.151420+1 6.840395-7 1.226142-3 4.317639-2 0.000000+07925 2151     \n"
    " 7.786000+4 2.151329+1 1.766406-6 1.226118-3 4.317683-2 0.000000+07925 2151     \n"
    " 7.789000+4 2.151239+1 3.143910-6 1.226094-3 4.317727-2 0.000000+07925 2151     \n"
    " 7.793000+4 2.151118+1 5.350312-6 1.226061-3 4.317786-2 0.000000+07925 2151     \n"
    " 7.798000+4 2.150967+1 8.607523-6 1.226021-3 4.317860-2 0.000000+07925 2151     \n"
    " 7.807000+4 2.150695+1 1.563689-5 1.225948-3 4.317992-2 0.000000+07925 2151     \n"
    " 7.820000+4 2.150303+1 2.798421-5 1.225843-3 4.318184-2 0.000000+07925 2151     \n"
    " 7.835000+4 2.149850+1 4.491526-5 1.225721-3 4.318405-2 0.000000+07925 2151     \n"
    " 7.862000+4 2.149035+1 8.140265-5 1.225502-3 4.318802-2 0.000000+07925 2151     \n"
    " 7.890000+4 2.148190+1 1.261262-4 1.225274-3 4.319215-2 0.000000+07925 2151     \n"
    " 7.945000+4 2.146532+1 2.308750-4 1.224826-3 4.320025-2 0.000000+07925 2151     \n"
    " 8.000000+4 2.144874+1 3.544968-4 1.224377-3 4.320835-2 0.000000+07925 2151     \n"
    " 8.100000+4 2.141865+1 6.193226-4 1.223558-3 4.322307-2 0.000000+07925 2151     \n"
    " 8.200000+4 2.138860+1 9.278012-4 1.222734-3 4.323780-2 0.000000+07925 2151     \n"
    " 8.300000+4 2.135859+1 1.273647-3 1.221907-3 4.325253-2 0.000000+07925 2151     \n"
    " 8.500000+4 2.129871+1 2.060958-3 1.220241-3 4.328198-2 0.000000+07925 2151     \n"
    " 8.750000+4 2.122411+1 3.195900-3 1.218137-3 4.331880-2 0.000000+07925 2151     \n"
    " 9.000000+4 2.114977+1 4.470771-3 1.216009-3 4.335562-2 0.000000+07925 2151     \n"
    " 9.200000+4 2.109050+1 5.577225-3 1.214291-3 4.338508-2 0.000000+07925 2151     \n"
    " 9.350000+4 2.104616+1 6.452076-3 1.212993-3 4.340717-2 0.000000+07925 2151     \n"
    " 9.500000+4 2.100192+1 7.362275-3 1.211687-3 4.342926-2 0.000000+07925 2151     \n"
    " 9.750000+4 2.092840+1 8.951117-3 1.209492-3 4.346608-2 0.000000+07925 2151     \n"
    " 1.000000+5 2.085515+1 1.062126-2 1.207276-3 4.350289-2 0.000000+07925 2151     \n"
    " 3.000000+0 0.000000+0          2          0        606        1007925 2151     \n"
    " 0.000000+0 0.000000+0 2.000000+0 1.000000+0 0.000000+0 0.000000+07925 2151     \n"
    " 2.000000+3 1.804535+1 0.000000+0 9.592847-4 4.205961-2 0.000000+07925 2151     \n"
    " 2.250000+3 1.803897+1 0.000000+0 9.592208-4 4.206329-2 0.000000+07925 2151     \n"
    " 2.500000+3 1.803259+1 0.000000+0 9.591561-4 4.206698-2 0.000000+07925 2151     \n"
    " 2.750000+3 1.802621+1 0.000000+0 9.590908-4 4.207066-2 0.000000+07925 2151     \n"
    " 3.000000+3 1.801984+1 0.000000+0 9.590248-4 4.207434-2 0.000000+07925 2151     \n"
    " 3.250000+3 1.801347+1 0.000000+0 9.589582-4 4.207802-2 0.000000+07925 2151     \n"
    " 3.500000+3 1.800710+1 0.000000+0 9.588910-4 4.208170-2 0.000000+07925 2151     \n"
    " 3.750000+3 1.800073+1 0.000000+0 9.588232-4 4.208538-2 0.000000+07925 2151     \n"
    " 4.000000+3 1.799436+1 0.000000+0 9.587548-4 4.208907-2 0.000000+07925 2151     \n"
    " 4.250000+3 1.798800+1 0.000000+0 9.586859-4 4.209275-2 0.000000+07925 2151     \n"
    " 4.500000+3 1.798164+1 0.000000+0 9.586164-4 4.209643-2 0.000000+07925 2151     \n"
    " 4.750000+3 1.797528+1 0.000000+0 9.585463-4 4.210011-2 0.000000+07925 2151     \n"
    " 5.000000+3 1.796893+1 0.000000+0 9.584757-4 4.210379-2 0.000000+07925 2151     \n"
    " 5.250000+3 1.796257+1 0.000000+0 9.584046-4 4.210748-2 0.000000+07925 2151     \n"
    " 5.500000+3 1.795622+1 0.000000+0 9.583330-4 4.211116-2 0.000000+07925 2151     \n"
    " 6.000000+3 1.794353+1 0.000000+0 9.581882-4 4.211852-2 0.000000+07925 2151     \n"
    " 6.500000+3 1.793084+1 0.000000+0 9.580414-4 4.212589-2 0.000000+07925 2151     \n"
    " 7.000000+3 1.791817+1 0.000000+0 9.578927-4 4.213325-2 0.000000+07925 2151     \n"
    " 7.500000+3 1.790550+1 0.000000+0 9.577420-4 4.214061-2 0.000000+07925 2151     \n"
    " 8.000000+3 1.789284+1 0.000000+0 9.575895-4 4.214798-2 0.000000+07925 2151     \n"
    " 8.500000+3 1.788019+1 0.000000+0 9.574352-4 4.215534-2 0.000000+07925 2151     \n"
    " 9.000000+3 1.786756+1 0.000000+0 9.572791-4 4.216270-2 0.000000+07925 2151     \n"
    " 9.500000+3 1.785493+1 0.000000+0 9.571213-4 4.217007-2 0.000000+07925 2151     \n"
    " 1.000000+4 1.784231+1 0.000000+0 9.569617-4 4.217743-2 0.000000+07925 2151     \n"
    " 1.100000+4 1.781709+1 0.000000+0 9.566375-4 4.219216-2 0.000000+07925 2151     \n"
    " 1.200000+4 1.779192+1 0.000000+0 9.563067-4 4.220689-2 0.000000+07925 2151     \n"
    " 1.300000+4 1.776678+1 0.000000+0 9.559695-4 4.222161-2 0.000000+07925 2151     \n"
    " 1.400000+4 1.774168+1 0.000000+0 9.556260-4 4.223634-2 0.000000+07925 2151     \n"
    " 1.500000+4 1.771661+1 0.000000+0 9.552764-4 4.225107-2 0.000000+07925 2151     \n"
    " 1.600000+4 1.769159+1 0.000000+0 9.549209-4 4.226580-2 0.000000+07925 2151     \n"
    " 1.700000+4 1.766660+1 0.000000+0 9.545594-4 4.228052-2 0.000000+07925 2151     \n"
    " 1.800000+4 1.764164+1 0.000000+0 9.541922-4 4.229525-2 0.000000+07925 2151     \n"
    " 1.900000+4 1.761673+1 0.000000+0 9.538194-4 4.230998-2 0.000000+07925 2151     \n"
    " 2.000000+4 1.759185+1 0.000000+0 9.534410-4 4.232470-2 0.000000+07925 2151     \n"
    " 2.100000+4 1.756700+1 0.000000+0 9.530571-4 4.233943-2 0.000000+07925 2151     \n"
    " 2.200000+4 1.754220+1 0.000000+0 9.526679-4 4.235416-2 0.000000+07925 2151     \n"
    " 2.300000+4 1.751743+1 0.000000+0 9.522733-4 4.236889-2 0.000000+07925 2151     \n"
    " 2.400000+4 1.749269+1 0.000000+0 9.518735-4 4.238361-2 0.000000+07925 2151     \n"
    " 2.500000+4 1.746799+1 0.000000+0 9.514686-4 4.239834-2 0.000000+07925 2151     \n"
    " 2.600000+4 1.744333+1 0.000000+0 9.510586-4 4.241307-2 0.000000+07925 2151     \n"
    " 2.700000+4 1.741871+1 0.000000+0 9.506436-4 4.242780-2 0.000000+07925 2151     \n"
    " 2.800000+4 1.739412+1 0.000000+0 9.502237-4 4.244252-2 0.000000+07925 2151     \n"
    " 2.900000+4 1.736957+1 0.000000+0 9.497989-4 4.245725-2 0.000000+07925 2151     \n"
    " 3.000000+4 1.734505+1 0.000000+0 9.493692-4 4.247198-2 0.000000+07925 2151     \n"
    " 3.100000+4 1.732057+1 0.000000+0 9.489348-4 4.248671-2 0.000000+07925 2151     \n"
    " 3.200000+4 1.729613+1 0.000000+0 9.484957-4 4.250143-2 0.000000+07925 2151     \n"
    " 3.300000+4 1.727172+1 0.000000+0 9.480519-4 4.251616-2 0.000000+07925 2151     \n"
    " 3.400000+4 1.724735+1 0.000000+0 9.476035-4 4.253089-2 0.000000+07925 2151     \n"
    " 3.500000+4 1.722301+1 0.000000+0 9.471506-4 4.254562-2 0.000000+07925 2151     \n"
    " 3.600000+4 1.719871+1 0.000000+0 9.466931-4 4.256034-2 0.000000+07925 2151     \n"
    " 3.700000+4 1.717444+1 0.000000+0 9.462312-4 4.257507-2 0.000000+07925 2151     \n"
    " 3.800000+4 1.715022+1 0.000000+0 9.457650-4 4.258980-2 0.000000+07925 2151     \n"
    " 3.900000+4 1.712602+1 0.000000+0 9.452943-4 4.260452-2 0.000000+07925 2151     \n"
    " 4.000000+4 1.710186+1 0.000000+0 9.448193-4 4.261925-2 0.000000+07925 2151     \n"
    " 4.100000+4 1.707774+1 0.000000+0 9.443401-4 4.263398-2 0.000000+07925 2151     \n"
    " 4.200000+4 1.705365+1 0.000000+0 9.438566-4 4.264871-2 0.000000+07925 2151     \n"
    " 4.300000+4 1.702960+1 0.000000+0 9.433690-4 4.266343-2 0.000000+07925 2151     \n"
    " 4.400000+4 1.700559+1 0.000000+0 9.428772-4 4.267816-2 0.000000+07925 2151     \n"
    " 4.500000+4 1.698161+1 0.000000+0 9.423813-4 4.269289-2 0.000000+07925 2151     \n"
    " 4.600000+4 1.695766+1 0.000000+0 9.418813-4 4.270762-2 0.000000+07925 2151     \n"
    " 4.700000+4 1.693375+1 0.000000+0 9.413773-4 4.272234-2 0.000000+07925 2151     \n"
    " 4.800000+4 1.690988+1 0.000000+0 9.408693-4 4.273707-2 0.000000+07925 2151     \n"
    " 5.000000+4 1.686223+1 0.000000+0 9.398415-4 4.276653-2 0.000000+07925 2151     \n"
    " 5.250000+4 1.680287+1 0.000000+0 9.385350-4 4.280334-2 0.000000+07925 2151     \n"
    " 5.500000+4 1.674373+1 0.000000+0 9.372049-4 4.284016-2 0.000000+07925 2151     \n"
    " 5.750000+4 1.668481+1 0.000000+0 9.358516-4 4.287698-2 0.000000+07925 2151     \n"
    " 6.000000+4 1.662610+1 0.000000+0 9.344757-4 4.291380-2 0.000000+07925 2151     \n"
    " 6.250000+4 1.656760+1 0.000000+0 9.330776-4 4.295062-2 0.000000+07925 2151     \n"
    " 6.500000+4 1.650933+1 0.000000+0 9.316579-4 4.298744-2 0.000000+07925 2151     \n"
    " 6.750000+4 1.645126+1 0.000000+0 9.302170-4 4.302425-2 0.000000+07925 2151     \n"
    " 7.000000+4 1.639341+1 0.000000+0 9.287554-4 4.306107-2 0.000000+07925 2151     \n"
    " 7.200000+4 1.634728+1 0.000000+0 9.275715-4 4.309053-2 0.000000+07925 2151     \n"
    " 7.400000+4 1.630129+1 0.000000+0 9.263748-4 4.311998-2 0.000000+07925 2151     \n"
    " 7.500000+4 1.627834+1 0.000000+0 9.257718-4 4.313471-2 0.000000+07925 2151     \n"
    " 7.775000+4 1.621541+1 0.000000+0 9.240974-4 4.317521-2 0.000000+07925 2151     \n"
    " 7.779600+4 1.621436+1 0.000000+0 9.240692-4 4.317589-2 0.000000+07925 2151     \n"
    " 7.781000+4 1.621404+1 0.000000+0 9.240606-4 4.317609-2 0.000000+07925 2151     \n"
    " 7.783000+4 1.621359+1 0.000000+0 9.240483-4 4.317639-2 0.000000+07925 2151     \n"
    " 7.786000+4 1.621290+1 0.000000+0 9.240299-4 4.317683-2 0.000000+07925 2151     \n"
    " 7.789000+4 1.621222+1 0.000000+0 9.240115-4 4.317727-2 0.000000+07925 2151     \n"
    " 7.793000+4 1.621130+1 0.000000+0 9.239870-4 4.317786-2 0.000000+07925 2151     \n"
    " 7.798000+4 1.621016+1 0.000000+0 9.239563-4 4.317860-2 0.000000+07925 2151     \n"
    " 7.807000+4 1.620811+1 0.000000+0 9.239011-4 4.317992-2 0.000000+07925 2151     \n"
    " 7.820000+4 1.620514+1 0.000000+0 9.238212-4 4.318184-2 0.000000+07925 2151     \n"
    " 7.835000+4 1.620172+1 0.000000+0 9.237290-4 4.318405-2 0.000000+07925 2151     \n"
    " 7.862000+4 1.619556+1 0.000000+0 9.235629-4 4.318802-2 0.000000+07925 2151     \n"
    " 7.890000+4 1.618917+1 0.000000+0 9.233904-4 4.319215-2 0.000000+07925 2151     \n"
    " 7.945000+4 1.617664+1 0.000000+0 9.230508-4 4.320025-2 0.000000+07925 2151     \n"
    " 8.000000+4 1.616412+1 0.000000+0 9.227104-4 4.320835-2 0.000000+07925 2151     \n"
    " 8.100000+4 1.614137+1 0.000000+0 9.220890-4 4.322307-2 0.000000+07925 2151     \n"
    " 8.200000+4 1.611866+1 0.000000+0 9.214648-4 4.323780-2 0.000000+07925 2151     \n"
    " 8.300000+4 1.609598+1 0.000000+0 9.208376-4 4.325253-2 0.000000+07925 2151     \n"
    " 8.500000+4 1.605073+1 0.000000+0 9.195744-4 4.328198-2 0.000000+07925 2151     \n"
    " 8.750000+4 1.599434+1 0.000000+0 9.179795-4 4.331880-2 0.000000+07925 2151     \n"
    " 9.000000+4 1.593816+1 0.000000+0 9.163670-4 4.335562-2 0.000000+07925 2151     \n"
    " 9.200000+4 1.589337+1 0.000000+0 9.150647-4 4.338508-2 0.000000+07925 2151     \n"
    " 9.350000+4 1.585986+1 0.000000+0 9.140809-4 4.340717-2 0.000000+07925 2151     \n"
    " 9.500000+4 1.582643+1 0.000000+0 9.130911-4 4.342926-2 0.000000+07925 2151     \n"
    " 9.750000+4 1.577086+1 0.000000+0 9.114283-4 4.346608-2 0.000000+07925 2151     \n"
    " 1.000000+5 1.571550+1 0.000000+0 9.097494-4 4.350289-2 0.000000+07925 2151     \n"
    " 1.952740+2 0.000000+0          2          0          5          07925 2151     \n"
    " 0.000000+0 0.000000+0          2          0        606        1007925 2151     \n"
    " 0.000000+0 0.000000+0 1.000000+0 1.000000+0 0.000000+0 0.000000+07925 2151     \n"
    " 2.000000+3 1.134047+2 0.000000+0 3.981167-2 1.292688-1 0.000000+07925 2151     \n"
    " 2.250000+3 1.133648+2 0.000000+0 3.979636-2 1.292801-1 0.000000+07925 2151     \n"
    " 2.500000+3 1.133250+2 0.000000+0 3.978105-2 1.292914-1 0.000000+07925 2151     \n"
    " 2.750000+3 1.132851+2 0.000000+0 3.976574-2 1.293027-1 0.000000+07925 2151     \n"
    " 3.000000+3 1.132453+2 0.000000+0 3.975044-2 1.293141-1 0.000000+07925 2151     \n"
    " 3.250000+3 1.132055+2 0.000000+0 3.973513-2 1.293254-1 0.000000+07925 2151     \n"
    " 3.500000+3 1.131657+2 0.000000+0 3.971983-2 1.293367-1 0.000000+07925 2151     \n"
    " 3.750000+3 1.131259+2 0.000000+0 3.970452-2 1.293480-1 0.000000+07925 2151     \n"
    " 4.000000+3 1.130861+2 0.000000+0 3.968921-2 1.293593-1 0.000000+07925 2151     \n"
    " 4.250000+3 1.130464+2 0.000000+0 3.967391-2 1.293706-1 0.000000+07925 2151     \n"
    " 4.500000+3 1.130067+2 0.000000+0 3.965860-2 1.293820-1 0.000000+07925 2151     \n"
    " 4.750000+3 1.129669+2 0.000000+0 3.964330-2 1.293933-1 0.000000+07925 2151     \n"
    " 5.000000+3 1.129272+2 0.000000+0 3.962800-2 1.294046-1 0.000000+07925 2151     \n"
    " 5.250000+3 1.128875+2 0.000000+0 3.961269-2 1.294159-1 0.000000+07925 2151     \n"
    " 5.500000+3 1.128478+2 0.000000+0 3.959739-2 1.294272-1 0.000000+07925 2151     \n"
    " 6.000000+3 1.127685+2 0.000000+0 3.956678-2 1.294498-1 0.000000+07925 2151     \n"
    " 6.500000+3 1.126893+2 0.000000+0 3.953617-2 1.294725-1 0.000000+07925 2151     \n"
    " 7.000000+3 1.126101+2 0.000000+0 3.950557-2 1.294951-1 0.000000+07925 2151     \n"
    " 7.500000+3 1.125309+2 0.000000+0 3.947496-2 1.295177-1 0.000000+07925 2151     \n"
    " 8.000000+3 1.124518+2 0.000000+0 3.944436-2 1.295404-1 0.000000+07925 2151     \n"
    " 8.500000+3 1.123728+2 0.000000+0 3.941375-2 1.295630-1 0.000000+07925 2151     \n"
    " 9.000000+3 1.122938+2 0.000000+0 3.938315-2 1.295856-1 0.000000+07925 2151     \n"
    " 9.500000+3 1.122149+2 0.000000+0 3.935255-2 1.296083-1 0.000000+07925 2151     \n"
    " 1.000000+4 1.121361+2 0.000000+0 3.932194-2 1.296309-1 0.000000+07925 2151     \n"
    " 1.100000+4 1.119785+2 0.000000+0 3.926073-2 1.296762-1 0.000000+07925 2151     \n"
    " 1.200000+4 1.118212+2 0.000000+0 3.919953-2 1.297214-1 0.000000+07925 2151     \n"
    " 1.300000+4 1.116641+2 0.000000+0 3.913832-2 1.297667-1 0.000000+07925 2151     \n"
    " 1.400000+4 1.115073+2 0.000000+0 3.907711-2 1.298120-1 0.000000+07925 2151     \n"
    " 1.500000+4 1.113507+2 0.000000+0 3.901590-2 1.298572-1 0.000000+07925 2151     \n"
    " 1.600000+4 1.111943+2 0.000000+0 3.895468-2 1.299025-1 0.000000+07925 2151     \n"
    " 1.700000+4 1.110381+2 0.000000+0 3.889347-2 1.299478-1 0.000000+07925 2151     \n"
    " 1.800000+4 1.108822+2 0.000000+0 3.883225-2 1.299930-1 0.000000+07925 2151     \n"
    " 1.900000+4 1.107265+2 0.000000+0 3.877103-2 1.300383-1 0.000000+07925 2151     \n"
    " 2.000000+4 1.105710+2 0.000000+0 3.870981-2 1.300835-1 0.000000+07925 2151     \n"
    " 2.100000+4 1.104158+2 0.000000+0 3.864858-2 1.301288-1 0.000000+07925 2151     \n"
    " 2.200000+4 1.102608+2 0.000000+0 3.858735-2 1.301741-1 0.000000+07925 2151     \n"
    " 2.300000+4 1.101060+2 0.000000+0 3.852612-2 1.302193-1 0.000000+07925 2151     \n"
    " 2.400000+4 1.099514+2 0.000000+0 3.846488-2 1.302646-1 0.000000+07925 2151     \n"
    " 2.500000+4 1.097971+2 0.000000+0 3.840364-2 1.303099-1 0.000000+07925 2151     \n"
    " 2.600000+4 1.096429+2 0.000000+0 3.834240-2 1.303551-1 0.000000+07925 2151     \n"
    " 2.700000+4 1.094891+2 0.000000+0 3.828115-2 1.304004-1 0.000000+07925 2151     \n"
    " 2.800000+4 1.093354+2 0.000000+0 3.821990-2 1.304457-1 0.000000+07925 2151     \n"
    " 2.900000+4 1.091820+2 0.000000+0 3.815864-2 1.304909-1 0.000000+07925 2151     \n"
    " 3.000000+4 1.090287+2 0.000000+0 3.809738-2 1.305362-1 0.000000+07925 2151     \n"
    " 3.100000+4 1.088757+2 0.000000+0 3.803612-2 1.305814-1 0.000000+07925 2151     \n"
    " 3.200000+4 1.087230+2 0.000000+0 3.797485-2 1.306267-1 0.000000+07925 2151     \n"
    " 3.300000+4 1.085704+2 0.000000+0 3.791358-2 1.306720-1 0.000000+07925 2151     \n"
    " 3.400000+4 1.084181+2 0.000000+0 3.785230-2 1.307172-1 0.000000+07925 2151     \n"
    " 3.500000+4 1.082660+2 0.000000+0 3.779102-2 1.307625-1 0.000000+07925 2151     \n"
    " 3.600000+4 1.081141+2 0.000000+0 3.772973-2 1.308078-1 0.000000+07925 2151     \n"
    " 3.700000+4 1.079625+2 0.000000+0 3.766844-2 1.308530-1 0.000000+07925 2151     \n"
    " 3.800000+4 1.078111+2 0.000000+0 3.760715-2 1.308983-1 0.000000+07925 2151     \n"
    " 3.900000+4 1.076599+2 0.000000+0 3.754585-2 1.309436-1 0.000000+07925 2151     \n"
    " 4.000000+4 1.075089+2 0.000000+0 3.748455-2 1.309888-1 0.000000+07925 2151     \n"
    " 4.100000+4 1.073581+2 0.000000+0 3.742324-2 1.310341-1 0.000000+07925 2151     \n"
    " 4.200000+4 1.072076+2 0.000000+0 3.736193-2 1.310794-1 0.000000+07925 2151     \n"
    " 4.300000+4 1.070572+2 0.000000+0 3.730061-2 1.311246-1 0.000000+07925 2151     \n"
    " 4.400000+4 1.069071+2 0.000000+0 3.723929-2 1.311699-1 0.000000+07925 2151     \n"
    " 4.500000+4 1.067572+2 0.000000+0 3.717797-2 1.312151-1 0.000000+07925 2151     \n"
    " 4.600000+4 1.066076+2 0.000000+0 3.711664-2 1.312604-1 0.000000+07925 2151     \n"
    " 4.700000+4 1.064581+2 0.000000+0 3.705531-2 1.313057-1 0.000000+07925 2151     \n"
    " 4.800000+4 1.063089+2 0.000000+0 3.699397-2 1.313509-1 0.000000+07925 2151     \n"
    " 5.000000+4 1.060111+2 0.000000+0 3.687129-2 1.314415-1 0.000000+07925 2151     \n"
    " 5.250000+4 1.056400+2 0.000000+0 3.671792-2 1.315546-1 0.000000+07925 2151     \n"
    " 5.500000+4 1.052704+2 0.000000+0 3.656452-2 1.316678-1 0.000000+07925 2151     \n"
    " 5.750000+4 1.049020+2 0.000000+0 3.641111-2 1.317809-1 0.000000+07925 2151     \n"
    " 6.000000+4 1.045350+2 0.000000+0 3.625768-2 1.318941-1 0.000000+07925 2151     \n"
    " 6.250000+4 1.041694+2 0.000000+0 3.610424-2 1.320073-1 0.000000+07925 2151     \n"
    " 6.500000+4 1.038050+2 0.000000+0 3.595079-2 1.321204-1 0.000000+07925 2151     \n"
    " 6.750000+4 1.034420+2 0.000000+0 3.579732-2 1.322336-1 0.000000+07925 2151     \n"
    " 7.000000+4 1.030804+2 0.000000+0 3.564386-2 1.323467-1 0.000000+07925 2151     \n"
    " 7.200000+4 1.027920+2 0.000000+0 3.552108-2 1.324373-1 0.000000+07925 2151     \n"
    " 7.400000+4 1.025044+2 0.000000+0 3.539831-2 1.325278-1 0.000000+07925 2151     \n"
    " 7.500000+4 1.023610+2 0.000000+0 3.533692-2 1.325731-1 0.000000+07925 2151     \n"
    " 7.775000+4 1.019675+2 0.000000+0 3.516811-2 1.326975-1 0.000000+07925 2151     \n"
    " 7.779600+4 1.019610+2 0.000000+0 3.516528-2 1.326996-1 0.000000+07925 2151     \n"
    " 7.781000+4 1.019590+2 6.529884-2 3.516442-2 1.327003-1 0.000000+07925 2151     \n"
    " 7.783000+4 1.019561+2 1.016447-1 3.516320-2 1.327012-1 0.000000+07925 2151     \n"
    " 7.786000+4 1.019518+2 1.392890-1 3.516135-2 1.327025-1 0.000000+07925 2151     \n"
    " 7.789000+4 1.019476+2 1.686488-1 3.515951-2 1.327039-1 0.000000+07925 2151     \n"
    " 7.793000+4 1.019418+2 2.011487-1 3.515706-2 1.327057-1 0.000000+07925 2151     \n"
    " 7.798000+4 1.019347+2 2.354455-1 3.515399-2 1.327080-1 0.000000+07925 2151     \n"
    " 7.807000+4 1.019219+2 2.868302-1 3.514846-2 1.327120-1 0.000000+07925 2151     \n"
    " 7.820000+4 1.019033+2 3.475866-1 3.514048-2 1.327179-1 0.000000+07925 2151     \n"
    " 7.835000+4 1.018819+2 4.062266-1 3.513128-2 1.327247-1 0.000000+07925 2151     \n"
    " 7.862000+4 1.018434+2 4.939346-1 3.511470-2 1.327369-1 0.000000+07925 2151     \n"
    " 7.890000+4 1.018035+2 5.702040-1 3.509751-2 1.327496-1 0.000000+07925 2151     \n"
    " 7.945000+4 1.017251+2 6.948175-1 3.506375-2 1.327745-1 0.000000+07925 2151     \n"
    " 8.000000+4 1.016468+2 7.989844-1 3.502999-2 1.327994-1 0.000000+07925 2151     \n"
    " 8.100000+4 1.015046+2 9.575282-1 3.496861-2 1.328447-1 0.000000+07925 2151     \n"
    " 8.200000+4 1.013626+2 1.091042+0 3.490724-2 1.328899-1 0.000000+07925 2151     \n"
    " 8.300000+4 1.012208+2 1.208096+0 3.484586-2 1.329352-1 0.000000+07925 2151     \n"
    " 8.500000+4 1.009378+2 1.409286+0 3.472311-2 1.330257-1 0.000000+07925 2151     \n"
    " 8.750000+4 1.005853+2 1.620208+0 3.456970-2 1.331389-1 0.000000+07925 2151     \n"
    " 9.000000+4 1.002340+2 1.801435+0 3.441630-2 1.332520-1 0.000000+07925 2151     \n"
    " 9.200000+4 9.995393+1 1.931041+0 3.429361-2 1.333426-1 0.000000+07925 2151     \n"
    " 9.350000+4 9.974439+1 2.021146+0 3.420160-2 1.334105-1 0.000000+07925 2151     \n"
    " 9.500000+4 9.953531+1 2.106109+0 3.410960-2 1.334783-1 0.000000+07925 2151     \n"
    " 9.750000+4 9.918787+1 2.237995+0 3.395629-2 1.335915-1 0.000000+07925 2151     \n"
    " 1.000000+5 9.884168+1 2.359643+0 3.380303-2 1.337047-1 0.000000+07925 2151     \n"
    " 1.000000+0 0.000000+0          2          0        606        1007925 2151     \n"
    " 0.000000+0 0.000000+0 1.000000+0 2.000000+0 0.000000+0 0.000000+07925 2151     \n"
    " 2.000000+3 3.848710+1 0.000000+0 1.351122-2 1.000000-9 0.000000+07925 2151     \n"
    " 2.250000+3 3.847356+1 0.000000+0 1.350602-2 1.000000-9 0.000000+07925 2151     \n"
    " 2.500000+3 3.846002+1 0.000000+0 1.350082-2 1.000000-9 0.000000+07925 2151     \n"
    " 2.750000+3 3.844648+1 0.000000+0 1.349562-2 1.000000-9 0.000000+07925 2151     \n"
    " 3.000000+3 3.843295+1 0.000000+0 1.349042-2 1.000000-9 0.000000+07925 2151     \n"
    " 3.250000+3 3.841943+1 0.000000+0 1.348522-2 1.000000-9 0.000000+07925 2151     \n"
    " 3.500000+3 3.840591+1 0.000000+0 1.348002-2 1.000000-9 0.000000+07925 2151     \n"
    " 3.750000+3 3.839239+1 0.000000+0 1.347482-2 1.000000-9 0.000000+07925 2151     \n"
    " 4.000000+3 3.837888+1 0.000000+0 1.346962-2 1.000000-9 0.000000+07925 2151     \n"
    " 4.250000+3 3.836538+1 0.000000+0 1.346442-2 1.000000-9 0.000000+07925 2151     \n"
    " 4.500000+3 3.835188+1 0.000000+0 1.345923-2 1.000000-9 0.000000+07925 2151     \n"
    " 4.750000+3 3.833838+1 0.000000+0 1.345403-2 1.000000-9 0.000000+07925 2151     \n"
    " 5.000000+3 3.832489+1 0.000000+0 1.344883-2 1.000000-9 0.000000+07925 2151     \n"
    " 5.250000+3 3.831141+1 0.000000+0 1.344363-2 1.000000-9 0.000000+07925 2151     \n"
    " 5.500000+3 3.829793+1 0.000000+0 1.343843-2 1.000000-9 0.000000+07925 2151     \n"
    " 6.000000+3 3.827098+1 0.000000+0 1.342804-2 1.000000-9 0.000000+07925 2151     \n"
    " 6.500000+3 3.824406+1 0.000000+0 1.341764-2 1.000000-9 0.000000+07925 2151     \n"
    " 7.000000+3 3.821715+1 0.000000+0 1.340724-2 1.000000-9 0.000000+07925 2151     \n"
    " 7.500000+3 3.819027+1 0.000000+0 1.339685-2 1.000000-9 0.000000+07925 2151     \n"
    " 8.000000+3 3.816340+1 0.000000+0 1.338645-2 1.000000-9 0.000000+07925 2151     \n"
    " 8.500000+3 3.813655+1 0.000000+0 1.337606-2 1.000000-9 0.000000+07925 2151     \n"
    " 9.000000+3 3.810973+1 0.000000+0 1.336566-2 1.000000-9 0.000000+07925 2151     \n"
    " 9.500000+3 3.808292+1 0.000000+0 1.335526-2 1.000000-9 0.000000+07925 2151     \n"
    " 1.000000+4 3.805613+1 0.000000+0 1.334487-2 1.000000-9 0.000000+07925 2151     \n"
    " 1.100000+4 3.800262+1 0.000000+0 1.332408-2 1.000000-9 0.000000+07925 2151     \n"
    " 1.200000+4 3.794918+1 0.000000+0 1.330329-2 1.000000-9 0.000000+07925 2151     \n"
    " 1.300000+4 3.789582+1 0.000000+0 1.328250-2 1.000000-9 0.000000+07925 2151     \n"
    " 1.400000+4 3.784254+1 0.000000+0 1.326171-2 1.000000-9 0.000000+07925 2151     \n"
    " 1.500000+4 3.778933+1 0.000000+0 1.324091-2 1.000000-9 0.000000+07925 2151     \n"
    " 1.600000+4 3.773621+1 0.000000+0 1.322012-2 1.000000-9 0.000000+07925 2151     \n"
    " 1.700000+4 3.768316+1 0.000000+0 1.319933-2 1.000000-9 0.000000+07925 2151     \n"
    " 1.800000+4 3.763019+1 0.000000+0 1.317854-2 1.000000-9 0.000000+07925 2151     \n"
    " 1.900000+4 3.757730+1 0.000000+0 1.315774-2 1.000000-9 0.000000+07925 2151     \n"
    " 2.000000+4 3.752449+1 0.000000+0 1.313695-2 1.000000-9 0.000000+07925 2151     \n"
    " 2.100000+4 3.747175+1 0.000000+0 1.311615-2 1.000000-9 0.000000+07925 2151     \n"
    " 2.200000+4 3.741909+1 0.000000+0 1.309535-2 1.000000-9 0.000000+07925 2151     \n"
    " 2.300000+4 3.736651+1 0.000000+0 1.307456-2 1.000000-9 0.000000+07925 2151     \n"
    " 2.400000+4 3.731400+1 0.000000+0 1.305376-2 1.000000-9 0.000000+07925 2151     \n"
    " 2.500000+4 3.726157+1 0.000000+0 1.303296-2 1.000000-9 0.000000+07925 2151     \n"
    " 2.600000+4 3.720922+1 0.000000+0 1.301215-2 1.000000-9 0.000000+07925 2151     \n"
    " 2.700000+4 3.715694+1 0.000000+0 1.299135-2 1.000000-9 0.000000+07925 2151     \n"
    " 2.800000+4 3.710475+1 0.000000+0 1.297055-2 1.000000-9 0.000000+07925 2151     \n"
    " 2.900000+4 3.705262+1 0.000000+0 1.294974-2 1.000000-9 0.000000+07925 2151     \n"
    " 3.000000+4 3.700058+1 0.000000+0 1.292893-2 1.000000-9 0.000000+07925 2151     \n"
    " 3.100000+4 3.694861+1 0.000000+0 1.290812-2 1.000000-9 0.000000+07925 2151     \n"
    " 3.200000+4 3.689671+1 0.000000+0 1.288731-2 1.000000-9 0.000000+07925 2151     \n"
    " 3.300000+4 3.684490+1 0.000000+0 1.286650-2 1.000000-9 0.000000+07925 2151     \n"
    " 3.400000+4 3.679315+1 0.000000+0 1.284569-2 1.000000-9 0.000000+07925 2151     \n"
    " 3.500000+4 3.674149+1 0.000000+0 1.282488-2 1.000000-9 0.000000+07925 2151     \n"
    " 3.600000+4 3.668990+1 0.000000+0 1.280406-2 1.000000-9 0.000000+07925 2151     \n"
    " 3.700000+4 3.663838+1 0.000000+0 1.278324-2 1.000000-9 0.000000+07925 2151     \n"
    " 3.800000+4 3.658694+1 0.000000+0 1.276243-2 1.000000-9 0.000000+07925 2151     \n"
    " 3.900000+4 3.653558+1 0.000000+0 1.274161-2 1.000000-9 0.000000+07925 2151     \n"
    " 4.000000+4 3.648429+1 0.000000+0 1.272078-2 1.000000-9 0.000000+07925 2151     \n"
    " 4.100000+4 3.643308+1 0.000000+0 1.269996-2 1.000000-9 0.000000+07925 2151     \n"
    " 4.200000+4 3.638194+1 0.000000+0 1.267914-2 1.000000-9 0.000000+07925 2151     \n"
    " 4.300000+4 3.633087+1 0.000000+0 1.265831-2 1.000000-9 0.000000+07925 2151     \n"
    " 4.400000+4 3.627988+1 0.000000+0 1.263749-2 1.000000-9 0.000000+07925 2151     \n"
    " 4.500000+4 3.622897+1 0.000000+0 1.261666-2 1.000000-9 0.000000+07925 2151     \n"
    " 4.600000+4 3.617813+1 0.000000+0 1.259583-2 1.000000-9 0.000000+07925 2151     \n"
    " 4.700000+4 3.612736+1 0.000000+0 1.257500-2 1.000000-9 0.000000+07925 2151     \n"
    " 4.800000+4 3.607667+1 0.000000+0 1.255417-2 1.000000-9 0.000000+07925 2151     \n"
    " 5.000000+4 3.597551+1 0.000000+0 1.251250-2 1.000000-9 0.000000+07925 2151     \n"
    " 5.250000+4 3.584947+1 0.000000+0 1.246041-2 1.000000-9 0.000000+07925 2151     \n"
    " 5.500000+4 3.572390+1 0.000000+0 1.240831-2 1.000000-9 0.000000+07925 2151     \n"
    " 5.750000+4 3.559878+1 0.000000+0 1.235621-2 1.000000-9 0.000000+07925 2151     \n"
    " 6.000000+4 3.547412+1 0.000000+0 1.230410-2 1.000000-9 0.000000+07925 2151     \n"
    " 6.250000+4 3.534992+1 0.000000+0 1.225199-2 1.000000-9 0.000000+07925 2151     \n"
    " 6.500000+4 3.522616+1 0.000000+0 1.219987-2 1.000000-9 0.000000+07925 2151     \n"
    " 6.750000+4 3.510286+1 0.000000+0 1.214775-2 1.000000-9 0.000000+07925 2151     \n"
    " 7.000000+4 3.498001+1 0.000000+0 1.209563-2 1.000000-9 0.000000+07925 2151     \n"
    " 7.200000+4 3.488205+1 0.000000+0 1.205394-2 1.000000-9 0.000000+07925 2151     \n"
    " 7.400000+4 3.478438+1 0.000000+0 1.201224-2 1.000000-9 0.000000+07925 2151     \n"
    " 7.500000+4 3.473565+1 0.000000+0 1.199140-2 1.000000-9 0.000000+07925 2151     \n"
    " 7.775000+4 3.460201+1 0.000000+0 1.193407-2 1.000000-9 0.000000+07925 2151     \n"
    " 7.779600+4 3.459978+1 0.000000+0 1.193311-2 1.000000-9 0.000000+07925 2151     \n"
    " 7.781000+4 3.459910+1 1.000000-9 1.193281-2 1.000000-9 0.000000+07925 2151     \n"
    " 7.783000+4 3.459813+1 1.000000-9 1.193240-2 1.000000-9 0.000000+07925 2151     \n"
    " 7.786000+4 3.459668+1 1.000000-9 1.193177-2 1.000000-9 0.000000+07925 2151     \n"
    " 7.789000+4 3.459523+1 1.000000-9 1.193115-2 1.000000-9 0.000000+07925 2151     \n"
    " 7.793000+4 3.459329+1 1.000000-9 1.193031-2 1.000000-9 0.000000+07925 2151     \n"
    " 7.798000+4 3.459086+1 1.000000-9 1.192927-2 1.000000-9 0.000000+07925 2151     \n"
    " 7.807000+4 3.458650+1 1.000000-9 1.192739-2 1.000000-9 0.000000+07925 2151     \n"
    " 7.820000+4 3.458020+1 1.000000-9 1.192468-2 1.000000-9 0.000000+07925 2151     \n"
    " 7.835000+4 3.457293+1 1.000000-9 1.192156-2 1.000000-9 0.000000+07925 2151     \n"
    " 7.862000+4 3.455985+1 1.000000-9 1.191593-2 1.000000-9 0.000000+07925 2151     \n"
    " 7.890000+4 3.454629+1 1.000000-9 1.191009-2 1.000000-9 0.000000+07925 2151     \n"
    " 7.945000+4 3.451967+1 1.000000-9 1.189863-2 1.000000-9 0.000000+07925 2151     \n"
    " 8.000000+4 3.449307+1 1.000000-9 1.188716-2 1.000000-9 0.000000+07925 2151     \n"
    " 8.100000+4 3.444477+1 1.000000-9 1.186632-2 1.000000-9 0.000000+07925 2151     \n"
    " 8.200000+4 3.439654+1 1.000000-9 1.184547-2 1.000000-9 0.000000+07925 2151     \n"
    " 8.300000+4 3.434837+1 1.000000-9 1.182463-2 1.000000-9 0.000000+07925 2151     \n"
    " 8.500000+4 3.425226+1 1.000000-9 1.178294-2 1.000000-9 0.000000+07925 2151     \n"
    " 8.750000+4 3.413251+1 1.000000-9 1.173085-2 1.000000-9 0.000000+07925 2151     \n"
    " 9.000000+4 3.401320+1 1.000000-9 1.167875-2 1.000000-9 0.000000+07925 2151     \n"
    " 9.200000+4 3.391806+1 1.000000-9 1.163709-2 1.000000-9 0.000000+07925 2151     \n"
    " 9.350000+4 3.384689+1 1.000000-9 1.160584-2 1.000000-9 0.000000+07925 2151     \n"
    " 9.500000+4 3.377587+1 1.000000-9 1.157460-2 1.000000-9 0.000000+07925 2151     \n"
    " 9.750000+4 3.365786+1 1.000000-9 1.152254-2 1.000000-9 0.000000+07925 2151     \n"
    " 1.000000+5 3.354027+1 1.000000-9 1.147049-2 1.000000-9 0.000000+07925 2151     \n"
    " 2.000000+0 0.000000+0          2          0        606        1007925 2151     \n"
    " 0.000000+0 0.000000+0 2.000000+0 2.000000+0 0.000000+0 0.000000+07925 2151     \n"
    " 2.000000+3 2.393741+1 0.000000+0 8.403428-3 1.000000-9 0.000000+07925 2151     \n"
    " 2.250000+3 2.392897+1 0.000000+0 8.400189-3 1.000000-9 0.000000+07925 2151     \n"
    " 2.500000+3 2.392054+1 0.000000+0 8.396949-3 1.000000-9 0.000000+07925 2151     \n"
    " 2.750000+3 2.391210+1 0.000000+0 8.393709-3 1.000000-9 0.000000+07925 2151     \n"
    " 3.000000+3 2.390367+1 0.000000+0 8.390470-3 1.000000-9 0.000000+07925 2151     \n"
    " 3.250000+3 2.389524+1 0.000000+0 8.387230-3 1.000000-9 0.000000+07925 2151     \n"
    " 3.500000+3 2.388682+1 0.000000+0 8.383991-3 1.000000-9 0.000000+07925 2151     \n"
    " 3.750000+3 2.387839+1 0.000000+0 8.380752-3 1.000000-9 0.000000+07925 2151     \n"
    " 4.000000+3 2.386998+1 0.000000+0 8.377512-3 1.000000-9 0.000000+07925 2151     \n"
    " 4.250000+3 2.386156+1 0.000000+0 8.374273-3 1.000000-9 0.000000+07925 2151     \n"
    " 4.500000+3 2.385315+1 0.000000+0 8.371034-3 1.000000-9 0.000000+07925 2151     \n"
    " 4.750000+3 2.384474+1 0.000000+0 8.367795-3 1.000000-9 0.000000+07925 2151     \n"
    " 5.000000+3 2.383633+1 0.000000+0 8.364556-3 1.000000-9 0.000000+07925 2151     \n"
    " 5.250000+3 2.382793+1 0.000000+0 8.361317-3 1.000000-9 0.000000+07925 2151     \n"
    " 5.500000+3 2.381953+1 0.000000+0 8.358078-3 1.000000-9 0.000000+07925 2151     \n"
    " 6.000000+3 2.380274+1 0.000000+0 8.351601-3 1.000000-9 0.000000+07925 2151     \n"
    " 6.500000+3 2.378596+1 0.000000+0 8.345123-3 1.000000-9 0.000000+07925 2151     \n"
    " 7.000000+3 2.376919+1 0.000000+0 8.338646-3 1.000000-9 0.000000+07925 2151     \n"
    " 7.500000+3 2.375244+1 0.000000+0 8.332169-3 1.000000-9 0.000000+07925 2151     \n"
    " 8.000000+3 2.373569+1 0.000000+0 8.325692-3 1.000000-9 0.000000+07925 2151     \n"
    " 8.500000+3 2.371896+1 0.000000+0 8.319215-3 1.000000-9 0.000000+07925 2151     \n"
    " 9.000000+3 2.370225+1 0.000000+0 8.312738-3 1.000000-9 0.000000+07925 2151     \n"
    " 9.500000+3 2.368554+1 0.000000+0 8.306261-3 1.000000-9 0.000000+07925 2151     \n"
    " 1.000000+4 2.366885+1 0.000000+0 8.299785-3 1.000000-9 0.000000+07925 2151     \n"
    " 1.100000+4 2.363550+1 0.000000+0 8.286831-3 1.000000-9 0.000000+07925 2151     \n"
    " 1.200000+4 2.360220+1 0.000000+0 8.273878-3 1.000000-9 0.000000+07925 2151     \n"
    " 1.300000+4 2.356895+1 0.000000+0 8.260925-3 1.000000-9 0.000000+07925 2151     \n"
    " 1.400000+4 2.353575+1 0.000000+0 8.247971-3 1.000000-9 0.000000+07925 2151     \n"
    " 1.500000+4 2.350259+1 0.000000+0 8.235018-3 1.000000-9 0.000000+07925 2151     \n"
    " 1.600000+4 2.346949+1 0.000000+0 8.222064-3 1.000000-9 0.000000+07925 2151     \n"
    " 1.700000+4 2.343643+1 0.000000+0 8.209110-3 1.000000-9 0.000000+07925 2151     \n"
    " 1.800000+4 2.340343+1 0.000000+0 8.196155-3 1.000000-9 0.000000+07925 2151     \n"
    " 1.900000+4 2.337047+1 0.000000+0 8.183200-3 1.000000-9 0.000000+07925 2151     \n"
    " 2.000000+4 2.333756+1 0.000000+0 8.170245-3 1.000000-9 0.000000+07925 2151     \n"
    " 2.100000+4 2.330469+1 0.000000+0 8.157289-3 1.000000-9 0.000000+07925 2151     \n"
    " 2.200000+4 2.327188+1 0.000000+0 8.144332-3 1.000000-9 0.000000+07925 2151     \n"
    " 2.300000+4 2.323911+1 0.000000+0 8.131375-3 1.000000-9 0.000000+07925 2151     \n"
    " 2.400000+4 2.320640+1 0.000000+0 8.118417-3 1.000000-9 0.000000+07925 2151     \n"
    " 2.500000+4 2.317373+1 0.000000+0 8.105459-3 1.000000-9 0.000000+07925 2151     \n"
    " 2.600000+4 2.314110+1 0.000000+0 8.092500-3 1.000000-9 0.000000+07925 2151     \n"
    " 2.700000+4 2.310853+1 0.000000+0 8.079540-3 1.000000-9 0.000000+07925 2151     \n"
    " 2.800000+4 2.307600+1 0.000000+0 8.066579-3 1.000000-9 0.000000+07925 2151     \n"
    " 2.900000+4 2.304353+1 0.000000+0 8.053618-3 1.000000-9 0.000000+07925 2151     \n"
    " 3.000000+4 2.301110+1 0.000000+0 8.040656-3 1.000000-9 0.000000+07925 2151     \n"
    " 3.100000+4 2.297871+1 0.000000+0 8.027693-3 1.000000-9 0.000000+07925 2151     \n"
    " 3.200000+4 2.294638+1 0.000000+0 8.014729-3 1.000000-9 0.000000+07925 2151     \n"
    " 3.300000+4 2.291409+1 0.000000+0 8.001765-3 1.000000-9 0.000000+07925 2151     \n"
    " 3.400000+4 2.288185+1 0.000000+0 7.988799-3 1.000000-9 0.000000+07925 2151     \n"
    " 3.500000+4 2.284965+1 0.000000+0 7.975833-3 1.000000-9 0.000000+07925 2151     \n"
    " 3.600000+4 2.281751+1 0.000000+0 7.962866-3 1.000000-9 0.000000+07925 2151     \n"
    " 3.700000+4 2.278541+1 0.000000+0 7.949899-3 1.000000-9 0.000000+07925 2151     \n"
    " 3.800000+4 2.275336+1 0.000000+0 7.936930-3 1.000000-9 0.000000+07925 2151     \n"
    " 3.900000+4 2.272135+1 0.000000+0 7.923961-3 1.000000-9 0.000000+07925 2151     \n"
    " 4.000000+4 2.268939+1 0.000000+0 7.910991-3 1.000000-9 0.000000+07925 2151     \n"
    " 4.100000+4 2.265748+1 0.000000+0 7.898020-3 1.000000-9 0.000000+07925 2151     \n"
    " 4.200000+4 2.262562+1 0.000000+0 7.885048-3 1.000000-9 0.000000+07925 2151     \n"
    " 4.300000+4 2.259380+1 0.000000+0 7.872075-3 1.000000-9 0.000000+07925 2151     \n"
    " 4.400000+4 2.256203+1 0.000000+0 7.859102-3 1.000000-9 0.000000+07925 2151     \n"
    " 4.500000+4 2.253030+1 0.000000+0 7.846128-3 1.000000-9 0.000000+07925 2151     \n"
    " 4.600000+4 2.249863+1 0.000000+0 7.833153-3 1.000000-9 0.000000+07925 2151     \n"
    " 4.700000+4 2.246699+1 0.000000+0 7.820178-3 1.000000-9 0.000000+07925 2151     \n"
    " 4.800000+4 2.243541+1 0.000000+0 7.807202-3 1.000000-9 0.000000+07925 2151     \n"
    " 5.000000+4 2.237238+1 0.000000+0 7.781247-3 1.000000-9 0.000000+07925 2151     \n"
    " 5.250000+4 2.229385+1 0.000000+0 7.748801-3 1.000000-9 0.000000+07925 2151     \n"
    " 5.500000+4 2.221561+1 0.000000+0 7.716351-3 1.000000-9 0.000000+07925 2151     \n"
    " 5.750000+4 2.213765+1 0.000000+0 7.683898-3 1.000000-9 0.000000+07925 2151     \n"
    " 6.000000+4 2.205998+1 0.000000+0 7.651442-3 1.000000-9 0.000000+07925 2151     \n"
    " 6.250000+4 2.198259+1 0.000000+0 7.618984-3 1.000000-9 0.000000+07925 2151     \n"
    " 6.500000+4 2.190549+1 0.000000+0 7.586524-3 1.000000-9 0.000000+07925 2151     \n"
    " 6.750000+4 2.182866+1 0.000000+0 7.554064-3 1.000000-9 0.000000+07925 2151     \n"
    " 7.000000+4 2.175212+1 0.000000+0 7.521602-3 1.000000-9 0.000000+07925 2151     \n"
    " 7.200000+4 2.169109+1 0.000000+0 7.495633-3 1.000000-9 0.000000+07925 2151     \n"
    " 7.400000+4 2.163024+1 0.000000+0 7.469665-3 1.000000-9 0.000000+07925 2151     \n"
    " 7.500000+4 2.159988+1 0.000000+0 7.456681-3 1.000000-9 0.000000+07925 2151     \n"
    " 7.775000+4 2.151662+1 0.000000+0 7.420976-3 1.000000-9 0.000000+07925 2151     \n"
    " 7.779600+4 2.151523+1 0.000000+0 7.420379-3 1.000000-9 0.000000+07925 2151     \n"
    " 7.781000+4 2.151481+1 0.000000+0 7.420197-3 1.000000-9 0.000000+07925 2151     \n"
    " 7.783000+4 2.151420+1 0.000000+0 7.419938-3 1.000000-9 0.000000+07925 2151     \n"
    " 7.786000+4 2.151329+1 0.000000+0 7.419548-3 1.000000-9 0.000000+07925 2151     \n"
    " 7.789000+4 2.151239+1 1.000000-9 7.419159-3 1.000000-9 0.000000+07925 2151     \n"
    " 7.793000+4 2.151118+1 1.000000-9 7.418639-3 1.000000-9 0.000000+07925 2151     \n"
    " 7.798000+4 2.150967+1 1.000000-9 7.417990-3 1.000000-9 0.000000+07925 2151     \n"
    " 7.807000+4 2.150695+1 1.000000-9 7.416822-3 1.000000-9 0.000000+07925 2151     \n"
    " 7.820000+4 2.150303+1 1.000000-9 7.415134-3 1.000000-9 0.000000+07925 2151     \n"
    " 7.835000+4 2.149850+1 1.000000-9 7.413187-3 1.000000-9 0.000000+07925 2151     \n"
    " 7.862000+4 2.149035+1 1.000000-9 7.409681-3 1.000000-9 0.000000+07925 2151     \n"
    " 7.890000+4 2.148190+1 1.000000-9 7.406046-3 1.000000-9 0.000000+07925 2151     \n"
    " 7.945000+4 2.146532+1 1.000000-9 7.398906-3 1.000000-9 0.000000+07925 2151     \n"
    " 8.000000+4 2.144874+1 1.000000-9 7.391765-3 1.000000-9 0.000000+07925 2151     \n"
    " 8.100000+4 2.141865+1 1.000000-9 7.378784-3 1.000000-9 0.000000+07925 2151     \n"
    " 8.200000+4 2.138860+1 1.000000-9 7.365802-3 1.000000-9 0.000000+07925 2151     \n"
    " 8.300000+4 2.135859+1 1.000000-9 7.352821-3 1.000000-9 0.000000+07925 2151     \n"
    " 8.500000+4 2.129871+1 1.000000-9 7.326862-3 1.000000-9 0.000000+07925 2151     \n"
    " 8.750000+4 2.122411+1 1.000000-9 7.294416-3 1.000000-9 0.000000+07925 2151     \n"
    " 9.000000+4 2.114977+1 1.000000-9 7.261976-3 1.000000-9 0.000000+07925 2151     \n"
    " 9.200000+4 2.109050+1 1.000000-9 7.236028-3 1.000000-9 0.000000+07925 2151     \n"
    " 9.350000+4 2.104616+1 1.000000-9 7.216570-3 1.000000-9 0.000000+07925 2151     \n"
    " 9.500000+4 2.100192+1 1.000000-9 7.197115-3 1.000000-9 0.000000+07925 2151     \n"
    " 9.750000+4 2.092840+1 1.000000-9 7.164696-3 1.000000-9 0.000000+07925 2151     \n"
    " 1.000000+5 2.085515+1 1.000000-9 7.132285-3 1.000000-9 0.000000+07925 2151     \n"
    " 3.000000+0 0.000000+0          2          0        606        1007925 2151     \n"
    " 0.000000+0 0.000000+0 1.000000+0 2.000000+0 0.000000+0 0.000000+07925 2151     \n"
    " 2.000000+3 1.804535+1 0.000000+0 6.334971-3 1.292688-1 0.000000+07925 2151     \n"
    " 2.250000+3 1.803897+1 0.000000+0 6.332522-3 1.292801-1 0.000000+07925 2151     \n"
    " 2.500000+3 1.803259+1 0.000000+0 6.330073-3 1.292914-1 0.000000+07925 2151     \n"
    " 2.750000+3 1.802621+1 0.000000+0 6.327624-3 1.293027-1 0.000000+07925 2151     \n"
    " 3.000000+3 1.801984+1 0.000000+0 6.325176-3 1.293141-1 0.000000+07925 2151     \n"
    " 3.250000+3 1.801347+1 0.000000+0 6.322727-3 1.293254-1 0.000000+07925 2151     \n"
    " 3.500000+3 1.800710+1 0.000000+0 6.320279-3 1.293367-1 0.000000+07925 2151     \n"
    " 3.750000+3 1.800073+1 0.000000+0 6.317830-3 1.293480-1 0.000000+07925 2151     \n"
    " 4.000000+3 1.799436+1 0.000000+0 6.315382-3 1.293593-1 0.000000+07925 2151     \n"
    " 4.250000+3 1.798800+1 0.000000+0 6.312934-3 1.293706-1 0.000000+07925 2151     \n"
    " 4.500000+3 1.798164+1 0.000000+0 6.310485-3 1.293820-1 0.000000+07925 2151     \n"
    " 4.750000+3 1.797528+1 0.000000+0 6.308037-3 1.293933-1 0.000000+07925 2151     \n"
    " 5.000000+3 1.796893+1 0.000000+0 6.305589-3 1.294046-1 0.000000+07925 2151     \n"
    " 5.250000+3 1.796257+1 0.000000+0 6.303141-3 1.294159-1 0.000000+07925 2151     \n"
    " 5.500000+3 1.795622+1 0.000000+0 6.300693-3 1.294272-1 0.000000+07925 2151     \n"
    " 6.000000+3 1.794353+1 0.000000+0 6.295797-3 1.294498-1 0.000000+07925 2151     \n"
    " 6.500000+3 1.793084+1 0.000000+0 6.290901-3 1.294725-1 0.000000+07925 2151     \n"
    " 7.000000+3 1.791817+1 0.000000+0 6.286005-3 1.294951-1 0.000000+07925 2151     \n"
    " 7.500000+3 1.790550+1 0.000000+0 6.281110-3 1.295177-1 0.000000+07925 2151     \n"
    " 8.000000+3 1.789284+1 0.000000+0 6.276214-3 1.295404-1 0.000000+07925 2151     \n"
    " 8.500000+3 1.788019+1 0.000000+0 6.271319-3 1.295630-1 0.000000+07925 2151     \n"
    " 9.000000+3 1.786756+1 0.000000+0 6.266423-3 1.295856-1 0.000000+07925 2151     \n"
    " 9.500000+3 1.785493+1 0.000000+0 6.261528-3 1.296083-1 0.000000+07925 2151     \n"
    " 1.000000+4 1.784231+1 0.000000+0 6.256633-3 1.296309-1 0.000000+07925 2151     \n"
    " 1.100000+4 1.781709+1 0.000000+0 6.246843-3 1.296762-1 0.000000+07925 2151     \n"
    " 1.200000+4 1.779192+1 0.000000+0 6.237052-3 1.297214-1 0.000000+07925 2151     \n"
    " 1.300000+4 1.776678+1 0.000000+0 6.227262-3 1.297667-1 0.000000+07925 2151     \n"
    " 1.400000+4 1.774168+1 0.000000+0 6.217472-3 1.298120-1 0.000000+07925 2151     \n"
    " 1.500000+4 1.771661+1 0.000000+0 6.207682-3 1.298572-1 0.000000+07925 2151     \n"
    " 1.600000+4 1.769159+1 0.000000+0 6.197892-3 1.299025-1 0.000000+07925 2151     \n"
    " 1.700000+4 1.766660+1 0.000000+0 6.188102-3 1.299478-1 0.000000+07925 2151     \n"
    " 1.800000+4 1.764164+1 0.000000+0 6.178311-3 1.299930-1 0.000000+07925 2151     \n"
    " 1.900000+4 1.761673+1 0.000000+0 6.168520-3 1.300383-1 0.000000+07925 2151     \n"
    " 2.000000+4 1.759185+1 0.000000+0 6.158729-3 1.300835-1 0.000000+07925 2151     \n"
    " 2.100000+4 1.756700+1 0.000000+0 6.148938-3 1.301288-1 0.000000+07925 2151     \n"
    " 2.200000+4 1.754220+1 0.000000+0 6.139146-3 1.301741-1 0.000000+07925 2151     \n"
    " 2.300000+4 1.751743+1 0.000000+0 6.129354-3 1.302193-1 0.000000+07925 2151     \n"
    " 2.400000+4 1.749269+1 0.000000+0 6.119562-3 1.302646-1 0.000000+07925 2151     \n"
    " 2.500000+4 1.746799+1 0.000000+0 6.109769-3 1.303099-1 0.000000+07925 2151     \n"
    " 2.600000+4 1.744333+1 0.000000+0 6.099975-3 1.303551-1 0.000000+07925 2151     \n"
    " 2.700000+4 1.741871+1 0.000000+0 6.090181-3 1.304004-1 0.000000+07925 2151     \n"
    " 2.800000+4 1.739412+1 0.000000+0 6.080387-3 1.304457-1 0.000000+07925 2151     \n"
    " 2.900000+4 1.736957+1 0.000000+0 6.070592-3 1.304909-1 0.000000+07925 2151     \n"
    " 3.000000+4 1.734505+1 0.000000+0 6.060797-3 1.305362-1 0.000000+07925 2151     \n"
    " 3.100000+4 1.732057+1 0.000000+0 6.051001-3 1.305814-1 0.000000+07925 2151     \n"
    " 3.200000+4 1.729613+1 0.000000+0 6.041205-3 1.306267-1 0.000000+07925 2151     \n"
    " 3.300000+4 1.727172+1 0.000000+0 6.031408-3 1.306720-1 0.000000+07925 2151     \n"
    " 3.400000+4 1.724735+1 0.000000+0 6.021611-3 1.307172-1 0.000000+07925 2151     \n"
    " 3.500000+4 1.722301+1 0.000000+0 6.011813-3 1.307625-1 0.000000+07925 2151     \n"
    " 3.600000+4 1.719871+1 0.000000+0 6.002015-3 1.308078-1 0.000000+07925 2151     \n"
    " 3.700000+4 1.717444+1 0.000000+0 5.992216-3 1.308530-1 0.000000+07925 2151     \n"
    " 3.800000+4 1.715022+1 0.000000+0 5.982417-3 1.308983-1 0.000000+07925 2151     \n"
    " 3.900000+4 1.712602+1 0.000000+0 5.972617-3 1.309436-1 0.000000+07925 2151     \n"
    " 4.000000+4 1.710186+1 0.000000+0 5.962816-3 1.309888-1 0.000000+07925 2151     \n"
    " 4.100000+4 1.707774+1 0.000000+0 5.953015-3 1.310341-1 0.000000+07925 2151     \n"
    " 4.200000+4 1.705365+1 0.000000+0 5.943214-3 1.310794-1 0.000000+07925 2151     \n"
    " 4.300000+4 1.702960+1 0.000000+0 5.933412-3 1.311246-1 0.000000+07925 2151     \n"
    " 4.400000+4 1.700559+1 0.000000+0 5.923609-3 1.311699-1 0.000000+07925 2151     \n"
    " 4.500000+4 1.698161+1 0.000000+0 5.913806-3 1.312151-1 0.000000+07925 2151     \n"
    " 4.600000+4 1.695766+1 0.000000+0 5.904003-3 1.312604-1 0.000000+07925 2151     \n"
    " 4.700000+4 1.693375+1 0.000000+0 5.894199-3 1.313057-1 0.000000+07925 2151     \n"
    " 4.800000+4 1.690988+1 0.000000+0 5.884395-3 1.313509-1 0.000000+07925 2151     \n"
    " 5.000000+4 1.686223+1 0.000000+0 5.864785-3 1.314415-1 0.000000+07925 2151     \n"
    " 5.250000+4 1.680287+1 0.000000+0 5.840271-3 1.315546-1 0.000000+07925 2151     \n"
    " 5.500000+4 1.674373+1 0.000000+0 5.815754-3 1.316678-1 0.000000+07925 2151     \n"
    " 5.750000+4 1.668481+1 0.000000+0 5.791236-3 1.317809-1 0.000000+07925 2151     \n"
    " 6.000000+4 1.662610+1 0.000000+0 5.766716-3 1.318941-1 0.000000+07925 2151     \n"
    " 6.250000+4 1.656760+1 0.000000+0 5.742195-3 1.320073-1 0.000000+07925 2151     \n"
    " 6.500000+4 1.650933+1 0.000000+0 5.717673-3 1.321204-1 0.000000+07925 2151     \n"
    " 6.750000+4 1.645126+1 0.000000+0 5.693151-3 1.322336-1 0.000000+07925 2151     \n"
    " 7.000000+4 1.639341+1 0.000000+0 5.668629-3 1.323467-1 0.000000+07925 2151     \n"
    " 7.200000+4 1.634728+1 0.000000+0 5.649012-3 1.324373-1 0.000000+07925 2151     \n"
    " 7.400000+4 1.630129+1 0.000000+0 5.629395-3 1.325278-1 0.000000+07925 2151     \n"
    " 7.500000+4 1.627834+1 0.000000+0 5.619587-3 1.325731-1 0.000000+07925 2151     \n"
    " 7.775000+4 1.621541+1 0.000000+0 5.592617-3 1.326975-1 0.000000+07925 2151     \n"
    " 7.779600+4 1.621436+1 0.000000+0 5.592166-3 1.326996-1 0.000000+07925 2151     \n"
    " 7.781000+4 1.621404+1 0.000000+0 5.592029-3 1.327003-1 0.000000+07925 2151     \n"
    " 7.783000+4 1.621359+1 0.000000+0 5.591832-3 1.327012-1 0.000000+07925 2151     \n"
    " 7.786000+4 1.621290+1 0.000000+0 5.591538-3 1.327025-1 0.000000+07925 2151     \n"
    " 7.789000+4 1.621222+1 0.000000+0 5.591244-3 1.327039-1 0.000000+07925 2151     \n"
    " 7.793000+4 1.621130+1 1.202129-9 5.590852-3 1.327057-1 0.000000+07925 2151     \n"
    " 7.798000+4 1.621016+1 2.655696-9 5.590361-3 1.327080-1 0.000000+07925 2151     \n"
    " 7.807000+4 1.620811+1 7.184735-9 5.589479-3 1.327120-1 0.000000+07925 2151     \n"
    " 7.820000+4 1.620514+1 1.896025-8 5.588204-3 1.327179-1 0.000000+07925 2151     \n"
    " 7.835000+4 1.620172+1 4.173490-8 5.586733-3 1.327247-1 0.000000+07925 2151     \n"
    " 7.862000+4 1.619556+1 1.125238-7 5.584085-3 1.327369-1 0.000000+07925 2151     \n"
    " 7.890000+4 1.618917+1 2.336358-7 5.581339-3 1.327496-1 0.000000+07925 2151     \n"
    " 7.945000+4 1.617664+1 6.409837-7 5.575946-3 1.327745-1 0.000000+07925 2151     \n"
    " 8.000000+4 1.616412+1 1.311983-6 5.570552-3 1.327994-1 0.000000+07925 2151     \n"
    " 8.100000+4 1.614137+1 3.334435-6 5.560747-3 1.328447-1 0.000000+07925 2151     \n"
    " 8.200000+4 1.611866+1 6.559006-6 5.550941-3 1.328899-1 0.000000+07925 2151     \n"
    " 8.300000+4 1.609598+1 1.115356-5 5.541136-3 1.329352-1 0.000000+07925 2151     \n"
    " 8.500000+4 1.605073+1 2.501946-5 5.521528-3 1.330257-1 0.000000+07925 2151     \n"
    " 8.750000+4 1.599434+1 5.235151-5 5.497022-3 1.331389-1 0.000000+07925 2151     \n"
    " 9.000000+4 1.593816+1 9.225868-5 5.472520-3 1.332520-1 0.000000+07925 2151     \n"
    " 9.200000+4 1.589337+1 1.341317-4 5.452922-3 1.333426-1 0.000000+07925 2151     \n"
    " 9.350000+4 1.585986+1 1.717280-4 5.438226-3 1.334105-1 0.000000+07925 2151     \n"
    " 9.500000+4 1.582643+1 2.148801-4 5.423533-3 1.334783-1 0.000000+07925 2151     \n"
    " 9.750000+4 1.577086+1 2.996953-4 5.399048-3 1.335915-1 0.000000+07925 2151     \n"
    " 1.000000+5 1.571550+1 4.013598-4 5.374571-3 1.337047-1 0.000000+07925 2151     \n"
    " 4.000000+0 0.000000+0          2          0        606        1007925 2151     \n"
    " 0.000000+0 0.000000+0 2.000000+0 1.000000+0 0.000000+0 0.000000+07925 2151     \n"
    " 2.000000+3 1.508143+1 0.000000+0 5.294460-3 1.292688-1 0.000000+07925 2151     \n"
    " 2.250000+3 1.507607+1 0.000000+0 5.292406-3 1.292801-1 0.000000+07925 2151     \n"
    " 2.500000+3 1.507072+1 0.000000+0 5.290352-3 1.292914-1 0.000000+07925 2151     \n"
    " 2.750000+3 1.506537+1 0.000000+0 5.288299-3 1.293027-1 0.000000+07925 2151     \n"
    " 3.000000+3 1.506002+1 0.000000+0 5.286245-3 1.293141-1 0.000000+07925 2151     \n"
    " 3.250000+3 1.505468+1 0.000000+0 5.284191-3 1.293254-1 0.000000+07925 2151     \n"
    " 3.500000+3 1.504933+1 0.000000+0 5.282138-3 1.293367-1 0.000000+07925 2151     \n"
    " 3.750000+3 1.504399+1 0.000000+0 5.280084-3 1.293480-1 0.000000+07925 2151     \n"
    " 4.000000+3 1.503865+1 0.000000+0 5.278031-3 1.293593-1 0.000000+07925 2151     \n"
    " 4.250000+3 1.503331+1 0.000000+0 5.275977-3 1.293706-1 0.000000+07925 2151     \n"
    " 4.500000+3 1.502797+1 0.000000+0 5.273924-3 1.293820-1 0.000000+07925 2151     \n"
    " 4.750000+3 1.502264+1 0.000000+0 5.271871-3 1.293933-1 0.000000+07925 2151     \n"
    " 5.000000+3 1.501731+1 0.000000+0 5.269818-3 1.294046-1 0.000000+07925 2151     \n"
    " 5.250000+3 1.501198+1 0.000000+0 5.267764-3 1.294159-1 0.000000+07925 2151     \n"
    " 5.500000+3 1.500665+1 0.000000+0 5.265711-3 1.294272-1 0.000000+07925 2151     \n"
    " 6.000000+3 1.499600+1 0.000000+0 5.261605-3 1.294498-1 0.000000+07925 2151     \n"
    " 6.500000+3 1.498536+1 0.000000+0 5.257499-3 1.294725-1 0.000000+07925 2151     \n"
    " 7.000000+3 1.497472+1 0.000000+0 5.253393-3 1.294951-1 0.000000+07925 2151     \n"
    " 7.500000+3 1.496409+1 0.000000+0 5.249287-3 1.295177-1 0.000000+07925 2151     \n"
    " 8.000000+3 1.495348+1 0.000000+0 5.245182-3 1.295404-1 0.000000+07925 2151     \n"
    " 8.500000+3 1.494286+1 0.000000+0 5.241076-3 1.295630-1 0.000000+07925 2151     \n"
    " 9.000000+3 1.493226+1 0.000000+0 5.236971-3 1.295856-1 0.000000+07925 2151     \n"
    " 9.500000+3 1.492167+1 0.000000+0 5.232865-3 1.296083-1 0.000000+07925 2151     \n"
    " 1.000000+4 1.491108+1 0.000000+0 5.228760-3 1.296309-1 0.000000+07925 2151     \n"
    " 1.100000+4 1.488993+1 0.000000+0 5.220550-3 1.296762-1 0.000000+07925 2151     \n"
    " 1.200000+4 1.486880+1 0.000000+0 5.212339-3 1.297214-1 0.000000+07925 2151     \n"
    " 1.300000+4 1.484772+1 0.000000+0 5.204129-3 1.297667-1 0.000000+07925 2151     \n"
    " 1.400000+4 1.482666+1 0.000000+0 5.195919-3 1.298120-1 0.000000+07925 2151     \n"
    " 1.500000+4 1.480563+1 0.000000+0 5.187709-3 1.298572-1 0.000000+07925 2151     \n"
    " 1.600000+4 1.478463+1 0.000000+0 5.179500-3 1.299025-1 0.000000+07925 2151     \n"
    " 1.700000+4 1.476367+1 0.000000+0 5.171290-3 1.299478-1 0.000000+07925 2151     \n"
    " 1.800000+4 1.474274+1 0.000000+0 5.163080-3 1.299930-1 0.000000+07925 2151     \n"
    " 1.900000+4 1.472183+1 0.000000+0 5.154870-3 1.300383-1 0.000000+07925 2151     \n"
    " 2.000000+4 1.470096+1 0.000000+0 5.146659-3 1.300835-1 0.000000+07925 2151     \n"
    " 2.100000+4 1.468012+1 0.000000+0 5.138449-3 1.301288-1 0.000000+07925 2151     \n"
    " 2.200000+4 1.465931+1 0.000000+0 5.130238-3 1.301741-1 0.000000+07925 2151     \n"
    " 2.300000+4 1.463853+1 0.000000+0 5.122027-3 1.302193-1 0.000000+07925 2151     \n"
    " 2.400000+4 1.461778+1 0.000000+0 5.113816-3 1.302646-1 0.000000+07925 2151     \n"
    " 2.500000+4 1.459706+1 0.000000+0 5.105605-3 1.303099-1 0.000000+07925 2151     \n"
    " 2.600000+4 1.457638+1 0.000000+0 5.097394-3 1.303551-1 0.000000+07925 2151     \n"
    " 2.700000+4 1.455572+1 0.000000+0 5.089182-3 1.304004-1 0.000000+07925 2151     \n"
    " 2.800000+4 1.453509+1 0.000000+0 5.080970-3 1.304457-1 0.000000+07925 2151     \n"
    " 2.900000+4 1.451450+1 0.000000+0 5.072757-3 1.304909-1 0.000000+07925 2151     \n"
    " 3.000000+4 1.449393+1 0.000000+0 5.064544-3 1.305362-1 0.000000+07925 2151     \n"
    " 3.100000+4 1.447340+1 0.000000+0 5.056331-3 1.305814-1 0.000000+07925 2151     \n"
    " 3.200000+4 1.445289+1 0.000000+0 5.048118-3 1.306267-1 0.000000+07925 2151     \n"
    " 3.300000+4 1.443242+1 0.000000+0 5.039904-3 1.306720-1 0.000000+07925 2151     \n"
    " 3.400000+4 1.441197+1 0.000000+0 5.031690-3 1.307172-1 0.000000+07925 2151     \n"
    " 3.500000+4 1.439156+1 0.000000+0 5.023475-3 1.307625-1 0.000000+07925 2151     \n"
    " 3.600000+4 1.437118+1 0.000000+0 5.015261-3 1.308078-1 0.000000+07925 2151     \n"
    " 3.700000+4 1.435082+1 0.000000+0 5.007046-3 1.308530-1 0.000000+07925 2151     \n"
    " 3.800000+4 1.433050+1 0.000000+0 4.998830-3 1.308983-1 0.000000+07925 2151     \n"
    " 3.900000+4 1.431020+1 0.000000+0 4.990614-3 1.309436-1 0.000000+07925 2151     \n"
    " 4.000000+4 1.428994+1 0.000000+0 4.982398-3 1.309888-1 0.000000+07925 2151     \n"
    " 4.100000+4 1.426971+1 0.000000+0 4.974182-3 1.310341-1 0.000000+07925 2151     \n"
    " 4.200000+4 1.424950+1 0.000000+0 4.965965-3 1.310794-1 0.000000+07925 2151     \n"
    " 4.300000+4 1.422933+1 0.000000+0 4.957748-3 1.311246-1 0.000000+07925 2151     \n"
    " 4.400000+4 1.420919+1 0.000000+0 4.949530-3 1.311699-1 0.000000+07925 2151     \n"
    " 4.500000+4 1.418907+1 0.000000+0 4.941312-3 1.312151-1 0.000000+07925 2151     \n"
    " 4.600000+4 1.416899+1 0.000000+0 4.933094-3 1.312604-1 0.000000+07925 2151     \n"
    " 4.700000+4 1.414893+1 0.000000+0 4.924876-3 1.313057-1 0.000000+07925 2151     \n"
    " 4.800000+4 1.412891+1 0.000000+0 4.916657-3 1.313509-1 0.000000+07925 2151     \n"
    " 5.000000+4 1.408894+1 0.000000+0 4.900219-3 1.314415-1 0.000000+07925 2151     \n"
    " 5.250000+4 1.403916+1 0.000000+0 4.879671-3 1.315546-1 0.000000+07925 2151     \n"
    " 5.500000+4 1.398955+1 0.000000+0 4.859121-3 1.316678-1 0.000000+07925 2151     \n"
    " 5.750000+4 1.394013+1 0.000000+0 4.838570-3 1.317809-1 0.000000+07925 2151     \n"
    " 6.000000+4 1.389089+1 0.000000+0 4.818018-3 1.318941-1 0.000000+07925 2151     \n"
    " 6.250000+4 1.384184+1 0.000000+0 4.797466-3 1.320073-1 0.000000+07925 2151     \n"
    " 6.500000+4 1.379296+1 0.000000+0 4.776914-3 1.321204-1 0.000000+07925 2151     \n"
    " 6.750000+4 1.374426+1 0.000000+0 4.756363-3 1.322336-1 0.000000+07925 2151     \n"
    " 7.000000+4 1.369575+1 0.000000+0 4.735812-3 1.323467-1 0.000000+07925 2151     \n"
    " 7.200000+4 1.365706+1 0.000000+0 4.719372-3 1.324373-1 0.000000+07925 2151     \n"
    " 7.400000+4 1.361849+1 0.000000+0 4.702933-3 1.325278-1 0.000000+07925 2151     \n"
    " 7.500000+4 1.359925+1 0.000000+0 4.694714-3 1.325731-1 0.000000+07925 2151     \n"
    " 7.775000+4 1.354648+1 0.000000+0 4.672113-3 1.326975-1 0.000000+07925 2151     \n"
    " 7.779600+4 1.354559+1 0.000000+0 4.671735-3 1.326996-1 0.000000+07925 2151     \n"
    " 7.781000+4 1.354533+1 0.000000+0 4.671620-3 1.327003-1 0.000000+07925 2151     \n"
    " 7.783000+4 1.354494+1 0.000000+0 4.671456-3 1.327012-1 0.000000+07925 2151     \n"
    " 7.786000+4 1.354437+1 0.000000+0 4.671209-3 1.327025-1 0.000000+07925 2151     \n"
    " 7.789000+4 1.354379+1 0.000000+0 4.670963-3 1.327039-1 0.000000+07925 2151     \n"
    " 7.793000+4 1.354303+1 0.000000+0 4.670634-3 1.327057-1 0.000000+07925 2151     \n"
    " 7.798000+4 1.354207+1 0.000000+0 4.670223-3 1.327080-1 0.000000+07925 2151     \n"
    " 7.807000+4 1.354035+1 0.000000+0 4.669483-3 1.327120-1 0.000000+07925 2151     \n"
    " 7.820000+4 1.353786+1 0.000000+0 4.668415-3 1.327179-1 0.000000+07925 2151     \n"
    " 7.835000+4 1.353499+1 0.000000+0 4.667183-3 1.327247-1 0.000000+07925 2151     \n"
    " 7.862000+4 1.352983+1 0.000000+0 4.664964-3 1.327369-1 0.000000+07925 2151     \n"
    " 7.890000+4 1.352447+1 0.000000+0 4.662663-3 1.327496-1 0.000000+07925 2151     \n"
    " 7.945000+4 1.351396+1 0.000000+0 4.658143-3 1.327745-1 0.000000+07925 2151     \n"
    " 8.000000+4 1.350346+1 0.000000+0 4.653624-3 1.327994-1 0.000000+07925 2151     \n"
    " 8.100000+4 1.348438+1 0.000000+0 4.645407-3 1.328447-1 0.000000+07925 2151     \n"
    " 8.200000+4 1.346534+1 0.000000+0 4.637191-3 1.328899-1 0.000000+07925 2151     \n"
    " 8.300000+4 1.344632+1 0.000000+0 4.628975-3 1.329352-1 0.000000+07925 2151     \n"
    " 8.500000+4 1.340837+1 0.000000+0 4.612545-3 1.330257-1 0.000000+07925 2151     \n"
    " 8.750000+4 1.336109+1 0.000000+0 4.592012-3 1.331389-1 0.000000+07925 2151     \n"
    " 9.000000+4 1.331398+1 0.000000+0 4.571482-3 1.332520-1 0.000000+07925 2151     \n"
    " 9.200000+4 1.327642+1 0.000000+0 4.555062-3 1.333426-1 0.000000+07925 2151     \n"
    " 9.350000+4 1.324832+1 0.000000+0 4.542749-3 1.334105-1 0.000000+07925 2151     \n"
    " 9.500000+4 1.322029+1 0.000000+0 4.530439-3 1.334783-1 0.000000+07925 2151     \n"
    " 9.750000+4 1.317370+1 0.000000+0 4.509926-3 1.335915-1 0.000000+07925 2151     \n"
    " 1.000000+5 1.312728+1 0.000000+0 4.489419-3 1.337047-1 0.000000+07925 2151     \n"
    "                                                                  7925 2  0     \n" );
}

std::string Ag107() {

  // Ag107 ENDF/B-VIII.0 LRU=2 resonance evaluation

  return
    " 4.710700+4 1.059870+2          0          0          1          04725 2151     \n"
    " 4.710700+4 1.000000+0          0          0          1          04725 2151     \n"
    " 6.500000+3 1.000000+5          2          2          0          04725 2151     \n"
    " 5.000000-1 6.600000-1          0          0          3          04725 2151     \n"
    " 1.059870+2 0.000000+0          0          0          2          04725 2151     \n"
    " 0.000000+0 0.000000+0          2          0         84         134725 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 1.000000+0 0.000000+0 0.000000+04725 2151     \n"
    " 6.500000+3 9.011400+1 0.000000+0 4.649900-3 1.400000-1 0.000000+04725 2151     \n"
    " 7.400000+3 9.126400+1 0.000000+0 4.709200-3 1.400000-1 0.000000+04725 2151     \n"
    " 1.067000+4 8.962300+1 0.000000+0 4.624500-3 1.400000-1 0.000000+04725 2151     \n"
    " 1.394000+4 7.685800+1 0.000000+0 3.965900-3 1.400000-1 0.000000+04725 2151     \n"
    " 1.721000+4 7.469800+1 0.000000+0 3.664100-3 1.400000-1 0.000000+04725 2151     \n"
    " 3.181000+4 8.210600+1 0.000000+0 3.582300-3 1.400000-1 0.000000+04725 2151     \n"
    " 4.641000+4 8.895200+1 0.000000+0 4.318800-3 1.400000-1 0.000000+04725 2151     \n"
    " 5.371000+4 8.746500+1 0.000000+0 4.077700-3 1.400000-1 0.000000+04725 2151     \n"
    " 6.100000+4 9.018500+1 0.000000+0 3.835200-3 1.400000-1 0.000000+04725 2151     \n"
    " 6.831000+4 9.117800+1 0.000000+0 3.136500-3 1.400000-1 0.000000+04725 2151     \n"
    " 7.561000+4 8.891500+1 0.000000+0 3.058700-3 1.400000-1 0.000000+04725 2151     \n"
    " 9.174000+4 7.875900+1 0.000000+0 2.709300-3 1.400000-1 0.000000+04725 2151     \n"
    " 1.000000+5 7.546000+1 0.000000+0 2.595800-3 1.400000-1 0.000000+04725 2151     \n"
    " 1.000000+0 0.000000+0          2          0         84         134725 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 1.000000+0 0.000000+0 0.000000+04725 2151     \n"
    " 6.500000+3 3.003800+1 0.000000+0 1.550000-3 1.400000-1 0.000000+04725 2151     \n"
    " 7.400000+3 3.042100+1 0.000000+0 1.569700-3 1.400000-1 0.000000+04725 2151     \n"
    " 1.067000+4 2.987400+1 0.000000+0 1.541500-3 1.400000-1 0.000000+04725 2151     \n"
    " 1.394000+4 2.561900+1 0.000000+0 1.322000-3 1.400000-1 0.000000+04725 2151     \n"
    " 1.721000+4 2.489900+1 0.000000+0 1.221400-3 1.400000-1 0.000000+04725 2151     \n"
    " 3.181000+4 2.736900+1 0.000000+0 1.194100-3 1.400000-1 0.000000+04725 2151     \n"
    " 4.641000+4 2.965100+1 0.000000+0 1.439600-3 1.400000-1 0.000000+04725 2151     \n"
    " 5.371000+4 2.915500+1 0.000000+0 1.359200-3 1.400000-1 0.000000+04725 2151     \n"
    " 6.100000+4 3.006200+1 0.000000+0 1.278400-3 1.400000-1 0.000000+04725 2151     \n"
    " 6.831000+4 3.039300+1 0.000000+0 1.045500-3 1.400000-1 0.000000+04725 2151     \n"
    " 7.561000+4 2.963800+1 0.000000+0 1.019600-3 1.400000-1 0.000000+04725 2151     \n"
    " 9.174000+4 2.625300+1 0.000000+0 9.031100-4 1.400000-1 0.000000+04725 2151     \n"
    " 1.000000+5 2.515300+1 0.000000+0 8.652800-4 1.400000-1 0.000000+04725 2151     \n"
    " 1.059870+2 0.000000+0          1          0          3          04725 2151     \n"
    " 0.000000+0 0.000000+0          2          0         84         134725 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 1.000000+0 0.000000+0 0.000000+04725 2151     \n"
    " 6.500000+3 9.011400+1 0.000000+0 4.006600-2 1.400000-1 0.000000+04725 2151     \n"
    " 7.400000+3 9.126400+1 0.000000+0 4.112000-2 1.400000-1 0.000000+04725 2151     \n"
    " 1.067000+4 8.962300+1 0.000000+0 3.871400-2 1.400000-1 0.000000+04725 2151     \n"
    " 1.394000+4 7.685800+1 0.000000+0 3.018300-2 1.400000-1 0.000000+04725 2151     \n"
    " 1.721000+4 7.469800+1 0.000000+0 2.741600-2 1.400000-1 0.000000+04725 2151     \n"
    " 3.181000+4 8.210600+1 0.000000+0 3.210900-2 1.400000-1 0.000000+04725 2151     \n"
    " 4.641000+4 8.895200+1 0.000000+0 3.301900-2 1.400000-1 0.000000+04725 2151     \n"
    " 5.371000+4 8.746500+1 0.000000+0 3.257900-2 1.400000-1 0.000000+04725 2151     \n"
    " 6.100000+4 9.018500+1 0.000000+0 3.369600-2 1.400000-1 0.000000+04725 2151     \n"
    " 6.831000+4 9.117800+1 0.000000+0 3.529600-2 1.400000-1 0.000000+04725 2151     \n"
    " 7.561000+4 8.891500+1 0.000000+0 3.346500-2 1.400000-1 0.000000+04725 2151     \n"
    " 9.174000+4 7.875900+1 0.000000+0 2.785700-2 1.400000-1 0.000000+04725 2151     \n"
    " 1.000000+5 7.546000+1 0.000000+0 2.638600-2 1.400000-1 0.000000+04725 2151     \n"
    " 1.000000+0 0.000000+0          2          0         84         134725 2151     \n"
    " 0.000000+0 0.000000+0 1.000000+0 2.000000+0 0.000000+0 0.000000+04725 2151     \n"
    " 6.500000+3 3.003800+1 0.000000+0 1.335500-2 1.400000-1 0.000000+04725 2151     \n"
    " 7.400000+3 3.042100+1 0.000000+0 1.370700-2 1.400000-1 0.000000+04725 2151     \n"
    " 1.067000+4 2.987400+1 0.000000+0 1.290500-2 1.400000-1 0.000000+04725 2151     \n"
    " 1.394000+4 2.561900+1 0.000000+0 1.006100-2 1.400000-1 0.000000+04725 2151     \n"
    " 1.721000+4 2.489900+1 0.000000+0 9.138500-3 1.400000-1 0.000000+04725 2151     \n"
    " 3.181000+4 2.736900+1 0.000000+0 1.070300-2 1.400000-1 0.000000+04725 2151     \n"
    " 4.641000+4 2.965100+1 0.000000+0 1.100600-2 1.400000-1 0.000000+04725 2151     \n"
    " 5.371000+4 2.915500+1 0.000000+0 1.086000-2 1.400000-1 0.000000+04725 2151     \n"
    " 6.100000+4 3.006200+1 0.000000+0 1.123200-2 1.400000-1 0.000000+04725 2151     \n"
    " 6.831000+4 3.039300+1 0.000000+0 1.176500-2 1.400000-1 0.000000+04725 2151     \n"
    " 7.561000+4 2.963800+1 0.000000+0 1.115500-2 1.400000-1 0.000000+04725 2151     \n"
    " 9.174000+4 2.625300+1 0.000000+0 9.285600-3 1.400000-1 0.000000+04725 2151     \n"
    " 1.000000+5 2.515300+1 1.788900-6 8.795300-3 1.400000-1 0.000000+04725 2151     \n"
    " 2.000000+0 0.000000+0          2          0         84         134725 2151     \n"
    " 0.000000+0 0.000000+0 2.000000+0 1.000000+0 0.000000+0 0.000000+04725 2151     \n"
    " 6.500000+3 1.802300+1 0.000000+0 8.013200-3 1.400000-1 0.000000+04725 2151     \n"
    " 7.400000+3 1.825300+1 0.000000+0 8.223900-3 1.400000-1 0.000000+04725 2151     \n"
    " 1.067000+4 1.792500+1 0.000000+0 7.742900-3 1.400000-1 0.000000+04725 2151     \n"
    " 1.394000+4 1.537200+1 0.000000+0 6.036500-3 1.400000-1 0.000000+04725 2151     \n"
    " 1.721000+4 1.494000+1 0.000000+0 5.483100-3 1.400000-1 0.000000+04725 2151     \n"
    " 3.181000+4 1.642100+1 0.000000+0 6.421800-3 1.400000-1 0.000000+04725 2151     \n"
    " 4.641000+4 1.779000+1 0.000000+0 6.603700-3 1.400000-1 0.000000+04725 2151     \n"
    " 5.371000+4 1.749300+1 0.000000+0 6.515700-3 1.400000-1 0.000000+04725 2151     \n"
    " 6.100000+4 1.803700+1 0.000000+0 6.739200-3 1.400000-1 0.000000+04725 2151     \n"
    " 6.831000+4 1.823600+1 0.000000+0 7.059100-3 1.400000-1 0.000000+04725 2151     \n"
    " 7.561000+4 1.778300+1 0.000000+0 6.693000-3 1.400000-1 0.000000+04725 2151     \n"
    " 9.174000+4 1.575200+1 0.000000+0 5.571400-3 1.400000-1 0.000000+04725 2151     \n"
    " 1.000000+5 1.509200+1 2.146700-6 5.277200-3 1.400000-1 0.000000+04725 2151     \n"
    " 1.059870+2 0.000000+0          2          0          3          04725 2151     \n"
    " 1.000000+0 0.000000+0          2          0         84         134725 2151     \n"
    " 0.000000+0 0.000000+0 0.000000+0 1.000000+0 0.000000+0 0.000000+04725 2151     \n"
    " 6.500000+3 3.003800+1 0.000000+0 1.592000-3 1.400000-1 0.000000+04725 2151     \n"
    " 7.400000+3 3.042100+1 0.000000+0 1.612300-3 1.400000-1 0.000000+04725 2151     \n"
    " 1.067000+4 2.987400+1 0.000000+0 1.583300-3 1.400000-1 0.000000+04725 2151     \n"
    " 1.394000+4 2.561900+1 0.000000+0 1.357800-3 1.400000-1 0.000000+04725 2151     \n"
    " 1.721000+4 2.489900+1 0.000000+0 1.319700-3 1.400000-1 0.000000+04725 2151     \n"
    " 3.181000+4 2.736900+1 0.000000+0 1.450500-3 1.400000-1 0.000000+04725 2151     \n"
    " 4.641000+4 2.965100+1 0.000000+0 1.571500-3 1.400000-1 0.000000+04725 2151     \n"
    " 5.371000+4 2.915500+1 0.000000+0 1.545200-3 1.400000-1 0.000000+04725 2151     \n"
    " 6.100000+4 3.006200+1 0.000000+0 1.593300-3 1.400000-1 0.000000+04725 2151     \n"
    " 6.831000+4 3.039300+1 0.000000+0 1.610800-3 1.400000-1 0.000000+04725 2151     \n"
    " 7.561000+4 2.963800+1 0.000000+0 1.570800-3 1.400000-1 0.000000+04725 2151     \n"
    " 9.174000+4 2.625300+1 0.000000+0 1.391400-3 1.400000-1 0.000000+04725 2151     \n"
    " 1.000000+5 2.515300+1 0.000000+0 1.333100-3 1.400000-1 0.000000+04725 2151     \n"
    " 2.000000+0 0.000000+0          2          0         84         134725 2151     \n"
    " 0.000000+0 0.000000+0 1.000000+0 2.000000+0 0.000000+0 0.000000+04725 2151     \n"
    " 6.500000+3 1.802300+1 0.000000+0 9.552100-4 1.400000-1 0.000000+04725 2151     \n"
    " 7.400000+3 1.825300+1 0.000000+0 9.674000-4 1.400000-1 0.000000+04725 2151     \n"
    " 1.067000+4 1.792500+1 0.000000+0 9.500000-4 1.400000-1 0.000000+04725 2151     \n"
    " 1.394000+4 1.537200+1 0.000000+0 8.146900-4 1.400000-1 0.000000+04725 2151     \n"
    " 1.721000+4 1.494000+1 0.000000+0 7.918000-4 1.400000-1 0.000000+04725 2151     \n"
    " 3.181000+4 1.642100+1 0.000000+0 8.703200-4 1.400000-1 0.000000+04725 2151     \n"
    " 4.641000+4 1.779000+1 0.000000+0 9.428900-4 1.400000-1 0.000000+04725 2151     \n"
    " 5.371000+4 1.749300+1 0.000000+0 9.271300-4 1.400000-1 0.000000+04725 2151     \n"
    " 6.100000+4 1.803700+1 0.000000+0 9.559600-4 1.400000-1 0.000000+04725 2151     \n"
    " 6.831000+4 1.823600+1 0.000000+0 9.664900-4 1.400000-1 0.000000+04725 2151     \n"
    " 7.561000+4 1.778300+1 0.000000+0 9.425000-4 1.400000-1 0.000000+04725 2151     \n"
    " 9.174000+4 1.575200+1 0.000000+0 8.348500-4 1.400000-1 0.000000+04725 2151     \n"
    " 1.000000+5 1.509200+1 5.056200-3 7.998800-4 1.400000-1 0.000000+04725 2151     \n"
    " 3.000000+0 0.000000+0          2          0         84         134725 2151     \n"
    " 0.000000+0 0.000000+0 2.000000+0 1.000000+0 0.000000+0 0.000000+04725 2151     \n"
    " 6.500000+3 1.287300+1 0.000000+0 6.822900-4 1.400000-1 0.000000+04725 2151     \n"
    " 7.400000+3 1.303800+1 0.000000+0 6.910000-4 1.400000-1 0.000000+04725 2151     \n"
    " 1.067000+4 1.280300+1 0.000000+0 6.785700-4 1.400000-1 0.000000+04725 2151     \n"
    " 1.394000+4 1.098000+1 0.000000+0 5.819200-4 1.400000-1 0.000000+04725 2151     \n"
    " 1.721000+4 1.067100+1 0.000000+0 5.655700-4 1.400000-1 0.000000+04725 2151     \n"
    " 3.181000+4 1.172900+1 0.000000+0 6.216600-4 1.400000-1 0.000000+04725 2151     \n"
    " 4.641000+4 1.270700+1 0.000000+0 6.734900-4 1.400000-1 0.000000+04725 2151     \n"
    " 5.371000+4 1.249500+1 0.000000+0 6.622400-4 1.400000-1 0.000000+04725 2151     \n"
    " 6.100000+4 1.288400+1 0.000000+0 6.828300-4 1.400000-1 0.000000+04725 2151     \n"
    " 6.831000+4 1.302500+1 0.000000+0 6.903500-4 1.400000-1 0.000000+04725 2151     \n"
    " 7.561000+4 1.270200+1 0.000000+0 6.732200-4 1.400000-1 0.000000+04725 2151     \n"
    " 9.174000+4 1.125100+1 0.000000+0 5.963200-4 1.400000-1 0.000000+04725 2151     \n"
    " 1.000000+5 1.078000+1 7.223100-3 5.713400-4 1.400000-1 0.000000+04725 2151     \n"
    " 0.000000+0 0.000000+0          0          0          0          04725 2  0     \n";
}
