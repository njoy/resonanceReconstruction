std::string Na22();
std::string Pu239();

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

      auto compoundsystem = std::get< legacy::unresolved::CompoundSystem >( resonances.compoundSystem() );

      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // content verification
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      // spin groups
      auto spingroups = compoundsystem.spinGroups();
      CHECK( 5 == spingroups.size() );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // spin group 0
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      auto spingroup0 = spingroups[0];
    } // THEN

    THEN( "cross sections can be reconstructed" ) {

      ReactionID elas( "n,Pu239->n,Pu239" );
      ReactionID capt( "n,Pu239->capture" );
      ReactionID fiss( "n,Pu239->fission" );
      std::map< ReactionID, CrossSection > xs;

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
