SCENARIO("channel radius"){
  using namespace dimwits;

  std::array< double, 4 > awri = {{ 9.991673E-1,    // hydrogen-1
                                    1.585751E+1,    // oxygen-16
                                    5.545400E+1,    // iron-56
                                    2.360058E+2 }}; // uranium-238

  std::array< double, 4 > reference = {{ 0.20331999182508331,
                                         0.38990773394347489,
                                         0.5503975406881233,
                                         0.84230419293656023 }};

  auto trial =
    awri
    | ranges::view::transform
      ( []( double d ){ return channelRadius(d)( 1.0 * electronVolts ).value; } );

  RANGES_FOR( const auto pair, ranges::view::zip( trial, reference ) ){
    const auto trial = std::get<0>(pair);
    const auto reference = std::get<1>(pair);
    REQUIRE( trial == Approx( reference ) );
  }
}
