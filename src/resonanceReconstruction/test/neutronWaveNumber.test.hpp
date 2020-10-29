SCENARIO("wave number"){
  using namespace dimwits;

  std::array< double, 4 > awri = {{ 9.991673E-1,    // hydrogen-1
                                    1.585751E+1,    // oxygen-16
                                    5.545400E+1,    // iron-56
                                    2.360058E+2 }}; // uranium-238

  auto energies =
    ranges::view::linear_distribute( -10., 30., 40 )
    | ranges::view::transform( [] ( double d ){ return std::pow( 2., d ); } )
    | ranges::view::transform( [] ( double d ){ return d * electronVolts; } );

  auto trial =
    ranges::view::cartesian_product( awri, energies )
    | ranges::view::transform
      ( []( auto&& tuple ){
          return neutronWaveNumber( std::get<0>(tuple) )
                                  ( std::get<1>(tuple) ).value; } );

  std::vector< double > reference =
    { 3.43108231e-05,   4.89559512e-05,   6.98521616e-05,
      9.96676473e-05,   1.42209485e-04,   2.02909752e-04,
      2.89519138e-04,   4.13096613e-04,   5.89421523e-04,
      8.41008426e-04,   1.19998192e-03,   1.71217858e-03,
      2.44299972e-03,   3.48576235e-03,   4.97361462e-03,
      7.09653727e-03,   1.01256018e-02,   1.44475831e-02,
      2.06143458e-02,   2.94133110e-02,   4.19679998e-02,
      5.98814941e-02,   8.54411303e-02,   1.21910564e-01,
      1.73946502e-01,   2.48193301e-01,   3.54131380e-01,
      5.05287749e-01,   7.20963246e-01,   1.02869702e+00,
      1.46778295e+00,   2.09428699e+00,   2.98820613e+00,
      4.26368301e+00,   6.08358060e+00,   8.68027777e+00,
      1.23853413e+01,   1.76718630e+01,   2.52148677e+01,
      3.59775056e+01,   6.45778574e-05,   9.21420750e-05,
      1.31471720e-04,   1.87588712e-04,   2.67658512e-04,
      3.81905062e-04,   5.44916266e-04,   7.77506680e-04,
      1.10937528e-03,   1.58289768e-03,   2.25853694e-03,
      3.22256403e-03,   4.59807354e-03,   6.56070136e-03,
      9.36105132e-03,   1.33566942e-02,   1.90578252e-02,
      2.71924098e-02,   3.87991359e-02,   5.53600419e-02,
      7.89897550e-02,   1.12705503e-01,   1.60812380e-01,
      2.29453051e-01,   3.27392098e-01,   4.67135153e-01,
      6.66525711e-01,   9.51023533e-01,   1.35695555e+00,
      1.93615435e+00,   2.76257663e+00,   3.94174650e+00,
      5.62422968e+00,   8.02485889e+00,   1.14501654e+01,
      1.63375195e+01,   2.33109769e+01,   3.32609639e+01,
      4.74579733e+01,   6.77147914e+01,   6.74342020e-05,
      9.62176132e-05,   1.37286849e-04,   1.95885953e-04,
      2.79497321e-04,   3.98797114e-04,   5.69018470e-04,
      8.11896595e-04,   1.15844409e-03,   1.65291086e-03,
      2.35843434e-03,   3.36510134e-03,   4.80145103e-03,
      6.85088789e-03,   9.77510019e-03,   1.39474744e-02,
      1.99007723e-02,   2.83951579e-02,   4.05152613e-02,
      5.78086731e-02,   8.24835526e-02,   1.17690583e-01,
      1.67925276e-01,   2.39601994e-01,   3.41872986e-01,
      4.87797019e-01,   6.96006825e-01,   9.93088277e-01,
      1.41697508e+00,   2.02179244e+00,   2.88476822e+00,
      4.11609398e+00,   5.87299510e+00,   8.37980660e+00,
      1.19566179e+01,   1.70601446e+01,   2.43420453e+01,
      3.47321304e+01,   4.95570880e+01,   7.07098858e+01,
      6.83605842e-05,   9.75394096e-05,   1.39172837e-04,
      1.98576950e-04,   2.83336936e-04,   4.04275618e-04,
      5.76835402e-04,   8.23050083e-04,   1.17435829e-03,
      1.67561784e-03,   2.39083350e-03,   3.41132966e-03,
      4.86741131e-03,   6.94500245e-03,   9.90938633e-03,
      1.41390789e-02,   2.01741606e-02,   2.87852384e-02,
      4.10718426e-02,   5.86028239e-02,   8.36166763e-02,
      1.19307366e-01,   1.70232162e-01,   2.42893544e-01,
      3.46569490e-01,   4.94498165e-01,   7.05568269e-01,
      1.00673090e+00,   1.43644087e+00,   2.04956696e+00,
      2.92439793e+00,   4.17263911e+00,   5.95367580e+00,
      8.49492479e+00,   1.21208728e+01,   1.72945096e+01,
      2.46764459e+01,   3.52092655e+01,   5.02378821e+01,
      7.16812680e+01 };

  for ( const auto& pair : ranges::view::zip( trial, reference ) ){
    const auto trial = std::get<0>(pair);
    const auto reference = std::get<1>(pair);
    REQUIRE( trial == Approx( reference ) );
  }

}
