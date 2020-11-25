SCENARIO( "evaluate" ) {

  GIVEN( "Pu239 unresolved resonance data" ) {

    // test based on Pu239 ENDF/B-VIII.0 LRU=2 unresolved resonance evaluation
    // cross section values extracted from NJOY2016.57

    // this test is essentially for full energy dependent parameters without
    // competition

    // scattering radius from equation D.14
    double a = 0.123 * std::pow( 236.9986 * 1.008664, 1. / 3. ) + 0.08;

    // particles
    Particle neutron( ParticleID( "n" ), neutronMass, 0.0 * coulombs, 0.5, +1);
    Particle pu239( ParticleID( "Pu239" ), 236.9986 * neutronMass,
                    94.0 * coulombs, 0.5, +1); // unknown parity

    // particle pairs
    ParticlePair in( neutron, pu239 );

    // channels with l=0 and l=1
    Channel< Neutron > elastic00( in, in, 0. * electronVolt, { 0, 0.5, 0.0, +1 }, { a * rootBarn, 0.946 * rootBarn } );
    Channel< Neutron > elastic10( in, in, 0. * electronVolt, { 1, 0.5, 0.0, +1 }, { a * rootBarn, 0.946 * rootBarn } );

    // corresponding unresolved resonance tables
    ResonanceTable table00(
      { Resonance( 2.500000e+3 * electronVolt, 8.917200e+0 * electronVolt, 9.508100e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.842000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.550000e+3 * electronVolt, 8.915500e+0 * electronVolt, 8.673000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 4.020000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 2.650000e+3 * electronVolt, 8.913900e+0 * electronVolt, 1.113200e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.841000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.750000e+3 * electronVolt, 8.912200e+0 * electronVolt, 8.952900e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.840000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.850000e+3 * electronVolt, 8.910500e+0 * electronVolt, 9.309800e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.840000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.950000e+3 * electronVolt, 8.908800e+0 * electronVolt, 1.450600e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 7.840000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.050000e+3 * electronVolt, 8.907200e+0 * electronVolt, 7.845300e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.839000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 3.150000e+3 * electronVolt, 8.905500e+0 * electronVolt, 1.079100e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.838000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 3.250000e+3 * electronVolt, 8.903800e+0 * electronVolt, 1.065500e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.838000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 3.350000e+3 * electronVolt, 8.902200e+0 * electronVolt, 6.778100e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 1.150000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 3.450000e+3 * electronVolt, 8.900500e+0 * electronVolt, 9.239600e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 9.960000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.550000e+3 * electronVolt, 8.898800e+0 * electronVolt, 7.025700e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.836000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 3.650000e+3 * electronVolt, 8.897100e+0 * electronVolt, 1.215900e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 3.930000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 3.750000e+3 * electronVolt, 8.895500e+0 * electronVolt, 7.980200e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.835000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 3.850000e+3 * electronVolt, 8.893800e+0 * electronVolt, 1.082000e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.835000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 3.950000e+3 * electronVolt, 8.892100e+0 * electronVolt, 1.093700e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.834000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 4.125000e+3 * electronVolt, 8.889100e+0 * electronVolt, 9.226000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.829000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 4.375000e+3 * electronVolt, 8.884900e+0 * electronVolt, 9.985000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.828000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 4.625000e+3 * electronVolt, 8.880800e+0 * electronVolt, 9.217000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.827000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 4.875000e+3 * electronVolt, 8.876600e+0 * electronVolt, 9.974000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.826000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 5.125000e+3 * electronVolt, 8.872400e+0 * electronVolt, 1.011400e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.820000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 5.375000e+3 * electronVolt, 8.868200e+0 * electronVolt, 1.011000e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.820000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 5.625000e+3 * electronVolt, 8.864100e+0 * electronVolt, 1.010500e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.820000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 5.875000e+3 * electronVolt, 8.859900e+0 * electronVolt, 1.010000e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.820000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 6.125000e+3 * electronVolt, 8.855700e+0 * electronVolt, 1.009500e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.820000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 6.375000e+3 * electronVolt, 8.851500e+0 * electronVolt, 1.009100e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.820000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 6.625000e+3 * electronVolt, 8.847300e+0 * electronVolt, 1.008600e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 9.320000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 6.875000e+3 * electronVolt, 8.843100e+0 * electronVolt, 1.008100e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 1.056500e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 7.125000e+3 * electronVolt, 8.839000e+0 * electronVolt, 1.019800e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 4.814000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 7.375000e+3 * electronVolt, 8.834800e+0 * electronVolt, 1.014200e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 1.130000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 7.625000e+3 * electronVolt, 8.830700e+0 * electronVolt, 1.013100e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 1.309500e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 7.875000e+3 * electronVolt, 8.826500e+0 * electronVolt, 1.013400e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.078000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 8.125000e+3 * electronVolt, 8.822300e+0 * electronVolt, 1.005700e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.140000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 8.375000e+3 * electronVolt, 8.818200e+0 * electronVolt, 1.005300e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.807000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 8.625000e+3 * electronVolt, 8.814000e+0 * electronVolt, 1.004800e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.806000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 8.875000e+3 * electronVolt, 8.809900e+0 * electronVolt, 1.004300e-3 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.804000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 9.125000e+3 * electronVolt, 8.805700e+0 * electronVolt, 8.810000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.803000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 9.375000e+3 * electronVolt, 8.801400e+0 * electronVolt, 8.622000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 1.249000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 9.625000e+3 * electronVolt, 8.797400e+0 * electronVolt, 9.643000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 1.250000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 9.875000e+3 * electronVolt, 8.793200e+0 * electronVolt, 9.609000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 1.250000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.025000e+4 * electronVolt, 8.787100e+0 * electronVolt, 9.823000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 1.182000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.075000e+4 * electronVolt, 8.778900e+0 * electronVolt, 9.815000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.794000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.125000e+4 * electronVolt, 8.770500e+0 * electronVolt, 9.805000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.792000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.175000e+4 * electronVolt, 8.762200e+0 * electronVolt, 9.797000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.789000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.225000e+4 * electronVolt, 8.754100e+0 * electronVolt, 9.787000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.786000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.275000e+4 * electronVolt, 8.745800e+0 * electronVolt, 9.778000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.784000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.325000e+4 * electronVolt, 8.737600e+0 * electronVolt, 9.768000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.781000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.375000e+4 * electronVolt, 8.729400e+0 * electronVolt, 9.759000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.779000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.425000e+4 * electronVolt, 8.721200e+0 * electronVolt, 9.750000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 6.620000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.475000e+4 * electronVolt, 8.712900e+0 * electronVolt, 9.740000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.773000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.525000e+4 * electronVolt, 8.704700e+0 * electronVolt, 9.732000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.771000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.575000e+4 * electronVolt, 8.696500e+0 * electronVolt, 9.723000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.768000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.625000e+4 * electronVolt, 8.688300e+0 * electronVolt, 9.714000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.766000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.675000e+4 * electronVolt, 8.680100e+0 * electronVolt, 9.694000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.763000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.725000e+4 * electronVolt, 8.672000e+0 * electronVolt, 9.697000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 9.465000e-1 * electronVolt, 0. * electronVolt ),
        Resonance( 1.775000e+4 * electronVolt, 8.663800e+0 * electronVolt, 9.687000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.758000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.825000e+4 * electronVolt, 8.655700e+0 * electronVolt, 9.678000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.755000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.875000e+4 * electronVolt, 8.647500e+0 * electronVolt, 9.669000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.753000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.925000e+4 * electronVolt, 8.639300e+0 * electronVolt, 9.660000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.750000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 1.975000e+4 * electronVolt, 8.631200e+0 * electronVolt, 9.651000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.747000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.050000e+4 * electronVolt, 8.619000e+0 * electronVolt, 9.636000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.743000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.150000e+4 * electronVolt, 8.602800e+0 * electronVolt, 9.618000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.738000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.250000e+4 * electronVolt, 8.586700e+0 * electronVolt, 9.600000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.733000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.350000e+4 * electronVolt, 8.570500e+0 * electronVolt, 9.582000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.728000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.450000e+4 * electronVolt, 8.554400e+0 * electronVolt, 9.564000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.723000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.550000e+4 * electronVolt, 8.538200e+0 * electronVolt, 9.546000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.718000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.650000e+4 * electronVolt, 8.522200e+0 * electronVolt, 9.528000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.713000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.750000e+4 * electronVolt, 8.506200e+0 * electronVolt, 9.510000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.708000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.850000e+4 * electronVolt, 8.490100e+0 * electronVolt, 9.492000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.703000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 2.950000e+4 * electronVolt, 8.474100e+0 * electronVolt, 9.474000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.697000e+0 * electronVolt, 0. * electronVolt ),
        Resonance( 3.000000e+4 * electronVolt, 8.465900e+0 * electronVolt, 9.465000e-4 * rootElectronVolt, 4.070000e-2 * electronVolt, 2.693000e+0 * electronVolt, 0. * electronVolt )
      },
      { 1, 0, 2, 0 } );
    ResonanceTable table10(
      {
        Resonance( 2.500000e+3 * electronVolt, 8.917200e+0 * electronVolt, 1.573800e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.550000e+3 * electronVolt, 8.915500e+0 * electronVolt, 1.449800e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.650000e+3 * electronVolt, 8.913900e+0 * electronVolt, 1.824400e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.750000e+3 * electronVolt, 8.912200e+0 * electronVolt, 1.497800e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.850000e+3 * electronVolt, 8.910500e+0 * electronVolt, 1.555000e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.950000e+3 * electronVolt, 8.908800e+0 * electronVolt, 2.338500e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 3.050000e+3 * electronVolt, 8.907200e+0 * electronVolt, 1.337900e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 3.150000e+3 * electronVolt, 8.905500e+0 * electronVolt, 1.785600e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 3.250000e+3 * electronVolt, 8.903800e+0 * electronVolt, 1.767600e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 3.350000e+3 * electronVolt, 8.902200e+0 * electronVolt, 1.182400e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 3.450000e+3 * electronVolt, 8.900500e+0 * electronVolt, 1.558300e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 3.550000e+3 * electronVolt, 8.898800e+0 * electronVolt, 1.224400e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 3.650000e+3 * electronVolt, 8.897100e+0 * electronVolt, 2.003600e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 3.750000e+3 * electronVolt, 8.895500e+0 * electronVolt, 1.373800e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 3.850000e+3 * electronVolt, 8.893800e+0 * electronVolt, 1.806300e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 3.950000e+3 * electronVolt, 8.892100e+0 * electronVolt, 1.826000e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 4.125000e+3 * electronVolt, 8.889100e+0 * electronVolt, 1.349000e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 4.375000e+3 * electronVolt, 8.884900e+0 * electronVolt, 1.458600e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 4.625000e+3 * electronVolt, 8.880800e+0 * electronVolt, 1.348000e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 4.875000e+3 * electronVolt, 8.876600e+0 * electronVolt, 1.459600e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 5.125000e+3 * electronVolt, 8.872400e+0 * electronVolt, 1.508300e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 5.375000e+3 * electronVolt, 8.868200e+0 * electronVolt, 1.507600e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 5.625000e+3 * electronVolt, 8.864100e+0 * electronVolt, 1.506900e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 5.875000e+3 * electronVolt, 8.859900e+0 * electronVolt, 1.506200e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 6.125000e+3 * electronVolt, 8.855700e+0 * electronVolt, 1.505500e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 6.375000e+3 * electronVolt, 8.851500e+0 * electronVolt, 1.504800e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 6.625000e+3 * electronVolt, 8.847300e+0 * electronVolt, 1.437500e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 6.875000e+3 * electronVolt, 8.843100e+0 * electronVolt, 1.437500e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 7.125000e+3 * electronVolt, 8.839000e+0 * electronVolt, 1.484100e-3 * rootElectronVolt, 1.170000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 7.375000e+3 * electronVolt, 8.834800e+0 * electronVolt, 1.483400e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 7.625000e+3 * electronVolt, 8.830700e+0 * electronVolt, 1.482800e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 7.875000e+3 * electronVolt, 8.826500e+0 * electronVolt, 1.481100e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 8.125000e+3 * electronVolt, 8.822300e+0 * electronVolt, 1.499800e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 8.375000e+3 * electronVolt, 8.818200e+0 * electronVolt, 1.499100e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 8.625000e+3 * electronVolt, 8.814000e+0 * electronVolt, 1.498400e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 8.875000e+3 * electronVolt, 8.809900e+0 * electronVolt, 1.497700e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 9.125000e+3 * electronVolt, 8.805700e+0 * electronVolt, 1.479500e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 9.375000e+3 * electronVolt, 8.801400e+0 * electronVolt, 1.466200e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 9.625000e+3 * electronVolt, 8.797400e+0 * electronVolt, 1.438700e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 9.875000e+3 * electronVolt, 8.793200e+0 * electronVolt, 1.424500e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.025000e+4 * electronVolt, 8.787100e+0 * electronVolt, 1.463500e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.075000e+4 * electronVolt, 8.778900e+0 * electronVolt, 1.462300e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.125000e+4 * electronVolt, 8.770500e+0 * electronVolt, 1.461100e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.175000e+4 * electronVolt, 8.762200e+0 * electronVolt, 1.459800e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.225000e+4 * electronVolt, 8.754100e+0 * electronVolt, 1.458600e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.275000e+4 * electronVolt, 8.745800e+0 * electronVolt, 1.457400e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.325000e+4 * electronVolt, 8.737600e+0 * electronVolt, 1.456200e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.375000e+4 * electronVolt, 8.729400e+0 * electronVolt, 1.455000e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.425000e+4 * electronVolt, 8.721200e+0 * electronVolt, 1.453700e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.475000e+4 * electronVolt, 8.712900e+0 * electronVolt, 1.452400e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.525000e+4 * electronVolt, 8.704700e+0 * electronVolt, 1.451200e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.575000e+4 * electronVolt, 8.696500e+0 * electronVolt, 1.449900e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.625000e+4 * electronVolt, 8.688300e+0 * electronVolt, 1.448500e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.675000e+4 * electronVolt, 8.680100e+0 * electronVolt, 1.447400e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.725000e+4 * electronVolt, 8.672000e+0 * electronVolt, 1.445700e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.775000e+4 * electronVolt, 8.663800e+0 * electronVolt, 1.444400e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.825000e+4 * electronVolt, 8.655700e+0 * electronVolt, 1.443000e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.875000e+4 * electronVolt, 8.647500e+0 * electronVolt, 1.441700e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.925000e+4 * electronVolt, 8.639300e+0 * electronVolt, 1.440200e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 1.975000e+4 * electronVolt, 8.631200e+0 * electronVolt, 1.438900e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.050000e+4 * electronVolt, 8.619000e+0 * electronVolt, 1.436800e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.150000e+4 * electronVolt, 8.602800e+0 * electronVolt, 1.434100e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.250000e+4 * electronVolt, 8.586700e+0 * electronVolt, 1.431400e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.350000e+4 * electronVolt, 8.570500e+0 * electronVolt, 1.428700e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.450000e+4 * electronVolt, 8.554400e+0 * electronVolt, 1.426000e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.550000e+4 * electronVolt, 8.538200e+0 * electronVolt, 1.423300e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.650000e+4 * electronVolt, 8.522200e+0 * electronVolt, 1.420600e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.750000e+4 * electronVolt, 8.506200e+0 * electronVolt, 1.418000e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.850000e+4 * electronVolt, 8.490100e+0 * electronVolt, 1.415300e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 2.950000e+4 * electronVolt, 8.474100e+0 * electronVolt, 1.412600e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt ),
        Resonance( 3.000000e+4 * electronVolt, 8.465900e+0 * electronVolt, 1.411300e-3 * rootElectronVolt, 1.150000e-2 * electronVolt, 0. * electronVolt, 0. * electronVolt )
      },
      { 1, 0, 0, 0 } );

    // spin groups
    SpinGroup group00( std::move( elastic00 ), std::move( table00 ) );
    SpinGroup group10( std::move( elastic10 ), std::move( table10 ) );

    ReactionID elas( "n,Pu239->n,Pu239" );
    ReactionID capt( "n,Pu239->capture" );
    ReactionID fiss( "n,Pu239->fission" );

    THEN( "cross sections can be calculated for l=0" ) {

      Map< ReactionID, CrossSection > xs;
      group00.evaluate( 2500. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.19150687014446557 == Approx( xs[ elas ].value ) );
      CHECK( 8.0066413975172654E-002 == Approx( xs[ capt ].value ) );
      CHECK( 1.8803733373232319 == Approx( xs[ fiss ].value ) );
      xs.clear();

        group00.evaluate( 2550. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.53240949386785641   == Approx( xs[ elas ].value ) );
      CHECK( 0.22921106466840691   == Approx( xs[ capt ].value ) );
      CHECK( 1.1815059631732758    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 2650. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.25213434053502726   == Approx( xs[ elas ].value ) );
      CHECK( 8.7206221031042350E-2 == Approx( xs[ capt ].value ) );
      CHECK( 2.1055104012706138    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 2750. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.16605057320659217   == Approx( xs[ elas ].value ) );
      CHECK( 7.2159372952314810E-2 == Approx( xs[ capt ].value ) );
      CHECK( 1.6906668483562362    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 2850. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.17675483864795602   == Approx( xs[ elas ].value ) );
      CHECK( 7.2775647174449468E-2 == Approx( xs[ capt ].value ) );
      CHECK( 1.7193973114159513    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 2950. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.84250183403028034   == Approx( xs[ elas ].value ) );
      CHECK( 0.20878087391704500   == Approx( xs[ capt ].value ) );
      CHECK( 1.9620917072271309    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 3050. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.12432775265158547   == Approx( xs[ elas ].value ) );
      CHECK( 6.1163013364684882E-2 == Approx( xs[ capt ].value ) );
      CHECK( 1.4161884802958304    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 3150. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.22561500125736328   == Approx( xs[ elas ].value ) );
      CHECK( 7.6615072973622106E-2 == Approx( xs[ capt ].value ) );
      CHECK( 1.8641092249919333    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 3250. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.21807395240201172   == Approx( xs[ elas ].value ) );
      CHECK( 7.4436332619761736E-2 == Approx( xs[ capt ].value ) );
      CHECK( 1.8119259941642061    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 3350. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.18137097686043865   == Approx( xs[ elas ].value ) );
      CHECK( 9.1682510866057609E-2 == Approx( xs[ capt ].value ) );
      CHECK( 1.0446218998888765    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 3450. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.33811092095257589   == Approx( xs[ elas ].value ) );
      CHECK( 0.12238618507556205   == Approx( xs[ capt ].value ) );
      CHECK( 1.3082674145369975    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 3550. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 9.4222884537941615E-2 == Approx( xs[ elas ].value ) );
      CHECK( 5.1232066882442959E-2 == Approx( xs[ capt ].value ) );
      CHECK( 1.1795122971392324    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 3650. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.81634590125451501   == Approx( xs[ elas ].value ) );
      CHECK( 0.22475774563898732   == Approx( xs[ capt ].value ) );
      CHECK( 1.2187705332254992    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 3750. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.11897521443851497   == Approx( xs[ elas ].value ) );
      CHECK( 5.4750023067088878E-2 == Approx( xs[ capt ].value ) );
      CHECK( 1.2885495823943456    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 3850. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.21157087224135340   == Approx( xs[ elas ].value ) );
      CHECK( 6.7833282397839756E-2 == Approx( xs[ capt ].value ) );
      CHECK( 1.6759726300286755    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 3950. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.21388843622900414   == Approx( xs[ elas ].value ) );
      CHECK( 6.7308701798889745E-2 == Approx( xs[ capt ].value ) );
      CHECK( 1.6688078116143830    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 4125. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.15237594699406456   == Approx( xs[ elas ].value ) );
      CHECK( 5.7777607420276038E-2 == Approx( xs[ capt ].value ) );
      CHECK( 1.3975930514284070    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 4375. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.17329753544874527   == Approx( xs[ elas ].value ) );
      CHECK( 5.9124776480665817E-2 == Approx( xs[ capt ].value ) );
      CHECK( 1.4542360468174687    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 4625. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.14502485444365804   == Approx( xs[ elas ].value ) );
      CHECK( 5.3828896024935111E-2 == Approx( xs[ capt ].value ) );
      CHECK( 1.3127833824654798    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 4875. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.16529825198479864   == Approx( xs[ elas ].value ) );
      CHECK( 5.5261784713079722E-2 == Approx( xs[ capt ].value ) );
      CHECK( 1.3699821983863170    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 5125. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.16639305934254739   == Approx( xs[ elas ].value ) );
      CHECK( 5.4194283697655776E-2 == Approx( xs[ capt ].value ) );
      CHECK( 1.3497290393616075    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 5375. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.16268953380668033   == Approx( xs[ elas ].value ) );
      CHECK( 5.2593025082403114E-2 == Approx( xs[ capt ].value ) );
      CHECK( 1.3148171569261782    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 5625. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.15908515734560352   == Approx( xs[ elas ].value ) );
      CHECK( 5.1102328845063341E-2 == Approx( xs[ capt ].value ) );
      CHECK( 1.2821714835929943    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 5875. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.15560491238110280   == Approx( xs[ elas ].value ) );
      CHECK( 4.9714536041432320E-2 == Approx( xs[ capt ].value ) );
      CHECK( 1.2516809644354141    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 6125. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.15223770752423660   == Approx( xs[ elas ].value ) );
      CHECK( 4.8418171339677213E-2 == Approx( xs[ capt ].value ) );
      CHECK( 1.2231086206453567    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 6375. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.14900465818238756   == Approx( xs[ elas ].value ) );
      CHECK( 4.7207275676135384E-2 == Approx( xs[ capt ].value ) );
      CHECK( 1.1963685107448991    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 6625. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.32504433781686132   == Approx( xs[ elas ].value ) );
      CHECK( 8.7852880947296541E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.95009113947593116   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 6875. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.29698267356801500   == Approx( xs[ elas ].value ) );
      CHECK( 8.0149267268396440E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.95783948856106582   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 7125. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.45120032091091256   == Approx( xs[ elas ].value ) );
      CHECK( 0.11634139860823513   == Approx( xs[ capt ].value ) );
      CHECK( 0.75667658507584223   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 7375. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.27945357105790858   == Approx( xs[ elas ].value ) );
      CHECK( 7.4074143679746379E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.93862979091754140   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 7625. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.24978775397206535   == Approx( xs[ elas ].value ) );
      CHECK( 6.6770420508196357E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.95061160280050871   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 7875. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.17471812450908425   == Approx( xs[ elas ].value ) );
      CHECK( 5.0042782419628028E-2 == Approx( xs[ capt ].value ) );
      CHECK( 1.0203015082084290    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 8125. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.16517775807741159   == Approx( xs[ elas ].value ) );
      CHECK( 4.7939138970887675E-2 == Approx( xs[ capt ].value ) );
      CHECK( 1.0011774519612162    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 8375. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.12656110246154728   == Approx( xs[ elas ].value ) );
      CHECK( 3.9783399640737488E-2 == Approx( xs[ capt ].value ) );
      CHECK( 1.0270825325205948    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 8625. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.12402553688241182   == Approx( xs[ elas ].value ) );
      CHECK( 3.9048832521618407E-2 == Approx( xs[ capt ].value ) );
      CHECK( 1.0102594436876879    == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 8875. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.12158851772384549   == Approx( xs[ elas ].value ) );
      CHECK( 3.8355976872688796E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.99410218912765091   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 9125. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 8.9767991253907886E-2 == Approx( xs[ elas ].value ) );
      CHECK( 3.4341375668305384E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.87251080101173106   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 9375. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.17469973200728195   == Approx( xs[ elas ].value ) );
      CHECK( 5.3831014710830731E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.73202418481106291   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 9625. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.21092794454508190   == Approx( xs[ elas ].value ) );
      CHECK( 5.6932375941261736E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.79047923048262947   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 9875. * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.20677640201394967   == Approx( xs[ elas ].value ) );
      CHECK( 5.5853963098592793E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.77667870699437702   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 10250 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.21919893022516124   == Approx( xs[ elas ].value ) );
      CHECK( 5.6984097658455336E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.76382692427893362   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 10750 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.10024737268829077   == Approx( xs[ elas ].value ) );
      CHECK( 3.3523197099629781E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.87726524713903720   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 11250 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 9.6193673136968458E-2 == Approx( xs[ elas ].value ) );
      CHECK( 3.2568234590924651E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.85497979604323049   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 11750 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 9.2357629141740866E-2 == Approx( xs[ elas ].value ) );
      CHECK( 3.1690897348028782E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.83425349013582895   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 12250 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 8.8599803807338812E-2 == Approx( xs[ elas ].value ) );
      CHECK( 3.0867228079327628E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.81467397261674990   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 12750 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 8.4941741376990385E-2 == Approx( xs[ elas ].value ) );
      CHECK( 3.0092697057058231E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.79641476410250800   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 13250 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 8.1408750435283991E-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.9370147489819625E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.77910002175173432   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 13750 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 7.7958264907679153E-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.8686446802197486E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.76286974987028933   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 14250 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.25838417963363508   == Approx( xs[ elas ].value ) );
      CHECK( 6.0769455120262474E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.53101203625685478   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 14750 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 7.1391262343816819E-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.7444101139721691E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.73283592905799144   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 15250 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.8225086163257231E-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.6870074614766724E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.71909744664567621   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 15750 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.5155415118072213E-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.6329957772667847E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.70593064356737856   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 16250 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 6.2122826939736989E-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.5811583726595026E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.69343255095364220   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 16750 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 5.9006945450954071E-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.5304786536288051E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.68081182722193678   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 17250 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 0.18029420265670587   == Approx( xs[ elas ].value ) );
      CHECK( 4.4837467813490010E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.52617932716597859   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 17750 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 5.3503124484172909E-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.4412126453887236E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.65915097816883272   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 18250 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 5.0765225821458287E-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.3989545957567117E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.64867861962495932   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 18750 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 4.8049583492932119E-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.3580731458462432E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.63868589570728851   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 19250 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 4.5423595720461757E-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.3194508045064011E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.62905565593094659   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 19750 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 4.2848890047287536E-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.2824218938674472E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.61979543333848852   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 20500 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 3.9040229962690548E-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.2292284901336470E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.60650484880031297   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 21500 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 3.4142605039863771E-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.1632269739881034E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.59001606956918151   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 22500 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.9416627690282526E-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.1020864229385282E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.57466203593246667   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 23500 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.4850017224460008E-2 == Approx( xs[ elas ].value ) );
      CHECK( 2.0453070931874793E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.56033208206683116   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 24500 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.0430638209501251E-2 == Approx( xs[ elas ].value ) );
      CHECK( 1.9923728405953654E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.54690680719353613   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 25500 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 1.6148804657479837E-2 == Approx( xs[ elas ].value ) );
      CHECK( 1.9429298330690313E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.53430780566088021   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 26500 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 1.1994829180871192E-2 == Approx( xs[ elas ].value ) );
      CHECK( 1.8965638230391178E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.52243721946785282   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 27500 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 7.9610080686992046E-3 == Approx( xs[ elas ].value ) );
      CHECK( 1.8530211055199686E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.51123896147655168   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 28500 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 4.0400218614800310E-3 == Approx( xs[ elas ].value ) );
      CHECK( 1.8120610891847969E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.50065865315291458   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 29500 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( 2.5628355549606185E-4 == Approx( xs[ elas ].value ) );
      CHECK( 1.7737946615604176E-2 == Approx( xs[ capt ].value ) );
      CHECK( 0.49059633787340712   == Approx( xs[ fiss ].value ) );
      xs.clear();

      group00.evaluate( 30000 * electronVolt, xs );
      CHECK( 3 == xs.size() );
      CHECK( -1.5668427469941182E-003 == Approx( xs[ elas ].value ) );
      CHECK( 1.7558869727382194E-002 == Approx( xs[ capt ].value ) );
      CHECK( 0.48573682117577405   == Approx( xs[ fiss ].value ) );
      xs.clear();
    } // THEN

    THEN( "cross sections can be calculated for l=1" ) {

      Map< ReactionID, CrossSection > xs;
      group10.evaluate( 2500. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.2023020985907092E-003 == Approx( xs[ elas ].value ) );
      CHECK( 2.6510443636947958E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 2550. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.7476510412513841E-003 == Approx( xs[ elas ].value ) );
      CHECK( 2.4827385712725372E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 2650. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 6.0300911305825097E-003 == Approx( xs[ elas ].value ) );
      CHECK( 3.0620671465621987E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 2750. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.5143076580705255E-003 == Approx( xs[ elas ].value ) );
      CHECK( 2.6133296257599359E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 2850. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 5.1215418939693259E-003 == Approx( xs[ elas ].value ) );
      CHECK( 2.7265074598598778E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 2950. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 1.1003145077357635E-002 == Approx( xs[ elas ].value ) );
      CHECK( 3.8541600606563454E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 3050. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.3852082200904006E-003 == Approx( xs[ elas ].value ) );
      CHECK( 2.4432168818782105E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 3150. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 7.6741516396198351E-003 == Approx( xs[ elas ].value ) );
      CHECK( 3.1406061367627705E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 3250. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 7.9222556184783723E-003 == Approx( xs[ elas ].value ) );
      CHECK( 3.1367542545700268E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 3350. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.1114861669344074E-003 == Approx( xs[ elas ].value ) );
      CHECK( 2.2567713080291964E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 3450. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 7.0021438787956948E-003 == Approx( xs[ elas ].value ) );
      CHECK( 2.8674448528858854E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 3550. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.8111807953741207E-003 == Approx( xs[ elas ].value ) );
      CHECK( 2.3620143700588472E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 3650. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 1.1710152544242994E-002 == Approx( xs[ elas ].value ) );
      CHECK( 3.5458553155074193E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 3750. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 6.4238612677885044E-003 == Approx( xs[ elas ].value ) );
      CHECK( 2.6353028517637503E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 3850. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 1.0663093087351881E-002 == Approx( xs[ elas ].value ) );
      CHECK( 3.2997122593268577E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 3950. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 1.1281639913675631E-002 == Approx( xs[ elas ].value ) );
      CHECK( 3.3417796718398278E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 4125. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 7.2270339667436278E-003 == Approx( xs[ elas ].value ) );
      CHECK( 2.6510858543463627E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 4375. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 9.0249990374670248E-003 == Approx( xs[ elas ].value ) );
      CHECK( 2.8529339448156586E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 4625. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 8.5937818567934574E-003 == Approx( xs[ elas ].value ) );
      CHECK( 2.7077358528070795E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 4875. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 1.0605151136233318E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.9034883890202201E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 5125. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 1.2028435263069327E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.9956119407914070E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 5375. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 1.2857916074084166E-002 == Approx( xs[ elas ].value ) );
      CHECK( 3.0102913170317849E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 5625. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 1.3689703865613531E-002 == Approx( xs[ elas ].value ) );
      CHECK( 3.0222060074310275E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 5875. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 1.4522944853138441E-002 == Approx( xs[ elas ].value ) );
      CHECK( 3.0316986822297144E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 6125. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 1.5356477392062660E-002 == Approx( xs[ elas ].value ) );
      CHECK( 3.0389796992134316E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 6375. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 1.6189405875248475E-002 == Approx( xs[ elas ].value ) );
      CHECK( 3.0442669583553121E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 6625. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 1.5859480598382138E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.9535829644726366E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 6875. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 1.6652468287478753E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.9574893135334200E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 7125. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 1.8153439390859307E-002 == Approx( xs[ elas ].value ) );
      CHECK( 3.0414611359299636E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 7375. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 1.9133107975699823E-002 == Approx( xs[ elas ].value ) );
      CHECK( 3.0238790218370653E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 7625. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 1.9941529106464018E-002 == Approx( xs[ elas ].value ) );
      CHECK( 3.0221524145615214E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 7875. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.0723157753384215E-002 == Approx( xs[ elas ].value ) );
      CHECK( 3.0179134552114441E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 8125. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.1944847636254566E-002 == Approx( xs[ elas ].value ) );
      CHECK( 3.0393415162355636E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 8375. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.2751021946641043E-002 == Approx( xs[ elas ].value ) );
      CHECK( 3.0342229616533642E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 8625. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.3552418169258117E-002 == Approx( xs[ elas ].value ) );
      CHECK( 3.0283356440788316E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 8875. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.4348282853588782E-002 == Approx( xs[ elas ].value ) );
      CHECK( 3.0216878340151986E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 9125. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.4709619283506994E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.9927264204408215E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 9375. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.5167414331653544E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.9696495341874293E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 9625. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.5240183496679103E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.9287552132184650E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 9875. * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.5624916368333736E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.9042310635306885E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 10250 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.7794948055394522E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.9394418530553138E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 10750 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 2.9273857703703336E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.9204075513201116E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 11250 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.0730186968881988E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.9002152415992615E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 11750 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.2159637699698852E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.8789017598357655E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 12250 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.3567764665301779E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.8568639162883575E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 12750 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.4953172171326277E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.8342947009609890E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 13250 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.6314877676278376E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.8112433209211442E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 13750 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.7653510509607642E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.7878666675839256E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 14250 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 3.8965764526431058E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.7641728406305605E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 14750 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.0255868877926267E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.7403835760061002E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 15250 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.1527117763205480E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.7166030640551286E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 15750 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.2772855193461599E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.6927391539082687E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 16250 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.3993234443098120E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.6688515168023185E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 16750 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.5204780634726258E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.6453364882882184E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 17250 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.6370637899436490E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.6214239419383906E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 17750 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.7533371970861243E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.5980661075231740E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 18250 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.8671662386102321E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.5748026741041340E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 18750 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 4.9795688562071963E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.5518787663501210E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 19250 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 5.0892062034895061E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.5290445568664981E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 19750 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 5.1978519385038967E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.5065991746301322E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 20500 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 5.3568268943189036E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.4733473486602910E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 21500 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 5.5632156459095979E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.4300734550035865E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 22500 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 5.7629840810745007E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.3879590858922682E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 23500 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 5.9565766017836107E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.3470807853849304E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 24500 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 6.1441524457448743E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.3073866702311579E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 25500 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 6.3261320352394679E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.2689218670092510E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 26500 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 6.5025560645448288E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.2315883019425296E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 27500 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 6.6743841415021252E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.1954721448188826E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 28500 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 6.8408050122278491E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.1604534085559356E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 29500 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 7.0024340691918596E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.1265074702546013E-002 == Approx( xs[ capt ].value ) );
      xs.clear();

      group10.evaluate( 30000 * electronVolt, xs );
      CHECK( 2 == xs.size() );
      CHECK( 7.0820124212988111E-002 == Approx( xs[ elas ].value ) );
      CHECK( 2.1100119938459196E-002 == Approx( xs[ capt ].value ) );
      xs.clear();
    } // THEN
  } // GIVEN
} // SCENARIO
