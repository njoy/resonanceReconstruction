/**
 *  @typedef
 *  @brief Channel types
 *
 *  Some of the data for a given channel will actually depend on the particle
 *  type of the channel. The penetrability, shift factor, phase shift and 
 *  coulomb phase shift will depend on whether the channel is a neutron, photon,
 *  charged particle or fission channel. This variant will allow us to capture
 *  that distinction.
 */
using ParticleChannel = std::variant< Channel< Neutron >,
                                      Channel< Photon >,
                                      Channel< ChargedParticle >,
                                      Channel< Fission > >;

/**
 *  @class
 *  @brief A spin group corresponding to a J,pi value
 */
class SpinGroup {

  /* fields */
  std::vector< ParticleChannel > channels_; // first channel is entrance
  ResonanceTable parameters_;

  mutable Matrix< std::complex< double > > matrix_;

  auto penetrabilities( const Energy& energy ) const {
    return this->channels_
             | ranges::view::transform(
                 [&] ( const auto& channel )
                     { return std::visit(
                         [&] ( const auto& channel )
                             { return channel.penetrability( energy ); },
                         channel ); } );
  }

  auto shiftFactors( const Energy& energy ) const {
    return this->channels_
             | ranges::view::transform(
                 [&] ( const auto& channel )
                     { return std::visit(
                         [&] ( const auto& channel )
                             { return channel.shiftFactor( energy ); },
                         channel ); } );
  }

  auto phaseShifts( const Energy& energy ) const {
    return this->channels_
             | ranges::view::transform(
                 [&] ( const auto& channel )
                     { return std::visit(
                         [&] ( const auto& channel )
                             { return channel.phaseShift( energy ); },
                         channel ); } );
  }

  auto coulombShifts( const Energy& energy ) const {
    return this->channels_
             | ranges::view::transform(
                 [&] ( const auto& channel )
                     { return std::visit(
                         [&] ( const auto& channel )
                             { return channel.coulombPhaseShift( energy ); },
                         channel ); } );
  }

  auto boundaryConditions() const {
    return this->channels_
             | ranges::view::transform(
                 [&] ( const auto& channel )
                     { return std::visit(
                         [&] ( const auto& channel )
                             { return channel.boundaryCondition(); },
                         channel ); } );
  }

  auto reactions() const {
    return this->channels_
             | ranges::view::transform(
                 [&] ( const auto& channel )
                     { return std::visit(
                         [&] ( const auto& channel )
                             { return channel.particlePair().reaction(); },
                         channel ); } );
  }

  auto factor( const Energy& energy ) const {

    auto factor = [&] ( const auto& channel ) {
      const auto waveNumber = channel.particlePair().waveNumber( energy );
      const auto squaredWaveNumber = waveNumber * waveNumber;
      const auto spinFactor = channel.statisticalSpinFactor();
      return pi / squaredWaveNumber * spinFactor;
    };
    return std::visit( factor, this->channels_.front() );
  }

public:

  /* constructor */
  SpinGroup( std::vector< ParticleChannel >&& channels,
             ResonanceTable&& table ) :
    channels_( std::move( channels ) ), parameters_( std::move( table ) ),
    matrix_( channels.size(), channels.size() ) {}

  auto resonanceTable() const { return this->parameters_; }
  auto channels() const { return ranges::view::all( this->channels_ ); }

  auto evaluate( const Energy& energy ) const {

    // penetrability, shift factor, phase shift and Coulomb phase shift
    // for each channel except the eliminated capture channel
    const auto penetrabilities = this->penetrabilities( energy );
    const auto shiftFactors = this->shiftFactors( energy );
    const auto phaseShifts = this->phaseShifts( energy );
    const auto coulombShifts = this->coulombShifts( energy );
    const auto boundary = this->boundaryConditions();

    // the diagonal of the L matrix = S - B, iP
    std::vector< std::complex< double > > diagonalLMatrix =
      ranges::view::zip_with(
          [] ( const double S, const double B, const double P )
             { return std::complex< double >( S - B, P ); },
          shiftFactors, boundary, penetrabilities );

    // the diagonal of the sqrt(P) matrix
    const auto diagonalSqrtPMatrix =
      penetrabilities
        | ranges::view::transform(
              [] ( const double penetrability ) -> double
                 { return std::sqrt( penetrability ); } );

    // the diagonal of the Omega matrix = exp( i ( w - phi ) )
    const auto diagonalOmegaMatrix =
      ranges::view::zip_with(
          [] ( const double w, const double phi )
             { return std::exp( std::complex< double >( 0.0, w - phi ) ); },
          coulombShifts, phaseShifts );

    // get the T = ( 1 - RL )^-1 R matrix
    this->matrix_ = this->resonanceTable().tmatrix( energy, diagonalLMatrix );

    // the index of the current incident channel
    const unsigned int c = 0;

    // lambda to derive a kronecker delta array for the current incident channel
    const unsigned int size = this->channels().size();
    auto delta = [=] ( const auto value ) {
      return ranges::view::concat(
                 ranges::view::repeat_n( 0., c ),
                 ranges::view::single( value ),
                 ranges::view::repeat_n( 0., size - c - 1 ) );
    };

    // the elements of the ( 1 - RL )^-1 R matrix for each channel (assumes
    // the incident channel is the first channel)
    const auto row =
      ranges::make_iterator_range(
          this->matrix_.row(c).data(),
          this->matrix_.row(c).data() + size );

    // the row of the U matrix corresponding with the incident channel
    const auto incidentSqrtP = diagonalSqrtPMatrix[c];
    const auto incidentOmega = diagonalOmegaMatrix[c];
    const auto uElements =
      ranges::view::zip_with( 
          [=] ( const auto delta, const auto tValue,
                const auto sqrtP, const auto omega )
              { return incidentOmega *
                       ( delta + std::complex< double >( 0., 2. ) *
                                 incidentSqrtP * tValue * sqrtP ) * omega; },
          delta( 1.0 ), row, diagonalSqrtPMatrix, diagonalOmegaMatrix );

    // the exponential of the coulomb phase shift for the incident channel
    const auto exponential =
      std::exp( std::complex< double >( 0., coulombShifts[c] ) );

    // the pi/k2 factor
    const auto factor = this->factor( energy );

    // lambda to calculate a norm squared
    auto normSquared = [] ( const auto value ) -> double
                          { return std::pow( std::abs( value ), 2. ); };

    // the cross section values
    const auto crossSections =
      ranges::view::concat(
          ranges::view::zip_with(
              [&] ( const auto delta, const auto uValue )
                  { return normSquared( delta - uValue ); },
              delta( exponential ),
              uElements ),
          ranges::view::single(
              ranges::accumulate(
                  uElements | ranges::view::transform( normSquared ), 1.,
                  ranges::minus() ) ) )
        | ranges::view::transform(
              [=] ( const auto value ) -> Quantity< Barn >
                  { return factor * value; } );

    // the reaction names
    const auto reactions =
        ranges::view::concat( this->reactions(),
                              ranges::view::single( "capture" ) );

    // return a range of pairs (name, xs)
    return ranges::view::zip( reactions, crossSections );
  }
};
