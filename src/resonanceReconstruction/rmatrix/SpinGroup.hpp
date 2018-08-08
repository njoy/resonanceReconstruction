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
template < typename BoundaryOption >
class SpinGroup {

  /* fields */
  std::vector< ParticleChannel > channels_;
  std::vector< unsigned int > incident_;
  ResonanceTable parameters_;

  Matrix< std::complex< double > > uMatrix_;
  std::vector< std::complex< double > > diagonalLMatrix_;

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

  auto channelIDs() const {
    return this->channels_
             | ranges::view::transform(
                 [&] ( const auto& channel )
                     { return std::visit(
                         [&] ( const auto& channel )
                             { return channel.channelID(); },
                         channel ); } );
  }

  auto reactionsIDs() const {
    return this->channels_
             | ranges::view::transform(
                 [&] ( const auto& channel )
                     { return std::visit(
                         [&] ( const auto& channel )
                             { return channel.particlePair().reaction(); },
                         channel ); } );
  }

  auto identifiers() const {
    return ranges::view::concat(
               this->reactionsIDs(),
               ranges::view::single( "capture" ) );
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

  auto incidentChannels() const { return ranges::view::all( this->incident_ ); }

public:

  /* constructor */
  SpinGroup( std::vector< ParticleChannel >&& channels,
             std::vector< unsigned int >&& incidentChannels,
             ResonanceTable&& table ) :
    channels_( std::move( channels ) ),
    incident_( std::move( incidentChannels ) ),
    parameters_( std::move( table ) ),
    uMatrix_( channels.size(), channels.size() ),
    diagonalLMatrix_( channels.size() ) {}

  const ResonanceTable& resonanceTable() const { return this->parameters_; }
  auto channels() const { return ranges::view::all( this->channels_ ); }

  void evaluate( const Energy& energy,
                 tsl::hopscotch_map< ReactionID, Quantity< Barn > >& result ) {

    // penetrability, shift factor, phase shift and Coulomb phase shift
    // for each channel except the eliminated capture channel
    const auto penetrabilities = this->penetrabilities( energy );
    const auto shiftFactors = this->shiftFactors( energy );
    const auto phaseShifts = this->phaseShifts( energy );
    const auto coulombShifts = this->coulombShifts( energy );
    const auto boundary = this->boundaryConditions();

    // the diagonal of the L matrix = S - B, iP
    this->diagonalLMatrix_ =
      ranges::view::zip_with(
          calculateLValue< BoundaryOption >,
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
    const unsigned int size = this->channels().size();
    this->uMatrix_ = Matrix< double >::Identity( size, size );
    this->resonanceTable().tmatrix( energy, this->diagonalLMatrix_, this->uMatrix_ );

    // the pi/k2 * gJ factor
    const auto factor = this->factor( energy );

    // the cross section identifiers
    const auto identifiers = this->identifiers();

    auto processIncidentChannel = [&] ( const unsigned int c ) {

      // lambda to derive a kronecker delta array for the current incident channel
      auto delta = [c,size] ( const auto value ) {
        return ranges::view::concat(
                   ranges::view::repeat_n( 0., c ),
                   ranges::view::single( value ),
                   ranges::view::repeat_n( 0., size - c - 1 ) );
      };

      // the elements of the ( 1 - RL )^-1 R matrix for each channel (assumes
      // the incident channel is the first channel)
      const auto row =
        ranges::make_iterator_range(
            this->uMatrix_.row(c).data(),
            this->uMatrix_.row(c).data() + size );

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

      // return the cross section values
      return ranges::view::zip( identifiers, crossSections );
    };

    // process the incident channels
    ranges::for_each(
        this->incidentChannels()
          | ranges::view::transform( processIncidentChannel )
          | ranges::view::join,
        [&] ( auto&& pair ) -> void
            { result[ std::get< 0 >( pair ) ] += std::get< 1 >( pair ); } );
  }
};
