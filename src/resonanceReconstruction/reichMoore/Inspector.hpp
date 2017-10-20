struct InspectorBase{
  Quantity< ElectronVolts > energy;
  Quantity< InvRootBarn > waveNumber;
  double channelRatio;
  double scatteringRatio;

  InspectorBase( Quantity< ElectronVolts > energy,
                 Quantity< InvRootBarn > waveNumber,
                 double channelRatio,
                 double scatteringRatio ) :
    energy( energy ),
    waveNumber( waveNumber ),
    channelRatio( channelRatio ),
    scatteringRatio( scatteringRatio ){}
};

template< typename Tag >
struct Inspector;

// NRO == 0, NAPS == 1
template<>
struct Inspector< Both > : InspectorBase {
  using InspectorBase::InspectorBase;
  auto operator()( const Lvalue& lValue ) const {
    decltype(auto) ap = lValue.ap();
    if ( ap ){
      const double ratio = waveNumber * (*ap)(energy);
      return lValue( this->energy, ratio, ratio );
    }
    return lValue( this->energy, this->channelRatio, this->scatteringRatio );
  }
};

// NRO == 1, NAPS == 0 or NAPS == 1
template< >
struct Inspector< Neither > : InspectorBase {
  using InspectorBase::InspectorBase;
  auto operator()( const Lvalue& lValue ) const {
    return lValue( this->energy, this->channelRatio, this->scatteringRatio );
  }
};

// NRO == 0, NAPS == 0
template<>
struct Inspector< Scattering > : InspectorBase {
  using InspectorBase::InspectorBase;
  auto operator()( const Lvalue& lValue ) const {
    decltype(auto) ap = lValue.ap();
    if ( ap ){
      const double scatteringRatio = waveNumber * (*ap)(energy);
      return lValue( this->energy, this->channelRatio, scatteringRatio );
    }
    return lValue( this->energy, this->channelRatio, this->scatteringRatio );
  }
};

// NRO == 1, NAPS == 2
template<>
struct Inspector< Channel > : InspectorBase {
  using InspectorBase::InspectorBase;
  auto operator()( const Lvalue& lValue ) const {
    decltype(auto) ap = lValue.ap();
    if ( ap ){
      const double channelRatio = waveNumber * (*ap)(energy);
      return lValue( this->energy, channelRatio, this->scatteringRatio );
    }
    return lValue( this->energy, this->channelRatio, this->scatteringRatio );
  }
};

