template< typename... >
class Type;

template< typename ChannelRadius, typename ScattingRadius >
class Type< ChannelRadius, ScattingRadius > : protected Base {
protected:
  mutable ChannelRadius channelRadius;
  mutable ScatteringRadius scatteringRadius;

public:
  template<typename... Args>
  auto operator()( Args&&... args ) const {
    return Base::operator()( std::forward<Args>(args)...,
                             this->channelRadius,
                             this->scatteringRadius );
  }

  template< typename... Args >
  Type( ChannelRadius&& channelRadius,
        ScatteringRadius&& scatteringRadius,
        Args&&... args ) :
    Base( std::forward< Args >( args )... ),
    channelRadius( std::move(channelRadius) ),
    scatteringRadius( std::move(scatteringRadius) ){}
};

template< typename Radius >
class Type< Radius > : protected Base {
protected:
  mutable Radius radius;

public:
  template<typename... Args>
  auto operator()( Args&&... args ) const {
    return Base::operator()( std::forward<Args>(args)..., this->radius );
  }

  template< typename... Args >
  Type( Radius&& radius, Args&&... args ) :
    Base( std::forward< Args >( args )... ),
    radius( std::move(radius) ){}
};
