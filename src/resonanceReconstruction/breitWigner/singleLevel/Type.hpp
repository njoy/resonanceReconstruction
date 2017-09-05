template< typename... >
class Type;

template< typename Radius >
class Type< Radius > : public Base, public CRTP< Type< Radius > > {
  Radius radius;
  
protected:
  friend CRTP< Type< Radius > >;
  using Parent = CRTP< Type< Radius > >;  
  #include "resonanceReconstruction/breitWigner/singleLevel/Type/src/evaluate1.hpp"

public:
  #include "resonanceReconstruction/breitWigner/singleLevel/Type/src/ctor1.hpp"
}; 

template< typename ChannelRadius, typename ScatteringRadius >
class Type< ChannelRadius, ScatteringRadius > :
  public Base,
  public CRTP< Type< ChannelRadius, ScatteringRadius > > {
  ChannelRadius channelRadius;
  ScatteringRadius scatteringRadius;
  
protected:
  friend CRTP< Type< ChannelRadius, ScatteringRadius > >;
  using Parent = CRTP< Type< ChannelRadius, ScatteringRadius > >;
  #include "resonanceReconstruction/breitWigner/singleLevel/Type/src/evaluate2.hpp"

public:
  #include "resonanceReconstruction/breitWigner/singleLevel/Type/src/ctor2.hpp"
};
