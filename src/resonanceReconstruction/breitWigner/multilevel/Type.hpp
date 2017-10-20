template< typename... >
class Type;

template< typename Radius >
class Type< Radius > : public Base< Type< Radius > > {
  Radius radius;
  
protected:
  friend Base< Type< Radius > >;
  using Parent = Base< Type< Radius > >;  
  #include "resonanceReconstruction/breitWigner/multilevel/Type/src/evaluate1.hpp"

public:
  #include "resonanceReconstruction/breitWigner/multilevel/Type/src/ctor1.hpp"
}; 

template< typename ChannelRadius, typename ScatteringRadius >
class Type< ChannelRadius, ScatteringRadius > :
  public Base< Type< ChannelRadius, ScatteringRadius > > {
  ChannelRadius channelRadius;
  ScatteringRadius scatteringRadius;
  
protected:
  friend Base< Type< ChannelRadius, ScatteringRadius > >;
  using Parent = Base< Type< ChannelRadius, ScatteringRadius > >;
  #include "resonanceReconstruction/breitWigner/multilevel/Type/src/evaluate2.hpp"

public:
  #include "resonanceReconstruction/breitWigner/multilevel/Type/src/ctor2.hpp"
};
