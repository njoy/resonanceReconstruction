Lvalue( std::vector< Resonance >&& resonances, int l, AP&& ap ) :
  ap_( std::move(ap) ),
  resonances( std::move( resonances ) ),
  orbitalAngularMomentum( l ){}

Lvalue( std::vector< Resonance >&& resonances, int l ) :
  ap_( std::nullopt ),
  resonances( std::move( resonances ) ),
  orbitalAngularMomentum( l ){}
