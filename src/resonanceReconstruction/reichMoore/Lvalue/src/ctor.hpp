Lvalue( std::vector< Resonance >&& resonances, int l, AP&& ap ) :
  ap_( std::move(ap) ),
  resonances( std::move( resonances ) ),
  orbitalAngularMomentum( l ){}

Lvalue( std::vector< Resonance >&& resonances, int l ) :
  ap_( std::nullopt ),
  resonances( std::move( resonances ) ),
  orbitalAngularMomentum( l ){}

Lvalue( const Lvalue& ) = default;

//working arouynd gcc bug regarding move ctors of lambdas with captures
Lvalue( Lvalue&& other ) :
  ap_( other.ap_ ),
  resonances( std::move( other.resonances ) ),
  orbitalAngularMomentum( other.orbitalAngularMomentum ){}

Lvalue& operator=( const Lvalue& ) = default;
Lvalue& operator=( Lvalue&& other ){
  (*this).~Lvalue();
  new(this) Lvalue( std::move(other) );
  return *this;  
}
