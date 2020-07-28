template< typename Derived >
class CRTP {
  const Derived& derived() const {
    return static_cast< const Derived& >( *this );
  }

protected:
  #include "resonanceReconstruction/breitWigner/CRTP/src/evaluate.hpp"

public:
  template< typename... Args >
  auto operator()( const Args&&... args ) const {
    return this->derived().evaluate( std::forward< Args >( args )... );
  }

};
