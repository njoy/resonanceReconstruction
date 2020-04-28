template< typename Derived >
class Base : public breitWigner::Type {
protected:
  using Parent = breitWigner::Type;

  const Derived& derived() const {
    return static_cast< const Derived& >( *this );
  }

  double targetSpin;

  #include "resonanceReconstruction/breitWigner/multiLevel/Base/src/lvalues.hpp"
  #include "resonanceReconstruction/breitWigner/multiLevel/Base/src/evaluate.hpp"
  #include "resonanceReconstruction/breitWigner/multiLevel/Base/src/order.hpp"
  #include "resonanceReconstruction/breitWigner/multiLevel/Base/src/ctor.hpp"

public:

  template< typename... Args >
  auto operator()( Args&&... args ) const {
    return this->derived().evaluate( std::forward< Args >( args )... );
  }
};
