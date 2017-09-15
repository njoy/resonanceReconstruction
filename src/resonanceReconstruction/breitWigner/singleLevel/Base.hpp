template< typename Derived >
class Base : public breitWigner::Type {
protected:
  using Parent = breitWigner::Type;

  const Derived& derived() const {
    return static_cast< const Derived& >( *this );
  }
  
  #include "resonanceReconstruction/breitWigner/singleLevel/Base/src/lvalues.hpp"
  #include "resonanceReconstruction/breitWigner/singleLevel/Base/src/evaluate.hpp"

  using breitWigner::Type::Type;

public:
  template< typename... Args >
  auto operator()( Args&&... args ) const {
    return this->derived().evaluate( std::forward< Args >( args )... );
  }
};
