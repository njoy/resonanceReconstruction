template< typename Derived >
class Base {
protected:
  #include "resonanceReconstruction/reichMoore/Base/src/order.hpp"

  const Derived& derived() const {
    return static_cast< const Derived& >( *this );
  }

  std::vector< Lvalue > lvalues;
  EnergyRange energyRange;
  double target2CompoundWeightRatio;
  std::variant< Both, Neither, Channel, Scattering > tag;

  #include "resonanceReconstruction/reichMoore/Base/src/ctor.hpp"
  #include "resonanceReconstruction/reichMoore/Base/src/neutronWaveNumber.hpp"
  #include "resonanceReconstruction/reichMoore/Base/src/evaluate.hpp"

public:
  template< typename... Args >
  auto operator()( Args&&... args ) const {
    return this->derived().evaluate( std::forward< Args >( args )... );
  }
};
