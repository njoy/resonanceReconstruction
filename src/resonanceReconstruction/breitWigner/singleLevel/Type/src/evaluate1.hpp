template< typename... Args >
auto evaluate( const Quantity<ElectronVolts> energy,
               Args&&... args ) const {
  const auto radius = this->radius( energy );
  const auto kernel = psiChi( energy, args... );
  return Parent::evaluate( energy,
                           kernel,
                           radius,
                           radius,
                           this->radius );
}
