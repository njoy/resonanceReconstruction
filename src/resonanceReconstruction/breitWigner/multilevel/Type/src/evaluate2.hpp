template< typename... Args >
auto evaluate( const Quantity<ElectronVolts> energy,
               Args&&... args ) const {
  const auto channelRadius = this->channelRadius( energy );
  const auto scatteringRadius = this->scatteringRadius( energy );
  const auto kernel = psiChi( energy, args... );
  return Parent::evaluate( energy,
                           kernel,
                           channelRadius,
                           scatteringRadius,
                           this->radius );
}
