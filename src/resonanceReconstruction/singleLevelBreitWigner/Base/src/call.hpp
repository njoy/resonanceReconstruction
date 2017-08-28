template< typename Arg, typename... Args >
auto operator()( const Quantity<ElectronVolts> energy,
                 Arg&& arg,
                 Args&&... args ) const
  -> std::enable_if_t< not isTemperature(arg), CrossSection > {
  const auto kernel = this->psiChi( energy );
  return this->evaluate( energy, kernel,
                         std::forward<Arg>(arg),
                         std::forward<Args>(args)... );
}

template< typename Arg, typename... Args >
auto operator()( const Quantity<ElectronVolts> energy,
                 Arg&& arg,
                 Args&&... args ) const
  -> std::enable_if_t< isTemperature(arg), CrossSection > {
  const auto kernel = this->psiChi( energy, arg );
  return this->evaluate( energy, kernel, std::forward<Args>(args)... );
}
