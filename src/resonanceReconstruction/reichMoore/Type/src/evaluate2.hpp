auto evaluate( const Quantity<ElectronVolts> energy ) const {
  const auto channelRadius = this->channelRadius( energy );
  const auto scatteringRadius = this->scatteringRadius( energy );
  return Parent::evaluate( energy, channelRadius, scatteringRadius );
}
