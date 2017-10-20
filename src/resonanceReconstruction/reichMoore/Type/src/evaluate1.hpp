auto evaluate( const Quantity<ElectronVolts> energy ) const {
  const auto radius = this->radius( energy );
  return Parent::evaluate( energy, radius, radius );
}
