NeutronWaveNumber
neutronWaveNumber() const
{
   return NeutronWaveNumber(target2CompoundWeightRatio,0);
}

Quantity<InvRootBarns>
neutronWaveNumber(const Quantity<ElectronVolts> energy) const
{
   return neutronWaveNumber()(energy);
}
