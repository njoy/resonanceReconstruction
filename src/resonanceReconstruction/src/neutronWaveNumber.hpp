class NeutronWaveNumber
{
   using units = decltype(pow(Barns()*ElectronVolts(),Ratio<-1,2>));
   static constexpr Quantity<units> neutronConstant =
      5.787793139E-14 * root(kilo(grams))/constant::dirac;
   const Quantity<units> constant;

public:

   explicit NeutronWaveNumber(const double atomicWeightRatio) :
      constant(neutronConstant * atomicWeightRatio/(atomicWeightRatio+1)) { }

   // unused int is to distinguish from the above; both ctors are used
   NeutronWaveNumber(const double weightFraction, int) :
      constant(neutronConstant * weightFraction) { }

   Quantity<InvRootBarns> operator()(const Quantity<ElectronVolts> energy) const
   {
      return constant * sqrt(std::abs(energy));
   }
};

inline NeutronWaveNumber neutronWaveNumber(const double atomicWeightRatio)
{
   return NeutronWaveNumber(atomicWeightRatio);
}
