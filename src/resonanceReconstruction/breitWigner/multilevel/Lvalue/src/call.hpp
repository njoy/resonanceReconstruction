template< typename PsiChi, typename ChannelRatio >
auto operator()( const Quantity<ElectronVolts> energy,
                 const PsiChi& kernel,
                 const double channelRatio,
                 const double scatteringRatio,
                 const double targetSpin,
                 ChannelRatio&& rho ) const {
  return ( this->weightedQValue ) ?
    this->evaluate( kernel, channelRatio, scatteringRatio, targetSpin,
                    this->competitiveWidth( energy, rho ) ) :
    this->evaluate( kernel, channelRatio, scatteringRatio, targetSpin );
}
