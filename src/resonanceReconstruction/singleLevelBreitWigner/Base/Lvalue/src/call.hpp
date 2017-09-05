template< typename PsiChi, typename ChannelRadius >
double operator()( const Quantity<ElectronVolts> energy,
                 const PsiChi& kernel,
                 const double channelRatio,
                 const double scatteringRatio,
                 ChannelRadius&& channelRadius ) const {
  return ( this->weightedQValue ) ?
    (*this)( kernel, channelRatio, scatteringRatio,
             this->withCompetitiveWidth( energy, channelRadius ) ) :
    (*this)( kernel, channelRatio, scatteringRatio,
             this->withoutCompetitiveWidth() );
}
