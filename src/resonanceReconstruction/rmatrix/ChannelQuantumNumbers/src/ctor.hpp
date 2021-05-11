//! @todo pybind11 variant needs default constructor workaround
#ifdef PYBIND11
/**
 *  @brief Default constructor - only enabled for pybind11
 */
ChannelQuantumNumbers() = default;
#endif

/**
 *  @brief Constructor
 *
 *  @param[in] l        the orbital angular momentum
 *  @param[in] s        the channel spin
 *  @param[in] J        the total angular momentum
 *  @param[in] parity   the parity
 */
ChannelQuantumNumbers( const OrbitalAngularMomentum& l,
                       const Spin& s,
                       const TotalAngularMomentum& J,
                       const Parity& parity ) :
  l_( l ), s_( s ), J_( J ), parity_( parity ) {}
