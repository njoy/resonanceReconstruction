/**
 *  @brief Constructor
 *
 *  @param[in] mass     the particle or nuclide mass (in amu or dalton)
 *  @param[in] charge   the charge of the particle or nuclide (in Coulomb)
 *  @param[in] spin     the particle spin
 *  @param[in] parity   the particle parity
 */
Particle( const AtomicMass& mass,
          const ElectricalCharge& charge,
          const Spin& spin,
          const Parity& parity ) :
  mass_( mass ), charge_( charge ), spin_( spin ), parity_( parity ) {}

