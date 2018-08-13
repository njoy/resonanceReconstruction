/**
 *  @brief Constructor
 *
 *  @param[in] id       the particle id or name (e.g. n, U235, U235_e1)
 *  @param[in] mass     the particle or nuclide mass (in amu or dalton)
 *  @param[in] charge   the charge of the particle or nuclide (in Coulomb)
 *  @param[in] spin     the particle spin
 *  @param[in] parity   the particle parity
 */
Particle( const ParticleID& id,
          const AtomicMass& mass,
          const ElectricalCharge& charge,
          const Spin& spin,
          const Parity& parity ) :
  id_( id ), mass_( mass ), charge_( charge ), spin_( spin ),
  parity_( parity ) {}

