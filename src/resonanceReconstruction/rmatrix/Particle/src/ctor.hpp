//! @todo pybind11 variant needs default constructor workaround
#ifdef PYBIND11
/**
 *  @brief Default constructor - only enabled for pybind11
 */
Particle() = default;
#endif

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
  parity_( parity ) {

  verifyNotNegative( this->mass_, "Mass" );
  verifyNotNegative( this->charge_, "Charge" );
}

/**
 *  @brief Constructor
 *
 *  @param[in] id       the particle id or name (e.g. n, U235, U235_e1)
 *  @param[in] mass     the particle or nuclide mass (in amu or dalton)
 *  @param[in] spin     the particle spin
 *  @param[in] parity   the particle parity
 */
Particle( const ParticleID& id,
          const AtomicMass& mass,
          const Spin& spin,
          const Parity& parity ) :
  Particle( id, mass,
            ( id.za() - id.za() % 1000 ) / 1000 * constant::elementaryCharge,
            spin, parity ) {}
