/**
 *  @class
 *  @brief Class used to reconstruct cross sections from ENDF resonances
 */
class Reconstructor {

  /* fields */
  Energy lower_;
  Energy upper_;
  CompoundSystemVariant system_;

public:

  /* constructor */
  Reconstructor( const Energy& lower, const Energy& upper,
                 const CompoundSystemVariant& reconstructor ) :
      lower_( lower ), upper_( upper ),
      system_( reconstructor ) {}

  /**
   *  @brief Return the compund system
   */
  CompoundSystemVariant& compoundSystem() { return this->system_; }

  /**
   *  @brief Return the lower energy boundary
   */
  const Energy& lowerEnergy() const { return this->lower_; }

  /**
   *  @brief Return the lower energy boundary
   */
  const Energy& upperEnergy() const { return this->upper_; }

  /**
   *  @brief Return the minimal energy grid derived from the resonance
   *         parameters
   */
  std::vector< Energy > grid() const {

    return std::visit( [&] ( auto& system )
                           { return system.grid(); },
                       this->system_ );
  }

  /**
   *  @brief Reconstruct the cross sections at the given energy
   */
  std::map< ReactionID, CrossSection > operator()( const Energy& energy ) {

    std::map< ReactionID, CrossSection > result;
    if ( ( energy >= this->lowerEnergy() ) and
         ( energy <= this->upperEnergy() ) ) {

      std::visit( [&] ( auto& system )
                      { system.evaluate( energy, result ); },
                  this->system_ );
    }
    return result;
  }
};
