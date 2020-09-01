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
   *  @brief Return the compound system
   */
  CompoundSystemVariant& compoundSystem() { return this->system_; }

  /**
   *  @brief Return whether or not the reconstructor is for resolved resonances
   */
  bool isResolved() {

    return this->system_.index() !=
           std::variant_size_v< CompoundSystemVariant > - 1;
  }

  /**
   *  @brief Return whether or not the reconstructor is for unresolved resonances
   */
  bool isUnresolved() {

    return not this->isResolved();
  }

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

    auto filter = [&] ( const auto& energy ) {

      return ( this->lowerEnergy() <= energy ) and
             ( energy <= this->upperEnergy() );
    };

    auto energies = std::visit( [&] ( auto& system ) { return system.grid(); },
                                this->system_ );
    return energies | ranges::view::filter( filter );
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
