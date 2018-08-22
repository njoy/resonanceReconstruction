/**
 *  @brief Produce a string identifier for the set of quantum numbers 
 *
 *  This function produces a string identifier for the set of quantum numbers
 *  formatted like this: {l,s,Jpi} in which l values are integers, s and J
 *  half integers and pi either + or -.
 */
std::string toString() const {

  auto toHalfIntegerString = 
      [] ( const double a )
         { double half;
           return std::modf( a, &half ) == 0. ?
                      std::to_string( static_cast<int>( half ) ) :
                      std::to_string( 2 * static_cast<int>( half ) + 1 ) 
                        + "/2"; };

  return "{" + std::to_string( this->l_ ) + ","
             + toHalfIntegerString( this->s_ ) + ","
             + toHalfIntegerString( this->J_ )
             + ( this->parity_ > 0 ? "+" : "-" ) + "}";
}

