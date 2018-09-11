template < typename Quantity >
static
void verifyNotNegative( const Quantity& value, const std::string& name ) {

  // the quantity cannot be negative

  if ( value.value < 0. ) {

    Log::error( "{} value cannot be negative.", name );
    throw std::exception();
  }
}
