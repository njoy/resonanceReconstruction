struct CrossSection {
  Quantity< Barn > elastic;
  Quantity< Barn > capture;
  Quantity< Barn > fission;

  template< typename T1, typename T2, typename T3 >
  CrossSection( const Pack< T1, T2, T3 >& pack ) : elastic( std::get<0>(pack) ),
                                                   capture( std::get<1>(pack) ),
                                                   fission( std::get<2>(pack) ){}

  CrossSection( Quantity< Barn > elastic,
                Quantity< Barn > capture,
                Quantity< Barn > fission ) : elastic( elastic ),
                                             capture( capture ),
                                             fission( fission ){}
  
  CrossSection operator+( const CrossSection& other ) const {
    return { this->elastic + other.elastic,
             this->capture + other.capture,
             this->fission + other.fission }; 
  }

  CrossSection operator-( const CrossSection& other ) const {
    return { this->elastic - other.elastic,
             this->capture - other.capture,
             this->fission - other.fission }; 
  }
};

