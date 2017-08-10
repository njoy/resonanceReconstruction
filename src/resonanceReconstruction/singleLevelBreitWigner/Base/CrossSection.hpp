struct CrossSection {
  Quantity< Barn > elastic;
  Quantity< Barn > capture;
  Quantity< Barn > fission;

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
