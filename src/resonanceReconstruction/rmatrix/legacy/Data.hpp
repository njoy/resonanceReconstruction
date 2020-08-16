/**
 *  @class
 *  @brief Generic data helper class for legacy unresolved data
 */
template < typename Quantity >
struct Data {

  Quantity elastic;
  Quantity capture;
  Quantity fission;
  Quantity competition;

  Data( Quantity elastic, Quantity capture,
        Quantity fission, Quantity competition ) :
    elastic( elastic ), capture( capture ),
    fission( fission ), competition( competition ) {}

  Data operator+( const Data& other ) const {
    return { this->elastic + other.elastic,
             this->capture + other.capture,
             this->fission + other.fission,
             this->competition + other.competition };
  }

  Data operator-( const Data& other ) const {
    return { this->elastic - other.elastic,
             this->capture - other.capture,
             this->fission - other.fission,
             this->competition - other.competition };
  }

  Data& operator+=( const Data& other ) {
    this->elastic += other.elastic;
    this->capture += other.capture;
    this->fission += other.fission;
    this->competition += other.competition;
    return *this;
  }

  Data& operator-=( const Data& other ) {
    this->elastic -= other.elastic;
    this->capture -= other.capture;
    this->fission -= other.fission;
    this->competition -= other.competition;
    return *this;
  }

  bool hasElastic() const { return this->elastic != Quantity(); }
  bool hasCapture() const { return this->capture != Quantity(); }
  bool hasFission() const { return this->fission != Quantity(); }
  bool hasCompetition() const { return this->competition != Quantity(); }
};
