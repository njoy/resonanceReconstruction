/**
 *  @class
 *  @brief A spin group corresponding to a J,pi value
 */
class SpinGroup {

  /* fields */
  std::vector< Channel > channels_;
  ResonanceTable parameters_;

public:

  /* constructor */

  auto resonanceTable() const { return this->parameters_; }

  auto generateCollisionMatrix( const Energy& energy ) const {

    // this code calculates the R-matrix using one eliminated channel

    // initialise the R-matrix
    Matrix< std::complex< double > > rMatrix =
        this->resonanceTable().rmatrix( energy );

  }
};
