/** @brief Horner evaluation of a polynomial function for different x and y
 *         types
 *
 *  A polynomial is defined as a sequence of coefficients c_i so that the
 *  polynomial of order n is given by:
 *
 *    y = c_n * x^n + c_(n - 1) * x^(n - 1) + c_(n - 2) * x^(n - 2) + ...
 *
 *  With the Horner scheme, a polynomial is evaluated as follows:
 *
 *    y = c_n * (x + c_(n-1) * (x + c_(n - 2) * (x + ...
 *
 *  The main reason for using the Horner scheme is computational efficiency.
 *
 *  Source: Numerical recipes - Third edition, p201-202
 *
 *  @param[in] first   an input iterator to the initial position in a sequence
 *                     (must be the highest order coefficient)
 *  @param[in] last    an input iterator to the final position in a sequence
 *  @param[in] x       the value of X
 *
 *  @return The evaluated Y value
 */
template < typename X, typename Y = X, typename Iter >
Y horner( Iter first, Iter last, const X& x ) noexcept {

  Y y = *first;
  ++first;
  for (; first != last; ++first) {

    y *= x;
    y += *first;
  }
  return y;
}
