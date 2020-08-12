class PenetrationShift {
   const int L; // an angular momentum
public:
   explicit PenetrationShift(const int ell) : L(ell) { }

   std::array<double,2> operator()(const double channelRatio) const
   {
      const double rho2 = channelRatio * channelRatio;
      const double rho4 = rho2 * rho2;
      double r, a, b; // r for reciprocal

      // use Horner's method
      if (L == 0) {
         r = 1;
         a = 1;
         b = 0;
      } else if (L == 1) {
         r = 1/(1 + rho2);
         a = rho2;
         b = 1;
      } else if (L == 2) {
         r = 1/(9 + rho2*(3 + rho2));
         a = rho4;
         b = 18 + rho2*3;
      } else if (L == 3) {
         r = 1/(225 + rho2*(45 + rho2*(6 + rho2)));
         a = rho4*rho2;
         b = 675 + rho2*(90 + 6*rho2);
      } else if (L == 4) {
         r = 1/(11025 + rho2*(1575 + rho2*(135 + rho2*(10 + rho2))));
         a = rho4*rho4;
         b = 44100 + rho2*(4725 + rho2*(270 + 10*rho2));
      } else {
         throw std::exception();
      }

      return {{ (channelRatio > 0)*channelRatio*r*a, -r*b }};
   }
};

template<int L>
inline PenetrationShift
penetrationShift(const Integer<L>)
{
   return PenetrationShift(L);
}

inline std::array<double,2>
penetrationShift(const int L, const double channelRatio)
{
   return PenetrationShift(L)(channelRatio);
}
