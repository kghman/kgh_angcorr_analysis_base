#ifndef __constants_h
#define __constants_h

#include <cmath>

const long double
  RADTODEG          = 180./M_PI,
  DEGTORAD          = M_PI/180.,
  MEVTOJ            = 1.60218E-13, //J per MeV
  UTOMEV            = 931.4940954, //MeV per u
  SQMMTOBARN        = 1E22, //barn per square mm
  SQFMTOBARN        = 0.01, //barn per square fm
  C                 = 2.9979E8, //m/s
  UNIT_CHARGE       = 1.602176634E-19, //Coulombs
  RESTMASS_ELECTRON = 0.000548579909, //amu
  AVOGADROS_NUMBER  = 6.02214076E23, //particles per mol
  HBAR_C            = 197.3269804, //MeV*fm
  FINE_STRUCTURE    = 1./(137.036);

#endif
