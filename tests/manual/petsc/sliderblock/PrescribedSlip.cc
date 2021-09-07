#include <portinfo>

#include "PrescribedSlip.hh"

#include <cmath>

namespace _PrescribedSlip {
    static const double finalSlip = 1.0;
    static const double riseTime = 25.0;
}

// --------------------------------------------------------------------------------------------------
double
PrescribedSlip::slip(const double t) {
    const double df = _PrescribedSlip::finalSlip;
    const double tr = _PrescribedSlip::riseTime;
    const double pi = M_PI;

    double d = 0.0;
    if (t < tr) {
      d = -df / (2.0*pi) * sin(2.0*pi*t/tr) + df*t/tr;
    } else {
      d = df;
    }
    return d;
}


// --------------------------------------------------------------------------------------------------
double
PrescribedSlip::slipRate(const double t) {
    const double df = _PrescribedSlip::finalSlip;
    const double tr = t / _PrescribedSlip::riseTime;
    const double pi = M_PI;

    double v = 0.0;
    if (tr < tr) {
      v = -df/tr * cos(2.0*pi*t/tr) + df/tr;
    }
    return v;
}


// --------------------------------------------------------------------------------------------------
double
PrescribedSlip::slipAcc(const double t) {
    const double df = _PrescribedSlip::finalSlip;
    const double tr = _PrescribedSlip::riseTime;
    const double pi = M_PI;
    
    double a = 0.0;
    if (t < tr) {
      a = 2.0*pi / (tr*tr) * df * sin(2.0*pi*t/tr);
    }
    return a;
}


// End of file
