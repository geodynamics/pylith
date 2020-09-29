#include <portinfo>

#include "PrescribedSlip.hh"

namespace _PrescribedSlip {
    static const double finalSlip = 1.0;
    static const double riseTime = 25.0;
}

// --------------------------------------------------------------------------------------------------
double
PrescribedSlip::slip(const double t) {
    const double df = _PrescribedSlip::finalSlip;
    const double tr = t / _PrescribedSlip::riseTime;

    double d = 0.0;
    if (tr < 0.5) {
        d = 2.0 * df * tr*tr;
    } else if (tr < 1.0) {
        d = 4.0 * df * (-0.5*(tr*tr) + tr - 0.25);
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

    double v = 0.0;
    if (tr < 0.5) {
        v = 4.0 * df / _PrescribedSlip::riseTime * tr;
    } else if (tr < 1.0) {
        v = 4.0 * df / _PrescribedSlip::riseTime * (1.0 - tr);
    }
    return v;
}


// --------------------------------------------------------------------------------------------------
double
PrescribedSlip::slipAcc(const double t) {
    const double df = _PrescribedSlip::finalSlip;
    const double tr = t / _PrescribedSlip::riseTime;

    double a = 0.0;
    if (tr < 0.5) {
        a = 4.0 * df / (_PrescribedSlip::riseTime * _PrescribedSlip::riseTime);
    } else if (tr < 1.0) {
        a = -4.0 * df / (_PrescribedSlip::riseTime * _PrescribedSlip::riseTime);
    }
    return a;
}


// End of file
