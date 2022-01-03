#include <portinfo>

#include "PrescribedSlip.hh"

#include <cmath>

namespace _PrescribedSlip {
    static const double finalSlip = 3.0;
    static const double riseTime = 4.0;

    // Duraction of acceleration impulse in smoothed ramp slip time function.
    static const double impulseDuration = 0.5;

    // Maximum value of acceleration impulse in smoothed ramp slip time function.
    double amax(void) {
        const double df = finalSlip;
        const double tr = riseTime;
        const double ta = impulseDuration;
    
        return 0.5 * df * ta / (1.0/6.0*pow(tr,3) -1.0/3.0*pow(tr-0.5*ta,3) +0.5*tr*pow(tr-0.5*ta,2) -0.5*pow(tr,2)*(tr-0.5*ta) +0.5*(tr-ta)*pow(tr-0.5*ta,2) -0.5*pow(tr-ta,2)*(tr-0.5*ta) +0.25*pow(ta,2)*(tr-0.5*ta) +1.0/6.0*pow(tr-ta,3) -0.125*pow(ta,3));
    } // amax
}

// -------------------------------------------------------------------------------------------------
double
SlipFnCos::slip(const double t) {
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
SlipFnCos::slipRate(const double t) {
    const double df = _PrescribedSlip::finalSlip;
    const double tr = t / _PrescribedSlip::riseTime;
    const double pi = M_PI;

    double v = 0.0;
    if (t < tr) {
      v = -df/tr * cos(2.0*pi*t/tr) + df/tr;
    }

    return v;
}


// --------------------------------------------------------------------------------------------------
double
SlipFnCos::slipAcc(const double t) {
    const double df = _PrescribedSlip::finalSlip;
    const double tr = _PrescribedSlip::riseTime;
    const double pi = M_PI;
    
    double a = 0.0;
    if (t < tr) {
      a = 2.0*pi / (tr*tr) * df * sin(2.0*pi*t/tr);
    }
    
    return a;
}


// -------------------------------------------------------------------------------------------------
double
SlipFnRamp::slip(const double t) {
    const double tr = _PrescribedSlip::riseTime;
    const double ta = _PrescribedSlip::impulseDuration;
    const double c = 2.0 * _PrescribedSlip::amax() / ta;

    double d = 0.0;
    if (t <= 0.5*ta) {
      d = c * 1/6.0*pow(t,3);
    } else if (t <= ta) {
      d = c * (0.5*ta*pow(t,2) -1/6.0*pow(t,3) -0.25*pow(ta,2)*t +1.0/24.0*pow(ta,3));
    } else if (t <= tr-ta) {
      d = c * (0.25*pow(ta,2)*t -0.125*pow(ta,3));
    } else if (t <= tr-0.5*ta) {
      d = -c * (1/6.0*pow(t,3) -0.5*(tr-ta)*pow(t,2) +0.5*pow(tr-ta,2)*t -0.25*pow(ta,2)*t -1/6.0*pow(tr-ta,3) +0.125*pow(ta,3));
    } else if (t <= tr) {
      const double ti = tr - 0.5*ta;
      d = c * (1/6.0*pow(t,3) -0.5*tr*pow(t,2) +0.5*pow(tr,2)*t -1/3.0*pow(ti,3) +0.5*tr*pow(ti,2) -0.5*pow(tr,2)*ti +0.5*(tr-ta)*pow(ti,2) -0.5*pow(tr-ta,2)*ti +0.25*pow(ta,2)*ti +1/6.0*pow(tr-ta,3) -0.125*pow(ta,3));
    } else {
      const double ti = tr - 0.5*ta;
      d = c * (1/6.0*pow(tr,3) -1/3.0*pow(ti,3) +0.5*tr*pow(ti,2) -0.5*pow(tr,2)*ti +0.5*(tr-ta)*pow(ti,2) -0.5*pow(tr-ta,2)*ti +0.25*pow(ta,2)*ti +1/6.0*pow(tr-ta,3) -0.125*pow(ta,3));
    } // if/else
    
    return d;
}


// --------------------------------------------------------------------------------------------------
double
SlipFnRamp::slipRate(const double t) {
    const double tr = t / _PrescribedSlip::riseTime;
    const double ta = _PrescribedSlip::impulseDuration;
    const double c = 2.0 * _PrescribedSlip::amax() / ta;

    double v = 0.0;
    if (t <= 0.5*ta) {
      v = c * 0.5*pow(t,2);
    } else if (t <= ta) {
      v = c * (ta*t -0.5*pow(t,2) -0.25*pow(ta,2));
    } else if (t <= tr-ta) {
      v = c * 0.25*pow(ta,2);
    } else if (t <= tr-0.5*ta) {
      v = -c * (0.5*pow(t,2) -(tr-ta)*t +0.5*pow(tr-ta,2) -0.25*pow(ta,2));
    } else if (t <= tr) {
      v = c * (0.5*pow(t,2) -tr*t +0.5*pow(tr,2));
    } else {
      v = 0.0;
    } // if/else
    
    return v;
}


// --------------------------------------------------------------------------------------------------
double
SlipFnRamp::slipAcc(const double t) {
    const double tr = _PrescribedSlip::riseTime;
    const double ta = _PrescribedSlip::impulseDuration;
    const double c = 2.0 * _PrescribedSlip::amax() / ta;
    
    double a = 0.0;
    if (t <= 0.5*ta) {
      a = c * t;
    } else if (t <= ta) {
      a = c * (ta-t);
    } else if (t <= tr-ta) {
      a = 0.0;
    } else if (t <= tr-0.5*ta) {
      a = -c * (t-(tr-ta));
    } else if (t <= tr) {
      a = c * (t-tr);
    } else {
      a = 0.0;
    } // if/else

    return a;
}


// End of file
