#include <portinfo>

#include "Friction.hh"

#include <cmath>
#include <strings.h> 


// --------------------------------------------------------------------------------------------------
namespace _Friction {
    static const double zeroTolerance = 1.0e-10;
}

Friction::Friction(void) {}


Friction::~Friction(void) {}


// --------------------------------------------------------------------------------------------------
namespace _StaticFriction {
    static const double friction = 1.0;
}

StaticFriction::StaticFriction(void) :
    _lockedSlip(0.0) {}


StaticFriction::~StaticFriction(void) {}


double
StaticFriction::traction(const double slip,
                         const double slipRate) {
    return _StaticFriction::friction;
}


double
StaticFriction::lockedSlip(void) const {
    return _lockedSlip;
}


double
StaticFriction::jacobianSlip(const double slip) {
    return 0.0;
}


double
StaticFriction::jacobianSlipRate(const double slipRate) {
    return 0.0;
}


void
StaticFriction::updateState(const double slip,
                            const double slipRate,
			    const double dt) {
    if (slipRate <= _Friction::zeroTolerance) {
        _lockedSlip = slip;
    } // if
}


// --------------------------------------------------------------------------------------------------
namespace _SlipWeakeningFriction {
    static const double frictionStatic = 1.0;
    static const double frictionDynamic = 0.2;
    static const double weakeningSlip = 1.0e-6;
}

SlipWeakeningFriction::SlipWeakeningFriction(const bool forceHealing) :
    _lockedSlip(0.0),
    _forceHealing(forceHealing) {}


SlipWeakeningFriction::~SlipWeakeningFriction(void) {}


double
SlipWeakeningFriction::traction(const double slip,
                        const double slipRate) {
    const double d0 = _SlipWeakeningFriction::weakeningSlip;
    const double fs = _SlipWeakeningFriction::frictionStatic;
    const double fd = _SlipWeakeningFriction::frictionDynamic;

    const double d = slip - _lockedSlip;
    
    double f = fs;
    if (d > d0) {
        f = fd;
    } else if (d > 0.0) {
        f = fs - (fs-fd) * d / d0;
    }
    return f;
}


double
SlipWeakeningFriction::lockedSlip(void) const {
    return _lockedSlip;
}


double
SlipWeakeningFriction::jacobianSlip(const double slip) {
    const double d0 = _SlipWeakeningFriction::weakeningSlip;
    const double fs = _SlipWeakeningFriction::frictionStatic;
    const double fd = _SlipWeakeningFriction::frictionDynamic;

    double jacobian = 0.0;
    if ((fabs(slip) > 0.0) && (fabs(slip) < d0)) {
        jacobian = -(fs - fd) / d0;
    }
    return jacobian;
}


double
SlipWeakeningFriction::jacobianSlipRate(const double slipRate) {
    return 0.0;
}


void
SlipWeakeningFriction::updateState(const double slip,
				   const double slipRate,
				   const double dt) {
    if (slipRate <= _Friction::zeroTolerance || _forceHealing) {
        _lockedSlip = slip;
    } // if
}


// --------------------------------------------------------------------------------------------------
namespace _ViscousFriction {
    static const double frictionStatic = 0.2;
    static const double viscosity = 0.1;
    static const double rateParameter = 1.0;
}

ViscousFriction::ViscousFriction(void) :
  _lockedSlip(0.0) {}


ViscousFriction::~ViscousFriction(void) {}


double
ViscousFriction::traction(const double slip,
                        const double slipRate) {
    const double f0 = _ViscousFriction::frictionStatic;
    const double v0 = _ViscousFriction::rateParameter;
    const double viscosity = _ViscousFriction::viscosity;

    const double f = f0 + viscosity * slipRate / v0;
    return f;
}


double
ViscousFriction::lockedSlip(void) const {
  return _lockedSlip;
}


double
ViscousFriction::jacobianSlip(const double slip) {
  return 0.0;
}


double
ViscousFriction::jacobianSlipRate(const double slipRate) {
    const double v0 = _ViscousFriction::rateParameter;
    const double viscosity = _ViscousFriction::viscosity;

    return viscosity / v0;
}


void
ViscousFriction::updateState(const double slip,
			     const double slipRate,
			     const double dt) {
    if (slipRate <= _Friction::zeroTolerance) {
        _lockedSlip = slip;
    } // if
}


// --------------------------------------------------------------------------------------------------
namespace _RateStateFriction {
    static const double aStable = 0.05;
    static const double aUnstable= 0.2;
    static const double b = 0.1;
    static const double l = 0.2;
    static const double v0 = 0.01;
    static const double frictionStatic = 1.0;
    static const double vLinear = 1.0e-8;
}

RateStateFriction::RateStateFriction(const char* paramset) :
  _lockedSlip(0.0),
  _theta(_RateStateFriction::l/_RateStateFriction::v0) {
  if (0 == strcasecmp(paramset, "stable")) {
    _a = _RateStateFriction::aStable;
  } else if (0 == strcasecmp(paramset, "unstable")) {
    _a = _RateStateFriction::aUnstable;
  } else {
    _a = _RateStateFriction::aStable;
  }
}


RateStateFriction::~RateStateFriction(void) {}


double
RateStateFriction::traction(const double slip,
                        const double slipRate) {
    const double f0 = _RateStateFriction::frictionStatic;
    const double v0 = _RateStateFriction::v0;
    const double l = _RateStateFriction::l;
    const double vlinear = _RateStateFriction::vLinear;
    const double b = _RateStateFriction::b;

    double f = f0;
    if (slipRate > vlinear) {
      f += _a*log(slipRate/v0) + b*log(v0*_theta/l);
    } else {
      f += _a*log(vlinear/v0) + b*log(v0*_theta/l) - _a*(1.0-slipRate/vlinear);
    }
    return f;
}


double
RateStateFriction::lockedSlip(void) const {
  return _lockedSlip;
}


double
RateStateFriction::jacobianSlip(const double slip) {
  return 0.0;
}


double
RateStateFriction::jacobianSlipRate(const double slipRate) {
    const double vlinear = _RateStateFriction::vLinear;

    return (slipRate > vlinear) ? _a/slipRate : _a/vlinear;
}


void
RateStateFriction::updateState(const double slip,
			       const double slipRate,
			       const double dt) {
  const double l = _RateStateFriction::l;
  
    if (slipRate <= _Friction::zeroTolerance) {
        _lockedSlip = slip;
    } // if

    if (slipRate*dt/l > 1.0e-5) {
      _theta = _theta*exp(-slipRate*dt/l) + l/slipRate*(1.0-exp(-slipRate*dt/l));
    } else {
      _theta = _theta*exp(-slipRate*dt/l) + dt-0.5*slipRate*dt*dt/l;
    } // if/else
}


// End of file
