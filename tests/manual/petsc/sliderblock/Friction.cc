#include <portinfo>

#include "Friction.hh"

#include <cmath>

namespace _Friction {
    static const double zeroTolerance = 1.0e-10;
}

namespace _StaticFriction {
    static const double friction = 1.0;
}

namespace _SlipWeakeningFriction {
    static const double frictionStatic = 1.0;
    static const double frictionDynamic = 0.2;
    static const double weakeningSlip = 1.0e-6;
}

namespace _ViscousFriction {
    static const double frictionStatic = 0.2;
    static const double viscosity = 0.1;
    static const double rateParameter = 1.0;
}

// --------------------------------------------------------------------------------------------------
Friction::Friction(void) {}


Friction::~Friction(void) {}


// --------------------------------------------------------------------------------------------------
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
                            const double slipRate) {
    if (slipRate <= _Friction::zeroTolerance) {
        _lockedSlip = slip;
    } // if
}


// --------------------------------------------------------------------------------------------------
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
                           const double slipRate) {
    if (slipRate <= _Friction::zeroTolerance || _forceHealing) {
        _lockedSlip = slip;
    } // if
}


// --------------------------------------------------------------------------------------------------
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
			     const double slipRate) {
    if (slipRate <= _Friction::zeroTolerance) {
        _lockedSlip = slip;
    } // if
}


// End of file
