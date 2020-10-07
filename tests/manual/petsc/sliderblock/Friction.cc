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
SlipWeakening::SlipWeakening(void) :
    _lockedSlip(0.0) {}


SlipWeakening::~SlipWeakening(void) {}


double
SlipWeakening::traction(const double slip,
                        const double slipRate) {
    const double d0 = _SlipWeakeningFriction::weakeningSlip;
    const double fs = _SlipWeakeningFriction::frictionStatic;
    const double fd = _SlipWeakeningFriction::frictionDynamic;

    double f = fs;
    if (slip > d0) {
        f = fd;
    } else if (slip > 0.0) {
        f = fs - (fs-fd) * fabs(slip) / d0;
    }
    return f;
}


double
SlipWeakening::lockedSlip(void) const {
    return _lockedSlip;
}


double
SlipWeakening::jacobianSlip(const double slip) {
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
SlipWeakening::jacobianSlipRate(const double slipRate) {
    return 0.0;
}


void
SlipWeakening::updateState(const double slip,
                           const double slipRate) {
    if (slipRate <= _Friction::zeroTolerance) {
        _lockedSlip = slip;
    } // if
}


// End of file
