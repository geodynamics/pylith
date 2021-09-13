#if !defined(friction_hh)
#define friction_hh

#include <portinfo>

class Friction {
public:

    Friction(void);
    virtual ~Friction(void);

    virtual
    double traction(const double slip,
                    const double slipRate) = 0;

    virtual
    double lockedSlip(void) const = 0;

    virtual
    double jacobianSlip(const double slip) = 0;

    virtual
    double jacobianSlipRate(const double slipRate) = 0;

    virtual
    void updateState(const double slip,
                     const double slipRate,
		     const double dt) = 0;

private:

    // Not implemented
    Friction(const Friction&);
    const Friction& operator=(const Friction&);

};

class StaticFriction : public Friction {
public:

    StaticFriction(void);
    ~StaticFriction(void);

    double traction(const double slip,
                    const double slipRate);

    double lockedSlip(void) const;

    double jacobianSlip(const double slip);

    double jacobianSlipRate(const double slipRate);

    void updateState(const double slip,
                     const double slipRate,
		     const double dt);

private:

    double _lockedSlip;

private:

    // Not implemented
    StaticFriction(const StaticFriction&);
    const StaticFriction& operator=(const StaticFriction&);

};

class SlipWeakeningFriction : public Friction {
public:

    SlipWeakeningFriction(const bool forceHealing);
    ~SlipWeakeningFriction(void);

    double traction(const double slip,
                    const double slipRate);

    double lockedSlip(void) const;

    double jacobianSlip(const double slip);

    double jacobianSlipRate(const double slipRate);

    void updateState(const double slip,
                     const double slipRate,
		     const double dt);

private:

    double _lockedSlip;
    bool _forceHealing;

private:

    // Not implemented
    SlipWeakeningFriction(void);
    SlipWeakeningFriction(const SlipWeakeningFriction&);
    const SlipWeakeningFriction& operator=(const SlipWeakeningFriction&);

};

class ViscousFriction : public Friction {
public:

    ViscousFriction(void);
    ~ViscousFriction(void);

    double traction(const double slip,
                    const double slipRate);

    double lockedSlip(void) const;

    double jacobianSlip(const double slip);

    double jacobianSlipRate(const double slipRate);

    void updateState(const double slip,
                     const double slipRate,
		     const double dt);

private:

    double _lockedSlip;

private:

    // Not implemented
    ViscousFriction(const ViscousFriction&);
    const ViscousFriction& operator=(const ViscousFriction&);

};

class RateStateFriction : public Friction {
public:

    RateStateFriction(const char* paramset);
    ~RateStateFriction(void);

    double traction(const double slip,
                    const double slipRate);

    double lockedSlip(void) const;

    double jacobianSlip(const double slip);

    double jacobianSlipRate(const double slipRate);

    void updateState(const double slip,
                     const double slipRate,
		     const double dt);

private:

    double _lockedSlip;
    double _theta;
    double _a;

private:

    // Not implemented
    RateStateFriction(void);
    RateStateFriction(const RateStateFriction&);
    const RateStateFriction& operator=(const RateStateFriction&);

};

#endif

// End of file
