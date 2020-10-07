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
                     const double slipRate) = 0;

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
                     const double slipRate);

private:

    double _lockedSlip;

private:

    // Not implemented
    StaticFriction(const StaticFriction&);
    const StaticFriction& operator=(const StaticFriction&);

};

class SlipWeakening : public Friction {
public:

    SlipWeakening(void);
    ~SlipWeakening(void);

    double traction(const double slip,
                    const double slipRate);

    double lockedSlip(void) const;

    double jacobianSlip(const double slip);

    double jacobianSlipRate(const double slipRate);

    void updateState(const double slip,
                     const double slipRate);

private:

    double _lockedSlip;

private:

    // Not implemented
    SlipWeakening(const SlipWeakening&);
    const SlipWeakening& operator=(const SlipWeakening&);

};

#endif

// End of file
