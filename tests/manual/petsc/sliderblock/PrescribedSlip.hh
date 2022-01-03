#if !defined(prescribedslip_hh)
#define prescribedslip_hh

#include <portinfo>

class SlipFnCos {
public:

    static
    double slip(const double t);

    static
    double slipRate(const double t);

    static
    double slipAcc(const double t);

};

class SlipFnRamp {
public:

    static
    double slip(const double t);

    static
    double slipRate(const double t);

    static
    double slipAcc(const double t);

};

#endif

// End of file
