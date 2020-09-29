#if !defined(prescribedslip_hh)
#define prescribedslip_hh

#include <portinfo>

class PrescribedSlip {
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
