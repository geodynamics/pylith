#if !defined(quasistaticprescribedslip_hh)
#define quasistaticprescribedslip_hh

#include <portinfo>

#include "Formulation.hh" // ISA Formulation

class QuasistaticPrescribedSlip : public Formulation {
public:

    QuasistaticPrescribedSlip(void);
    ~QuasistaticPrescribedSlip(void);

protected:

    void _computeLHSResidual(const PetscReal t,
                             const PetscVec solution,
                             const PetscVec solutionDot,
                             PetscVec residual);

    void _computeLHSJacobian(const PetscReal t,
                             const PetscVec solution,
                             const PetscVec solutionDot,
                             const PetscReal stshift,
                             PetscMat jacobian,
                             PetscMat preconditioner);

private:

    // Not implemented
    QuasistaticPrescribedSlip(const QuasistaticPrescribedSlip&);
    const QuasistaticPrescribedSlip& operator=(const QuasistaticPrescribedSlip&);

};

#endif

// End of file
