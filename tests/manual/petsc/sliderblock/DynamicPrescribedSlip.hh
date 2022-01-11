#if !defined(dynamicprescribedslip_hh)
#define dynamicprescribedslip_hh

#include <portinfo>

#include "Formulation.hh" // ISA Formulation

class DynamicPrescribedSlip : public Formulation {
public:

    DynamicPrescribedSlip(void);
    ~DynamicPrescribedSlip(void);

protected:

    void _computeLHSResidual(const PetscReal t,
                             const PetscReal dt,
                             const PetscVec solution,
                             const PetscVec solutionDot,
                             PetscVec residual);

    void _computeRHSResidual(const PetscReal t,
                             const PetscReal dt,
                             const PetscVec solution,
                             PetscVec residual);

    void _computeLHSJacobian(const PetscReal t,
			     const PetscReal dt,
                             const PetscVec solution,
                             const PetscVec solutionDot,
                             const PetscReal shift,
                             PetscMat jacobian,
                             PetscMat preconditioner);

private:

    // Not implemented
    DynamicPrescribedSlip(const DynamicPrescribedSlip&);
    const DynamicPrescribedSlip& operator=(const DynamicPrescribedSlip&);

};

#endif

// End of file
