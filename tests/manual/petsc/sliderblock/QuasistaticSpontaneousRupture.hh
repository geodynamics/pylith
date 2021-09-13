#if !defined(quasistaticspontaneousrupture_hh)
#define quasistaticspontaneousrupture_hh

#include <portinfo>

#include "Formulation.hh" // ISA Formulation

class Friction;

class QuasistaticSpontaneousRupture : public Formulation {
public:

    QuasistaticSpontaneousRupture(const char* friction);
    ~QuasistaticSpontaneousRupture(void);

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

    void _updateState(const double dt);

    void _setSolutionBounds(PetscTS ts);

private:

    Friction* _friction;

private:

    // Not implemented
    QuasistaticSpontaneousRupture(void);
    QuasistaticSpontaneousRupture(const QuasistaticSpontaneousRupture&);
    const QuasistaticSpontaneousRupture& operator=(const QuasistaticSpontaneousRupture&);

};

#endif

// End of file
