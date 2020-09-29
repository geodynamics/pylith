#if !defined(formulation_hh)
#define formulation_hh

#include <portinfo>

#include "petsc.h"

typedef struct _p_Mat* PetscMat;
typedef struct _p_Vec* PetscVec;
typedef struct _p_TS* PetscTS;

#define CHECK_ERROR(err) do {if (PetscUnlikely(err)) {PetscError(PETSC_COMM_SELF,__LINE__,PETSC_FUNCTION_NAME,__FILE__,err,PETSC_ERROR_REPEAT,0);throw std::runtime_error("Error detected while in PETSc function.");}} while (0)

class Formulation {
public:

    Formulation(void);

    virtual ~Formulation(void);

    void initialize(const char* const outputFilename);

    void solve(void);

    static PetscErrorCode computeLHSResidual(PetscTS ts,
                                             PetscReal t,
                                             PetscVec solution,
                                             PetscVec solutionDot,
                                             PetscVec residual,
                                             void *context);

    static PetscErrorCode computeLHSJacobian(PetscTS ts,
                                             PetscReal t,
                                             PetscVec solution,
                                             PetscVec solutionDot,
                                             PetscReal shift,
                                             PetscMat jacobian,
                                             PetscMat preconditioner,
                                             void *context);

    static PetscErrorCode computeRHSResidual(PetscTS ts,
                                             PetscReal t,
                                             PetscVec solution,
                                             PetscVec residual,
                                             void *context);

    static PetscErrorCode computeRHSJacobian(PetscTS ts,
                                             PetscReal t,
                                             PetscVec solution,
                                             PetscMat jacobian,
                                             PetscMat preconditioner,
                                             void *context);

    static PetscErrorCode poststep(PetscTS ts);

protected:

    void _createMatsAndVecs(void);

    void _poststep(void);

    virtual
    void _computeLHSResidual(const PetscReal,
                             const PetscVec,
                             const PetscVec,
                             PetscVec);

    virtual
    void _computeRHSResidual(const PetscReal,
                             const PetscVec,
                             PetscVec);

    virtual
    void _computeLHSJacobian(const PetscReal,
                             const PetscVec,
                             const PetscVec,
                             const PetscReal,
                             PetscMat,
                             PetscMat);

    virtual
    void _computeRHSJacobian(const PetscReal,
                             const PetscVec,
                             PetscMat,
                             PetscMat);

    bool _hasLHSResidual;
    bool _hasRHSResidual;
    bool _hasLHSJacobian;
    bool _hasRHSJacobian;
    bool _isDynamic;
    TSType _tsAlgorithm;

private:

    PetscTS _ts;
    PetscVec _tstamp;
    PetscViewer _viewer;

    static const double _duration;
    static const double _dtInitial;

protected:

    PetscVec _solution;
    PetscMat _jacobianLHS;
    PetscMat _jacobianRHS;
    size_t _numDOFAll;

    static const double _ka;
    static const double _ma;
    static const double _kb;
    static const double _mb;

    static const size_t _numDOFDisp;

private:

    // Not implemented
    Formulation(const Formulation&);
    const Formulation& operator=(const Formulation&);

};

#endif

// End of file
