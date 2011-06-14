// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "SolverNonlinear.hh" // implementation of class methods

#include "Formulation.hh" // USES Formulation
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include <petscsnes.h> // USES PetscSNES

#include "pylith/utils/petscerror.h" // USES CHECK_PETSC_ERROR

// KLUDGE, Fixes issue with PetscIsInfOrNanReal and include cmath
// instead of math.h.
#define isnan std::isnan // TEMPORARY
#define isinf std::isinf // TEMPORARY

#include <private/snesimpl.h>

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

namespace pylith {
  namespace problems {
    namespace _SolverNonlinear {

      typedef struct {
	PetscErrorCode           (*LineSearch)(PetscSNES,void*,PetscVec,PetscVec,PetscVec,PetscVec,PetscVec,PetscReal,PetscReal,PetscReal*,PetscReal*,PetscBool *);
	void                     *lsP;                              /* user-defined line-search context (optional) */
	/* --------------- Parameters used by line search method ----------------- */
	PetscReal                alpha;                                                                  /* used to determine sufficient reduction */
	PetscReal                maxstep;                                                          /* maximum step size */
	PetscReal                minlambda;                                                        /* determines smallest line search lambda used */
	PetscErrorCode           (*precheckstep)(PetscSNES,PetscVec,PetscVec,void*,PetscBool *);                  /* step-checking routine (optional) */
	void                     *precheck;                                                        /* user-defined step-checking context (optional) */
	PetscErrorCode           (*postcheckstep)(PetscSNES,PetscVec,PetscVec,PetscVec,void*,PetscBool *,PetscBool *); /* step-checking routine (optional) */
	void                     *postcheck;                                                       /* user-defined step-checking context (optional) */
	PetscViewerASCIIMonitor  monitor;
      } SNES_LS;
    } // _SolverNonlinear
  } // problems
} // pylith

// ----------------------------------------------------------------------
// Constructor
pylith::problems::SolverNonlinear::SolverNonlinear(void) :
  _snes(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::problems::SolverNonlinear::~SolverNonlinear(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate data structures.
void
pylith::problems::SolverNonlinear::deallocate(void)
{ // deallocate
  Solver::deallocate();

  if (0 != _snes) {
    PetscErrorCode err = SNESDestroy(&_snes); _snes = 0;
    CHECK_PETSC_ERROR(err);
  } // if
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize solver.
void
pylith::problems::SolverNonlinear::initialize(
			             const topology::SolutionFields& fields,
				     const topology::Jacobian& jacobian,
				     Formulation* formulation)
{ // initialize
  assert(0 != formulation);

  _initializeLogger();
  Solver::initialize(fields, jacobian, formulation);

  PetscErrorCode err = 0;
  if (0 != _snes) {
    err = SNESDestroy(&_snes); _snes = 0;
    CHECK_PETSC_ERROR(err);
  } // if    
  err = SNESCreate(fields.mesh().comm(), &_snes); CHECK_PETSC_ERROR(err);

  const topology::Field<topology::Mesh>& residual = fields.get("residual");
  const PetscVec residualVec = residual.vector();
  err = SNESSetFunction(_snes, residualVec, reformResidual,
			(void*) formulation);
  CHECK_PETSC_ERROR(err);

  err = SNESSetJacobian(_snes, jacobian.matrix(), _jacobianPC, reformJacobian, (void*) formulation);
  CHECK_PETSC_ERROR(err);

  err = SNESSetFromOptions(_snes); CHECK_PETSC_ERROR(err);
  err = SNESLineSearchSet(_snes, lineSearch, (void*) formulation); CHECK_PETSC_ERROR(err);

  if (formulation->splitFields()) {
    PetscKSP ksp = 0;
    PetscPC pc = 0;
    err = SNESGetKSP(_snes, &ksp); CHECK_PETSC_ERROR(err);
    err = KSPGetPC(ksp, &pc); CHECK_PETSC_ERROR(err);
    _setupFieldSplit(&pc, formulation, jacobian, fields);
  } // if
} // initialize

// ----------------------------------------------------------------------
// Solve the system.
void
pylith::problems::SolverNonlinear::solve(
			      topology::Field<topology::Mesh>* solution,
			      topology::Jacobian* jacobian,
			      const topology::Field<topology::Mesh>& residual)
{ // solve
  assert(0 != solution);

#if 0
  // Update KSP operators with custom preconditioner if necessary.
  const ALE::Obj<RealSection>& solutionSection = solution->section();
  assert(!solutionSection.isNull());
  const ALE::Obj<SieveMesh>& sieveMesh = solution->mesh().sieveMesh();
   assert(!sieveMesh.isNull());
  if (solutionSection->getNumSpaces() > sieveMesh->getDimension() &&
      0 != _jacobianPCFault) {
    PetscKSP ksp = 0;
    PetscPC pc = 0;
    PetscKSP *ksps = 0;
    PetscMat A = 0;
    PetscInt num = 0;

    PetscErrorCode err = 0;
    err = SNESGetKSP(_snes, &ksp); CHECK_PETSC_ERROR(err);
    err = KSPSetUp(ksp); CHECK_PETSC_ERROR(err);
    err = KSPGetPC(ksp, &pc); CHECK_PETSC_ERROR(err);
    err = PCFieldSplitGetSubKSP(pc, &num, &ksps); CHECK_PETSC_ERROR(err);
    assert(solutionSection->getNumSpaces() == num);

#if 0 // debugging
    std::cout << "Preconditioner Matrix" << std::endl;
    MatView(_jacobianPCFault, PETSC_VIEWER_STDOUT_WORLD);
#endif


    MatStructure flag;
    err = KSPGetOperators(ksps[num-1], &A, 
			  PETSC_NULL, &flag); CHECK_PETSC_ERROR(err);
    err = PetscObjectReference((PetscObject) A); CHECK_PETSC_ERROR(err);
    err = KSPSetOperators(ksps[num-1], A, _jacobianPCFault, 
			  flag); CHECK_PETSC_ERROR(err);
    err = PetscFree(ksps); CHECK_PETSC_ERROR(err);
  } // if
#endif

  const int solveEvent = _logger->eventId("SoNl solve");
  const int scatterEvent = _logger->eventId("SoNl scatter");
  _logger->eventBegin(solveEvent);

  PetscErrorCode err = 0;
  const PetscVec solutionVec = solution->vector();

  err = SNESSolve(_snes, PETSC_NULL, solutionVec); CHECK_PETSC_ERROR(err);
  
  _logger->eventEnd(solveEvent);
  _logger->eventBegin(scatterEvent);

  // Update section view of field.
  solution->scatterVectorToSection();

  _logger->eventEnd(scatterEvent);

  // Update rate fields to be consistent with current solution.
  _formulation->calcRateFields();
} // solve

// ----------------------------------------------------------------------
// Generic C interface for reformResidual for integration with
// PETSc SNES solvers.
PetscErrorCode
pylith::problems::SolverNonlinear::reformResidual(PetscSNES snes,
						  PetscVec tmpSolutionVec,
						  PetscVec tmpResidualVec,
						  void* context)
{ // reformResidual
  assert(0 != context);
  Formulation* formulation = (Formulation*) context;
  assert(0 != formulation);

  // Make sure we have an admissible Lagrange force (\lambda)
  formulation->constrainSolnSpace(&tmpSolutionVec);

  // Reform residual
  formulation->reformResidual(&tmpResidualVec, &tmpSolutionVec);

  return 0;
} // reformResidual

// ----------------------------------------------------------------------
// Generic C interface for reformJacobian for integration with
// PETSc SNES solvers.
PetscErrorCode
pylith::problems::SolverNonlinear::reformJacobian(PetscSNES snes,
						  PetscVec tmpSolutionVec,
						  PetscMat* jacobianMat,
						  PetscMat* preconditionerMat,
						  MatStructure* preconditionerLayout,
						  void* context)
{ // reformJacobian
  assert(0 != context);
  Formulation* formulation = (Formulation*) context;
  assert(0 != formulation);

  formulation->reformJacobian(&tmpSolutionVec);

  return 0;
} // reformJacobian

// ----------------------------------------------------------------------
// Generic C interface for customized PETSc line search.
#undef __FUNCT__
#define __FUNCT__ "lineSearch"
PetscErrorCode
pylith::problems::SolverNonlinear::lineSearch(PetscSNES snes,
					      void *lsctx,
					      PetscVec x,
					      PetscVec f,
					      PetscVec g,
					      PetscVec y,
					      PetscVec w,
					      PetscReal fnorm,
					      PetscReal xnorm,
					      PetscReal *ynorm,
					      PetscReal *gnorm,
					      PetscBool *flag)
{ // lineSearch
  // Note that for line search purposes we work with with the related
  // minimization problem:
  // min  z(x):  R^n -> R,
  // where z(x) = .5 * fnorm*fnorm, and fnorm = || f ||_2.
  typedef pylith::problems::_SolverNonlinear::SNES_LS SNES_LS;
       
  PetscReal      initslope,lambdaprev,gnormprev,a,b,d,t1,t2,rellength;
  PetscReal      minlambda,lambda,lambdatemp;
#if defined(PETSC_USE_COMPLEX)
  PetscScalar    cinitslope;
#endif
  PetscErrorCode ierr;
  PetscInt       count;
  SNES_LS        *neP = (SNES_LS*)snes->data;
  PetscBool     changed_w = PETSC_FALSE,changed_y = PETSC_FALSE;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(SNES_LineSearch,snes,x,f,g);CHKERRQ(ierr);
  *flag   = PETSC_TRUE;

  ierr = VecNorm(y,NORM_2,ynorm);CHKERRQ(ierr);
  if (!*ynorm) {
    ierr = PetscInfo(snes,"Search direction and size is 0\n");CHKERRQ(ierr);
    *gnorm = fnorm;
    ierr   = VecCopy(x,w);CHKERRQ(ierr);
    ierr   = VecCopy(f,g);CHKERRQ(ierr);
    *flag  = PETSC_FALSE;
    goto theend1;
  }
  if (*ynorm > neP->maxstep) {	/* Step too big, so scale back */
    ierr = PetscInfo2(snes,"Scaling step by %G old ynorm %G\n",
		      neP->maxstep/(*ynorm),*ynorm);CHKERRQ(ierr);
    ierr = VecScale(y,neP->maxstep/(*ynorm));CHKERRQ(ierr);
    *ynorm = neP->maxstep;
  }
  ierr      = VecMaxPointwiseDivide(y,x,&rellength);CHKERRQ(ierr);
  minlambda = neP->minlambda/rellength;
  ierr      = MatMult(snes->jacobian,y,w);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
  ierr      = VecDot(f,w,&cinitslope);CHKERRQ(ierr);
  initslope = PetscRealPart(cinitslope);
#else
  ierr      = VecDot(f,w,&initslope);CHKERRQ(ierr);
#endif
  if (initslope > 0.0)  initslope = -initslope;
  if (initslope == 0.0) initslope = -1.0;

  ierr = VecWAXPY(w,-1.0,y,x);CHKERRQ(ierr);
  if (snes->nfuncs >= snes->max_funcs) {
    ierr  = PetscInfo(snes,"Exceeded maximum function evaluations, "
		      "while checking full step length!\n");CHKERRQ(ierr);
    *flag = PETSC_FALSE;
    snes->reason = SNES_DIVERGED_FUNCTION_COUNT;
    goto theend1;
  } // if
  // TEMPORARY: update w?

  ierr = SNESComputeFunction(snes,w,g);CHKERRQ(ierr);
  if (snes->domainerror) {
    ierr = PetscLogEventEnd(SNES_LineSearch,snes,x,f,g);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  ierr = VecNorm(g,NORM_2,gnorm);CHKERRQ(ierr);
  if (PetscIsInfOrNanReal(*gnorm))
    SETERRQ(PETSC_COMM_SELF,
	    PETSC_ERR_FP, "User provided compute function generated a "
	    "Not-a-Number");
  ierr = PetscInfo2(snes,"Initial fnorm %G gnorm %G\n",fnorm,*gnorm);CHKERRQ(ierr);
  if (.5*(*gnorm)*(*gnorm) <= .5*fnorm*fnorm + neP->alpha*initslope) {
    // Sufficient reduction
    ierr = PetscInfo2(snes,"Using full step: fnorm %G gnorm %G\n",
		      fnorm,*gnorm);CHKERRQ(ierr);
    goto theend1;
  } // if

  // Fit points with quadratic
  lambda     = 1.0;
  lambdatemp = -initslope/((*gnorm)*(*gnorm) - fnorm*fnorm - 2.0*initslope);
  lambdaprev = lambda;
  gnormprev  = *gnorm;
  if (lambdatemp > .5*lambda)
    lambdatemp = .5*lambda;
  if (lambdatemp <= .1*lambda)
    lambda = .1*lambda; 
  else
    lambda = lambdatemp;

  ierr  = VecWAXPY(w,-lambda,y,x);CHKERRQ(ierr);
  if (snes->nfuncs >= snes->max_funcs) {
    ierr  = PetscInfo1(snes,"Exceeded maximum function evaluations, while attempting quadratic backtracking! %D \n",snes->nfuncs);CHKERRQ(ierr);
    *flag = PETSC_FALSE;
    snes->reason = SNES_DIVERGED_FUNCTION_COUNT;
    goto theend1;
  } // if
  // TEMPORARY: update w?

  ierr = SNESComputeFunction(snes,w,g);CHKERRQ(ierr);
  if (snes->domainerror) {
    ierr = PetscLogEventEnd(SNES_LineSearch,snes,x,f,g);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  } // if
  ierr = VecNorm(g,NORM_2,gnorm);CHKERRQ(ierr);
  if (PetscIsInfOrNanReal(*gnorm)) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FP,"User provided compute function generated a Not-a-Number");
  ierr = PetscInfo1(snes,"gnorm after quadratic fit %G\n",*gnorm);CHKERRQ(ierr);
  if (.5*(*gnorm)*(*gnorm) < .5*fnorm*fnorm + lambda*neP->alpha*initslope) { /* sufficient reduction */
    ierr = PetscInfo1(snes,"Quadratically determined step, lambda=%18.16e\n",lambda);CHKERRQ(ierr);
    goto theend1;
  } // if

  // Fit points with cubic
  count = 1;
  while (PETSC_TRUE) {
    if (lambda <= minlambda) { 
      ierr = PetscInfo1(snes,"Unable to find good step length! After %D "
			"tries \n",count);CHKERRQ(ierr);
      ierr = PetscInfo6(snes,"fnorm=%18.16e, gnorm=%18.16e, ynorm=%18.16e, "
			"minlambda=%18.16e, lambda=%18.16e, "
			"initial slope=%18.16e\n",
			fnorm,*gnorm,*ynorm,minlambda,lambda,initslope);
      CHKERRQ(ierr);
      ierr = PetscInfo1(snes,"Using last lambda tried %g\n",lambda);CHKERRQ(ierr);
      break;
    }
    t1 = .5*((*gnorm)*(*gnorm) - fnorm*fnorm) - lambda*initslope;
    t2 = .5*(gnormprev*gnormprev  - fnorm*fnorm) - lambdaprev*initslope;
    a  = (t1/(lambda*lambda) - t2/(lambdaprev*lambdaprev))/(lambda-lambdaprev);
    b  = (-lambdaprev*t1/(lambda*lambda) + 
	  lambda*t2/(lambdaprev*lambdaprev))/(lambda-lambdaprev);
    d  = b*b - 3*a*initslope;
    if (d < 0.0) d = 0.0;
    if (a == 0.0) {
      lambdatemp = -initslope/(2.0*b);
    } else {
      lambdatemp = (-b + sqrt(d))/(3.0*a);
    } // if/else
    lambdaprev = lambda;
    gnormprev  = *gnorm;
    if (lambdatemp > .5*lambda)
      lambdatemp = .5*lambda;
    if (lambdatemp <= .1*lambda)
      lambda = .1*lambda;
    else
      lambda = lambdatemp;
    ierr  = VecWAXPY(w,-lambda,y,x);CHKERRQ(ierr);
    if (snes->nfuncs >= snes->max_funcs) {
      ierr = PetscInfo1(snes,"Exceeded maximum function evaluations, "
			"while looking for good step length! %D \n",count);
      CHKERRQ(ierr);
      ierr = PetscInfo5(snes,"fnorm=%18.16e, gnorm=%18.16e, ynorm=%18.16e, "
			"lambda=%18.16e, initial slope=%18.16e\n",
			fnorm,*gnorm,*ynorm,lambda,initslope);CHKERRQ(ierr);
      *flag = PETSC_FALSE;
      snes->reason = SNES_DIVERGED_FUNCTION_COUNT;
      break;
    } // if
    // TEMPORARY: update w?

    ierr = SNESComputeFunction(snes,w,g);CHKERRQ(ierr);
    if (snes->domainerror) {
      ierr = PetscLogEventEnd(SNES_LineSearch,snes,x,f,g);CHKERRQ(ierr);
      PetscFunctionReturn(0);
    } // if
    ierr = VecNorm(g,NORM_2,gnorm);CHKERRQ(ierr);
    if (PetscIsInfOrNanReal(*gnorm))
      SETERRQ(PETSC_COMM_SELF, 
	      PETSC_ERR_FP,"User provided compute function generated a "
	      "Not-a-Number");
    if (.5*(*gnorm)*(*gnorm) < .5*fnorm*fnorm + lambda*neP->alpha*initslope) {
      // is reduction enough?
      ierr = PetscInfo2(snes,"Cubically determined step, current gnorm %G "
			"lambda=%18.16e\n",*gnorm,lambda);CHKERRQ(ierr);
      break;
    } else {
      ierr = PetscInfo2(snes,"Cubic step no good, shrinking lambda, "
			"current gnorem %G lambda=%18.16e\n",
			*gnorm,lambda);CHKERRQ(ierr);
    } // if/else
    count++;
  } // while

  theend1:
  /* Optional user-defined check for line search step validity */
  if (neP->postcheckstep && *flag) {
    ierr = (*neP->postcheckstep)(snes,x,y,w,neP->postcheck,&changed_y,
				 &changed_w);CHKERRQ(ierr);
    if (changed_y)
      ierr = VecWAXPY(w,-1.0,y,x);CHKERRQ(ierr);
    // TEMPORARY: update w?

    if (changed_y || changed_w) {
      // recompute the function if the step has changed 
      ierr = SNESComputeFunction(snes,w,g);CHKERRQ(ierr);
      if (snes->domainerror) {
        ierr = PetscLogEventEnd(SNES_LineSearch,snes,x,f,g);CHKERRQ(ierr);
        PetscFunctionReturn(0);
      }
      ierr = VecNormBegin(g,NORM_2,gnorm);CHKERRQ(ierr);
      if (PetscIsInfOrNanReal(*gnorm))
	SETERRQ(PETSC_COMM_SELF, 
		PETSC_ERR_FP, "User provided compute function generated a "
		"Not-a-Number");
      ierr = VecNormBegin(y,NORM_2,ynorm);CHKERRQ(ierr);
      ierr = VecNormEnd(g,NORM_2,gnorm);CHKERRQ(ierr);
      ierr = VecNormEnd(y,NORM_2,ynorm);CHKERRQ(ierr);
    } // if
  } // if
  ierr = PetscLogEventEnd(SNES_LineSearch,snes,x,f,g);CHKERRQ(ierr);

  // ======================================================================
  // Code to constrain solution space.
  assert(0 != lsctx);
  Formulation* formulation = (Formulation*) lsctx;
  assert(0 != formulation);
  formulation->constrainSolnSpace(&w);

  ierr = SNESComputeFunction(snes,w,g);CHKERRQ(ierr);
  if (snes->domainerror) {
    ierr = PetscLogEventEnd(SNES_LineSearch,snes,x,f,g);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  ierr = VecNormBegin(g,NORM_2,gnorm);CHKERRQ(ierr);
  if (PetscIsInfOrNanReal(*gnorm))
    SETERRQ(PETSC_COMM_SELF,
	    PETSC_ERR_FP, "User provided compute function generated a "
	    "Not-a-Number");
  ierr = VecNormBegin(y,NORM_2,ynorm);CHKERRQ(ierr);
  ierr = VecNormEnd(g,NORM_2,gnorm);CHKERRQ(ierr);
  ierr = VecNormEnd(y,NORM_2,ynorm);CHKERRQ(ierr);
  // ======================================================================

  PetscFunctionReturn(0);
} // lineSearch

// ----------------------------------------------------------------------
// Initialize logger.
void
pylith::problems::SolverNonlinear::_initializeLogger(void)
{ // initializeLogger
  delete _logger; _logger = new utils::EventLogger;
  assert(0 != _logger);
  _logger->className("SolverNonlinear");
  _logger->initialize();
  _logger->registerEvent("SoNl setup");
  _logger->registerEvent("SoNl solve");
  _logger->registerEvent("SoNl scatter");
} // initializeLogger


// End of file
