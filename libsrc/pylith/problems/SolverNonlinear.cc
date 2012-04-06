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
// Copyright (c) 2010-2012 University of California, Davis
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
#include "pylith/utils/constdefs.h" // USES PYLITH_MAXSCALAR
#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include <petscsnes.h> // USES PetscSNES

#include "pylith/utils/petscerror.h" // USES CHECK_PETSC_ERROR

// KLUDGE, Fixes issue with PetscIsInfOrNanReal and include cmath
// instead of math.h.
#define isnan std::isnan // TEMPORARY
#define isinf std::isinf // TEMPORARY

// Customized line search based on PETSc code in
// src/snes/linesearch/bt/linesearchbt.c.
#include <petsc-private/snesimpl.h>
#include <petsc-private/linesearchimpl.h>

typedef enum {
  SNES_LINESEARCH_BT_QUADRATIC, 
  SNES_LINESEARCH_BT_CUBIC,
} PetscSNESLineSearchBTOrder;

struct PetscSNESLineSearch_BT {
  PetscReal        alpha; /* sufficient decrease parameter */
  PetscSNESLineSearchBTOrder order;
};

#define PYLITH_CUSTOM_LINESEARCH 1

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

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

  PetscErrorCode err = SNESDestroy(&_snes);CHECK_PETSC_ERROR(err);
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
  err = SNESSetComputeInitialGuess(_snes, initialGuess, (void*) formulation); CHECK_PETSC_ERROR(err);
  PetscSNESLineSearch ls;

  err = SNESGetSNESLineSearch(_snes, &ls); CHECK_PETSC_ERROR(err);
  err = SNESLineSearchSetType(ls, SNESLINESEARCHSHELL); CHECK_PETSC_ERROR(err);
  err = SNESLineSearchShellSetUserFunc(ls, lineSearch, (void*) formulation); CHECK_PETSC_ERROR(err);

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
pylith::problems::SolverNonlinear::lineSearch(PetscSNESLineSearch linesearch,
					      void* lsctx)
{ // lineSearch
  // Note that for line search purposes we work with with the related
  // minimization problem:
  // min  z(x):  R^n -> R,
  // where z(x) = .5 * fnorm*fnorm, and fnorm = || f ||_2.
  //
  // Customized version of src/snes/linesearch/impls/bt/linesearchbt.c.

  PetscBool      changed_y =PETSC_FALSE, changed_w =PETSC_FALSE;
  PetscErrorCode ierr;
  Vec            X, F, Y, W, G;
  SNES           snes;
  PetscReal      fnorm, xnorm, ynorm, gnorm, gnormprev;
  PetscReal      lambda, lambdatemp, lambdaprev, minlambda, maxstep, rellength, initslope, alpha;
  PetscReal      t1, t2, a, b, d, steptol;
#if defined(PETSC_USE_COMPLEX)
  PetscScalar    cinitslope;
#endif
  PetscBool      domainerror;
  PetscViewer    monitor;
  PetscInt       max_its, count;
  PetscSNESLineSearch_BT  *bt;
  Mat            jac;


  PetscFunctionBegin;

  ierr = SNESLineSearchGetVecs(linesearch, &X, &F, &Y, &W, &G);CHKERRQ(ierr);
  ierr = SNESLineSearchGetNorms(linesearch, &xnorm, &fnorm, &ynorm);CHKERRQ(ierr);
  ierr = SNESLineSearchGetLambda(linesearch, &lambda);CHKERRQ(ierr);
  ierr = SNESLineSearchGetSNES(linesearch, &snes);CHKERRQ(ierr);
  ierr = SNESLineSearchGetMonitor(linesearch, &monitor);CHKERRQ(ierr);
  ierr = SNESLineSearchGetTolerances(linesearch, &steptol, &maxstep, PETSC_NULL, PETSC_NULL, PETSC_NULL, &max_its);
  bt = (PetscSNESLineSearch_BT *)linesearch->data;

  alpha = bt->alpha;

  ierr = SNESGetJacobian(snes, &jac, PETSC_NULL, PETSC_NULL, PETSC_NULL);CHKERRQ(ierr);
  if (!jac) {
    SETERRQ(((PetscObject)linesearch)->comm, PETSC_ERR_USER, "SNESLineSearchBT requires a Jacobian matrix");
  }
  /* precheck */
  ierr = SNESLineSearchPreCheck(linesearch, &changed_y);CHKERRQ(ierr);
  ierr = SNESLineSearchSetSuccess(linesearch, PETSC_TRUE);CHKERRQ(ierr);

  ierr = VecNorm(Y, NORM_2, &ynorm);CHKERRQ(ierr);
  if (ynorm == 0.0) {
    if (monitor) {
      ierr = PetscViewerASCIIAddTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(monitor,"    Line search: Initial direction and size is 0\n");CHKERRQ(ierr);
      ierr = PetscViewerASCIISubtractTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
    }
    ierr   = VecCopy(X,W);CHKERRQ(ierr);
    ierr   = VecCopy(F,G);CHKERRQ(ierr);
    ierr = SNESLineSearchSetSuccess(linesearch, PETSC_FALSE);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  if (ynorm > maxstep) {	/* Step too big, so scale back */
    if (monitor) {
      ierr = PetscViewerASCIIAddTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(monitor,"    Line search: Scaling step by %14.12e old ynorm %14.12e\n", (maxstep/ynorm),ynorm);CHKERRQ(ierr);
      ierr = PetscViewerASCIISubtractTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
    }
    ierr = VecScale(Y,maxstep/(ynorm));CHKERRQ(ierr);
    ynorm = maxstep;
  }
  ierr      = VecMaxPointwiseDivide(Y,X,&rellength);CHKERRQ(ierr);

#if defined(PYLITH_CUSTOM_LINESEARCH)
  // Place reasonable absolute limit on minimum lambda
  minlambda = std::max(steptol/rellength, 1.0/PYLITH_MAXSCALAR);
#else // ORIGINAL
  minlambda = steptol/rellength;
#endif

  ierr      = MatMult(jac,Y,W);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
  ierr      = VecDot(F,W,&cinitslope);CHKERRQ(ierr);
  initslope = PetscRealPart(cinitslope);
#else
  ierr      = VecDot(F,W,&initslope);CHKERRQ(ierr);
#endif
  if (initslope > 0.0)  initslope = -initslope;
  if (initslope == 0.0) initslope = -1.0;

  ierr = VecWAXPY(W,-lambda,Y,X);CHKERRQ(ierr);
  if (linesearch->ops->viproject) {
    ierr = (*linesearch->ops->viproject)(snes, W);CHKERRQ(ierr);
  }
  if (snes->nfuncs >= snes->max_funcs) {
    ierr  = PetscInfo(snes,"Exceeded maximum function evaluations, while checking full step length!\n");CHKERRQ(ierr);
    snes->reason = SNES_DIVERGED_FUNCTION_COUNT;
    ierr = SNESLineSearchSetSuccess(linesearch, PETSC_FALSE);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  ierr = SNESComputeFunction(snes,W,G);CHKERRQ(ierr);
  ierr = SNESGetFunctionDomainError(snes, &domainerror);CHKERRQ(ierr);
  if (domainerror) {
    ierr = SNESLineSearchSetSuccess(linesearch, PETSC_FALSE);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  if (linesearch->ops->vinorm) {
    gnorm = fnorm;
    ierr = (*linesearch->ops->vinorm)(snes, G, W, &gnorm);CHKERRQ(ierr);
  } else {
    ierr = VecNorm(G,NORM_2,&gnorm);CHKERRQ(ierr);
  }

  if (PetscIsInfOrNanReal(gnorm)) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FP,"User provided compute function generated a Not-a-Number");
  ierr = PetscInfo2(snes,"Initial fnorm %14.12e gnorm %14.12e\n", fnorm, gnorm);CHKERRQ(ierr);
  if (.5*gnorm*gnorm <= .5*fnorm*fnorm + lambda*alpha*initslope) { /* Sufficient reduction */
    if (monitor) {
      ierr = PetscViewerASCIIAddTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(monitor,"    Line search: Using full step: fnorm %14.12e gnorm %14.12e\n", fnorm, gnorm);CHKERRQ(ierr);
      ierr = PetscViewerASCIISubtractTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
    }
  } else {
    /* Fit points with quadratic */
    lambdatemp = -initslope/(gnorm*gnorm - fnorm*fnorm - 2.0*lambda*initslope);
    lambdaprev = lambda;
    gnormprev  = gnorm;
    if (lambdatemp > .5*lambda)  lambdatemp = .5*lambda;
    if (lambdatemp <= .1*lambda) lambda = .1*lambda;
    else                         lambda = lambdatemp;

    ierr  = VecWAXPY(W,-lambda,Y,X);CHKERRQ(ierr);
    if (linesearch->ops->viproject) {
      ierr = (*linesearch->ops->viproject)(snes, W);CHKERRQ(ierr);
    }
    if (snes->nfuncs >= snes->max_funcs) {
      ierr  = PetscInfo1(snes,"Exceeded maximum function evaluations, while attempting quadratic backtracking! %D \n",snes->nfuncs);CHKERRQ(ierr);
      snes->reason = SNES_DIVERGED_FUNCTION_COUNT;
      ierr = SNESLineSearchSetSuccess(linesearch, PETSC_FALSE);CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
    ierr = SNESComputeFunction(snes,W,G);CHKERRQ(ierr);
    ierr = SNESGetFunctionDomainError(snes, &domainerror);CHKERRQ(ierr);
    if (domainerror) {
      PetscFunctionReturn(0);
    }
    if (linesearch->ops->vinorm) {
      gnorm = fnorm;
      ierr = (*linesearch->ops->vinorm)(snes, G, W, &gnorm);CHKERRQ(ierr);
    } else {
      ierr = VecNorm(G,NORM_2,&gnorm);CHKERRQ(ierr);
    }
    if (PetscIsInfOrNanReal(gnorm)) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FP,"User provided compute function generated a Not-a-Number");
    if (monitor) {
      ierr = PetscViewerASCIIAddTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(monitor,"    Line search: gnorm after quadratic fit %14.12e\n",gnorm);CHKERRQ(ierr);
      ierr = PetscViewerASCIISubtractTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
    }
    if (.5*gnorm*gnorm < .5*fnorm*fnorm + lambda*alpha*initslope) { /* sufficient reduction */
      if (monitor) {
        ierr = PetscViewerASCIIAddTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(monitor,"    Line search: Quadratically determined step, lambda=%18.16e\n",(double)lambda);CHKERRQ(ierr);
        ierr = PetscViewerASCIISubtractTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
      }
    } else {
      /* Fit points with cubic */
      for (count = 0; count < max_its; count++) {
        if (lambda <= minlambda) {
          if (monitor) {
            ierr = PetscViewerASCIIAddTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
            ierr = PetscViewerASCIIPrintf(monitor,"    Line search: unable to find good step length! After %D tries \n",count);CHKERRQ(ierr);
            ierr = PetscViewerASCIIPrintf(monitor,
                                          "    Line search: fnorm=%18.16e, gnorm=%18.16e, ynorm=%18.16e, minlambda=%18.16e, lambda=%18.16e, initial slope=%18.16e\n",
                                          fnorm, gnorm, ynorm, minlambda, lambda, initslope);CHKERRQ(ierr);
            ierr = PetscViewerASCIISubtractTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
          }
#if defined(PYLITH_CUSTOM_LINESEARCH)
#if 0 // DEBUGGING
	  assert(lsctx);
	  Formulation* formulation = (Formulation*) lsctx;
	  assert(formulation);
	  formulation->printState(&w, &g, &x, &y);
	  std::cerr << "WARNING: Line search diverged ... continuing nonlinear iterations anyway in hopes that solution will converge anyway."
		    << std::endl;
#endif
	  break;
#else // ORIGINAL
          ierr = SNESLineSearchSetSuccess(linesearch, PETSC_FALSE);CHKERRQ(ierr);
          PetscFunctionReturn(0);
#endif
        }
        if (bt->order == SNES_LINESEARCH_BT_CUBIC) {
          t1 = .5*(gnorm*gnorm - fnorm*fnorm) - lambda*initslope;
          t2 = .5*(gnormprev*gnormprev  - fnorm*fnorm) - lambdaprev*initslope;
          a  = (t1/(lambda*lambda) - t2/(lambdaprev*lambdaprev))/(lambda-lambdaprev);
          b  = (-lambdaprev*t1/(lambda*lambda) + lambda*t2/(lambdaprev*lambdaprev))/(lambda-lambdaprev);
          d  = b*b - 3*a*initslope;
          if (d < 0.0) d = 0.0;
          if (a == 0.0) {
            lambdatemp = -initslope/(2.0*b);
          } else {
            lambdatemp = (-b + PetscSqrtReal(d))/(3.0*a);
          }
        } else if (bt->order == SNES_LINESEARCH_BT_QUADRATIC) {
          lambdatemp = -initslope/(gnorm*gnorm - fnorm*fnorm - 2.0*initslope);
        }
          lambdaprev = lambda;
          gnormprev  = gnorm;
        if (lambdatemp > .5*lambda)  lambdatemp = .5*lambda;
        if (lambdatemp <= .1*lambda) lambda     = .1*lambda;
        else                         lambda     = lambdatemp;
        ierr  = VecWAXPY(W,-lambda,Y,X);CHKERRQ(ierr);
        if (snes->nfuncs >= snes->max_funcs) {
          ierr = PetscInfo1(snes,"Exceeded maximum function evaluations, while looking for good step length! %D \n",count);CHKERRQ(ierr);
          ierr = PetscInfo5(snes,"fnorm=%18.16e, gnorm=%18.16e, ynorm=%18.16e, lambda=%18.16e, initial slope=%18.16e\n",
                            fnorm,gnorm,ynorm,lambda,initslope);CHKERRQ(ierr);
          ierr = SNESLineSearchSetSuccess(linesearch, PETSC_FALSE);CHKERRQ(ierr);
          snes->reason = SNES_DIVERGED_FUNCTION_COUNT;
          PetscFunctionReturn(0);
        }
        ierr = SNESComputeFunction(snes,W,G);CHKERRQ(ierr);
        ierr = SNESGetFunctionDomainError(snes, &domainerror);CHKERRQ(ierr);
        if (domainerror) {
          ierr = SNESLineSearchSetSuccess(linesearch, PETSC_FALSE);CHKERRQ(ierr);
          PetscFunctionReturn(0);
        }
        if (linesearch->ops->vinorm) {
          gnorm = fnorm;
          ierr = (*linesearch->ops->vinorm)(snes, G, W, &gnorm);CHKERRQ(ierr);
        } else {
          ierr = VecNorm(G,NORM_2,&gnorm);CHKERRQ(ierr);
        }
        if (PetscIsInfOrNanReal(gnorm)) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FP,"User provided compute function generated a Not-a-Number");
        if (.5*gnorm*gnorm < .5*fnorm*fnorm + lambda*alpha*initslope) { /* is reduction enough? */
          if (monitor) {
            ierr = PetscViewerASCIIAddTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
            if (bt->order == SNES_LINESEARCH_BT_CUBIC) {
              ierr = PetscViewerASCIIPrintf(monitor,"    Line search: Cubically determined step, current gnorm %14.12e lambda=%18.16e\n",gnorm,lambda);CHKERRQ(ierr);
            } else {
              ierr = PetscViewerASCIIPrintf(monitor,"    Line search: Quadratically determined step, current gnorm %14.12e lambda=%18.16e\n",gnorm,lambda);CHKERRQ(ierr);
            }
            ierr = PetscViewerASCIISubtractTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
          }
          break;
        } else {
          if (monitor) {
            ierr = PetscViewerASCIIAddTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
            if (bt->order == SNES_LINESEARCH_BT_CUBIC) {
              ierr = PetscViewerASCIIPrintf(monitor,"    Line search: Cubic step no good, shrinking lambda, current gnorm %12.12e lambda=%18.16e\n",gnorm,lambda);CHKERRQ(ierr);
            } else {
              ierr = PetscViewerASCIIPrintf(monitor,"    Line search: Quadratic step no good, shrinking lambda, current gnorm %12.12e lambda=%18.16e\n",gnorm,lambda);CHKERRQ(ierr);
            }
            ierr = PetscViewerASCIISubtractTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
          }
        }
      }
    }
  }

  /* postcheck */
  ierr = SNESLineSearchPostCheck(linesearch, &changed_y, &changed_w);CHKERRQ(ierr);
  if (changed_y) {
    ierr = VecWAXPY(W,-lambda,Y,X);CHKERRQ(ierr);
    if (linesearch->ops->viproject) {
      ierr = (*linesearch->ops->viproject)(snes, W);CHKERRQ(ierr);
    }
  }
#if defined(PYLITH_CUSTOM_LINESEARCH)
  if (true) {
#else // ORIGINAL
  if (changed_y || changed_w) { /* recompute the function if the step has changed */
#endif
    ierr = SNESComputeFunction(snes,W,G);CHKERRQ(ierr);
    ierr = SNESGetFunctionDomainError(snes, &domainerror);CHKERRQ(ierr);
    if (domainerror) {
      ierr = SNESLineSearchSetSuccess(linesearch, PETSC_FALSE);CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
    if (linesearch->ops->vinorm) {
      gnorm = fnorm;
      ierr = (*linesearch->ops->vinorm)(snes, G, W, &gnorm);CHKERRQ(ierr);
    } else {
      ierr = VecNorm(G,NORM_2,&gnorm);CHKERRQ(ierr);
    }
    ierr = VecNorm(Y,NORM_2,&ynorm);CHKERRQ(ierr);
    if (PetscIsInfOrNanReal(gnorm)) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FP,"User provided compute function generated a Not-a-Number");

  }

  /* copy the solution over */
  ierr = VecCopy(W, X);CHKERRQ(ierr);
  ierr = VecCopy(G, F);CHKERRQ(ierr);
  ierr = VecNorm(X, NORM_2, &xnorm);CHKERRQ(ierr);
  ierr = SNESLineSearchSetLambda(linesearch, lambda);CHKERRQ(ierr);
  ierr = SNESLineSearchSetNorms(linesearch, xnorm, gnorm, ynorm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
} // lineSearch

// ----------------------------------------------------------------------
// Generic C interface for customized PETSc initial guess.
#undef __FUNCT__
#define __FUNCT__ "initialGuess"
PetscErrorCode
pylith::problems::SolverNonlinear::initialGuess(PetscSNES snes,
						PetscVec initialGuessVec,
						void *lsctx)
{ // initialGuess
  VecSet(initialGuessVec, 0.0);
} // initialGuess

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
