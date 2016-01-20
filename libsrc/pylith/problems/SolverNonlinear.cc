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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "SolverNonlinear.hh" // implementation of class methods

#include "Formulation.hh" // USES Formulation

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian

#include "pylith/utils/constdefs.h" // USES PYLITH_MAXSCALAR
#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/error.h" // USES PYLITH_CHECK_ERROR

#include <petscsnes.h> // USES PetscSNES

// KLUDGE, Fixes issue with PetscIsInfOrNanReal and include cmath
// instead of math.h.
#define isnan std::isnan // TEMPORARY
#define isinf std::isinf // TEMPORARY

// Customized line search based on PETSc code in
// src/snes/linesearch/impls/bt/linesearchbt.c.
#include <petsc/private/snesimpl.h>
#include <petsc/private/linesearchimpl.h>

struct PetscSNESLineSearch_BT {
  PetscReal        alpha; /* sufficient decrease parameter */
};

#define PYLITH_CUSTOM_LINESEARCH 1

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
  PYLITH_METHOD_BEGIN;

  Solver::deallocate();

  PetscErrorCode err = SNESDestroy(&_snes);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Initialize solver.
void
pylith::problems::SolverNonlinear::initialize(const topology::SolutionFields& fields,
					      const topology::Jacobian& jacobian,
					      Formulation* formulation)
{ // initialize
  PYLITH_METHOD_BEGIN;

  assert(formulation);

  _initializeLogger();
  Solver::initialize(fields, jacobian, formulation);

  PetscErrorCode err = 0;
  if (_snes) {
    err = SNESDestroy(&_snes); _snes = 0;
    PYLITH_CHECK_ERROR(err);
  } // if    
  err = SNESCreate(fields.mesh().comm(), &_snes); PYLITH_CHECK_ERROR(err);

  const topology::Field& residual = fields.get("residual");
  const PetscVec residualVec = residual.globalVector();
  err = SNESSetFunction(_snes, residualVec, reformResidual, (void*) formulation);
  PYLITH_CHECK_ERROR(err);

  err = SNESSetJacobian(_snes, jacobian.matrix(), _jacobianPC, reformJacobian, (void*) formulation);PYLITH_CHECK_ERROR(err);

  // Set default line search type to SNESSHELL and use our custom line search
  PetscSNESLineSearch ls;
  err = SNESGetLineSearch(_snes, &ls);PYLITH_CHECK_ERROR(err);
  err = SNESLineSearchSetType(ls, SNESSHELL);PYLITH_CHECK_ERROR(err);
  err = SNESLineSearchSetOrder(ls, SNES_LINESEARCH_ORDER_CUBIC);PYLITH_CHECK_ERROR(err);
  err = SNESLineSearchShellSetUserFunc(ls, lineSearch, (void*) formulation);PYLITH_CHECK_ERROR(err);

  // Get SNES options and allow the user to override the line search type
  err = SNESSetFromOptions(_snes);PYLITH_CHECK_ERROR(err);
  err = SNESSetComputeInitialGuess(_snes, initialGuess, (void*) formulation);PYLITH_CHECK_ERROR(err);

  if (formulation->splitFields()) {
    PetscKSP ksp = 0;
    PetscPC pc = 0;
    err = SNESGetKSP(_snes, &ksp); PYLITH_CHECK_ERROR(err);
    err = KSPGetPC(ksp, &pc); PYLITH_CHECK_ERROR(err);
    _setupFieldSplit(&pc, formulation, jacobian, fields);
  } // if

  _createNullSpace(fields);

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Solve the system.
void
pylith::problems::SolverNonlinear::solve(topology::Field* solution,
					 topology::Jacobian* jacobian,
					 const topology::Field& residual)
{ // solve
  PYLITH_METHOD_BEGIN;

  assert(solution);

  const int solveEvent = _logger->eventId("SoNl solve");
  const int scatterEvent = _logger->eventId("SoNl scatter");
  _logger->eventBegin(solveEvent);

  PetscErrorCode err = 0;
  const PetscVec solutionVec = solution->globalVector();

  err = SNESSolve(_snes, PETSC_NULL, solutionVec); PYLITH_CHECK_ERROR(err);
  
  _logger->eventEnd(solveEvent);
  _logger->eventBegin(scatterEvent);

  // Update section view of field.
  solution->scatterGlobalToLocal();

  _logger->eventEnd(scatterEvent);

  // Update rate fields to be consistent with current solution.
  _formulation->calcRateFields();

  PYLITH_METHOD_END;
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
  PYLITH_METHOD_BEGIN;

  assert(context);
  Formulation* formulation = (Formulation*) context;
  assert(formulation);

  // Make sure we have an admissible Lagrange multiplier (\lambda)
  VecLockPop(tmpSolutionVec); // :KLUDGE: TEMPORARY
  formulation->constrainSolnSpace(&tmpSolutionVec);
  VecLockPush(tmpSolutionVec); // :KLUDGE: TEMPORARY

  // Reform residual
  formulation->reformResidual(&tmpResidualVec, &tmpSolutionVec);

  PYLITH_METHOD_RETURN(0);
} // reformResidual

// ----------------------------------------------------------------------
// Generic C interface for reformJacobian for integration with
// PETSc SNES solvers.
PetscErrorCode
pylith::problems::SolverNonlinear::reformJacobian(PetscSNES snes,
						  PetscVec tmpSolutionVec,
						  PetscMat jacobianMat,
						  PetscMat preconditionerMat,
						  void* context)
{ // reformJacobian
  PYLITH_METHOD_BEGIN;

  assert(context);
  Formulation* formulation = (Formulation*) context;
  assert(formulation);

  formulation->reformJacobian(&tmpSolutionVec);

  PYLITH_METHOD_RETURN(0);
} // reformJacobian

// ----------------------------------------------------------------------
// Generic C interface for customized PETSc line search.
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

  PetscBool         changed_y,changed_w;
  PetscErrorCode    ierr;
  Vec               X,F,Y,W,G;
  SNES              snes;
  PetscReal         fnorm, xnorm, ynorm, gnorm;
  PetscReal         lambda,lambdatemp,lambdaprev,minlambda,maxstep,initslope,alpha,stol;
  PetscReal         t1,t2,a,b,d;
  PetscReal         f;
  PetscReal         g,gprev;
  PetscViewer       monitor;
  PetscInt          max_its,count;
  PetscSNESLineSearch_BT *bt;
  Mat               jac;
  PetscErrorCode    (*objective)(SNES,Vec,PetscReal*,void*);

  PetscFunctionBegin;
  ierr = SNESLineSearchGetVecs(linesearch, &X, &F, &Y, &W, &G);CHKERRQ(ierr);
  ierr = SNESLineSearchGetNorms(linesearch, &xnorm, &fnorm, &ynorm);CHKERRQ(ierr);
  ierr = SNESLineSearchGetLambda(linesearch, &lambda);CHKERRQ(ierr);
  ierr = SNESLineSearchGetSNES(linesearch, &snes);CHKERRQ(ierr);
  ierr = SNESLineSearchGetDefaultMonitor(linesearch, &monitor);CHKERRQ(ierr);
  ierr = SNESLineSearchGetTolerances(linesearch,&minlambda,&maxstep,NULL,NULL,NULL,&max_its);CHKERRQ(ierr);
  ierr = SNESGetTolerances(snes,NULL,NULL,&stol,NULL,NULL);CHKERRQ(ierr);
  ierr = SNESGetObjective(snes,&objective,NULL);CHKERRQ(ierr);
  bt   = (PetscSNESLineSearch_BT*)linesearch->data;

  alpha = bt->alpha;

  ierr = SNESGetJacobian(snes, &jac, NULL, NULL, NULL);CHKERRQ(ierr);

  if (!jac && !objective) SETERRQ(PetscObjectComm((PetscObject)linesearch), PETSC_ERR_USER, "SNESLineSearchBT requires a Jacobian matrix");

  /* precheck */
  ierr = SNESLineSearchPreCheck(linesearch,X,Y,&changed_y);CHKERRQ(ierr);
  ierr = SNESLineSearchSetReason(linesearch, SNES_LINESEARCH_SUCCEEDED);CHKERRQ(ierr);

  ierr = VecNormBegin(Y, NORM_2, &ynorm);CHKERRQ(ierr);
  ierr = VecNormBegin(X, NORM_2, &xnorm);CHKERRQ(ierr);
  ierr = VecNormEnd(Y, NORM_2, &ynorm);CHKERRQ(ierr);
  ierr = VecNormEnd(X, NORM_2, &xnorm);CHKERRQ(ierr);

  if (ynorm == 0.0) {
    if (monitor) {
      ierr = PetscViewerASCIIAddTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(monitor,"    Line search: Initial direction and size is 0\n");CHKERRQ(ierr);
      ierr = PetscViewerASCIISubtractTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
    }
    ierr = VecCopy(X,W);CHKERRQ(ierr);
    ierr = VecCopy(F,G);CHKERRQ(ierr);
    ierr = SNESLineSearchSetNorms(linesearch,xnorm,fnorm,ynorm);CHKERRQ(ierr);
    ierr = SNESLineSearchSetReason(linesearch, SNES_LINESEARCH_FAILED_REDUCT);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  if (ynorm > maxstep) {        /* Step too big, so scale back */
    if (monitor) {
      ierr = PetscViewerASCIIAddTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(monitor,"    Line search: Scaling step by %14.12e old ynorm %14.12e\n", (double)(maxstep/ynorm),(double)ynorm);CHKERRQ(ierr);
      ierr = PetscViewerASCIISubtractTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
    }
    ierr  = VecScale(Y,maxstep/(ynorm));CHKERRQ(ierr);
    ynorm = maxstep;
  }

  /* if the SNES has an objective set, use that instead of the function value */
  if (objective) {
    ierr = SNESComputeObjective(snes,X,&f);CHKERRQ(ierr);
  } else {
    f = fnorm*fnorm;
  }

  /* compute the initial slope */
  if (objective) {
    /* slope comes from the function (assumed to be the gradient of the objective */
    ierr = VecDotRealPart(Y,F,&initslope);CHKERRQ(ierr);
  } else {
    /* slope comes from the normal equations */
    ierr = MatMult(jac,Y,W);CHKERRQ(ierr);
    ierr = VecDotRealPart(F,W,&initslope);CHKERRQ(ierr);
    if (initslope > 0.0)  initslope = -initslope;
    if (initslope == 0.0) initslope = -1.0;
  }

  ierr = VecWAXPY(W,-lambda,Y,X);CHKERRQ(ierr);
  if (linesearch->ops->viproject) {
    ierr = (*linesearch->ops->viproject)(snes, W);CHKERRQ(ierr);
  }
  if (snes->nfuncs >= snes->max_funcs) {
    ierr         = PetscInfo(snes,"Exceeded maximum function evaluations, while checking full step length!\n");CHKERRQ(ierr);
    snes->reason = SNES_DIVERGED_FUNCTION_COUNT;
    ierr         = SNESLineSearchSetReason(linesearch, SNES_LINESEARCH_FAILED_FUNCTION);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  if (objective) {
    ierr = SNESComputeObjective(snes,W,&g);CHKERRQ(ierr);
  } else {
    ierr = (*linesearch->ops->snesfunc)(snes,W,G);CHKERRQ(ierr);
    if (linesearch->ops->vinorm) {
      gnorm = fnorm;
      ierr  = (*linesearch->ops->vinorm)(snes, G, W, &gnorm);CHKERRQ(ierr);
    } else {
      ierr = VecNorm(G,NORM_2,&gnorm);CHKERRQ(ierr);
    }
    g = PetscSqr(gnorm);
  }

  if (PetscIsInfOrNanReal(g)) {
    ierr = SNESLineSearchSetReason(linesearch, SNES_LINESEARCH_FAILED_NANORINF);CHKERRQ(ierr);
    ierr = PetscInfo(snes,"Aborted due to Nan or Inf in function evaluation\n");CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  if (!objective) {
    ierr = PetscInfo2(snes,"Initial fnorm %14.12e gnorm %14.12e\n", (double)fnorm, (double)gnorm);CHKERRQ(ierr);
  }
  if (.5*g <= .5*f + lambda*alpha*initslope) { /* Sufficient reduction or step tolerance convergence */
    if (monitor) {
      ierr = PetscViewerASCIIAddTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
      if (!objective) {
        ierr = PetscViewerASCIIPrintf(monitor,"    Line search: Using full step: fnorm %14.12e gnorm %14.12e\n", (double)fnorm, (double)gnorm);CHKERRQ(ierr);
      } else {
        ierr = PetscViewerASCIIPrintf(monitor,"    Line search: Using full step: obj0 %14.12e obj %14.12e\n", (double)f, (double)g);CHKERRQ(ierr);
      }
      ierr = PetscViewerASCIISubtractTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
    }
  } else {
    /* Since the full step didn't work and the step is tiny, quit */
    if (stol*xnorm > ynorm) {
      ierr = SNESLineSearchSetNorms(linesearch,xnorm,fnorm,ynorm);CHKERRQ(ierr);
      ierr = SNESLineSearchSetReason(linesearch, SNES_LINESEARCH_FAILED_REDUCT);CHKERRQ(ierr);
      if (monitor) {
        ierr = PetscViewerASCIIAddTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(monitor,"    Line search: Aborted due to ynorm < stol*xnorm (%14.12e < %14.12e) and inadequate full step.\n",(double)ynorm,(double)stol*xnorm);CHKERRQ(ierr);
        ierr = PetscViewerASCIISubtractTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }
    /* Fit points with quadratic */
    lambdatemp = -initslope/(g - f - 2.0*lambda*initslope);
    lambdaprev = lambda;
    gprev      = g;
    if (lambdatemp > .5*lambda)  lambdatemp = .5*lambda;
    if (lambdatemp <= .1*lambda) lambda = .1*lambda;
    else                         lambda = lambdatemp;

    ierr  = VecWAXPY(W,-lambda,Y,X);CHKERRQ(ierr);
    if (linesearch->ops->viproject) {
      ierr = (*linesearch->ops->viproject)(snes, W);CHKERRQ(ierr);
    }
    if (snes->nfuncs >= snes->max_funcs) {
      ierr         = PetscInfo1(snes,"Exceeded maximum function evaluations, while attempting quadratic backtracking! %D \n",snes->nfuncs);CHKERRQ(ierr);
      snes->reason = SNES_DIVERGED_FUNCTION_COUNT;
      ierr         = SNESLineSearchSetReason(linesearch, SNES_LINESEARCH_FAILED_FUNCTION);CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
    if (objective) {
      ierr = SNESComputeObjective(snes,W,&g);CHKERRQ(ierr);
    } else {
      ierr = (*linesearch->ops->snesfunc)(snes,W,G);CHKERRQ(ierr);
      if (linesearch->ops->vinorm) {
        gnorm = fnorm;
        ierr = (*linesearch->ops->vinorm)(snes, G, W, &gnorm);CHKERRQ(ierr);
      } else {
        ierr = VecNorm(G,NORM_2,&gnorm);CHKERRQ(ierr);
      }
      g = PetscSqr(gnorm);
    }
    if (PetscIsInfOrNanReal(g)) {
      ierr = SNESLineSearchSetReason(linesearch, SNES_LINESEARCH_FAILED_NANORINF);CHKERRQ(ierr);
      ierr = PetscInfo(snes,"Aborted due to Nan or Inf in function evaluation\n");CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
    if (monitor) {
      ierr = PetscViewerASCIIAddTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
      if (!objective) {
        ierr = PetscViewerASCIIPrintf(monitor,"    Line search: gnorm after quadratic fit %14.12e\n",(double)gnorm);CHKERRQ(ierr);
      } else {
        ierr = PetscViewerASCIIPrintf(monitor,"    Line search: obj after quadratic fit %14.12e\n",(double)g);CHKERRQ(ierr);
      }
      ierr = PetscViewerASCIISubtractTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
    }
    if (.5*g < .5*f + lambda*alpha*initslope) { /* sufficient reduction */
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
            if (!objective) {
              ierr = PetscViewerASCIIPrintf(monitor,
                                            "    Line search: fnorm=%18.16e, gnorm=%18.16e, ynorm=%18.16e, minlambda=%18.16e, lambda=%18.16e, initial slope=%18.16e\n",
                                            (double)fnorm, (double)gnorm, (double)ynorm, (double)minlambda, (double)lambda, (double)initslope);CHKERRQ(ierr);
            } else {
              ierr = PetscViewerASCIIPrintf(monitor,
                                            "    Line search: obj(0)=%18.16e, obj=%18.16e, ynorm=%18.16e, minlambda=%18.16e, lambda=%18.16e, initial slope=%18.16e\n",
                                            (double)f, (double)g, (double)ynorm, (double)minlambda, (double)lambda, (double)initslope);CHKERRQ(ierr);
            }
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
          ierr = SNESLineSearchSetReason(linesearch, SNES_LINESEARCH_FAILED_REDUCT);CHKERRQ(ierr);
          PetscFunctionReturn(0);
#endif
        }
        if (linesearch->order == SNES_LINESEARCH_ORDER_CUBIC) {
          t1 = .5*(g - f) - lambda*initslope;
          t2 = .5*(gprev  - f) - lambdaprev*initslope;
          a  = (t1/(lambda*lambda) - t2/(lambdaprev*lambdaprev))/(lambda-lambdaprev);
          b  = (-lambdaprev*t1/(lambda*lambda) + lambda*t2/(lambdaprev*lambdaprev))/(lambda-lambdaprev);
          d  = b*b - 3*a*initslope;
          if (d < 0.0) d = 0.0;
          if (a == 0.0) lambdatemp = -initslope/(2.0*b);
          else lambdatemp = (-b + PetscSqrtReal(d))/(3.0*a);

        } else if (linesearch->order == SNES_LINESEARCH_ORDER_QUADRATIC) {
          lambdatemp = -initslope/(g - f - 2.0*initslope);
        } else SETERRQ(PetscObjectComm((PetscObject)linesearch), PETSC_ERR_SUP, "unsupported line search order for type bt");
        lambdaprev = lambda;
        gprev      = g;
        if (lambdatemp > .5*lambda)  lambdatemp = .5*lambda;
        if (lambdatemp <= .1*lambda) lambda     = .1*lambda;
        else                         lambda     = lambdatemp;
        ierr = VecWAXPY(W,-lambda,Y,X);CHKERRQ(ierr);
        if (linesearch->ops->viproject) {
          ierr = (*linesearch->ops->viproject)(snes,W);CHKERRQ(ierr);
        }
        if (snes->nfuncs >= snes->max_funcs) {
          ierr = PetscInfo1(snes,"Exceeded maximum function evaluations, while looking for good step length! %D \n",count);CHKERRQ(ierr);
          if (!objective) {
            ierr = PetscInfo5(snes,"fnorm=%18.16e, gnorm=%18.16e, ynorm=%18.16e, lambda=%18.16e, initial slope=%18.16e\n",
                              (double)fnorm,(double)gnorm,(double)ynorm,(double)lambda,(double)initslope);CHKERRQ(ierr);
          }
          ierr         = SNESLineSearchSetReason(linesearch, SNES_LINESEARCH_FAILED_FUNCTION);CHKERRQ(ierr);
          snes->reason = SNES_DIVERGED_FUNCTION_COUNT;
          PetscFunctionReturn(0);
        }
        if (objective) {
          ierr = SNESComputeObjective(snes,W,&g);CHKERRQ(ierr);
        } else {
          ierr = (*linesearch->ops->snesfunc)(snes,W,G);CHKERRQ(ierr);
          if (linesearch->ops->vinorm) {
            gnorm = fnorm;
            ierr  = (*linesearch->ops->vinorm)(snes, G, W, &gnorm);CHKERRQ(ierr);
          } else {
            ierr = VecNorm(G,NORM_2,&gnorm);CHKERRQ(ierr);
          }
          g = PetscSqr(gnorm);
        }
        if (PetscIsInfOrNanReal(gnorm)) {
          ierr = SNESLineSearchSetReason(linesearch, SNES_LINESEARCH_FAILED_NANORINF);CHKERRQ(ierr);
          ierr = PetscInfo(snes,"Aborted due to Nan or Inf in function evaluation\n");CHKERRQ(ierr);
          PetscFunctionReturn(0);
        }
        if (.5*g < .5*f + lambda*alpha*initslope) { /* is reduction enough? */
          if (monitor) {
            ierr = PetscViewerASCIIAddTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
            if (!objective) {
              if (linesearch->order == SNES_LINESEARCH_ORDER_CUBIC) {
                ierr = PetscViewerASCIIPrintf(monitor,"    Line search: Cubically determined step, current gnorm %14.12e lambda=%18.16e\n",(double)gnorm,(double)lambda);CHKERRQ(ierr);
              } else {
                ierr = PetscViewerASCIIPrintf(monitor,"    Line search: Quadratically determined step, current gnorm %14.12e lambda=%18.16e\n",(double)gnorm,(double)lambda);CHKERRQ(ierr);
              }
              ierr = PetscViewerASCIISubtractTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
            } else {
              if (linesearch->order == SNES_LINESEARCH_ORDER_CUBIC) {
                ierr = PetscViewerASCIIPrintf(monitor,"    Line search: Cubically determined step, obj %14.12e lambda=%18.16e\n",(double)g,(double)lambda);CHKERRQ(ierr);
              } else {
                ierr = PetscViewerASCIIPrintf(monitor,"    Line search: Quadratically determined step, obj %14.12e lambda=%18.16e\n",(double)g,(double)lambda);CHKERRQ(ierr);
              }
              ierr = PetscViewerASCIISubtractTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
            }
          }
          break;
        } else if (monitor) {
          ierr = PetscViewerASCIIAddTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
          if (!objective) {
            if (linesearch->order == SNES_LINESEARCH_ORDER_CUBIC) {
              ierr = PetscViewerASCIIPrintf(monitor,"    Line search: Cubic step no good, shrinking lambda, current gnorm %12.12e lambda=%18.16e\n",(double)gnorm,(double)lambda);CHKERRQ(ierr);
            } else {
              ierr = PetscViewerASCIIPrintf(monitor,"    Line search: Quadratic step no good, shrinking lambda, current gnorm %12.12e lambda=%18.16e\n",(double)gnorm,(double)lambda);CHKERRQ(ierr);
            }
            ierr = PetscViewerASCIISubtractTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
          } else {
            if (linesearch->order == SNES_LINESEARCH_ORDER_CUBIC) {
              ierr = PetscViewerASCIIPrintf(monitor,"    Line search: Cubic step no good, shrinking lambda, obj %12.12e lambda=%18.16e\n",(double)g,(double)lambda);CHKERRQ(ierr);
            } else {
              ierr = PetscViewerASCIIPrintf(monitor,"    Line search: Quadratic step no good, shrinking lambda, obj %12.12e lambda=%18.16e\n",(double)g,(double)lambda);CHKERRQ(ierr);
            }
            ierr = PetscViewerASCIISubtractTab(monitor,((PetscObject)linesearch)->tablevel);CHKERRQ(ierr);
          }
        }
      }
    }
  }

  /* postcheck */
  ierr = SNESLineSearchPostCheck(linesearch,X,Y,W,&changed_y,&changed_w);CHKERRQ(ierr);
  if (changed_y) {
    ierr = VecWAXPY(W,-lambda,Y,X);CHKERRQ(ierr);
    if (linesearch->ops->viproject) {
      ierr = (*linesearch->ops->viproject)(snes, W);CHKERRQ(ierr);
    }
  }
#if defined(PYLITH_CUSTOM_LINESEARCH)
  if (true) {
#else // ORIGINAL
  if (changed_y || changed_w || objective) { /* recompute the function norm if the step has changed or the objective isn't the norm */
#endif
    ierr = (*linesearch->ops->snesfunc)(snes,W,G);CHKERRQ(ierr);
    if (linesearch->ops->vinorm) {
      gnorm = fnorm;
      ierr  = (*linesearch->ops->vinorm)(snes, G, W, &gnorm);CHKERRQ(ierr);
    } else {
      ierr = VecNorm(G,NORM_2,&gnorm);CHKERRQ(ierr);
    }
    ierr = VecNorm(Y,NORM_2,&ynorm);CHKERRQ(ierr);
    if (PetscIsInfOrNanReal(gnorm)) {
      ierr = SNESLineSearchSetReason(linesearch,SNES_LINESEARCH_FAILED_NANORINF);CHKERRQ(ierr);
      ierr = PetscInfo(snes,"Aborted due to Nan or Inf in function evaluation\n");CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
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
PetscErrorCode
pylith::problems::SolverNonlinear::initialGuess(PetscSNES snes,
						PetscVec initialGuessVec,
						void *lsctx)
{ // initialGuess
  PYLITH_METHOD_BEGIN;

  PetscErrorCode err = VecSet(initialGuessVec, 0.0);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_RETURN(0);
} // initialGuess

// ----------------------------------------------------------------------
// Initialize logger.
void
pylith::problems::SolverNonlinear::_initializeLogger(void)
{ // initializeLogger
  PYLITH_METHOD_BEGIN;

  delete _logger; _logger = new utils::EventLogger;assert(_logger);
  _logger->className("SolverNonlinear");
  _logger->initialize();
  _logger->registerEvent("SoNl setup");
  _logger->registerEvent("SoNl solve");
  _logger->registerEvent("SoNl scatter");

  PYLITH_METHOD_END;
} // initializeLogger


// End of file
