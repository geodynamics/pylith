// -*- C++ -*-
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
//
//  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
//
//  Permission is hereby granted, free of charge, to any person obtaining
//  a copy of this software and associated documentation files (the
//  "Software"), to deal in the Software without restriction, including
//  without limitation the rights to use, copy, modify, merge, publish,
//  distribute, sublicense, and/or sell copies of the Software, and to
//  permit persons to whom the Software is furnished to do so, subject to
//  the following conditions:
//
//  The above copyright notice and this permission notice shall be
//  included in all copies or substantial portions of the Software.
//
//  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
//  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
//  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
//  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
//  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
//  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 

#include <portinfo>
#include "journal/debug.h"

#include <Python.h>

#include "libpylith3d.h"
#include "exceptionhandler.h"
#include "pylith3d_externs.h"
#include <stdio.h>
#include <string.h>


// Autocompute gravitational prestresses.

char pypylith3d_autoprestr__doc__[] = "";
char pypylith3d_autoprestr__name__[] = "autoprestr";

PyObject * pypylith3d_autoprestr(PyObject *, PyObject *args)
{
  PyObject* pyA;                             // Sparse matrix array
  PyObject* pyRhs;
  PyObject* pySol;
  PyObject* py_bextern;               // Force vectors
  PyObject* py_btraction;
  PyObject* py_bgravity;
  PyObject* py_bconcForce;
  PyObject* py_bintern;
  PyObject* py_bresid;
  PyObject* py_bwink;
  PyObject* py_bwinkx;
  PyObject* py_dispVec;
  PyObject* py_dprev;
  PyObject* py_nforce;
  PyObject* py_grav;
  PyObject* py_x;                             // Global arrays
  PyObject* py_d;
  PyObject* py_deld;
  PyObject* py_dcur;
  PyObject* py_id;
  PyObject* py_iwink;
  PyObject* py_wink;
  PyObject* py_nsysdat;
  PyObject* py_iddmat;
  PyObject* py_ibond;                 // Boundary condition arrays
  PyObject* py_bond;
  PyObject* py_dx;                    // Slippery node arrays
  PyObject* py_deldx;
  PyObject* py_dxcur;
  PyObject* py_diforc;
  PyObject* py_idx;
  PyObject* py_iwinkx;
  PyObject* py_winkx;
  PyObject* py_idslp;
  PyObject* py_ipslp;
  PyObject* py_idhist;
  PyObject* py_fault;                 // Split node arrays
  PyObject* py_nfault;
  PyObject* py_dfault;
  PyObject* py_tfault;
  PyObject* py_s;                     // Local stiffness matrix arrays
  PyObject* py_stemp;
  PyObject* py_state;                 // Element arrays
  PyObject* py_dstate;
  PyObject* py_state0;
  PyObject* py_dmat;
  PyObject* py_iens;
  PyObject* py_lm;
  PyObject* py_lmx;
  PyObject* py_lmf;
  PyObject* py_ivfamily;
  PyObject* py_npar;
  PyObject* py_ielindx;
  PyObject* py_tractionverts;         // Traction BC arrays
  PyObject* py_tractionvals;
  PyObject* py_gauss2d;
  PyObject* py_sh2d;
  PyObject* py_infetype2d;
  PyObject* py_prop; // Material property arrays
  PyObject* py_infmatmod;
  PyObject* py_gauss;                 // Element type arrays
  PyObject* py_sh;
  PyObject* py_shj;
  PyObject* py_infetype;
  PyObject* py_histry;                // Time information
  PyObject* py_rtimdat;
  PyObject* py_ntimdat;
  PyObject* py_nvisdat;
  PyObject* py_maxstp;
  PyObject* py_delt;
  PyObject* py_alfa;
  PyObject* py_maxit;
  PyObject* py_ntdinit;
  PyObject* py_lgdef;
  PyObject* py_utol;
  PyObject* py_ftol;
  PyObject* py_etol;
  PyObject* py_itmax;
  PyObject* py_rgiter;       // Iterative solution information
  PyObject* py_skew;                  // Skew rotation information
  PyObject* py_ncodat;       // Input/output information
  PyObject* py_nunits;
  PyObject* py_nprint;
  PyObject* py_istatout;
  PyObject* py_nstatout;
  char* ofile;                      // Output file names
  char* pfile;
  char* ucdroot;
  int autoprestrStage, iterateEvent;          // PETSc logging

  int ok = PyArg_ParseTuple(args, "OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOsssii:autoprestr",
                            &pyA,                             // Sparse matrix arrays
                            &pyRhs,
                            &pySol,
                            &py_bextern,               // Force vectors
                            &py_btraction,
                            &py_bgravity,
                            &py_bconcForce,
                            &py_bintern,
                            &py_bresid,
                            &py_bwink,
                            &py_bwinkx,
                            &py_dispVec,
                            &py_dprev,
                            &py_nforce,
                            &py_grav,
                            &py_x,                             // Global arrays
                            &py_d,
                            &py_deld,
                            &py_dcur,
                            &py_id,
                            &py_iwink,
                            &py_wink,
                            &py_nsysdat,
                            &py_iddmat,
                            &py_ibond,                 // Boundary condition arrays
                            &py_bond,
                            &py_dx,                    // Slippery node arrays
                            &py_deldx,
                            &py_dxcur,
                            &py_diforc,
                            &py_idx,
                            &py_iwinkx,
                            &py_winkx,
                            &py_idslp,
                            &py_ipslp,
                            &py_idhist,
                            &py_fault,                 // Split node arrays
                            &py_nfault,
                            &py_dfault,
                            &py_tfault,
                            &py_s,                     // Local stiffness matrix arrays
                            &py_stemp,
                            &py_state,                 // Element arrays
                            &py_dstate,
                            &py_state0,
                            &py_dmat,
                            &py_iens,
                            &py_lm,
                            &py_lmx,
                            &py_lmf,
                            &py_ivfamily,
                            &py_npar,
                            &py_ielindx,
                            &py_tractionverts,         // Traction BC arrays
                            &py_tractionvals,
                            &py_gauss2d,
                            &py_sh2d,
                            &py_infetype2d,
                            &py_prop, // Material property arrays
                            &py_infmatmod,
                            &py_gauss,                 // Element type arrays
                            &py_sh,
                            &py_shj,
                            &py_infetype,
                            &py_histry,                // Time information
                            &py_rtimdat,
                            &py_ntimdat,
                            &py_nvisdat,
                            &py_maxstp,
                            &py_delt,
                            &py_alfa,
                            &py_maxit,
                            &py_ntdinit,
                            &py_lgdef,
                            &py_utol,
                            &py_ftol,
                            &py_etol,
                            &py_itmax,
                            &py_rgiter,       // Iterative solution information
                            &py_skew,                  // Skew rotation information
                            &py_ncodat,       // Input/output information
                            &py_nunits,
                            &py_nprint,
                            &py_istatout,
                            &py_nstatout,
                            &ofile,                  // Output file names
                            &pfile,
                            &ucdroot,
                            &autoprestrStage,                  // PETSc logging
                            &iterateEvent);


  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  Mat      A = (Mat) PyCObject_AsVoidPtr(pyA);
  Vec      rhs = (Vec) PyCObject_AsVoidPtr(pyRhs);
  Vec      sol = (Vec) PyCObject_AsVoidPtr(pySol);
  double*  bextern = (double*) PyCObject_AsVoidPtr(py_bextern);
  double*  btraction = (double*) PyCObject_AsVoidPtr(py_btraction);
  double*  bgravity = (double*) PyCObject_AsVoidPtr(py_bgravity);
  double*  bconcForce = (double*) PyCObject_AsVoidPtr(py_bconcForce);
  double*  bintern = (double*) PyCObject_AsVoidPtr(py_bintern);
  double*  bresid = (double*) PyCObject_AsVoidPtr(py_bresid);
  double*  bwink = (double*) PyCObject_AsVoidPtr(py_bwink);
  double*  bwinkx = (double*) PyCObject_AsVoidPtr(py_bwinkx);
  double*  dispVec = (double*) PyCObject_AsVoidPtr(py_dispVec);
  double*  dprev = (double*) PyCObject_AsVoidPtr(py_dprev);
  int*  nforce = (int*) PyCObject_AsVoidPtr(py_nforce);
  double*  grav = (double*) PyCObject_AsVoidPtr(py_grav);
  double*  x = (double*) PyCObject_AsVoidPtr(py_x);
  double*  d = (double*) PyCObject_AsVoidPtr(py_d);
  double*  deld = (double*) PyCObject_AsVoidPtr(py_deld);
  double*  dcur = (double*) PyCObject_AsVoidPtr(py_dcur);
  int*  id = (int*) PyCObject_AsVoidPtr(py_id);
  int*  iwink = (int*) PyCObject_AsVoidPtr(py_iwink);
  double*  wink = (double*) PyCObject_AsVoidPtr(py_wink);
  int*  nsysdat = (int*) PyCObject_AsVoidPtr(py_nsysdat);
  int*  iddmat = (int*) PyCObject_AsVoidPtr(py_iddmat);
  int*  ibond = (int*) PyCObject_AsVoidPtr(py_ibond);
  double*  bond = (double*) PyCObject_AsVoidPtr(py_bond);
  double*  dx = (double*) PyCObject_AsVoidPtr(py_dx);
  double*  deldx = (double*) PyCObject_AsVoidPtr(py_deldx);
  double*  dxcur = (double*) PyCObject_AsVoidPtr(py_dxcur);
  double*  diforc = (double*) PyCObject_AsVoidPtr(py_diforc);
  int*  idx = (int*) PyCObject_AsVoidPtr(py_idx);
  int*  iwinkx = (int*) PyCObject_AsVoidPtr(py_iwinkx);
  double*  winkx = (double*) PyCObject_AsVoidPtr(py_winkx);
  int*  idslp = (int*) PyCObject_AsVoidPtr(py_idslp);
  int*  ipslp = (int*) PyCObject_AsVoidPtr(py_ipslp);
  int*  idhist = (int*) PyCObject_AsVoidPtr(py_idhist);
  double*  fault = (double*) PyCObject_AsVoidPtr(py_fault);
  int*  nfault = (int*) PyCObject_AsVoidPtr(py_nfault);
  double*  dfault = (double*) PyCObject_AsVoidPtr(py_dfault);
  double*  tfault = (double*) PyCObject_AsVoidPtr(py_tfault);
  double*  s = (double*) PyCObject_AsVoidPtr(py_s);
  double*  stemp = (double*) PyCObject_AsVoidPtr(py_stemp);
  double*  state = (double*) PyCObject_AsVoidPtr(py_state);
  double*  dstate = (double*) PyCObject_AsVoidPtr(py_dstate);
  double*  state0 = (double*) PyCObject_AsVoidPtr(py_state0);
  double*  dmat = (double*) PyCObject_AsVoidPtr(py_dmat);
  int*  iens = (int*) PyCObject_AsVoidPtr(py_iens);
  int*  lm = (int*) PyCObject_AsVoidPtr(py_lm);
  int*  lmx = (int*) PyCObject_AsVoidPtr(py_lmx);
  int*  lmf = (int*) PyCObject_AsVoidPtr(py_lmf);
  int*  ivfamily = (int*) PyCObject_AsVoidPtr(py_ivfamily);
  int*  npar = (int*) PyCObject_AsVoidPtr(py_npar);
  int*  ielindx = (int*) PyCObject_AsVoidPtr(py_ielindx);
  int*  tractionverts = (int*) PyCObject_AsVoidPtr(py_tractionverts);
  double*  tractionvals = (double*) PyCObject_AsVoidPtr(py_tractionvals);
  double*  gauss2d = (double*) PyCObject_AsVoidPtr(py_gauss2d);
  double*  sh2d = (double*) PyCObject_AsVoidPtr(py_sh2d);
  int*  infetype2d = (int*) PyCObject_AsVoidPtr(py_infetype2d);
  double*  prop = (double*) PyCObject_AsVoidPtr(py_prop);
  int*  infmatmod = (int*) PyCObject_AsVoidPtr(py_infmatmod);
  double*  gauss = (double*) PyCObject_AsVoidPtr(py_gauss);
  double*  sh = (double*) PyCObject_AsVoidPtr(py_sh);
  double*  shj = (double*) PyCObject_AsVoidPtr(py_shj);
  int*  infetype = (int*) PyCObject_AsVoidPtr(py_infetype);
  double*  histry = (double*) PyCObject_AsVoidPtr(py_histry);
  double*  rtimdat = (double*) PyCObject_AsVoidPtr(py_rtimdat);
  int*  ntimdat = (int*) PyCObject_AsVoidPtr(py_ntimdat);
  int*  nvisdat = (int*) PyCObject_AsVoidPtr(py_nvisdat);
  int*  maxstp = (int*) PyCObject_AsVoidPtr(py_maxstp);
  double*  delt = (double*) PyCObject_AsVoidPtr(py_delt);
  double*  alfa = (double*) PyCObject_AsVoidPtr(py_alfa);
  int*  maxit = (int*) PyCObject_AsVoidPtr(py_maxit);
  int*  ntdinit = (int*) PyCObject_AsVoidPtr(py_ntdinit);
  int*  lgdef = (int*) PyCObject_AsVoidPtr(py_lgdef);
  double*  utol = (double*) PyCObject_AsVoidPtr(py_utol);
  double*  ftol = (double*) PyCObject_AsVoidPtr(py_ftol);
  double*  etol = (double*) PyCObject_AsVoidPtr(py_etol);
  int*  itmax = (int*) PyCObject_AsVoidPtr(py_itmax);
  double*  rgiter = (double*) PyCObject_AsVoidPtr(py_rgiter);
  double*  skew = (double*) PyCObject_AsVoidPtr(py_skew);
  int*  ncodat = (int*) PyCObject_AsVoidPtr(py_ncodat);
  int*  nunits = (int*) PyCObject_AsVoidPtr(py_nunits);
  int*  nprint = (int*) PyCObject_AsVoidPtr(py_nprint);
  int*  istatout = (int*) PyCObject_AsVoidPtr(py_istatout);
  int*  nstatout = (int*) PyCObject_AsVoidPtr(py_nstatout);

  autoprestr_f(&A,                                // Sparse matrix arrays
               &rhs,
               &sol,
               bextern,                  // Force vectors
               btraction,
               bgravity,
               bconcForce,
               bintern,
               bresid,
               bwink,
               bwinkx,
               dispVec,
               dprev,
               nforce,
               grav,
               x,                                 // Global arrays
               d,
               deld,
               dcur,
               id,
               iwink,
               wink,
               nsysdat,
               iddmat,
               ibond,                    // Boundary condition arrays
               bond,
               dx,                       // Slippery node arrays
               deldx,
               dxcur,
               diforc,
               idx,
               iwinkx,
               winkx,
               idslp,
               ipslp,
               idhist,
               fault,                    // Split node arrays
               nfault,
               dfault,
               tfault,
               s,                        // Local stiffness matrix arrays
               stemp,
               state,                    // Element arrays
               dstate,
               state0,
               dmat,
               iens,
               lm,
               lmx,
               lmf,
               ivfamily,
               npar,
               ielindx,
               tractionverts,            // Traction BC arrays
               tractionvals,
               gauss2d,
               sh2d,
               infetype2d,
               prop,    // Material property arrays
               infmatmod,
               gauss,                    // Element type arrays
               sh,
               shj,
               infetype,
               histry,                   // Time information
               rtimdat,
               ntimdat,
               nvisdat,
               maxstp,
               delt,
               alfa,
               maxit,
               ntdinit,
               lgdef,
               utol,
               ftol,
               etol,
               itmax,
               rgiter,          // Iterative solution information
               skew,                     // Skew rotation information
               ncodat,          // Input/output information
               nunits,
               nprint,
               istatout,
               nstatout,
               ofile,                   // Output file names
               pfile,
               ucdroot,
               &autoprestrStage,
               &iterateEvent,
               &errorcode,                        // Error codes
               errorstring,
               strlen(ofile),           // String lengths
               strlen(pfile),
               strlen(ucdroot),
               sizeof(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "ofile:" << ofile
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}



// Perform the elastic solution.

char pypylith3d_elastc__doc__[] = "";
char pypylith3d_elastc__name__[] = "elastc";

PyObject * pypylith3d_elastc(PyObject *, PyObject *args)
{
  PyObject* pyA;                              // Sparse matrix arrays
  PyObject* pyRhs;
  PyObject* pySol;
  PyObject* py_bextern;               // Force vectors
  PyObject* py_btraction;
  PyObject* py_bgravity;
  PyObject* py_bconcForce;
  PyObject* py_bintern;
  PyObject* py_bresid;
  PyObject* py_bwink;
  PyObject* py_bwinkx;
  PyObject* py_dispVec;
  PyObject* py_dprev;
  PyObject* py_nforce;
  PyObject* py_grav;
  PyObject* py_x;                             // Global arrays
  PyObject* py_d;
  PyObject* py_deld;
  PyObject* py_dcur;
  PyObject* py_id;
  PyObject* py_iwink;
  PyObject* py_wink;
  PyObject* py_nsysdat;
  PyObject* py_iddmat;
  PyObject* py_ibond;                 // Boundary condition arrays
  PyObject* py_bond;
  PyObject* py_dx;                    // Slippery node arrays
  PyObject* py_deldx;
  PyObject* py_dxcur;
  PyObject* py_diforc;
  PyObject* py_idx;
  PyObject* py_iwinkx;
  PyObject* py_winkx;
  PyObject* py_idslp;
  PyObject* py_ipslp;
  PyObject* py_idhist;
  PyObject* py_fault;                 // Split node arrays
  PyObject* py_nfault;
  PyObject* py_dfault;
  PyObject* py_tfault;
  PyObject* py_s;                     // Local stiffness matrix arrays
  PyObject* py_stemp;
  PyObject* py_state;                 // Element arrays
  PyObject* py_dstate;
  PyObject* py_state0;
  PyObject* py_dmat;
  PyObject* py_ien;
  PyObject* py_lm;
  PyObject* py_lmx;
  PyObject* py_lmf;
  PyObject* py_ivfamily;
  PyObject* py_npar;
  PyObject* py_ielindx;
  PyObject* py_tractionverts;         // Traction BC arrays
  PyObject* py_tractionvals;
  PyObject* py_gauss2d;
  PyObject* py_sh2d;
  PyObject* py_infetype2d;
  PyObject* py_prop; // Material property arrays
  PyObject* py_infmatmod;
  PyObject* py_gauss;                 // Element type arrays
  PyObject* py_sh;
  PyObject* py_shj;
  PyObject* py_infetype;
  PyObject* py_histry;                // Time information
  PyObject* py_rtimdat;
  PyObject* py_ntimdat;
  PyObject* py_nvisdat;
  PyObject* py_maxstp;
  PyObject* py_delt;
  PyObject* py_alfa;
  PyObject* py_maxit;
  PyObject* py_ntdinit;
  PyObject* py_lgdef;
  PyObject* py_utol;
  PyObject* py_ftol;
  PyObject* py_etol;
  PyObject* py_itmax;
  PyObject* py_rgiter;       // Iterative solution information
  PyObject* py_skew;                  // Skew rotation information
  PyObject* py_ncodat;       // Input/output information
  PyObject* py_nunits;
  PyObject* py_nprint;
  PyObject* py_istatout;
  PyObject* py_nstatout;
  char* ofile;                      // Output file names
  char* pfile;
  char* ucdroot;
  int elasticStage, iterateEvent;             // PETSc logging

  int ok = PyArg_ParseTuple(args, "OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOsssii:elastc",
			    &pyA,                              // Sparse matrix arrays
			    &pyRhs,
			    &pySol,
			    &py_bextern,               // Force vectors
			    &py_btraction,
			    &py_bgravity,
			    &py_bconcForce,
			    &py_bintern,
			    &py_bresid,
			    &py_bwink,
			    &py_bwinkx,
			    &py_dispVec,
			    &py_dprev,
			    &py_nforce,
			    &py_grav,
			    &py_x,                             // Global arrays
			    &py_d,
			    &py_deld,
			    &py_dcur,
			    &py_id,
			    &py_iwink,
			    &py_wink,
			    &py_nsysdat,
			    &py_iddmat,
			    &py_ibond,                 // Boundary condition arrays
			    &py_bond,
			    &py_dx,                    // Slippery node arrays
			    &py_deldx,
			    &py_dxcur,
			    &py_diforc,
			    &py_idx,
			    &py_iwinkx,
			    &py_winkx,
			    &py_idslp,
			    &py_ipslp,
			    &py_idhist,
			    &py_fault,                 // Split node arrays
			    &py_nfault,
			    &py_dfault,
			    &py_tfault,
			    &py_s,                     // Local stiffness matrix arrays
			    &py_stemp,
			    &py_state,                 // Element arrays
			    &py_dstate,
			    &py_state0,
			    &py_dmat,
			    &py_ien,
			    &py_lm,
			    &py_lmx,
			    &py_lmf,
			    &py_ivfamily,
			    &py_npar,
                            &py_ielindx,
                            &py_tractionverts,         // Traction BC arrays
                            &py_tractionvals,
                            &py_gauss2d,
                            &py_sh2d,
                            &py_infetype2d,
			    &py_prop, // Material property arrays
			    &py_infmatmod,
			    &py_gauss,                 // Element type arrays
			    &py_sh,
			    &py_shj,
			    &py_infetype,
			    &py_histry,                // Time information
			    &py_rtimdat,
			    &py_ntimdat,
			    &py_nvisdat,
			    &py_maxstp,
			    &py_delt,
			    &py_alfa,
			    &py_maxit,
			    &py_ntdinit,
			    &py_lgdef,
			    &py_utol,
			    &py_ftol,
			    &py_etol,
			    &py_itmax,
			    &py_rgiter,       // Iterative solution information
			    &py_skew,                  // Skew rotation information
			    &py_ncodat,       // Input/output information
			    &py_nunits,
			    &py_nprint,
			    &py_istatout,
			    &py_nstatout,
			    &ofile,                  // Output file names
			    &pfile,
			    &ucdroot,
                            &elasticStage,                    // PETSc logging
                            &iterateEvent);


  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  Mat      A = (Mat) PyCObject_AsVoidPtr(pyA);
  Vec      rhs = (Vec) PyCObject_AsVoidPtr(pyRhs);
  Vec      sol = (Vec) PyCObject_AsVoidPtr(pySol);
  double*  bextern = (double*) PyCObject_AsVoidPtr(py_bextern);
  double*  btraction = (double*) PyCObject_AsVoidPtr(py_btraction);
  double*  bgravity = (double*) PyCObject_AsVoidPtr(py_bgravity);
  double*  bconcForce = (double*) PyCObject_AsVoidPtr(py_bconcForce);
  double*  bintern = (double*) PyCObject_AsVoidPtr(py_bintern);
  double*  bresid = (double*) PyCObject_AsVoidPtr(py_bresid);
  double*  bwink = (double*) PyCObject_AsVoidPtr(py_bwink);
  double*  bwinkx = (double*) PyCObject_AsVoidPtr(py_bwinkx);
  double*  dispVec = (double*) PyCObject_AsVoidPtr(py_dispVec);
  double*  dprev = (double*) PyCObject_AsVoidPtr(py_dprev);
  int*  nforce = (int*) PyCObject_AsVoidPtr(py_nforce);
  double*  grav = (double*) PyCObject_AsVoidPtr(py_grav);
  double*  x = (double*) PyCObject_AsVoidPtr(py_x);
  double*  d = (double*) PyCObject_AsVoidPtr(py_d);
  double*  deld = (double*) PyCObject_AsVoidPtr(py_deld);
  double*  dcur = (double*) PyCObject_AsVoidPtr(py_dcur);
  int*  id = (int*) PyCObject_AsVoidPtr(py_id);
  int*  iwink = (int*) PyCObject_AsVoidPtr(py_iwink);
  double*  wink = (double*) PyCObject_AsVoidPtr(py_wink);
  int*  nsysdat = (int*) PyCObject_AsVoidPtr(py_nsysdat);
  int*  iddmat = (int*) PyCObject_AsVoidPtr(py_iddmat);
  int*  ibond = (int*) PyCObject_AsVoidPtr(py_ibond);
  double*  bond = (double*) PyCObject_AsVoidPtr(py_bond);
  double*  dx = (double*) PyCObject_AsVoidPtr(py_dx);
  double*  deldx = (double*) PyCObject_AsVoidPtr(py_deldx);
  double*  dxcur = (double*) PyCObject_AsVoidPtr(py_dxcur);
  double*  diforc = (double*) PyCObject_AsVoidPtr(py_diforc);
  int*  idx = (int*) PyCObject_AsVoidPtr(py_idx);
  int*  iwinkx = (int*) PyCObject_AsVoidPtr(py_iwinkx);
  double*  winkx = (double*) PyCObject_AsVoidPtr(py_winkx);
  int*  idslp = (int*) PyCObject_AsVoidPtr(py_idslp);
  int*  ipslp = (int*) PyCObject_AsVoidPtr(py_ipslp);
  int*  idhist = (int*) PyCObject_AsVoidPtr(py_idhist);
  double*  fault = (double*) PyCObject_AsVoidPtr(py_fault);
  int*  nfault = (int*) PyCObject_AsVoidPtr(py_nfault);
  double*  dfault = (double*) PyCObject_AsVoidPtr(py_dfault);
  double*  tfault = (double*) PyCObject_AsVoidPtr(py_tfault);
  double*  s = (double*) PyCObject_AsVoidPtr(py_s);
  double*  stemp = (double*) PyCObject_AsVoidPtr(py_stemp);
  double*  state = (double*) PyCObject_AsVoidPtr(py_state);
  double*  dstate = (double*) PyCObject_AsVoidPtr(py_dstate);
  double*  state0 = (double*) PyCObject_AsVoidPtr(py_state0);
  double*  dmat = (double*) PyCObject_AsVoidPtr(py_dmat);
  int*  ien = (int*) PyCObject_AsVoidPtr(py_ien);
  int*  lm = (int*) PyCObject_AsVoidPtr(py_lm);
  int*  lmx = (int*) PyCObject_AsVoidPtr(py_lmx);
  int*  lmf = (int*) PyCObject_AsVoidPtr(py_lmf);
  int*  ivfamily = (int*) PyCObject_AsVoidPtr(py_ivfamily);
  int*  npar = (int*) PyCObject_AsVoidPtr(py_npar);
  int*  ielindx = (int*) PyCObject_AsVoidPtr(py_ielindx);
  int*  tractionverts = (int*) PyCObject_AsVoidPtr(py_tractionverts);
  double*  tractionvals = (double*) PyCObject_AsVoidPtr(py_tractionvals);
  double*  gauss2d = (double*) PyCObject_AsVoidPtr(py_gauss2d);
  double*  sh2d = (double*) PyCObject_AsVoidPtr(py_sh2d);
  int*  infetype2d = (int*) PyCObject_AsVoidPtr(py_infetype2d);
  double*  prop = (double*) PyCObject_AsVoidPtr(py_prop);
  int*  infmatmod = (int*) PyCObject_AsVoidPtr(py_infmatmod);
  double*  gauss = (double*) PyCObject_AsVoidPtr(py_gauss);
  double*  sh = (double*) PyCObject_AsVoidPtr(py_sh);
  double*  shj = (double*) PyCObject_AsVoidPtr(py_shj);
  int*  infetype = (int*) PyCObject_AsVoidPtr(py_infetype);
  double*  histry = (double*) PyCObject_AsVoidPtr(py_histry);
  double*  rtimdat = (double*) PyCObject_AsVoidPtr(py_rtimdat);
  int*  ntimdat = (int*) PyCObject_AsVoidPtr(py_ntimdat);
  int*  nvisdat = (int*) PyCObject_AsVoidPtr(py_nvisdat);
  int*  maxstp = (int*) PyCObject_AsVoidPtr(py_maxstp);
  double*  delt = (double*) PyCObject_AsVoidPtr(py_delt);
  double*  alfa = (double*) PyCObject_AsVoidPtr(py_alfa);
  int*  maxit = (int*) PyCObject_AsVoidPtr(py_maxit);
  int*  ntdinit = (int*) PyCObject_AsVoidPtr(py_ntdinit);
  int*  lgdef = (int*) PyCObject_AsVoidPtr(py_lgdef);
  double*  utol = (double*) PyCObject_AsVoidPtr(py_utol);
  double*  ftol = (double*) PyCObject_AsVoidPtr(py_ftol);
  double*  etol = (double*) PyCObject_AsVoidPtr(py_etol);
  int*  itmax = (int*) PyCObject_AsVoidPtr(py_itmax);
  double*  rgiter = (double*) PyCObject_AsVoidPtr(py_rgiter);
  double*  skew = (double*) PyCObject_AsVoidPtr(py_skew);
  int*  ncodat = (int*) PyCObject_AsVoidPtr(py_ncodat);
  int*  nunits = (int*) PyCObject_AsVoidPtr(py_nunits);
  int*  nprint = (int*) PyCObject_AsVoidPtr(py_nprint);
  int*  istatout = (int*) PyCObject_AsVoidPtr(py_istatout);
  int*  nstatout = (int*) PyCObject_AsVoidPtr(py_nstatout);

  elastc_f(&A,                                // Sparse matrix arrays
           &rhs,
	   &sol,
	   bextern,                  // Force vectors
	   btraction,
	   bgravity,
	   bconcForce,
	   bintern,
	   bresid,
	   bwink,
	   bwinkx,
	   dispVec,
	   dprev,
	   nforce,
	   grav,
	   x,                                 // Global arrays
	   d,
	   deld,
	   dcur,
	   id,
	   iwink,
	   wink,
	   nsysdat,
	   iddmat,
	   ibond,                    // Boundary condition arrays
	   bond,
	   dx,                       // Slippery node arrays
	   deldx,
	   dxcur,
	   diforc,
	   idx,
	   iwinkx,
	   winkx,
	   idslp,
	   ipslp,
	   idhist,
	   fault,                    // Split node arrays
	   nfault,
	   dfault,
	   tfault,
	   s,                        // Local stiffness matrix arrays
	   stemp,
	   state,                    // Element arrays
	   dstate,
	   state0,
	   dmat,
	   ien,
	   lm,
	   lmx,
	   lmf,
	   ivfamily,
	   npar,
	   ielindx,
	   tractionverts,            // Traction BC arrays
	   tractionvals,
	   gauss2d,
	   sh2d,
	   infetype2d,
	   prop,    // Material property arrays
	   infmatmod,
	   gauss,                    // Element type arrays
	   sh,
	   shj,
	   infetype,
	   histry,                   // Time information
	   rtimdat,
	   ntimdat,
	   nvisdat,
	   maxstp,
	   delt,
	   alfa,
	   maxit,
	   ntdinit,
	   lgdef,
	   utol,
	   ftol,
	   etol,
	   itmax,
	   rgiter,          // Iterative solution information
	   skew,                     // Skew rotation information
	   ncodat,          // Input/output information
	   nunits,
	   nprint,
	   istatout,
	   nstatout,
	   ofile,                   // Output file names
	   pfile,
	   ucdroot,
           &elasticStage,
           &iterateEvent,
	   &errorcode,                        // Error codes
	   errorstring,
	   strlen(ofile),           // String lengths
	   strlen(pfile),
	   strlen(ucdroot),
	   sizeof(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    //int *fail = NULL;
    //int i = fail[0];
    return NULL;
  }
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "ofile:" << ofile
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}




// Assign equation numbers to Iwink and Iwinkx arrays, and compact
// Wink and Winkx arrays to correspond to assigned BC.

char pypylith3d_assign_wink__doc__[] = "";
char pypylith3d_assign_wink__name__[] = "assign_wink";

PyObject * pypylith3d_assign_wink(PyObject *, PyObject *args)
{
  PyObject* py_winkdef;
  PyObject* py_wink;
  PyObject* py_iwinkdef;
  PyObject* py_iwinkid;
  PyObject* py_iwink;
  PyObject* py_id;
  int numnp;
  int nwink;
  int nwinke;

  int ok = PyArg_ParseTuple(args, "OOOOOOiii:assign_wink",
			    &py_winkdef,
			    &py_wink,
			    &py_iwinkdef,
			    &py_iwinkid,
			    &py_iwink,
			    &py_id,
			    &numnp,
			    &nwink,
			    &nwinke);

  if (!ok) {
    return 0;
  }

  double* winkdef = (double*) PyCObject_AsVoidPtr(py_winkdef);
  double* wink = (double*) PyCObject_AsVoidPtr(py_wink);
  int* iwinkdef = (int*) PyCObject_AsVoidPtr(py_iwinkdef);
  int* iwinkid = (int*) PyCObject_AsVoidPtr(py_iwinkid);
  int* iwink = (int*) PyCObject_AsVoidPtr(py_iwink);
  int* id = (int*) PyCObject_AsVoidPtr(py_id);

  assign_wink_f(winkdef,
		wink,
		iwinkdef,
		iwinkid,
		iwink,
		id,
		&numnp,
		&nwink,
		&nwinke);

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "nwink:" << nwink
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}
    

// Create Id and Idx arrays, and determine the number of equations.

char pypylith3d_create_id__doc__[] = "";
char pypylith3d_create_id__name__[] = "create_id";

PyObject * pypylith3d_create_id(PyObject *, PyObject *args)
{
  PyObject* py_id;
  PyObject* py_idx;
  PyObject* py_ibond;
  PyObject* py_nslip;
  PyObject* py_idslp;
  int numslp;
  int numnp;
  int totalNumberSlipperyNodes;

  int ok = PyArg_ParseTuple(args, "OOOOOiii:create_id",
			    &py_id,
			    &py_idx,
			    &py_ibond,
			    &py_nslip,
			    &py_idslp,
			    &numslp,
			    &numnp,
			    &totalNumberSlipperyNodes);

  if (!ok) {
    return 0;
  }

  int* id = (int*) PyCObject_AsVoidPtr(py_id);
  int* idx = (int*) PyCObject_AsVoidPtr(py_idx);
  int* ibond = (int*) PyCObject_AsVoidPtr(py_ibond);
  int* nslip = (int*) PyCObject_AsVoidPtr(py_nslip);
  int* idslp = (int*) PyCObject_AsVoidPtr(py_idslp);
  int neq = 0;

  create_id_f(id,
	  idx,
	  ibond,
	  nslip,
	  idslp,
	  &numslp,
	  &numnp,
	  &totalNumberSlipperyNodes,
	  &neq);

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "neq:" << neq
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i", neq);
}


// Form id array for split nodes

char pypylith3d_id_split__doc__[] = "";
char pypylith3d_id_split__name__[] = "id_split";

PyObject * pypylith3d_id_split(PyObject *, PyObject *args)
{
  PyObject* py_nfault;
  PyObject* py_idftn;
  int numnp;
  int numfn;
  int totalNumberSplitNodes;

  int ok = PyArg_ParseTuple(args, "OOiii:id_split",
			    &py_nfault,
			    &py_idftn,
			    &numnp,
			    &numfn,
			    &totalNumberSplitNodes);

  if (!ok) {
    return 0;
  }

  int* nfault = (int*) PyCObject_AsVoidPtr(py_nfault);
  int* idftn = (int*) PyCObject_AsVoidPtr(py_idftn);

  id_split_f(nfault,
	     idftn,
	     &numnp,
	     &numfn,
	     &totalNumberSplitNodes);

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numfn:" << numfn
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}
    

// Find closest fault neighbors for slippery nodes

char pypylith3d_nfind__doc__[] = "";
char pypylith3d_nfind__name__[] = "nfind";

PyObject * pypylith3d_nfind(PyObject *, PyObject *args)
{
  PyObject* py_x;
  PyObject* py_xtmp;
  PyObject* py_idslp;
  PyObject* py_ipslp;
  PyObject* py_itmp;
  PyObject* py_itmp1;
  PyObject* py_itmp2;
  PyObject* py_nslip;
  int numslp;
  int totalNumberSlipperyNodes;
  int numnp;

  int ok = PyArg_ParseTuple(args, "OOOOOOOOiii:nfind",
			    &py_x,
			    &py_xtmp,
			    &py_idslp,
			    &py_ipslp,
			    &py_itmp,
			    &py_itmp1,
			    &py_itmp2,
			    &py_nslip,
			    &numslp,
			    &totalNumberSlipperyNodes,
			    &numnp);

  if (!ok) {
    return 0;
  }

  double* x = (double*) PyCObject_AsVoidPtr(py_x);
  double* xtmp = (double*) PyCObject_AsVoidPtr(py_xtmp);
  int* idslp = (int*) PyCObject_AsVoidPtr(py_idslp);
  int* ipslp = (int*) PyCObject_AsVoidPtr(py_ipslp);
  int* itmp = (int*) PyCObject_AsVoidPtr(py_itmp);
  int* itmp1 = (int*) PyCObject_AsVoidPtr(py_itmp1);
  int* itmp2 = (int*) PyCObject_AsVoidPtr(py_itmp2);
  int* nslip = (int*) PyCObject_AsVoidPtr(py_nslip);

  nfind_f(x,
	  xtmp,
	  idslp,
	  ipslp,
	  itmp,
	  itmp1,
	  itmp2,
	  nslip,
	  &numslp,
	  &totalNumberSlipperyNodes,
	  &numnp);

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "x:" << x
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Read boundary conditions

char pypylith3d_read_bc__doc__[] = "";
char pypylith3d_read_bc__name__[] = "read_bc";

PyObject * pypylith3d_read_bc(PyObject *, PyObject *args)
{
  PyObject* py_bond;
  double dscale;
  double vscale;
  double fscale;
  PyObject* py_ibond;
  int numnp;
  int numbc;
  int kr;
  char* bcInputFile;

  int ok = PyArg_ParseTuple(args, "OdddOiiis:read_bc",
			    &py_bond,
			    &dscale,
			    &vscale,
			    &fscale,
			    &py_ibond,
			    &numnp,
			    &numbc,
			    &kr,
			    &bcInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* bond = (double*) PyCObject_AsVoidPtr(py_bond);
  int* ibond = (int*) PyCObject_AsVoidPtr(py_ibond);
  int numberConcForces = 0;

  read_bc_f(bond,
	    &dscale,
	    &vscale,
	    &fscale,
	    ibond,
	    &numnp,
	    &numbc,
	    &numberConcForces,
	    &kr,
	    bcInputFile,
	    &errorcode,
	    errorstring,
	    strlen(bcInputFile),
	    sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numbc:" << numbc
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i", numberConcForces);
}


// Read connectivities

char pypylith3d_read_connect__doc__[] = "";
char pypylith3d_read_connect__name__[] = "read_connect";

PyObject * pypylith3d_read_connect(PyObject *, PyObject *args)
{
  PyObject* py_ien;
  PyObject* py_mat;
  int nen;
  int numelv;
  int numnp;
  int nvfamilies;
  int kr;
  char* connectivityInputFile;

  int ok = PyArg_ParseTuple(args, "OOiiiiis:read_connect",
                            &py_ien,
                            &py_mat,
                            &nen,
                            &numelv,
                            &numnp,
                            &nvfamilies,
                            &kr,
                            &connectivityInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* ien = (int*) PyCObject_AsVoidPtr(py_ien);
  int* mat = (int*) PyCObject_AsVoidPtr(py_mat);

  read_connect_f(ien,
		 mat,
		 &nen,
                 &numelv,
                 &numnp,
                 &nvfamilies,
		 &kr,
		 connectivityInputFile,
		 &errorcode,
		 errorstring,
		 strlen(connectivityInputFile),
		 sizeof(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "nvfamilies:" << nvfamilies
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Read coordinates

char pypylith3d_read_coords__doc__[] = "";
char pypylith3d_read_coords__name__[] = "read_coords";

PyObject * pypylith3d_read_coords(PyObject *, PyObject *args)
{
  PyObject* py_x;
  double cscale;
  int numnp;
  int kr;
  char* coordinateInputFile;

  int ok = PyArg_ParseTuple(args, "Odiis:read_coords",
			    &py_x,
			    &cscale,
			    &numnp,
			    &kr,
			    &coordinateInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* x = (double*) PyCObject_AsVoidPtr(py_x);

  read_coords_f(x,
		&cscale,
		&numnp,
		&kr,
		coordinateInputFile,
		&errorcode,
		errorstring,
		strlen(coordinateInputFile),
		sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "cscale:" << cscale
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Read differential forces

char pypylith3d_read_diff__doc__[] = "";
char pypylith3d_read_diff__name__[] = "read_diff";

PyObject * pypylith3d_read_diff(PyObject *, PyObject *args)
{
  PyObject* py_diforc;
  PyObject* py_nslip;
  PyObject* py_idhist;
  int numslp;
  int numdif;
  int numnp;
  int kr;
  char* difile;

  int ok = PyArg_ParseTuple(args, "OOOiiiis:read_diff",
			    &py_diforc,
			    &py_nslip,
			    &py_idhist,
			    &numslp,
			    &numdif,
			    &numnp,
			    &kr,
			    &difile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* diforc = (double*) PyCObject_AsVoidPtr(py_diforc);
  int* nslip = (int*) PyCObject_AsVoidPtr(py_nslip);
  int* idhist = (int*) PyCObject_AsVoidPtr(py_idhist);

  read_diff_f(diforc,
	      nslip,
	      idhist,
	      &numslp,
	      &numdif,
	      &numnp,
	      &kr,
	      difile,
	      &errorcode,
	      errorstring,
	      strlen(difile),
	      sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numdif:" << numdif
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Read time steps at which full output is desired

char pypylith3d_read_fuldat__doc__[] = "";
char pypylith3d_read_fuldat__name__[] = "read_fuldat";

PyObject * pypylith3d_read_fuldat(PyObject *, PyObject *args)
{
  PyObject* py_iprint;
  int icontr;
  int icode;
  int ncycle;
  int lastep;
  int kr;
  char* fofile;

  int ok = PyArg_ParseTuple(args, "Oiiiiis:read_fuldat",
			    &py_iprint,
			    &icontr,
			    &icode,
			    &ncycle,
			    &lastep,
			    &kr,
			    &fofile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* iprint = (int*) PyCObject_AsVoidPtr(py_iprint);

  read_fuldat_f(iprint,
		&icontr,
		&icode,
		&ncycle,
		&lastep,
		&kr,
		fofile,
		&errorcode,
		errorstring,
		strlen(fofile),
		sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "icontr:" << icontr
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Read load histories

char pypylith3d_read_hist__doc__[] = "";
char pypylith3d_read_hist__name__[] = "read_hist";

PyObject * pypylith3d_read_hist(PyObject *, PyObject *args)
{
  PyObject* py_histry;
  PyObject* py_times;
  int nhist;
  int lastep;
  int kr;
  char* hfile;

  int ok = PyArg_ParseTuple(args, "OOiiis:read_hist",
			    &py_histry,
			    &py_times,
			    &nhist,
			    &lastep,
			    &kr,
			    &hfile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* histry = (double*) PyCObject_AsVoidPtr(py_histry);
  double* times = (double*) PyCObject_AsVoidPtr(py_times);

  read_hist_f(histry,
	      times,
	      &nhist,
	      &lastep,
	      &kr,
	      hfile,
	      &errorcode,
	      errorstring,
	      strlen(hfile),
	      sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "nhist:" << nhist
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Read element prestresses

// char pypylith3d_read_prestr__doc__[] = "";
// char pypylith3d_read_prestr__name__[] = "read_prestr";

// PyObject * pypylith3d_read_prestr(PyObject *, PyObject *args)
// {
  // PyObject* py_stn;
  // PyObject* py_st0;
  // PyObject* py_prscal;
  // int numberStressComponents;
  // int numberGaussPoints;
  // int numberPrestressGaussPoints;
  // int numberElements;
  // int numberPrestressEntries;
  // int ipstrs;
  // int idout;
  // int kr;
  // int kw;
  // char* prestressInputFile;
  // char* ofile;

  // int ok = PyArg_ParseTuple(args, "OOOiiiiiiiiiss:read_prestr",
			    // &py_stn,
			    // &py_st0,
			    // &py_prscal,
			    // &numberStressComponents,
			    // &numberGaussPoints,
			    // &numberPrestressGaussPoints,
			    // &numberElements,
			    // &numberPrestressEntries,
			    // &ipstrs,
			    // &idout,
			    // &kr,
			    // &kw,
			    // &prestressInputFile,
			    // &ofile);

  // if (!ok) {
    // return 0;
  // }

  // int errorcode = 0;
  // const int maxsize = 4096;
  // char errorstring[maxsize];
  // double* stn = (double*) PyCObject_AsVoidPtr(py_stn);
  // double* st0 = (double*) PyCObject_AsVoidPtr(py_st0);
  // double* prscal = (double*) PyCObject_AsVoidPtr(py_prscal);

  // read_prestr_f(stn,
		// st0,
		// prscal,
		// &numberStressComponents,
		// &numberGaussPoints,
		// &numberPrestressGaussPoints,
		// &numberElements,
		// &numberPrestressEntries,
		// &ipstrs,
		// &idout,
		// &kr,
		// &kw,
		// &errorcode,
		// prestressInputFile,
		// ofile,
		// strlen(prestressInputFile),
		// strlen(ofile));
    
// if(0 != exceptionhandler(errorcode, errorstring)) {
    // return 0;
  // }

  // journal::debug_t debug("pylith3d");
  // debug
    // << journal::at(__HERE__)
    // << "numberPrestressEntries:" << numberPrestressEntries
    // << journal::endl;

  // return
  // Py_INCREF(Py_None);
  // return Py_None;
// }


// Read local coordinate rotations

char pypylith3d_read_skew__doc__[] = "";
char pypylith3d_read_skew__name__[] = "read_skew";

PyObject * pypylith3d_read_skew(PyObject *, PyObject *args)
{
  PyObject* py_skew;
  double runits;
  int numrot;
  int numnp;
  int iskopt;
  int kr;
  char* skfile;

  int ok = PyArg_ParseTuple(args, "Odiiiis:read_skew",
			    &py_skew,
			    &runits,
			    &numrot,
			    &numnp,
			    &iskopt,
			    &kr,
			    &skfile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* skew = (double*) PyCObject_AsVoidPtr(py_skew);

  read_skew_f(skew,
	      &runits,
	      &numrot,
	      &numnp,
	      &iskopt,
	      &kr,
	      skfile,
	      &errorcode,
	      errorstring,
	      strlen(skfile),
	      sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numrot:" << numrot
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Read slippery node entries

char pypylith3d_read_slip__doc__[] = "";
char pypylith3d_read_slip__name__[] = "read_slip";

PyObject * pypylith3d_read_slip(PyObject *, PyObject *args)
{
  PyObject* py_nslip;
  int numslp;
  int numnp;
  int iskopt;
  int kr;
  char* slfile;

  int ok = PyArg_ParseTuple(args, "Oiiiis:read_slip",
			    &py_nslip,
			    &numslp,
			    &numnp,
			    &iskopt,
			    &kr,
			    &slfile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* nslip = (int*) PyCObject_AsVoidPtr(py_nslip);
  int totalNumberSlipperyNodes = 0;

  read_slip_f(nslip,
	      &numslp,
	      &totalNumberSlipperyNodes,
	      &numnp,
	      &iskopt,
	      &kr,
	      slfile,
	      &errorcode,
	      errorstring,
	      strlen(slfile),
	      sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "totalNumberSlipperyNodes:" << totalNumberSlipperyNodes
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i", totalNumberSlipperyNodes);
}


// Read split node entries

char pypylith3d_read_split__doc__[] = "";
char pypylith3d_read_split__name__[] = "read_split";

PyObject * pypylith3d_read_split(PyObject *, PyObject *args)
{
  PyObject* py_fault;
  PyObject* py_nfault;
  int numfn;
  int numnp;
  int numelv;
  int kr;
  char* spfile;

  int ok = PyArg_ParseTuple(args, "OOiiiis:read_split",
			    &py_fault,
			    &py_nfault,
			    &numfn,
			    &numnp,
			    &numelv,
			    &kr,
			    &spfile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* fault = (double*) PyCObject_AsVoidPtr(py_fault);
  int* nfault = (int*) PyCObject_AsVoidPtr(py_nfault);
  int totalNumberSplitNodes = 0;

  read_split_f(fault,
	       nfault,
	       &numfn,
	       &totalNumberSplitNodes,
	       &numnp,
	       &numelv,
	       &kr,
	       spfile,
	       &errorcode,
	       errorstring,
	       strlen(spfile),
	       sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "totalNumberSplitNodes:" << totalNumberSplitNodes
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i", totalNumberSplitNodes);
}


// Read state output information

char pypylith3d_read_stateout__doc__[] = "";
char pypylith3d_read_stateout__name__[] = "read_stateout";

PyObject * pypylith3d_read_stateout(PyObject *, PyObject *args)
{
  PyObject* py_istatout;
  PyObject* py_nstatout;
  int kr;
  char* stfile;

  int ok = PyArg_ParseTuple(args, "OOis:read_stateout",
			    &py_istatout,
			    &py_nstatout,
			    &kr,
			    &stfile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* istatout = (int*) PyCObject_AsVoidPtr(py_istatout);
  int* nstatout = (int*) PyCObject_AsVoidPtr(py_nstatout);

  read_stateout_f(istatout,
	          nstatout,
	          &kr,
	          stfile,
	          &errorcode,
	          errorstring,
	          strlen(stfile),
	          sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "istatout:" << istatout
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Read time step data

char pypylith3d_read_timdat__doc__[] = "";
char pypylith3d_read_timdat__name__[] = "read_timdat";

PyObject * pypylith3d_read_timdat(PyObject *, PyObject *args)
{
  PyObject* py_delt;
  PyObject* py_alfa;
  PyObject* py_utol;
  PyObject* py_ftol;
  PyObject* py_etol;
  PyObject* py_times;
  double tunits;
  PyObject* py_maxstp;
  PyObject* py_maxit;
  PyObject* py_ntdinit;
  PyObject* py_lgdef;
  PyObject* py_itmax;
  int nintg;
  int lastep;
  int kr;
  char* timeStepInputFile;

  int ok = PyArg_ParseTuple(args, "OOOOOOdOOOOOiiis:read_timdat",
			    &py_delt,
			    &py_alfa,
			    &py_utol,
			    &py_ftol,
			    &py_etol,
			    &py_times,
			    &tunits,
			    &py_maxstp,
			    &py_maxit,
			    &py_ntdinit,
			    &py_lgdef,
			    &py_itmax,
			    &nintg,
			    &lastep,
			    &kr,
			    &timeStepInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* delt = (double*) PyCObject_AsVoidPtr(py_delt);
  double* alfa = (double*) PyCObject_AsVoidPtr(py_alfa);
  double* utol = (double*) PyCObject_AsVoidPtr(py_utol);
  double* ftol = (double*) PyCObject_AsVoidPtr(py_ftol);
  double* etol = (double*) PyCObject_AsVoidPtr(py_etol);
  double* times = (double*) PyCObject_AsVoidPtr(py_times);
  int* maxstp = (int*) PyCObject_AsVoidPtr(py_maxstp);
  int* maxit = (int*) PyCObject_AsVoidPtr(py_maxit);
  int* ntdinit = (int*) PyCObject_AsVoidPtr(py_ntdinit);
  int* lgdef = (int*) PyCObject_AsVoidPtr(py_lgdef);
  int* itmax = (int*) PyCObject_AsVoidPtr(py_itmax);

  read_timdat_f(delt,
		alfa,
		utol,
		ftol,
		etol,
		times,
		&tunits,
		maxstp,
		maxit,
		ntdinit,
		lgdef,
		itmax,
		&nintg,
		&lastep,
		&kr,
		timeStepInputFile,
		&errorcode,
		errorstring,
		strlen(timeStepInputFile),
		sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "lastep:" << lastep
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Read traction BC

char pypylith3d_read_tractions__doc__[] = "";
char pypylith3d_read_tractions__name__[] = "read_tractions";

PyObject * pypylith3d_read_tractions(PyObject *, PyObject *args)
{
  PyObject* py_tractionverts;
  PyObject* py_tractionvals;
  double tscale;
  int numtractions;
  int numberSurfaceElementNodes;
  int kr;
  char* tractionInputFile;

  int ok = PyArg_ParseTuple(args, "OOdiiis:read_tractions",
			    &py_tractionverts,
			    &py_tractionvals,
			    &tscale,
			    &numtractions,
			    &numberSurfaceElementNodes,
			    &kr,
			    &tractionInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* tractionverts = (int*) PyCObject_AsVoidPtr(py_tractionverts);
  double* tractionvals = (double*) PyCObject_AsVoidPtr(py_tractionvals);

  read_tractions_f(tractionverts,
		   tractionvals,
		   &tscale,
		   &numtractions,
		   &numberSurfaceElementNodes,
		   &kr,
		   tractionInputFile,
		   &errorcode,
		   errorstring,
		   strlen(tractionInputFile),
		   sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "tscale:" << tscale
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Read winkler BC

char pypylith3d_read_wink__doc__[] = "";
char pypylith3d_read_wink__name__[] = "read_wink";

PyObject * pypylith3d_read_wink(PyObject *, PyObject *args)
{
  PyObject* py_winkdef;
  PyObject* py_wscal;
  PyObject* py_iwinkdef;
  PyObject* py_iwinkid;
  int nwink;
  int nwinke;
  int kr;
  char* wfile;

  int ok = PyArg_ParseTuple(args, "OOOOiiis:read_wink",
			    &py_winkdef,
			    &py_wscal,
			    &py_iwinkdef,
			    &py_iwinkid,
			    &nwink,
			    &nwinke,
			    &kr,
			    &wfile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* winkdef = (double*) PyCObject_AsVoidPtr(py_winkdef);
  double* wscal = (double*) PyCObject_AsVoidPtr(py_wscal);
  int* iwinkdef = (int*) PyCObject_AsVoidPtr(py_iwinkdef);
  int* iwinkid = (int*) PyCObject_AsVoidPtr(py_iwinkid);

  read_wink_f(winkdef,
	      wscal,
	      iwinkdef,
	      iwinkid,
	      &nwink,
	      &nwinke,
	      &kr,
	      wfile,
	      &errorcode,
	      errorstring,
	      strlen(wfile),
	      sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "nwink:" << nwink
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Sort elements into element families

char pypylith3d_sort_elements__doc__[] = "";
char pypylith3d_sort_elements__name__[] = "sort_elements";

PyObject * pypylith3d_sort_elements(PyObject *, PyObject *args)
{
  PyObject* py_ien;
  PyObject* py_mat;
  PyObject* py_infmatmod;
  PyObject* py_volumeElementFamilyList;
  PyObject* py_ivfamily;
  PyObject* py_iens;
  PyObject* py_ivftmp;
  PyObject* py_indxiel;
  PyObject* py_ielindx;
  int nen;
  int ngauss;
  int maxvfamilies;
  int nvfamilies;
  int prestressFlag;
  int numelv;
  int numnp;

  int ok = PyArg_ParseTuple(args, "OOOOOOOOOiiiiiii:sort_elements",
			    &py_ien,
			    &py_mat,
			    &py_infmatmod,
			    &py_volumeElementFamilyList,
			    &py_ivfamily,
			    &py_iens,
			    &py_ivftmp,
			    &py_indxiel,
			    &py_ielindx,
			    &nen,
			    &ngauss,
			    &maxvfamilies,
			    &nvfamilies,
			    &prestressFlag,
			    &numelv,
			    &numnp);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* ien = (int*) PyCObject_AsVoidPtr(py_ien);
  int* mat = (int*) PyCObject_AsVoidPtr(py_mat);
  int* infmatmod = (int*) PyCObject_AsVoidPtr(py_infmatmod);
  int* ivflist = (int*) PyCObject_AsVoidPtr(py_volumeElementFamilyList);
  int* ivfamily = (int*) PyCObject_AsVoidPtr(py_ivfamily);
  int* iens = (int*) PyCObject_AsVoidPtr(py_iens);
  int* ivftmp = (int*) PyCObject_AsVoidPtr(py_ivftmp);
  int* indxiel = (int*) PyCObject_AsVoidPtr(py_indxiel);
  int* ielindx = (int*) PyCObject_AsVoidPtr(py_ielindx);
  int stateSize = 0;
  int state0Size = 0;
  int propertySize = 0;

  sort_elements_f(ien,
		  mat,
		  infmatmod,
		  ivflist,
		  ivfamily,
		  iens,
		  ivftmp,
		  indxiel,
		  ielindx,
		  &nen,
		  &ngauss,
		  &maxvfamilies,
		  &nvfamilies,
		  &prestressFlag,
		  &numelv,
		  &numnp,
		  &stateSize,
		  &state0Size,
		  &propertySize,
		  &errorcode,
		  errorstring,
		  sizeof(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "stateSize:" << stateSize
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("iii",
		       stateSize,
		       state0Size,
		       propertySize);
}

// Sort slippery nodes according to reordered elements

char pypylith3d_sort_slip_nodes__doc__[] = "";
char pypylith3d_sort_slip_nodes__name__[] = "sort_slip_nodes";

PyObject * pypylith3d_sort_slip_nodes(PyObject *, PyObject *args)
{
  PyObject* py_nslip;
  PyObject* py_indxiel;
  int numslp;
  int numelv;

  int ok = PyArg_ParseTuple(args, "OOii:sort_slip_nodes",
			    &py_nslip,
			    &py_indxiel,
			    &numslp,
			    &numelv);

  if (!ok) {
    return 0;
  }

  int* nslip = (int*) PyCObject_AsVoidPtr(py_nslip);
  int* indxiel = (int*) PyCObject_AsVoidPtr(py_indxiel);

  sort_slip_nodes_f(nslip,
		    indxiel,
		    &numslp,
		    &numelv);

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numslp:" << numslp
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}
    


// Sort split nodes according to reordered elements

char pypylith3d_sort_split_nodes__doc__[] = "";
char pypylith3d_sort_split_nodes__name__[] = "sort_split_nodes";

PyObject * pypylith3d_sort_split_nodes(PyObject *, PyObject *args)
{
  PyObject* py_nfault;
  PyObject* py_indxiel;
  int numfn;
  int numelv;

  int ok = PyArg_ParseTuple(args, "OOii:sort_split_nodes",
			    &py_nfault,
			    &py_indxiel,
			    &numfn,
			    &numelv);

  if (!ok) {
    return 0;
  }

  int* nfault = (int*) PyCObject_AsVoidPtr(py_nfault);
  int* indxiel = (int*) PyCObject_AsVoidPtr(py_indxiel);

  sort_split_nodes_f(nfault,
		    indxiel,
		    &numfn,
		    &numelv);

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numfn:" << numfn
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Compute number of nonzero entries in stiffness matrix

char pypylith3d_cmp_stiffsz__doc__[] = "";
char pypylith3d_cmp_stiffsz__name__[] = "cmp_stiffsz";

PyObject * pypylith3d_cmp_stiffsz(PyObject *, PyObject *args)
{
  int neq;
  PyObject* py_lm;
  PyObject* py_lmx;
  int numelv;
  int totalNumberSlipperyNodes;
  int nen;

  int ok = PyArg_ParseTuple(args, "iOOiii:cmp_stiffsz",
			    &neq,
			    &py_lm,
			    &py_lmx,
			    &numelv,
			    &totalNumberSlipperyNodes,
			    &nen);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* lm = (int*) PyCObject_AsVoidPtr(py_lm);
  int* lmx = (int*) PyCObject_AsVoidPtr(py_lmx);
  int iwork = 0;

  cmp_stiffsz_f(&neq,
	   	lm,
	   	lmx,
	   	&numelv,
	   	&iwork,
	   	&totalNumberSlipperyNodes,
	   	&nen,
	   	&errorcode,
	   	errorstring,
	   	sizeof(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }
		  
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "iwork:" << iwork
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i", iwork);
}



// Create linked list of nonzero row and column entries in the stiffness matrix

char pypylith3d_lnklst__doc__[] = "";
char pypylith3d_lnklst__name__[] = "lnklst";

PyObject * pypylith3d_lnklst(PyObject *, PyObject *args)
{
  int neq;
  PyObject* py_lm;
  PyObject* py_lmx;
  int numelv;
  int nen;
  int numberVolumeElementEquations;
  PyObject* py_indx;
  PyObject* py_link;
  PyObject* py_nbrs;
  int iwork;
  int totalNumberSlipperyNodes;

  int ok = PyArg_ParseTuple(args, "iOOiiiOOOii:lnklst",
			    &neq,
			    &py_lm,
			    &py_lmx,
			    &numelv,
			    &nen,
			    &numberVolumeElementEquations,
			    &py_indx,
			    &py_link,
			    &py_nbrs,
			    &iwork,
			    &totalNumberSlipperyNodes);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* lm = (int*) PyCObject_AsVoidPtr(py_lm);
  int* lmx = (int*) PyCObject_AsVoidPtr(py_lmx);
  int* indx = (int*) PyCObject_AsVoidPtr(py_indx);
  int* link = (int*) PyCObject_AsVoidPtr(py_link);
  int* nbrs = (int*) PyCObject_AsVoidPtr(py_nbrs);
  int stiffnessMatrixSize =0;
  int stiffnessOffDiagonalSize =0;

  lnklst_f(&neq,
	   lm,
	   lmx,
	   &numelv,
	   &nen,
	   &numberVolumeElementEquations,
	   indx,
	   link,
	   nbrs,
	   &iwork,
	   &stiffnessOffDiagonalSize,
	   &stiffnessMatrixSize,
	   &totalNumberSlipperyNodes,
	   &errorcode,
	   errorstring,
	   sizeof(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }
		  
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "stiffnessMatrixSize:" << stiffnessMatrixSize
    << journal::endl;

  // return
  return Py_BuildValue("ii", stiffnessMatrixSize,
   		  stiffnessOffDiagonalSize);
}


// Localize id array for reference by element

char pypylith3d_local__doc__[] = "";
char pypylith3d_local__name__[] = "local";

PyObject * pypylith3d_local(PyObject *, PyObject *args)
{
  PyObject* py_id;
  int numnp;
  PyObject* py_iens;
  PyObject* py_lm;
  int numelv;
  int nen;

  int ok = PyArg_ParseTuple(args, "OiOOii:local",
			    &py_id,
			    &numnp,
			    &py_iens,
			    &py_lm,
			    &numelv,
			    &nen);

  if (!ok) {
    return 0;
  }

  int* id = (int*) PyCObject_AsVoidPtr(py_id);
  int* iens = (int*) PyCObject_AsVoidPtr(py_iens);
  int* lm = (int*) PyCObject_AsVoidPtr(py_lm);

  local_f(id,
	  &numnp,
	  iens,
	  lm,
	  &numelv,
	  &nen);

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numelv:" << numelv
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}
    

// Localize nfault array for reference by element

char pypylith3d_localf__doc__[] = "";
char pypylith3d_localf__name__[] = "localf";

PyObject * pypylith3d_localf(PyObject *, PyObject *args)
{
  PyObject* py_iens;
  PyObject* py_lmf;
  int numelv;
  PyObject* py_nfault;
  int numfn;
  int nen;

  int ok = PyArg_ParseTuple(args, "OOiOii:localf",
			    &py_iens,
			    &py_lmf,
			    &numelv,
			    &py_nfault,
			    &numfn,
			    &nen);

  if (!ok) {
    return 0;
  }

  int* iens = (int*) PyCObject_AsVoidPtr(py_iens);
  int* lmf = (int*) PyCObject_AsVoidPtr(py_lmf);
  int* nfault = (int*) PyCObject_AsVoidPtr(py_nfault);

  localf_f(iens,
	   lmf,
	   &numelv,
	   nfault,
	   &numfn,
	   &nen);
		  
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numfn:" << numfn
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}
    

// Localize idx array for reference by element

char pypylith3d_localx__doc__[] = "";
char pypylith3d_localx__name__[] = "localx";

PyObject * pypylith3d_localx(PyObject *, PyObject *args)
{
  PyObject* py_idx;
  int numnp;
  PyObject* py_iens;
  PyObject* py_lmx;
  int numelv;
  PyObject* py_nslip;
  int numslp;
  int nen;

  int ok = PyArg_ParseTuple(args, "OiOOiOii:localx",
			    &py_idx,
  			    &numnp,
  			    &py_iens,
  			    &py_lmx,
  			    &numelv,
  			    &py_nslip,
			    &numslp,
  			    &nen);

  if (!ok) {
    return 0;
  }

  int* idx = (int*) PyCObject_AsVoidPtr(py_idx);
  int* iens = (int*) PyCObject_AsVoidPtr(py_iens);
  int* lmx = (int*) PyCObject_AsVoidPtr(py_lmx);
  int* nslip = (int*) PyCObject_AsVoidPtr(py_nslip);

  localx_f(idx,
	   &numnp,
	   iens,
	   lmx,
	   &numelv,
	   nslip,
	   &numslp,
	   &nen);

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numelv:" << numelv
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}

// Transform linked list into index array for modified sparse row format

char pypylith3d_makemsr__doc__[] = "";
char pypylith3d_makemsr__name__[] = "makemsr";

PyObject * pypylith3d_makemsr(PyObject *, PyObject *args)
{
  PyObject* pyA;
  PyObject* py_indx;
  PyObject* py_link;
  PyObject* py_nbrs;
  int neq;
  int stiffnessMatrixSize;
  int iwork;

  int ok = PyArg_ParseTuple(args, "OOOOiii:makemsr",
                            &pyA,
			    &py_indx,
			    &py_link,
			    &py_nbrs,
			    &neq,
			    &stiffnessMatrixSize,
			    &iwork);

  if (!ok) {
    return 0;
  }
  Mat A = (Mat) PyCObject_AsVoidPtr(pyA);
  int* indx = (int*) PyCObject_AsVoidPtr(py_indx);
  int* link = (int*) PyCObject_AsVoidPtr(py_link);
  int* nbrs = (int*) PyCObject_AsVoidPtr(py_nbrs);
  int nmin = 0;
  int nmax = 0;
  double wavg = 0.0;

  makemsr_f(&A,
	    indx,
	    link,
	    nbrs,
	    &neq,
	    &stiffnessMatrixSize,
	    &iwork,
	    &nmin,
	    &nmax,
	    &wavg);

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "iwork:" << iwork
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  // return Py_None;
  return Py_BuildValue("iid",
		       nmin,
		       nmax,
		       wavg);
}



// Perform the time-dependent solution.

char pypylith3d_viscos__doc__[] = "";
char pypylith3d_viscos__name__[] = "viscos";

PyObject * pypylith3d_viscos(PyObject *, PyObject *args)
{
  PyObject* pyA;                              // Sparse matrix arrays
  PyObject* pyRhs;
  PyObject* pySol;
  PyObject* py_bextern;               // Force vectors
  PyObject* py_btraction;
  PyObject* py_bgravity;
  PyObject* py_bconcForce;
  PyObject* py_bintern;
  PyObject* py_bresid;
  PyObject* py_bwink;
  PyObject* py_bwinkx;
  PyObject* py_dispVec;
  PyObject* py_dprev;
  PyObject* py_nforce;
  PyObject* py_grav;
  PyObject* py_x;                             // Global arrays
  PyObject* py_d;
  PyObject* py_deld;
  PyObject* py_dcur;
  PyObject* py_id;
  PyObject* py_iwink;
  PyObject* py_wink;
  PyObject* py_nsysdat;
  PyObject* py_iddmat;
  PyObject* py_ibond;                 // Boundary condition arrays
  PyObject* py_bond;
  PyObject* py_dx;                    // Slippery node arrays
  PyObject* py_deldx;
  PyObject* py_dxcur;
  PyObject* py_diforc;
  PyObject* py_idx;
  PyObject* py_iwinkx;
  PyObject* py_winkx;
  PyObject* py_idslp;
  PyObject* py_ipslp;
  PyObject* py_idhist;
  PyObject* py_fault;                 // Split node arrays
  PyObject* py_nfault;
  PyObject* py_dfault;
  PyObject* py_tfault;
  PyObject* py_s;                     // Local stiffness matrix arrays
  PyObject* py_stemp;
  PyObject* py_state;                 // Element arrays
  PyObject* py_dstate;
  PyObject* py_state0;
  PyObject* py_dmat;
  PyObject* py_ien;
  PyObject* py_lm;
  PyObject* py_lmx;
  PyObject* py_lmf;
  PyObject* py_ivfamily;
  PyObject* py_npar;
  PyObject* py_ielindx;
  PyObject* py_tractionverts;         // Traction BC arrays
  PyObject* py_tractionvals;
  PyObject* py_gauss2d;
  PyObject* py_sh2d;
  PyObject* py_infetype2d;
  PyObject* py_prop; // Material property arrays
  PyObject* py_infmatmod;
  PyObject* py_gauss;                 // Element type arrays
  PyObject* py_sh;
  PyObject* py_shj;
  PyObject* py_infetype;
  PyObject* py_histry;                // Time information
  PyObject* py_rtimdat;
  PyObject* py_ntimdat;
  PyObject* py_nvisdat;
  PyObject* py_maxstp;
  PyObject* py_delt;
  PyObject* py_alfa;
  PyObject* py_maxit;
  PyObject* py_ntdinit;
  PyObject* py_lgdef;
  PyObject* py_utol;
  PyObject* py_ftol;
  PyObject* py_etol;
  PyObject* py_itmax;
  PyObject* py_rgiter;       // Iterative solution information
  PyObject* py_skew;                  // Skew rotation information
  PyObject* py_iprint;                // Input/output information
  PyObject* py_ncodat;
  PyObject* py_nunits;
  PyObject* py_nprint;
  PyObject* py_istatout;
  PyObject* py_nstatout;
  char* ofile;                      // Output file names
  char* pfile;
  char* ucdroot;
  int viscousStage;
  int iterateEvent;

  int ok = PyArg_ParseTuple(args, "OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOsssii:viscos",
			    &pyA,                              // Sparse matrix arrays
			    &pyRhs,
			    &pySol,
			    &py_bextern,               // Force vectors
			    &py_btraction,
			    &py_bgravity,
			    &py_bconcForce,
			    &py_bintern,
			    &py_bresid,
			    &py_bwink,
			    &py_bwinkx,
			    &py_dispVec,
			    &py_dprev,
			    &py_nforce,
			    &py_grav,
			    &py_x,                             // Global arrays
			    &py_d,
			    &py_deld,
			    &py_dcur,
			    &py_id,
			    &py_iwink,
			    &py_wink,
			    &py_nsysdat,
			    &py_iddmat,
			    &py_ibond,                 // Boundary condition arrays
			    &py_bond,
			    &py_dx,                    // Slippery node arrays
			    &py_deldx,
			    &py_dxcur,
			    &py_diforc,
			    &py_idx,
			    &py_iwinkx,
			    &py_winkx,
			    &py_idslp,
			    &py_ipslp,
			    &py_idhist,
			    &py_fault,                 // Split node arrays
			    &py_nfault,
			    &py_dfault,
			    &py_tfault,
			    &py_s,                     // Local stiffness matrix arrays
			    &py_stemp,
			    &py_state,                 // Element arrays
			    &py_dstate,
			    &py_state0,
			    &py_dmat,
			    &py_ien,
			    &py_lm,
			    &py_lmx,
			    &py_lmf,
			    &py_ivfamily,
			    &py_npar,
			    &py_ielindx,
                            &py_tractionverts,         // Traction BC arrays
                            &py_tractionvals,
                            &py_gauss2d,
                            &py_sh2d,
                            &py_infetype2d,
			    &py_prop, // Material property arrays
			    &py_infmatmod,
			    &py_gauss,                 // Element type arrays
			    &py_sh,
			    &py_shj,
			    &py_infetype,
			    &py_histry,                // Time information
			    &py_rtimdat,
			    &py_ntimdat,
			    &py_nvisdat,
			    &py_maxstp,
			    &py_delt,
			    &py_alfa,
			    &py_maxit,
			    &py_ntdinit,
			    &py_lgdef,
			    &py_utol,
			    &py_ftol,
			    &py_etol,
			    &py_itmax,
			    &py_rgiter,       // Iterative solution information
			    &py_skew,                  // Skew rotation information
			    &py_iprint,                // Input/output information
			    &py_ncodat,
			    &py_nunits,
			    &py_nprint,
			    &py_istatout,
			    &py_nstatout,
			    &ofile,                  // Output file names
			    &pfile,
			    &ucdroot,
                            &viscousStage,
                            &iterateEvent);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  Mat      A = (Mat) PyCObject_AsVoidPtr(pyA);
  Vec      rhs = (Vec) PyCObject_AsVoidPtr(pyRhs);
  Vec      sol = (Vec) PyCObject_AsVoidPtr(pySol);
  double*  bextern = (double*) PyCObject_AsVoidPtr(py_bextern);
  double*  btraction = (double*) PyCObject_AsVoidPtr(py_btraction);
  double*  bgravity = (double*) PyCObject_AsVoidPtr(py_bgravity);
  double*  bconcForce = (double*) PyCObject_AsVoidPtr(py_bconcForce);
  double*  bintern = (double*) PyCObject_AsVoidPtr(py_bintern);
  double*  bresid = (double*) PyCObject_AsVoidPtr(py_bresid);
  double*  bwink = (double*) PyCObject_AsVoidPtr(py_bwink);
  double*  bwinkx = (double*) PyCObject_AsVoidPtr(py_bwinkx);
  double*  dispVec = (double*) PyCObject_AsVoidPtr(py_dispVec);
  double*  dprev = (double*) PyCObject_AsVoidPtr(py_dprev);
  int*  nforce = (int*) PyCObject_AsVoidPtr(py_nforce);
  double*  grav = (double*) PyCObject_AsVoidPtr(py_grav);
  double*  x = (double*) PyCObject_AsVoidPtr(py_x);
  double*  d = (double*) PyCObject_AsVoidPtr(py_d);
  double*  deld = (double*) PyCObject_AsVoidPtr(py_deld);
  double*  dcur = (double*) PyCObject_AsVoidPtr(py_dcur);
  int*  id = (int*) PyCObject_AsVoidPtr(py_id);
  int*  iwink = (int*) PyCObject_AsVoidPtr(py_iwink);
  double*  wink = (double*) PyCObject_AsVoidPtr(py_wink);
  int*  nsysdat = (int*) PyCObject_AsVoidPtr(py_nsysdat);
  int*  iddmat = (int*) PyCObject_AsVoidPtr(py_iddmat);
  int*  ibond = (int*) PyCObject_AsVoidPtr(py_ibond);
  double*  bond = (double*) PyCObject_AsVoidPtr(py_bond);
  double*  dx = (double*) PyCObject_AsVoidPtr(py_dx);
  double*  deldx = (double*) PyCObject_AsVoidPtr(py_deldx);
  double*  dxcur = (double*) PyCObject_AsVoidPtr(py_dxcur);
  double*  diforc = (double*) PyCObject_AsVoidPtr(py_diforc);
  int*  idx = (int*) PyCObject_AsVoidPtr(py_idx);
  int*  iwinkx = (int*) PyCObject_AsVoidPtr(py_iwinkx);
  double*  winkx = (double*) PyCObject_AsVoidPtr(py_winkx);
  int*  idslp = (int*) PyCObject_AsVoidPtr(py_idslp);
  int*  ipslp = (int*) PyCObject_AsVoidPtr(py_ipslp);
  int*  idhist = (int*) PyCObject_AsVoidPtr(py_idhist);
  double*  fault = (double*) PyCObject_AsVoidPtr(py_fault);
  int*  nfault = (int*) PyCObject_AsVoidPtr(py_nfault);
  double*  dfault = (double*) PyCObject_AsVoidPtr(py_dfault);
  double*  tfault = (double*) PyCObject_AsVoidPtr(py_tfault);
  double*  s = (double*) PyCObject_AsVoidPtr(py_s);
  double*  stemp = (double*) PyCObject_AsVoidPtr(py_stemp);
  double*  state = (double*) PyCObject_AsVoidPtr(py_state);
  double*  dstate = (double*) PyCObject_AsVoidPtr(py_dstate);
  double*  state0 = (double*) PyCObject_AsVoidPtr(py_state0);
  double*  dmat = (double*) PyCObject_AsVoidPtr(py_dmat);
  int*  ien = (int*) PyCObject_AsVoidPtr(py_ien);
  int*  lm = (int*) PyCObject_AsVoidPtr(py_lm);
  int*  lmx = (int*) PyCObject_AsVoidPtr(py_lmx);
  int*  lmf = (int*) PyCObject_AsVoidPtr(py_lmf);
  int*  ivfamily = (int*) PyCObject_AsVoidPtr(py_ivfamily);
  int*  npar = (int*) PyCObject_AsVoidPtr(py_npar);
  int*  ielindx = (int*) PyCObject_AsVoidPtr(py_ielindx);
  int*  tractionverts = (int*) PyCObject_AsVoidPtr(py_tractionverts);
  double*  tractionvals = (double*) PyCObject_AsVoidPtr(py_tractionvals);
  double*  gauss2d = (double*) PyCObject_AsVoidPtr(py_gauss2d);
  double*  sh2d = (double*) PyCObject_AsVoidPtr(py_sh2d);
  int*  infetype2d = (int*) PyCObject_AsVoidPtr(py_infetype2d);
  double*  prop = (double*) PyCObject_AsVoidPtr(py_prop);
  int*  infmatmod = (int*) PyCObject_AsVoidPtr(py_infmatmod);
  double*  gauss = (double*) PyCObject_AsVoidPtr(py_gauss);
  double*  sh = (double*) PyCObject_AsVoidPtr(py_sh);
  double*  shj = (double*) PyCObject_AsVoidPtr(py_shj);
  int*  infetype = (int*) PyCObject_AsVoidPtr(py_infetype);
  double*  histry = (double*) PyCObject_AsVoidPtr(py_histry);
  double*  rtimdat = (double*) PyCObject_AsVoidPtr(py_rtimdat);
  int*  ntimdat = (int*) PyCObject_AsVoidPtr(py_ntimdat);
  int*  nvisdat = (int*) PyCObject_AsVoidPtr(py_nvisdat);
  int*  maxstp = (int*) PyCObject_AsVoidPtr(py_maxstp);
  double*  delt = (double*) PyCObject_AsVoidPtr(py_delt);
  double*  alfa = (double*) PyCObject_AsVoidPtr(py_alfa);
  int*  maxit = (int*) PyCObject_AsVoidPtr(py_maxit);
  int*  ntdinit = (int*) PyCObject_AsVoidPtr(py_ntdinit);
  int*  lgdef = (int*) PyCObject_AsVoidPtr(py_lgdef);
  double*  utol = (double*) PyCObject_AsVoidPtr(py_utol);
  double*  ftol = (double*) PyCObject_AsVoidPtr(py_ftol);
  double*  etol = (double*) PyCObject_AsVoidPtr(py_etol);
  int*  itmax = (int*) PyCObject_AsVoidPtr(py_itmax);
  double*  rgiter = (double*) PyCObject_AsVoidPtr(py_rgiter);
  double*  skew = (double*) PyCObject_AsVoidPtr(py_skew);
  int*  iprint = (int*) PyCObject_AsVoidPtr(py_iprint);
  int*  ncodat = (int*) PyCObject_AsVoidPtr(py_ncodat);
  int*  nunits = (int*) PyCObject_AsVoidPtr(py_nunits);
  int*  nprint = (int*) PyCObject_AsVoidPtr(py_nprint);
  int*  istatout = (int*) PyCObject_AsVoidPtr(py_istatout);
  int*  nstatout = (int*) PyCObject_AsVoidPtr(py_nstatout);


  viscos_f(&A,                            // 1: Sparse matrix arrays
	   &rhs,
	   &sol,
	   bextern,                  // Force vectors
	   btraction,
	   bgravity,
	   bconcForce,
	   bintern,
	   bresid,
	   bwink,                    // 10
	   bwinkx,
	   dispVec,
	   dprev,
	   nforce,
	   grav,
	   x,                                 // Global arrays
	   d,
	   deld,
	   dcur,
	   id,                       // 20
	   iwink,
	   wink,
	   nsysdat,
	   iddmat,
	   ibond,                    // Boundary condition arrays
	   bond,
	   dx,                       // Slippery node arrays
	   deldx,
	   dxcur,
	   diforc,                   // 30
	   idx,
	   iwinkx,
	   winkx,
	   idslp,
	   ipslp,
	   idhist,
	   fault,                    // Split node arrays
	   nfault,
	   dfault,
	   tfault,                   // 40
	   s,                        // Local stiffness matrix arrays
	   stemp,
	   state,                    // Element arrays
	   dstate,
	   state0,
	   dmat,
	   ien,
	   lm,
	   lmx,
	   lmf,                      // 50
	   ivfamily,
	   npar,
	   ielindx,
	   tractionverts,            // Traction BC arrays
	   tractionvals,
	   gauss2d,
	   sh2d,
	   infetype2d,
	   prop,    // Material property arrays
	   infmatmod,                         // 60
	   gauss,                    // Element type arrays
	   sh,
	   shj,
	   infetype,
	   histry,                   // Time information
	   rtimdat,
	   ntimdat,
	   nvisdat,
	   maxstp,
	   delt,                     // 70
	   alfa,
	   maxit,
	   ntdinit,
	   lgdef,
	   utol,
	   ftol,
	   etol,
	   itmax,
	   rgiter,          // Iterative solution information
	   skew,                     // 80: Skew rotation information
	   iprint,                   // Input/output information
	   ncodat,
	   nunits,
	   nprint,
	   istatout,
	   nstatout,
	   ofile,                   // Output file names
	   pfile,
	   ucdroot,
	   &viscousStage,                     // 90: PETSC logging
	   &iterateEvent,
	   &errorcode,                        // Error codes
	   errorstring,
	   strlen(ofile),           // String lengths
	   strlen(pfile),
	   strlen(ucdroot),
	   sizeof(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }


  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "ofile:" << ofile
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}

// Perform the time-dependent solution.

char pypylith3d_viscos_step__doc__[] = "";
char pypylith3d_viscos_step__name__[] = "viscos_step";

PyObject * pypylith3d_viscos_step(PyObject *, PyObject *args)
{
  PyObject* pyA;                              // Sparse matrix arrays
  PyObject* pyRhs;
  PyObject* pySol;
  PyObject* py_bextern;               // Force vectors
  PyObject* py_btraction;
  PyObject* py_bgravity;
  PyObject* py_bconcForce;
  PyObject* py_bintern;
  PyObject* py_bresid;
  PyObject* py_bwink;
  PyObject* py_bwinkx;
  PyObject* py_dispVec;
  PyObject* py_dprev;
  PyObject* py_nforce;
  PyObject* py_grav;
  PyObject* py_x;                             // Global arrays
  PyObject* py_d;
  PyObject* py_deld;
  PyObject* py_dcur;
  PyObject* py_id;
  PyObject* py_iwink;
  PyObject* py_wink;
  PyObject* py_nsysdat;
  PyObject* py_iddmat;
  PyObject* py_ibond;                 // Boundary condition arrays
  PyObject* py_bond;
  PyObject* py_dx;                    // Slippery node arrays
  PyObject* py_deldx;
  PyObject* py_dxcur;
  PyObject* py_diforc;
  PyObject* py_idx;
  PyObject* py_iwinkx;
  PyObject* py_winkx;
  PyObject* py_idslp;
  PyObject* py_ipslp;
  PyObject* py_idhist;
  PyObject* py_fault;                 // Split node arrays
  PyObject* py_nfault;
  PyObject* py_dfault;
  PyObject* py_tfault;
  PyObject* py_s;                     // Local stiffness matrix arrays
  PyObject* py_stemp;
  PyObject* py_state;                 // Element arrays
  PyObject* py_dstate;
  PyObject* py_state0;
  PyObject* py_dmat;
  PyObject* py_ien;
  PyObject* py_lm;
  PyObject* py_lmx;
  PyObject* py_lmf;
  PyObject* py_ivfamily;
  PyObject* py_npar;
  PyObject* py_ielindx;
  PyObject* py_tractionverts;         // Traction BC arrays
  PyObject* py_tractionvals;
  PyObject* py_gauss2d;
  PyObject* py_sh2d;
  PyObject* py_infetype2d;
  PyObject* py_prop; // Material property arrays
  PyObject* py_infmatmod;
  PyObject* py_gauss;                 // Element type arrays
  PyObject* py_sh;
  PyObject* py_shj;
  PyObject* py_infetype;
  PyObject* py_histry;                // Time information
  PyObject* py_rtimdat;
  PyObject* py_ntimdat;
  PyObject* py_nvisdat;
  PyObject* py_maxstp;
  PyObject* py_delt;
  PyObject* py_alfa;
  PyObject* py_maxit;
  PyObject* py_ntdinit;
  PyObject* py_lgdef;
  PyObject* py_utol;
  PyObject* py_ftol;
  PyObject* py_etol;
  PyObject* py_itmax;
  PyObject* py_rgiter;       // Iterative solution information
  PyObject* py_skew;                  // Skew rotation information
  PyObject* py_iprint;                // Input/output information
  PyObject* py_ncodat;
  PyObject* py_nunits;
  PyObject* py_nprint;
  PyObject* py_istatout;
  PyObject* py_nstatout;
  char* ofile;                      // Output file names
  char* pfile;
  char* ucdroot;
  int viscousStage;
  int iterateEvent;
  int ntot, ltim, indexx, jcyc, i, j, skc, nfirst, nstep, lgdefp;
  double time, deltp;
  PyObject *pyGtol;

  int ok = PyArg_ParseTuple(args, "OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOsssiiiiiiiiiiiddiO!:viscos_step",
			    &pyA,                              // Sparse matrix arrays
			    &pyRhs,
			    &pySol,
			    &py_bextern,               // Force vectors
			    &py_btraction,
			    &py_bgravity,
			    &py_bconcForce,
			    &py_bintern,
			    &py_bresid,
			    &py_bwink,
			    &py_bwinkx,
			    &py_dispVec,
			    &py_dprev,
			    &py_nforce,
			    &py_grav,
			    &py_x,                             // Global arrays
			    &py_d,
			    &py_deld,
			    &py_dcur,
			    &py_id,
			    &py_iwink,
			    &py_wink,
			    &py_nsysdat,
			    &py_iddmat,
			    &py_ibond,                 // Boundary condition arrays
			    &py_bond,
			    &py_dx,                    // Slippery node arrays
			    &py_deldx,
			    &py_dxcur,
			    &py_diforc,
			    &py_idx,
			    &py_iwinkx,
			    &py_winkx,
			    &py_idslp,
			    &py_ipslp,
			    &py_idhist,
			    &py_fault,                 // Split node arrays
			    &py_nfault,
			    &py_dfault,
			    &py_tfault,
			    &py_s,                     // Local stiffness matrix arrays
			    &py_stemp,
			    &py_state,                 // Element arrays
			    &py_dstate,
			    &py_state0,
			    &py_dmat,
			    &py_ien,
			    &py_lm,
			    &py_lmx,
			    &py_lmf,
			    &py_ivfamily,
			    &py_npar,
			    &py_ielindx,
                            &py_tractionverts,         // Traction BC arrays
                            &py_tractionvals,
                            &py_gauss2d,
                            &py_sh2d,
                            &py_infetype2d,
			    &py_prop, // Material property arrays
			    &py_infmatmod,
			    &py_gauss,                 // Element type arrays
			    &py_sh,
			    &py_shj,
			    &py_infetype,
			    &py_histry,                // Time information
			    &py_rtimdat,
			    &py_ntimdat,
			    &py_nvisdat,
			    &py_maxstp,
			    &py_delt,
			    &py_alfa,
			    &py_maxit,
			    &py_ntdinit,
			    &py_lgdef,
			    &py_utol,
			    &py_ftol,
			    &py_etol,
			    &py_itmax,
			    &py_rgiter,       // Iterative solution information
			    &py_skew,                  // Skew rotation information
			    &py_iprint,                // Input/output information
			    &py_ncodat,
			    &py_nunits,
			    &py_nprint,
			    &py_istatout,
			    &py_nstatout,
			    &ofile,                  // Output file names
			    &pfile,
			    &ucdroot,
                &viscousStage,
                &iterateEvent,
                            &ntot,
                            &ltim,
                            &indexx,
                            &jcyc,
                            &i,
                            &j,
                            &skc,
                            &nfirst,
                            &nstep,
                            &time,
                            &deltp,
                            &lgdefp,
                            &PyList_Type, &pyGtol);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  Mat      A = (Mat) PyCObject_AsVoidPtr(pyA);
  Vec      rhs = (Vec) PyCObject_AsVoidPtr(pyRhs);
  Vec      sol = (Vec) PyCObject_AsVoidPtr(pySol);
  double*  bextern = (double*) PyCObject_AsVoidPtr(py_bextern);
  double*  btraction = (double*) PyCObject_AsVoidPtr(py_btraction);
  double*  bgravity = (double*) PyCObject_AsVoidPtr(py_bgravity);
  double*  bconcForce = (double*) PyCObject_AsVoidPtr(py_bconcForce);
  double*  bintern = (double*) PyCObject_AsVoidPtr(py_bintern);
  double*  bresid = (double*) PyCObject_AsVoidPtr(py_bresid);
  double*  bwink = (double*) PyCObject_AsVoidPtr(py_bwink);
  double*  bwinkx = (double*) PyCObject_AsVoidPtr(py_bwinkx);
  double*  dispVec = (double*) PyCObject_AsVoidPtr(py_dispVec);
  double*  dprev = (double*) PyCObject_AsVoidPtr(py_dprev);
  int*  nforce = (int*) PyCObject_AsVoidPtr(py_nforce);
  double*  grav = (double*) PyCObject_AsVoidPtr(py_grav);
  double*  x = (double*) PyCObject_AsVoidPtr(py_x);
  double*  d = (double*) PyCObject_AsVoidPtr(py_d);
  double*  deld = (double*) PyCObject_AsVoidPtr(py_deld);
  double*  dcur = (double*) PyCObject_AsVoidPtr(py_dcur);
  int*  id = (int*) PyCObject_AsVoidPtr(py_id);
  int*  iwink = (int*) PyCObject_AsVoidPtr(py_iwink);
  double*  wink = (double*) PyCObject_AsVoidPtr(py_wink);
  int*  nsysdat = (int*) PyCObject_AsVoidPtr(py_nsysdat);
  int*  iddmat = (int*) PyCObject_AsVoidPtr(py_iddmat);
  int*  ibond = (int*) PyCObject_AsVoidPtr(py_ibond);
  double*  bond = (double*) PyCObject_AsVoidPtr(py_bond);
  double*  dx = (double*) PyCObject_AsVoidPtr(py_dx);
  double*  deldx = (double*) PyCObject_AsVoidPtr(py_deldx);
  double*  dxcur = (double*) PyCObject_AsVoidPtr(py_dxcur);
  double*  diforc = (double*) PyCObject_AsVoidPtr(py_diforc);
  int*  idx = (int*) PyCObject_AsVoidPtr(py_idx);
  int*  iwinkx = (int*) PyCObject_AsVoidPtr(py_iwinkx);
  double*  winkx = (double*) PyCObject_AsVoidPtr(py_winkx);
  int*  idslp = (int*) PyCObject_AsVoidPtr(py_idslp);
  int*  ipslp = (int*) PyCObject_AsVoidPtr(py_ipslp);
  int*  idhist = (int*) PyCObject_AsVoidPtr(py_idhist);
  double*  fault = (double*) PyCObject_AsVoidPtr(py_fault);
  int*  nfault = (int*) PyCObject_AsVoidPtr(py_nfault);
  double*  dfault = (double*) PyCObject_AsVoidPtr(py_dfault);
  double*  tfault = (double*) PyCObject_AsVoidPtr(py_tfault);
  double*  s = (double*) PyCObject_AsVoidPtr(py_s);
  double*  stemp = (double*) PyCObject_AsVoidPtr(py_stemp);
  double*  state = (double*) PyCObject_AsVoidPtr(py_state);
  double*  dstate = (double*) PyCObject_AsVoidPtr(py_dstate);
  double*  state0 = (double*) PyCObject_AsVoidPtr(py_state0);
  double*  dmat = (double*) PyCObject_AsVoidPtr(py_dmat);
  int*  ien = (int*) PyCObject_AsVoidPtr(py_ien);
  int*  lm = (int*) PyCObject_AsVoidPtr(py_lm);
  int*  lmx = (int*) PyCObject_AsVoidPtr(py_lmx);
  int*  lmf = (int*) PyCObject_AsVoidPtr(py_lmf);
  int*  ivfamily = (int*) PyCObject_AsVoidPtr(py_ivfamily);
  int*  npar = (int*) PyCObject_AsVoidPtr(py_npar);
  int*  ielindx = (int*) PyCObject_AsVoidPtr(py_ielindx);
  int*  tractionverts = (int*) PyCObject_AsVoidPtr(py_tractionverts);
  double*  tractionvals = (double*) PyCObject_AsVoidPtr(py_tractionvals);
  double*  gauss2d = (double*) PyCObject_AsVoidPtr(py_gauss2d);
  double*  sh2d = (double*) PyCObject_AsVoidPtr(py_sh2d);
  int*  infetype2d = (int*) PyCObject_AsVoidPtr(py_infetype2d);
  double*  prop = (double*) PyCObject_AsVoidPtr(py_prop);
  int*  infmatmod = (int*) PyCObject_AsVoidPtr(py_infmatmod);
  double*  gauss = (double*) PyCObject_AsVoidPtr(py_gauss);
  double*  sh = (double*) PyCObject_AsVoidPtr(py_sh);
  double*  shj = (double*) PyCObject_AsVoidPtr(py_shj);
  int*  infetype = (int*) PyCObject_AsVoidPtr(py_infetype);
  double*  histry = (double*) PyCObject_AsVoidPtr(py_histry);
  double*  rtimdat = (double*) PyCObject_AsVoidPtr(py_rtimdat);
  int*  ntimdat = (int*) PyCObject_AsVoidPtr(py_ntimdat);
  int*  nvisdat = (int*) PyCObject_AsVoidPtr(py_nvisdat);
  int*  maxstp = (int*) PyCObject_AsVoidPtr(py_maxstp);
  double*  delt = (double*) PyCObject_AsVoidPtr(py_delt);
  double*  alfa = (double*) PyCObject_AsVoidPtr(py_alfa);
  int*  maxit = (int*) PyCObject_AsVoidPtr(py_maxit);
  int*  ntdinit = (int*) PyCObject_AsVoidPtr(py_ntdinit);
  int*  lgdef = (int*) PyCObject_AsVoidPtr(py_lgdef);
  double*  utol = (double*) PyCObject_AsVoidPtr(py_utol);
  double*  ftol = (double*) PyCObject_AsVoidPtr(py_ftol);
  double*  etol = (double*) PyCObject_AsVoidPtr(py_etol);
  int*  itmax = (int*) PyCObject_AsVoidPtr(py_itmax);
  double*  rgiter = (double*) PyCObject_AsVoidPtr(py_rgiter);
  double*  skew = (double*) PyCObject_AsVoidPtr(py_skew);
  int*  iprint = (int*) PyCObject_AsVoidPtr(py_iprint);
  int*  ncodat = (int*) PyCObject_AsVoidPtr(py_ncodat);
  int*  nunits = (int*) PyCObject_AsVoidPtr(py_nunits);
  int*  nprint = (int*) PyCObject_AsVoidPtr(py_nprint);
  int*  istatout = (int*) PyCObject_AsVoidPtr(py_istatout);
  int*  nstatout = (int*) PyCObject_AsVoidPtr(py_nstatout);

  double gtol[4];
  if (PyList_Size(pyGtol) != 3) {
    return 0;
  }
  for(int v = 0; v < 3; v++) {
    PyObject *listitem = PyList_GetItem(pyGtol, v);

    if(!PyFloat_Check(listitem)) {
      return 0;
    }
    gtol[v] = PyFloat_AsDouble(listitem);
  }

  viscos_step_f(&A,                            // 1: Sparse matrix arrays
	   &rhs,
	   &sol,
	   bextern,                  // Force vectors
	   btraction,
	   bgravity,
	   bconcForce,
	   bintern,
	   bresid,
	   bwink,                    // 10
	   bwinkx,
	   dispVec,
	   dprev,
	   nforce,
	   grav,
	   x,                                 // Global arrays
	   d,
	   deld,
	   dcur,
	   id,                       // 20
	   iwink,
	   wink,
	   nsysdat,
	   iddmat,
	   ibond,                    // Boundary condition arrays
	   bond,
	   dx,                       // Slippery node arrays
	   deldx,
	   dxcur,
	   diforc,                   // 30
	   idx,
	   iwinkx,
	   winkx,
	   idslp,
	   ipslp,
	   idhist,
	   fault,                    // Split node arrays
	   nfault,
	   dfault,
	   tfault,                   // 40
	   s,                        // Local stiffness matrix arrays
	   stemp,
	   state,                    // Element arrays
	   dstate,
	   state0,
	   dmat,
	   ien,
	   lm,
	   lmx,
	   lmf,                      // 50
	   ivfamily,
	   npar,
		ielindx,
		tractionverts,            // Traction BC arrays
	   tractionvals,
	   gauss2d,
	   sh2d,
	   infetype2d,
	   prop,    // Material property arrays
	   infmatmod,                         // 60
	   gauss,                    // Element type arrays
	   sh,
	   shj,
	   infetype,
	   histry,                   // Time information
	   rtimdat,
	   ntimdat,
	   nvisdat,
	   maxstp,
	   delt,                     // 70
	   alfa,
	   maxit,
	   ntdinit,
	   lgdef,
	   utol,
	   ftol,
	   etol,
	   itmax,
	   rgiter,          // Iterative solution information
	   skew,                     // 80: Skew rotation information
	   iprint,                   // Input/output information
	   ncodat,
	   nunits,
	   nprint,
	   istatout,
	   nstatout,
	   ofile,                   // Output file names
	   pfile,
	   ucdroot,
	   &viscousStage,                     // 90: PETSC logging
	   &iterateEvent,
                &ntot,
                &ltim,
                &indexx,
                &jcyc,
                &i,
                &j,
                &skc,
                &nfirst,
                &nstep,
                &time,
                &deltp,
                &lgdefp,
                gtol,
	   &errorcode,                        // Error codes
	   errorstring,
	   strlen(ofile),           // String lengths
	   strlen(pfile),
	   strlen(ucdroot),
	   sizeof(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  // return
  Py_INCREF(Py_None);
  return Py_None;
}

// Setup the time-dependent solution.

char pypylith3d_viscos_setup__doc__[] = "";
char pypylith3d_viscos_setup__name__[] = "viscos_setup";

PyObject * pypylith3d_viscos_setup(PyObject *, PyObject *args)
{
  PyObject* py_nprint;
  PyObject* py_nunits;
  char* ofile;                      // Output file names
  char* pfile;
  int viscousStage;

  int ok = PyArg_ParseTuple(args, "OOssi:viscos_setup",
			    &py_nprint,
			    &py_nunits,
			    &ofile,                  // Output file names
                &pfile,
                &viscousStage);

  if (!ok) {
    return 0;
  }

  int  errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int *nprint = (int *) PyCObject_AsVoidPtr(py_nprint);
  int *nunits = (int *) PyCObject_AsVoidPtr(py_nunits);

  viscos_setup_f(nprint,
                 nunits,
                 ofile,                   // Output file names
                 pfile,
                 &viscousStage,                     // PETSC logging
                 &errorcode,                        // Error codes
                 errorstring,
                 strlen(ofile),           // String lengths
                 strlen(pfile),
                 sizeof(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }


  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "ofile:" << ofile
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}

// Cleanup the time-dependent solution.

char pypylith3d_viscos_cleanup__doc__[] = "";
char pypylith3d_viscos_cleanup__name__[] = "viscos_cleanup";

PyObject * pypylith3d_viscos_cleanup(PyObject *, PyObject *args)
{
  PyObject* py_ntimdat;
  PyObject* py_nprint;
  PyObject* py_nunits;

  int ok = PyArg_ParseTuple(args, "OOO:viscos_cleanup",&py_ntimdat, &py_nprint, &py_nunits);

  if (!ok) {
    return 0;
  }

  int  errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int *ntimdat = (int *) PyCObject_AsVoidPtr(py_ntimdat);
  int *nprint = (int *) PyCObject_AsVoidPtr(py_nprint);
  int *nunits  = (int *) PyCObject_AsVoidPtr(py_nunits);

  viscos_cleanup_f(ntimdat,
                   nprint,
                   nunits,
                   &errorcode,                        // Error codes
                   errorstring,
                   sizeof(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  // return
  Py_INCREF(Py_None);
  return Py_None;
}



// Write boundary conditions

char pypylith3d_write_bc__doc__[] = "";
char pypylith3d_write_bc__name__[] = "write_bc";

PyObject * pypylith3d_write_bc(PyObject *, PyObject *args)
{
  PyObject* py_bond;
  PyObject* py_ibond;
  int numnp;
  int kw;
  int idout;
  char* ofile;

  int ok = PyArg_ParseTuple(args, "OOiiis:write_bc",
			    &py_bond,
			    &py_ibond,
			    &numnp,
			    &kw,
			    &idout,
			    &ofile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* bond = (double*) PyCObject_AsVoidPtr(py_bond);
  int* ibond = (int*) PyCObject_AsVoidPtr(py_ibond);

  write_bc_f(bond,
	     ibond,
	     &numnp,
	     &kw,
	     &idout,
	     ofile,
	     &errorcode,
	     errorstring,
	     strlen(ofile),
	     sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numnp:" << numnp
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write connectivities

char pypylith3d_write_connect__doc__[] = "";
char pypylith3d_write_connect__name__[] = "write_connect";

PyObject * pypylith3d_write_connect(PyObject *, PyObject *args)
{
  PyObject* py_iens;
  PyObject* py_ivfamily;
  PyObject* py_indxiel;
  int nen;
  int ngauss;
  int numelv;
  int ietypev;
  int nvfamilies;
  int kw;
  int kp;
  int idout;
  int idsk;
  char* ofile;
  char* pfile;

  int ok = PyArg_ParseTuple(args, "OOOiiiiiiiiiss:write_connect",
                            &py_iens,
                            &py_ivfamily,
                            &py_indxiel,
                            &nen,
                            &ngauss,
                            &numelv,
                            &ietypev,
                            &nvfamilies,
                            &kw,
                            &kp,
                            &idout,
                            &idsk,
                            &ofile,
                            &pfile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* iens = (int*) PyCObject_AsVoidPtr(py_iens);
  int* ivfamily = (int*) PyCObject_AsVoidPtr(py_ivfamily);
  int* indxiel = (int*) PyCObject_AsVoidPtr(py_indxiel);

  write_connect_f(iens,
		  ivfamily,
		  indxiel,
		  &nen,
		  &ngauss,
		  &numelv,
		  &ietypev,
		  &nvfamilies,
		  &kw,
		  &kp,
		  &idout,
		  &idsk,
		  ofile,
		  pfile,
		  &errorcode,
		  errorstring,
		  strlen(ofile),
		  strlen(pfile),
		  sizeof(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numelv:" << numelv
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write coordinates

char pypylith3d_write_coords__doc__[] = "";
char pypylith3d_write_coords__name__[] = "write_coords";

PyObject * pypylith3d_write_coords(PyObject *, PyObject *args)
{
  PyObject* py_x;
  int numnp;
  int kw;
  int kp;
  int idout;
  int idsk;
  char* ofile;
  char* pfile;

  int ok = PyArg_ParseTuple(args, "Oiiiiiss:write_coords",
			    &py_x,
			    &numnp,
			    &kw,
			    &kp,
			    &idout,
			    &idsk,
			    &ofile,
			    &pfile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* x = (double*) PyCObject_AsVoidPtr(py_x);

  write_coords_f(x,
		 &numnp,
		 &kw,
		 &kp,
		 &idout,
		 &idsk,
		 ofile,
		 pfile,
		 &errorcode,
		 errorstring,
		 strlen(ofile),
		 strlen(pfile),
		 sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numnp:" << numnp
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write differential forces

char pypylith3d_write_diff__doc__[] = "";
char pypylith3d_write_diff__name__[] = "write_diff";

PyObject * pypylith3d_write_diff(PyObject *, PyObject *args)
{
  PyObject* py_diforc;
  PyObject* py_nslip;
  PyObject* py_idhist;
  int numslp;
  int numdif;
  int numnp;
  int kw;
  int idout;
  char* ofile;

  int ok = PyArg_ParseTuple(args, "OOOiiiiis:write_diff",
			    &py_diforc,
			    &py_nslip,
			    &py_idhist,
			    &numslp,
			    &numdif,
			    &numnp,
			    &kw,
			    &idout,
			    &ofile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* diforc = (double*) PyCObject_AsVoidPtr(py_diforc);
  int* nslip = (int*) PyCObject_AsVoidPtr(py_nslip);
  int* idhist = (int*) PyCObject_AsVoidPtr(py_idhist);

  write_diff_f(diforc,
	       nslip,
	       idhist,
	       &numslp,
	       &numdif,
	       &numnp,
	       &kw,
	       &idout,
	       ofile,
	       &errorcode,
	       errorstring,
	       strlen(ofile),
	       sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numdif:" << numdif
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}

// Write out element info

char pypylith3d_write_element_info__doc__[] = "";
char pypylith3d_write_element_info__name__[] = "write_element_info";

PyObject * pypylith3d_write_element_info(PyObject *, PyObject *args)
{
  int numelv;
  int nen;
  int ngauss;
  int ietypev;
  int intord;
  int ipstrs;
  int ipauto;
  double tpois;
  double tyoungs;
  int kw;
  int idout;
  char* ofile;

  int ok = PyArg_ParseTuple(args, "iiiiiiiddiis:write_element_info",
			    &numelv,
			    &nen,
			    &ngauss,
			    &ietypev,
			    &intord,
			    &ipstrs,
			    &ipauto,
			    &tpois,
			    &tyoungs,
			    &kw,
			    &idout,
			    &ofile);

  if (!ok) {
    return 0;
  }

  write_element_info_f(&numelv,
		       &nen,
		       &ngauss,
		       &ietypev,
		       &intord,
		       &ipstrs,
		       &ipauto,
		       &tpois,
		       &tyoungs,
		       &kw,
		       &idout,
		       ofile,strlen(ofile));

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numelv:" << numelv
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write time steps at which full output is desired

char pypylith3d_write_fuldat__doc__[] = "";
char pypylith3d_write_fuldat__name__[] = "write_fuldat";

PyObject * pypylith3d_write_fuldat(PyObject *, PyObject *args)
{
  PyObject* py_iprint;
  int icontr;
  int icode;
  int ncycle;
  int lastep;
  int kw;
  int kp;
  int idout;
  int idsk;
  char* ofile;
  char* pfile;

  int ok = PyArg_ParseTuple(args, "Oiiiiiiiiss:write_fuldat",
			    &py_iprint,
			    &icontr,
			    &icode,
			    &ncycle,
			    &lastep,
			    &kw,
			    &kp,
			    &idout,
			    &idsk,
			    &ofile,
			    &pfile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* iprint = (int*) PyCObject_AsVoidPtr(py_iprint);

  write_fuldat_f(iprint,
		 &icontr,
		 &icode,
		 &ncycle,
		 &lastep,
		 &kw,
		 &kp,
		 &idout,
		 &idsk,
		 ofile,
		 pfile,
		 &errorcode,
		 errorstring,
		 strlen(ofile),
		 strlen(pfile),
		 sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "icontr:" << icontr
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write out global info

char pypylith3d_write_global_info__doc__[] = "";
char pypylith3d_write_global_info__name__[] = "write_global_info";

PyObject * pypylith3d_write_global_info(PyObject *, PyObject *args)
{
  char* title;
  int idout;
  int idsk;
  int numnp;
  int icode;
  int idebug;
  int kw;
  int kp;
  char* ofile;
  char* pfile;

  int ok = PyArg_ParseTuple(args, "siiiiiiiss:write_global_info",
			    &title,
			    &idout,
			    &idsk,
			    &numnp,
			    &icode,
			    &idebug,
			    &kw,
			    &kp,
			    &ofile,
			    &pfile);

  if (!ok) {
    return 0;
  }

  write_global_info_f(title,
		      &idout,
		      &idsk,
		      &numnp,
		      &icode,
		      &idebug,
		      &kw,
		      &kp,
		      ofile,
		      pfile,
		      strlen(title),
		      strlen(ofile),
		      strlen(pfile));

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "idebug:" << idebug
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write load histories

char pypylith3d_write_hist__doc__[] = "";
char pypylith3d_write_hist__name__[] = "write_hist";

PyObject * pypylith3d_write_hist(PyObject *, PyObject *args)
{
  PyObject* py_histry;
  PyObject* py_times;
  int nhist;
  int lastep;
  int kw;
  int idout;
  char* ofile;

  int ok = PyArg_ParseTuple(args, "OOiiiis:write_hist",
			    &py_histry,
			    &py_times,
			    &nhist,
			    &lastep,
			    &kw,
			    &idout,
			    &ofile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* histry = (double*) PyCObject_AsVoidPtr(py_histry);
  double* times = (double*) PyCObject_AsVoidPtr(py_times);

  write_hist_f(histry,
	       times,
	       &nhist,
	       &lastep,
	       &kw,
	       &idout,
	       ofile,
	       &errorcode,
	       errorstring,
	       strlen(ofile),
	       sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "nhist:" << nhist
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write out material property info

char pypylith3d_write_props__doc__[] = "";
char pypylith3d_write_props__name__[] = "write_props";

PyObject * pypylith3d_write_props(PyObject *, PyObject *args)
{
  PyObject* py_prop;
  PyObject* py_grav;
  PyObject* py_ivfamily;
  PyObject* py_infmatmod;
  int nvfamilies;
  int propertySize;
  int idout;
  int idsk;
  int kw;
  int kp;
  char* ofile;
  char* pfile;

  int ok = PyArg_ParseTuple(args, "OOOOiiiiiiss:write_props",
			    &py_prop,
			    &py_grav,
			    &py_ivfamily,
			    &py_infmatmod,
			    &nvfamilies,
			    &propertySize,
			    &idout,
			    &idsk,
			    &kw,
			    &kp,
			    &ofile,
			    &pfile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* prop = (double*) PyCObject_AsVoidPtr(py_prop);
  double* grav = (double*) PyCObject_AsVoidPtr(py_grav);
  int* ivfamily = (int*) PyCObject_AsVoidPtr(py_ivfamily);
  int* infmatmod = (int*) PyCObject_AsVoidPtr(py_infmatmod);

  write_props_f(prop,
		grav,
		ivfamily,
		infmatmod,
		&nvfamilies,
		&propertySize,
		&idout,
		&idsk,
		&kw,
		&kp,
		ofile,
		pfile,
		&errorcode,
		errorstring,
		strlen(ofile),
		strlen(pfile),
		sizeof(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "nvfamilies:" << nvfamilies
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write local coordinate rotations

char pypylith3d_write_skew__doc__[] = "";
char pypylith3d_write_skew__name__[] = "write_skew";

PyObject * pypylith3d_write_skew(PyObject *, PyObject *args)
{
  PyObject* py_skew;
  int numrot;
  int iskopt;
  int numnp;
  int kw;
  int idout;
  char* ofile;

  int ok = PyArg_ParseTuple(args, "Oiiiiis:write_skew",
			    &py_skew,
			    &numrot,
			    &iskopt,
			    &numnp,
			    &kw,
			    &idout,
			    &ofile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* skew = (double*) PyCObject_AsVoidPtr(py_skew);

  write_skew_f(skew,
	       &numrot,
	       &iskopt,
	       &numnp,
	       &kw,
	       &idout,
	       ofile,
	       &errorcode,
	       errorstring,
	       strlen(ofile),
	       sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numrot:" << numrot
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write slippery node entries

char pypylith3d_write_slip__doc__[] = "";
char pypylith3d_write_slip__name__[] = "write_slip";

PyObject * pypylith3d_write_slip(PyObject *, PyObject *args)
{
  PyObject* py_nslip;
  int numslp;
  int totalNumberSlipperyNodes;
  int kw;
  int kp;
  int idout;
  int idsk;
  char* ofile;
  char* pfile;

  int ok = PyArg_ParseTuple(args, "Oiiiiiiss:write_slip",
			    &py_nslip,
			    &numslp,
			    &totalNumberSlipperyNodes,
			    &kw,
			    &kp,
			    &idout,
			    &idsk,
			    &ofile,
			    &pfile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* nslip = (int*) PyCObject_AsVoidPtr(py_nslip);

  write_slip_f(nslip,
	       &numslp,
	       &totalNumberSlipperyNodes,
	       &kw,
	       &kp,
	       &idout,
	       &idsk,
	       ofile,
	       pfile,
	       &errorcode,
	       errorstring,
	       strlen(ofile),
	       strlen(pfile),
	       sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "totalNumberSlipperyNodes:" << totalNumberSlipperyNodes
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write out sparse matrix info

char pypylith3d_write_sparse_info__doc__[] = "";
char pypylith3d_write_sparse_info__name__[] = "write_sparse_info";

PyObject * pypylith3d_write_sparse_info(PyObject *, PyObject *args)
{
  int neq;
  int stiffnessMatrixSize;
  int nmin;
  int nmax;
  double wavg;
  int idout;
  int kw;
  char* ofile;

  int ok = PyArg_ParseTuple(args, "iiiidiis:write_sparse_info",
			    &neq,
			    &stiffnessMatrixSize,
			    &nmin,
			    &nmax,
			    &wavg,
			    &idout,
			    &kw,
			    &ofile);

  if (!ok) {
    return 0;
  }

  write_sparse_info_f(&neq,
		      &stiffnessMatrixSize,
		      &nmin,
		      &nmax,
		      &wavg,
		      &idout,
		      &kw,
		      ofile,strlen(ofile));
		  
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "nmin:" << nmin
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write split node entries

char pypylith3d_write_split__doc__[] = "";
char pypylith3d_write_split__name__[] = "write_split";

PyObject * pypylith3d_write_split(PyObject *, PyObject *args)
{
  PyObject* py_fault;
  PyObject* py_nfault;
  int numfn;
  int kw;
  int kp;
  int idout;
  int idsk;
  char* ofile;
  char* pfile;

  int ok = PyArg_ParseTuple(args, "OOiiiiiss:write_split",
			    &py_fault,
			    &py_nfault,
			    &numfn,
			    &kw,
			    &kp,
			    &idout,
			    &idsk,
			    &ofile,
			    &pfile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* fault = (double*) PyCObject_AsVoidPtr(py_fault);
  int* nfault = (int*) PyCObject_AsVoidPtr(py_nfault);

  write_split_f(fault,
		nfault,
		&numfn,
		&kw,
		&kp,
		&idout,
		&idsk,
		ofile,
		pfile,
		&errorcode,
		errorstring,
		strlen(ofile),
		strlen(pfile),
		sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numfn:" << numfn
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write split node entries for plot output

char pypylith3d_write_split_plot__doc__[] = "";
char pypylith3d_write_split_plot__name__[] = "write_split_plot";

PyObject * pypylith3d_write_split_plot(PyObject *, PyObject *args)
{
  PyObject* py_idftn;
  int totalNumberSplitNodes;
  int kp;
  int idsk;
  char* pfile;

  int ok = PyArg_ParseTuple(args, "Oiiis:write_split_plot",
			    &py_idftn,
			    &totalNumberSplitNodes,
			    &kp,
			    &idsk,
			    &pfile);

  if (!ok) {
    return 0;
  }

  int* idftn = (int*) PyCObject_AsVoidPtr(py_idftn);

  write_split_plot_f(idftn,
		     &totalNumberSplitNodes,
		     &kp,
		     &idsk,
		     pfile,
		     strlen(pfile));
    

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "totalNumberSplitNodes:" << totalNumberSplitNodes
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write state output information

char pypylith3d_write_stateout__doc__[] = "";
char pypylith3d_write_stateout__name__[] = "write_stateout";

PyObject * pypylith3d_write_stateout(PyObject *, PyObject *args)
{
  PyObject* py_istatout;
  PyObject* py_nstatout;
  int kw;
  int kp;
  int idout;
  int idsk;
  char* ofile;
  char* pfile;

  int ok = PyArg_ParseTuple(args, "OOiiiiss:write_stateout",
			    &py_istatout,
			    &py_nstatout,
			    &kw,
			    &kp,
			    &idout,
			    &idsk,
			    &ofile,
			    &pfile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* istatout = (int*) PyCObject_AsVoidPtr(py_istatout);
  int* nstatout = (int*) PyCObject_AsVoidPtr(py_nstatout);

  write_stateout_f(istatout,
		   nstatout,
		   &kw,
		   &kp,
		   &idout,
		   &idsk,
		   ofile,
		   pfile,
		   &errorcode,
		   errorstring,
		   strlen(ofile),
		   strlen(pfile),
		   sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "istatout:" << istatout
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write out stress integration parameters

char pypylith3d_write_strscomp__doc__[] = "";
char pypylith3d_write_strscomp__name__[] = "write_strscomp";

PyObject * pypylith3d_write_strscomp(PyObject *, PyObject *args)
{
  double stol;
  double dtol;
  double epert;
  int kw;
  int idout;
  char* ofile;

  int ok = PyArg_ParseTuple(args, "dddiis:write_strscomp",
			    &stol,
			    &dtol,
			    &epert,
			    &kw,
			    &idout,
			    &ofile);

  if (!ok) {
    return 0;
  }

  write_strscomp_f(&stol,
		   &dtol,
		   &epert,
		   &kw,
		   &idout,
		   ofile,strlen(ofile));

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "epert:" << epert
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write out subiteration convergence parameters

char pypylith3d_write_subiter__doc__[] = "";
char pypylith3d_write_subiter__name__[] = "write_subiter";

PyObject * pypylith3d_write_subiter(PyObject *, PyObject *args)
{
  int nprevdflag;
  int kw;
  int idout;
  char* ofile;

  int ok = PyArg_ParseTuple(args, "iiis:write_subiter",
			    &nprevdflag,
			    &kw,
			    &idout,
			    &ofile);

  if (!ok) {
    return 0;
  }


  write_subiter_f(&nprevdflag,
		  &kw,
		  &idout,
		  ofile,strlen(ofile));
		  
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "nprevdflag:" << nprevdflag
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write time step data

char pypylith3d_write_timdat__doc__[] = "";
char pypylith3d_write_timdat__name__[] = "write_timdat";

PyObject * pypylith3d_write_timdat(PyObject *, PyObject *args)
{
  PyObject* py_delt;
  PyObject* py_alfa;
  PyObject* py_utol;
  PyObject* py_ftol;
  PyObject* py_etol;
  PyObject* py_times;
  PyObject* py_maxstp;
  PyObject* py_maxit;
  PyObject* py_ntdinit;
  PyObject* py_lgdef;
  PyObject* py_itmax;
  int nintg;
  int lastep;
  int kw;
  int idout;
  char* ofile;

  int ok = PyArg_ParseTuple(args, "OOOOOOOOOOOiiiis:write_timdat",
			    &py_delt,
			    &py_alfa,
			    &py_utol,
			    &py_ftol,
			    &py_etol,
			    &py_times,
			    &py_maxstp,
			    &py_maxit,
			    &py_ntdinit,
			    &py_lgdef,
			    &py_itmax,
			    &nintg,
			    &lastep,
			    &kw,
			    &idout,
			    &ofile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* delt = (double*) PyCObject_AsVoidPtr(py_delt);
  double* alfa = (double*) PyCObject_AsVoidPtr(py_alfa);
  double* utol = (double*) PyCObject_AsVoidPtr(py_utol);
  double* ftol = (double*) PyCObject_AsVoidPtr(py_ftol);
  double* etol = (double*) PyCObject_AsVoidPtr(py_etol);
  double* times = (double*) PyCObject_AsVoidPtr(py_times);
  int* maxstp = (int*) PyCObject_AsVoidPtr(py_maxstp);
  int* maxit = (int*) PyCObject_AsVoidPtr(py_maxit);
  int* ntdinit = (int*) PyCObject_AsVoidPtr(py_ntdinit);
  int* lgdef = (int*) PyCObject_AsVoidPtr(py_lgdef);
  int* itmax = (int*) PyCObject_AsVoidPtr(py_itmax);

  write_timdat_f(delt,
		 alfa,
		 utol,
		 ftol,
		 etol,
		 times,
		 maxstp,
		 maxit,
		 ntdinit,
		 lgdef,
		 itmax,
		 &nintg,
		 &lastep,
		 &kw,
		 &idout,
		 ofile,
		 &errorcode,
		 errorstring,
		 strlen(ofile),
		 sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "lastep:" << lastep
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write traction BC

char pypylith3d_write_tractions__doc__[] = "";
char pypylith3d_write_tractions__name__[] = "write_tractions";

PyObject * pypylith3d_write_tractions(PyObject *, PyObject *args)
{
  PyObject* py_tractionverts;
  PyObject* py_tractionvals;
  int numtractions;
  int numberSurfaceElementNodes;
  int kw;
  int idout;
  char* ofile;

  int ok = PyArg_ParseTuple(args, "OOiiiis:write_tractions",
			    &py_tractionverts,
			    &py_tractionvals,
			    &numtractions,
			    &numberSurfaceElementNodes,
			    &kw,
			    &idout,
			    &ofile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* tractionverts = (int*) PyCObject_AsVoidPtr(py_tractionverts);
  double* tractionvals = (double*) PyCObject_AsVoidPtr(py_tractionvals);

  write_tractions_f(tractionverts,
		    tractionvals,
		    &numtractions,
		    &numberSurfaceElementNodes,
		    &kw,
		    &idout,
		    ofile,
		    &errorcode,
		    errorstring,
		    strlen(ofile),
		    sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numtractions:" << numtractions
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write mesh info to UCD file

char pypylith3d_write_ucd_mesh__doc__[] = "";
char pypylith3d_write_ucd_mesh__name__[] = "write_ucd_mesh";

PyObject * pypylith3d_write_ucd_mesh(PyObject *, PyObject *args)
{
  PyObject* py_x;
  int numnp;
  PyObject* py_iens;
  PyObject* py_ivfamily;
  int numelv;
  int nvfamilies;
  PyObject* py_sh;
  int nen;
  int ngauss;
  int ietypev;
  PyObject* py_istatout;
  PyObject* py_nstatout;
  int kucd;
  int iucd;
  char* ucdroot;

  int ok = PyArg_ParseTuple(args, "OiOOiiOiiiOOiis:write_ucd_mesh",
			    &py_x,
			    &numnp,
			    &py_iens,
			    &py_ivfamily,
			    &numelv,
			    &nvfamilies,
			    &py_sh,
			    &nen,
			    &ngauss,
			    &ietypev,
			    &py_istatout,
			    &py_nstatout,
			    &kucd,
			    &iucd,
			    &ucdroot);

  if (!ok) {
    return 0;
  }

  double* x = (double*) PyCObject_AsVoidPtr(py_x);
  int* iens = (int*) PyCObject_AsVoidPtr(py_iens);
  int* ivfamily = (int*) PyCObject_AsVoidPtr(py_ivfamily);
  double* sh = (double*) PyCObject_AsVoidPtr(py_sh);
  int* istatout = (int*) PyCObject_AsVoidPtr(py_istatout);
  int* nstatout = (int*) PyCObject_AsVoidPtr(py_nstatout);

  write_ucd_mesh_f(x,
		   &numnp,
		   iens,
		   ivfamily,
		   &numelv,
		   &nvfamilies,
		   sh,
		   &nen,
		   &ngauss,
		   &ietypev,
		   istatout,
		   nstatout,
		   &kucd,
		   &iucd,
		   ucdroot,strlen(ucdroot));
		  
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numelv:" << numelv
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write winkler BC

char pypylith3d_write_wink__doc__[] = "";
char pypylith3d_write_wink__name__[] = "write_wink";

PyObject * pypylith3d_write_wink(PyObject *, PyObject *args)
{
  PyObject* py_winkdef;
  PyObject* py_iwinkdef;
  PyObject* py_iwinkid;
  int nwinke;
  int kw;
  int idout;
  char* ofile;

  int ok = PyArg_ParseTuple(args, "OOOiiis:write_wink",
			    &py_winkdef,
			    &py_iwinkdef,
			    &py_iwinkid,
			    &nwinke,
			    &kw,
			    &idout,
			    &ofile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* winkdef = (double*) PyCObject_AsVoidPtr(py_winkdef);
  int* iwinkdef = (int*) PyCObject_AsVoidPtr(py_iwinkdef);
  int* iwinkid = (int*) PyCObject_AsVoidPtr(py_iwinkid);

  write_wink_f(winkdef,
	       iwinkdef,
	       iwinkid,
	       &nwinke,
	       &kw,
	       &idout,
	       ofile,
	       &errorcode,
	       errorstring,
	       strlen(ofile),
	       sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "nwinke:" << nwinke
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Write slippery winkler BC

char pypylith3d_write_winkx__doc__[] = "";
char pypylith3d_write_winkx__name__[] = "write_winkx";

PyObject * pypylith3d_write_winkx(PyObject *, PyObject *args)
{
  PyObject* py_winkxdef;
  PyObject* py_iwinkxdef;
  PyObject* py_iwinkxid;
  int nwinkxe;
  int kw;
  int idout;
  char* ofile;

  int ok = PyArg_ParseTuple(args, "OOOiiis:write_winkx",
			    &py_winkxdef,
			    &py_iwinkxdef,
			    &py_iwinkxid,
			    &nwinkxe,
			    &kw,
			    &idout,
			    &ofile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* winkxdef = (double*) PyCObject_AsVoidPtr(py_winkxdef);
  int* iwinkxdef = (int*) PyCObject_AsVoidPtr(py_iwinkxdef);
  int* iwinkxid = (int*) PyCObject_AsVoidPtr(py_iwinkxid);

  write_winkx_f(winkxdef,
		iwinkxdef,
		iwinkxid,
		&nwinkxe,
		&kw,
		&idout,
		ofile,
		&errorcode,
		errorstring,
		strlen(ofile),
		sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "nwinkxe:" << nwinkxe
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// Scan boundary conditions

char pypylith3d_scan_bc__doc__[] = "";
char pypylith3d_scan_bc__name__[] = "scan_bc";

PyObject * pypylith3d_scan_bc(PyObject *, PyObject *args)
{
  int kr;
  char* displacement_units;
  char* velocity_units;
  char* force_units;
  char* bcInputFile;

  int ok = PyArg_ParseTuple(args, (char *) "issss:scan_bc",
			    &kr,
			    &displacement_units,
			    &velocity_units,
			    &force_units,
			    &bcInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int numbc = 0;

  scan_bc_f(&numbc,
	    &kr,
	    displacement_units,
	    velocity_units,
	    force_units,
	    bcInputFile,
	    &errorcode,
	    errorstring,
	    strlen(displacement_units),
	    strlen(velocity_units),
	    strlen(force_units),
	    strlen(bcInputFile),
	    sizeof(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numbc:" << numbc
    << journal::endl;

  // return
  return Py_BuildValue((char *) "i", numbc);
}


// Scan connectivities

char pypylith3d_scan_connect__doc__[] = "";
char pypylith3d_scan_connect__name__[] = "scan_connect";

PyObject * pypylith3d_scan_connect(PyObject *, PyObject *args)
{
  PyObject* py_numberElementNodesBase;
  PyObject* py_infmatmod;
  PyObject* py_infmat;
  PyObject* py_volumeElementFamilyList;
  int maxvfamilies;
  int numat;
  int kr;
  char* connectivityInputFile;

  int ok = PyArg_ParseTuple(args, (char *) "OOOOiiis:scan_connect",
			    &py_numberElementNodesBase,
			    &py_infmatmod,
                            &py_infmat,
			    &py_volumeElementFamilyList,
			    &maxvfamilies,
			    &numat,
			    &kr,
			    &connectivityInputFile);

  if (!ok) {
    return 0;
  }

  int* neni = (int*) PyCObject_AsVoidPtr(py_numberElementNodesBase);
  int* infmatmod = (int*) PyCObject_AsVoidPtr(py_infmatmod);
  int* infmat = (int*) PyCObject_AsVoidPtr(py_infmat);
  int* ivflist = (int*) PyCObject_AsVoidPtr(py_volumeElementFamilyList);
  int numelv = 0;
  int nvfamilies = 0;
  int ietypev = 0;
  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];

  scan_connect_f(neni,
		 infmatmod,
		 infmat,
		 ivflist,
		 &maxvfamilies,
		 &numat,
		 &numelv,
		 &nvfamilies,
		 &ietypev,
		 &kr,
		 connectivityInputFile,
		 &errorcode,
		 errorstring,
		 strlen(connectivityInputFile),
		 sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "nvfamilies:" << nvfamilies
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue((char *) "iii", numelv,
		       nvfamilies,
		       ietypev);
}


// Scan coordinates

char pypylith3d_scan_coords__doc__[] = "";
char pypylith3d_scan_coords__name__[] = "scan_coords";

PyObject * pypylith3d_scan_coords(PyObject *, PyObject *args)
{
  int kr;
  char *coord_units;
  char *coordinateInputFile;

  int ok = PyArg_ParseTuple(args, (char *) "iss:scan_coords",
			    &kr,
			    &coord_units,
			    &coordinateInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int numnp = 0;

  scan_coords_f(&numnp,
		&kr,
		coord_units,
		coordinateInputFile,
		&errorcode,
		errorstring,
		strlen(coord_units),
		strlen(coordinateInputFile),
		sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numnp:" << numnp
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue((char *) "i", numnp);
}


// Scan differential forces

char pypylith3d_scan_diff__doc__[] = "";
char pypylith3d_scan_diff__name__[] = "scan_diff";

PyObject * pypylith3d_scan_diff(PyObject *, PyObject *args)
{
  int numslp;
  int kr;
  char* difile;

  int ok = PyArg_ParseTuple(args, (char *) "iis:scan_diff",
			    &numslp,
			    &kr,
			    &difile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int numdif = 0;

  scan_diff_f(&numslp,
	      &numdif,
	      &kr,
	      difile,
	      &errorcode,
	      errorstring,
	      strlen(difile),
	      sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numdif:" << numdif
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue((char *) "i", numdif);
}


// Scan time steps at which full output is desired

char pypylith3d_scan_fuldat__doc__[] = "";
char pypylith3d_scan_fuldat__name__[] = "scan_fuldat";

PyObject * pypylith3d_scan_fuldat(PyObject *, PyObject *args)
{
  int icode;
  int lastep;
  int kr;
  char* fofile;

  int ok = PyArg_ParseTuple(args, (char *) "iiis:scan_fuldat",
			    &icode,
			    &lastep,
			    &kr,
			    &fofile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int icontr = 0;

  scan_fuldat_f(&icode,
		&lastep,
		&icontr,
		&kr,
		fofile,
		&errorcode,
		errorstring,
		strlen(fofile),
		sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "icontr:" << icontr
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue((char *) "i", icontr);
}


// Scan load histories

char pypylith3d_scan_hist__doc__[] = "";
char pypylith3d_scan_hist__name__[] = "scan_hist";

PyObject * pypylith3d_scan_hist(PyObject *, PyObject *args)
{
  int kr;
  char* hfile;

  int ok = PyArg_ParseTuple(args, (char *) "is:scan_hist",
			    &kr,
			    &hfile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int nhist = 0;

  scan_hist_f(&nhist,
	      &kr,
	      hfile,
	      &errorcode,
	      errorstring,
	      strlen(hfile),
	      sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "nhist:" << nhist
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue((char *) "i", nhist);
}


// Scan element prestresses

  // char pypylith3d_scan_prestr__doc__[] = "";
  // char pypylith3d_scan_prestr__name__[] = "scan_prestr";

  // PyObject * pypylith3d_scan_prestr(PyObject *, PyObject *args)
  // {
  //   int numberStressComponents;
  //   int numberPrestressGaussPoints;
  //   int numberElements;
  //   int ipstrs;
  //   int kr;
  //   char* prestressInputFile;

  //   int ok = PyArg_ParseTuple(args, "iiiiis:scan_prestr",
  // 			    &numberStressComponents,
  // 			    &numberPrestressGaussPoints,
  // 			    &numberElements,
  // 			    &ipstrs,
  // 			    &kr,
  // 			    &prestressInputFile);

  //   if (!ok) {
  //     return 0;
  //   }

  //   int errorcode = 0;
  //   const int maxsize = 4096;
  //   char errorstring[maxsize];
  //   int numberPrestressEntries = 0;

  //   scan_prestr_f(&numberStressComponents,
  // 		&numberPrestressGaussPoints,
  // 		&numberPrestressEntries,
  // 		&numberElements,
  // 		&ipstrs,
  // 		&kr,
  // 		&errorcode,
  // 		prestressInputFile,strlen(prestressInputFile));
    
  //   if(0 != exceptionhandler(errorcode, errorstring)) {
  //     return 0;
  //   }

  //   journal::debug_t debug("pylith3d");
  //   debug
  //     << journal::at(__HERE__)
  //     << "numberPrestressEntries:" << numberPrestressEntries
  //     << journal::endl;

  // return
  //   Py_INCREF(Py_None);
  //   return Py_BuildValue("i",numberPrestressEntries);
  // }


// Scan local coordinate rotations

char pypylith3d_scan_skew__doc__[] = "";
char pypylith3d_scan_skew__name__[] = "scan_skew";

PyObject * pypylith3d_scan_skew(PyObject *, PyObject *args)
{
  int kr;
  char* rotation_units;
  char* skfile;

  int ok = PyArg_ParseTuple(args, (char *) "iss:scan_skew",
			    &kr,
			    &rotation_units,
			    &skfile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int numrot = 0;

  scan_skew_f(&numrot,
	      &kr,
	      rotation_units,
	      skfile,
	      &errorcode,
	      errorstring,
	      strlen(rotation_units),
	      strlen(skfile),
	      sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numrot:" << numrot
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue((char *) "i", numrot);
}


// Scan slippery node entries

char pypylith3d_scan_slip__doc__[] = "";
char pypylith3d_scan_slip__name__[] = "scan_slip";

PyObject * pypylith3d_scan_slip(PyObject *, PyObject *args)
{
  int kr;
  char* slfile;

  int ok = PyArg_ParseTuple(args, (char *) "is:scan_slip",
			    &kr,
			    &slfile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int numslp = 0;

  scan_slip_f(&numslp,
	      &kr,
	      slfile,
	      &errorcode,
	      errorstring,
	      strlen(slfile),
	      sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numslp:" << numslp
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue((char *) "i", numslp);
}


// Scan split node entries

char pypylith3d_scan_split__doc__[] = "";
char pypylith3d_scan_split__name__[] = "scan_split";

PyObject * pypylith3d_scan_split(PyObject *, PyObject *args)
{
  int kr;
  char* spfile;

  int ok = PyArg_ParseTuple(args, (char *) "is:scan_split",
			    &kr,
			    &spfile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int numfn = 0;

  scan_split_f(&numfn,
	       &kr,
	       spfile,
	       &errorcode,
	       errorstring,
	       strlen(spfile),
	       sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numfn:" << numfn
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue((char *) "i", numfn);
}


// Scan time step data

char pypylith3d_scan_timdat__doc__[] = "";
char pypylith3d_scan_timdat__name__[] = "scan_timdat";

PyObject * pypylith3d_scan_timdat(PyObject *, PyObject *args)
{
  int kr;
  char* time_units;
  char* timeStepInputFile;

  int ok = PyArg_ParseTuple(args, (char *) "iss:scan_timdat",
			    &kr,
			    &time_units,
			    &timeStepInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int nintg = 0;
  int lastep = 0;

  scan_timdat_f(&lastep,
		&nintg,
		&kr,
		time_units,
		timeStepInputFile,
		&errorcode,
		errorstring,
		strlen(time_units),
		strlen(timeStepInputFile),
		sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "lastep:" << lastep
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue((char *) "ii", nintg,
		       lastep);
}



// Scan traction BC

char pypylith3d_scan_tractions__doc__[] = "";
char pypylith3d_scan_tractions__name__[] = "scan_tractions";

PyObject * pypylith3d_scan_tractions(PyObject *, PyObject *args)
{
  int nsnodesmax;
  int kr;
  char *traction_units;
  char *tractionInputFile;

  int ok = PyArg_ParseTuple(args, (char *) "iiss:scan_tractions",
                            &nsnodesmax,
			    &kr,
			    &traction_units,
			    &tractionInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int numtractions = 0;

  scan_tractions_f(&numtractions,
		   &nsnodesmax,
		   &kr,
		   traction_units,
		   tractionInputFile,
		   &errorcode,
		   errorstring,
		   strlen(traction_units),
		   strlen(tractionInputFile),
		   sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numtractions:" << numtractions
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue((char *) "i", numtractions);
}

// Read winkler BC

char pypylith3d_scan_wink__doc__[] = "";
char pypylith3d_scan_wink__name__[] = "scan_wink";

PyObject * pypylith3d_scan_wink(PyObject *, PyObject *args)
{
  int kr;
  char* wfile;

  int ok = PyArg_ParseTuple(args, (char *) "is:scan_wink",
			    &kr,
			    &wfile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int nwinke = 0;
  int nwink = 0;

  scan_wink_f(&nwinke,
	      &nwink,
	      &kr,
	      wfile,
	      &errorcode,
	      errorstring,
	      strlen(wfile),
	      sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "nwink:" << nwink
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue((char *) "ii", nwinke,
	               nwink);
}


// Read winkler BC for slippery nodes

char pypylith3d_scan_winkx__doc__[] = "";
char pypylith3d_scan_winkx__name__[] = "scan_winkx";

PyObject * pypylith3d_scan_winkx(PyObject *, PyObject *args)
{
  int numslp;
  int kr;
  char* wxfile;

  int ok = PyArg_ParseTuple(args, (char *) "iis:scan_winkx",
			    &numslp,
			    &kr,
			    &wxfile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int nwinkxe = 0;
  int nwinkx = 0;

  scan_winkx_f(&numslp,
	       &nwinkxe,
	       &nwinkx,
	       &kr,
	       wxfile,
	       &errorcode,
	       errorstring,
	       strlen(wxfile),
	       sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "nwinkx:" << nwinkx
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue((char *) "ii", nwinkxe,
		       nwinkx);
}
    

// Initialize material model info

char pypylith3d_matmod_def__doc__[] = "";
char pypylith3d_matmod_def__name__[] = "matmod_def";

PyObject * pypylith3d_matmod_def(PyObject *, PyObject *args)
{
  PyObject* py_infmatmod;

  int ok = PyArg_ParseTuple(args, "O:matmod_def",
			    &py_infmatmod);

  if (!ok) {
    return 0;
  }

  int* infmatmod = (int*) PyCObject_AsVoidPtr(py_infmatmod);

  matmod_def_f(infmatmod);
    
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "infmatmod:" << infmatmod
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}

// Precompute shape function info

char pypylith3d_preshape__doc__[] = "";
char pypylith3d_preshape__name__[] = "preshape";

PyObject * pypylith3d_preshape(PyObject *, PyObject *args)
{
  PyObject* py_sh;
  PyObject* py_shj;
  PyObject* py_gauss;
  int intord;
  int elementType;
  int nen;
  int ngauss;

  int ok = PyArg_ParseTuple(args, "OOOiiii:preshape",
			    &py_sh,
			    &py_shj,
			    &py_gauss,
			    &intord,
			    &elementType,
			    &nen,
			    &ngauss);

  if (!ok) {
    return 0;
  }

  double* sh = (double*) PyCObject_AsVoidPtr(py_sh);
  double* shj = (double*) PyCObject_AsVoidPtr(py_shj);
  double* gauss = (double*) PyCObject_AsVoidPtr(py_gauss);
  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];

  preshape_f(sh,
	     shj,
	     gauss,
	     &intord,
	     &elementType,
	     &nen,
	     &ngauss,
	     &errorcode,
	     errorstring,
	     sizeof(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }
    
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "intord:" << intord
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}

// Precompute shape function info for element faces

char pypylith3d_preshape2d__doc__[] = "";
char pypylith3d_preshape2d__name__[] = "preshape2d";

PyObject * pypylith3d_preshape2d(PyObject *, PyObject *args)
{
  PyObject* py_sh2d;
  PyObject* py_gauss2d;
  int intord;
  int elementType;
  int numberSurfaceElementNodes;
  int numberSurfaceElementGaussPoints;

  int ok = PyArg_ParseTuple(args, "OOiiii:preshape2d",
			    &py_sh2d,
			    &py_gauss2d,
			    &intord,
			    &elementType,
			    &numberSurfaceElementNodes,
			    &numberSurfaceElementGaussPoints);

  if (!ok) {
    return 0;
  }

  double* sh2d = (double*) PyCObject_AsVoidPtr(py_sh2d);
  double* gauss2d = (double*) PyCObject_AsVoidPtr(py_gauss2d);
  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];

  preshape2d_f(sh2d,
	       gauss2d,
	       &intord,
	       &elementType,
	       &numberSurfaceElementNodes,
	       &numberSurfaceElementGaussPoints,
	       &errorcode,
	       errorstring,
	       sizeof(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }
    
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "intord:" << intord
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}

// try_binio

char pypylith3d_try_binio__doc__[] = "";
char pypylith3d_try_binio__name__[] = "try_binio";

PyObject * pypylith3d_try_binio(PyObject *, PyObject *args)
{
  int unit;

  int ok = PyArg_ParseTuple(args, "i:try_binio",
			    &unit);
  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];

  try_binio_f(&unit,
	      &errorcode,
	      errorstring,
	      sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  // return
  Py_INCREF(Py_None);
  return Py_None;
}
    
// end of file
