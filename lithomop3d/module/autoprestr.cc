// -*- C++ -*-
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//                               Charles A. Williams
//                        Rensselaer Polytechnic Institute
//                        (C) 2005 All Rights Reserved
// 
//  All worldwide rights reserved.  A license to use, copy, modify and
//  distribute this software for non-commercial research purposes only
//  is hereby granted, provided that this copyright notice and
//  accompanying disclaimer is not modified or removed from the software.
//
//  DISCLAIMER:  The software is distributed "AS IS" without any express
//  or implied warranty, including but not limited to, any implied
//  warranties of merchantability or fitness for a particular purpose
//  or any warranty of non-infringement of any current or pending patent
//  rights.  The authors of the software make no representations about
//  the suitability of this software for any particular purpose.  The
//  entire risk as to the quality and performance of the software is with
//  the user.  Should the software prove defective, the user assumes the
//  cost of all necessary servicing, repair or correction.  In
//  particular, neither Rensselaer Polytechnic Institute, nor the authors
//  of the software are liable for any indirect, special, consequential,
//  or incidental damages related to the software, to the maximum extent
//  the law permits.
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 

#include <portinfo>
#include "journal/debug.h"

#include <Python.h>

#include "autoprestr.h"
#include "exceptionhandler.h"
#include "lithomop3d_externs.h"
#include <stdio.h>
#include <string.h>


// Autocompute gravitational prestresses.

char pylithomop3d_autoprestr__doc__[] = "";
char pylithomop3d_autoprestr__name__[] = "autoprestr";

PyObject * pylithomop3d_autoprestr(PyObject *, PyObject *args)
{
  PyObject* pyA;                             // Sparse matrix array
  PyObject* pyPointerToBextern;               // Force vectors
  PyObject* pyPointerToBtraction;
  PyObject* pyPointerToBgravity;
  PyObject* pyPointerToBconcForce;
  PyObject* pyPointerToBintern;
  PyObject* pyPointerToBresid;
  PyObject* pyPointerToBwink;
  PyObject* pyPointerToBwinkx;
  PyObject* pyPointerToDispVec;
  PyObject* pyPointerToDprev;
  PyObject* pyPointerToListArrayNforce;
  PyObject* pyPointerToListArrayGrav;
  PyObject* pyPointerToX;                     // Global arrays
  PyObject* pyPointerToD;
  PyObject* pyPointerToDeld;
  PyObject* pyPointerToDcur;
  PyObject* pyPointerToId;
  PyObject* pyPointerToIwink;
  PyObject* pyPointerToWink;
  PyObject* pyPointerToListArrayNsysdat;
  PyObject* pyPointerToListArrayIddmat;
  PyObject* pyPointerToIbond;                 // Boundary condition arrays
  PyObject* pyPointerToBond;
  PyObject* pyPointerToDx;                    // Slippery node arrays
  PyObject* pyPointerToDeldx;
  PyObject* pyPointerToDxcur;
  PyObject* pyPointerToDiforc;
  PyObject* pyPointerToIdx;
  PyObject* pyPointerToIwinkx;
  PyObject* pyPointerToWinkx;
  PyObject* pyPointerToIdslp;
  PyObject* pyPointerToIpslp;
  PyObject* pyPointerToIdhist;
  PyObject* pyPointerToFault;                 // Split node arrays
  PyObject* pyPointerToNfault;
  PyObject* pyPointerToDfault;
  PyObject* pyPointerToTfault;
  PyObject* pyPointerToS;                     // Local stiffness matrix arrays
  PyObject* pyPointerToStemp;
  PyObject* pyPointerToState;                 // Element arrays
  PyObject* pyPointerToDstate;
  PyObject* pyPointerToState0;
  PyObject* pyPointerToDmat;
  PyObject* pyPointerToIens;
  PyObject* pyPointerToLm;
  PyObject* pyPointerToLmx;
  PyObject* pyPointerToLmf;
  PyObject* pyPointerToIvfamily;
  PyObject* pyPointerToListArrayNpar;
  PyObject* pyPointerToIelno;                 // Traction BC arrays
  PyObject* pyPointerToIside;
  PyObject* pyPointerToIhistry;
  PyObject* pyPointerToPres;
  PyObject* pyPointerToPdir;
  PyObject* pyPointerToListArrayPropertyList; // Material property arrays
  PyObject* pyPointerToMaterialModelInfo;
  PyObject* pyPointerToGauss;                 // Element type arrays
  PyObject* pyPointerToSh;
  PyObject* pyPointerToShj;
  PyObject* pyPointerToListArrayElementTypeInfo;
  PyObject* pyPointerToHistry;                // Time information
  PyObject* pyPointerToListArrayRtimdat;
  PyObject* pyPointerToListArrayNtimdat;
  PyObject* pyPointerToListArrayNvisdat;
  PyObject* pyPointerToMaxstp;
  PyObject* pyPointerToDelt;
  PyObject* pyPointerToAlfa;
  PyObject* pyPointerToMaxit;
  PyObject* pyPointerToNtdinit;
  PyObject* pyPointerToLgdef;
  PyObject* pyPointerToUtol;
  PyObject* pyPointerToFtol;
  PyObject* pyPointerToEtol;
  PyObject* pyPointerToItmax;
  PyObject* pyPointerToListArrayRgiter;       // Iterative solution information
  PyObject* pyPointerToSkew;                  // Skew rotation information
  PyObject* pyPointerToListArrayNcodat;       // Input/output information
  PyObject* pyPointerToListArrayNunits;
  PyObject* pyPointerToListArrayNprint;
  PyObject* pyPointerToIstatout;
  PyObject* pyPointerToNstatout;
  char* asciiOutputFile;                      // Output file names
  char* plotOutputFile;
  char* ucdOutputRoot;
  int autoprestrStage, iterateEvent;          // PETSc logging

  int ok = PyArg_ParseTuple(args, "OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOsssii:autoprestr",
			    &pyA,                             // Sparse matrix arrays
			    &pyPointerToBextern,               // Force vectors
			    &pyPointerToBtraction,
			    &pyPointerToBgravity,
			    &pyPointerToBconcForce,
			    &pyPointerToBintern,
			    &pyPointerToBresid,
			    &pyPointerToBwink,
			    &pyPointerToBwinkx,
			    &pyPointerToDispVec,
			    &pyPointerToDprev,
			    &pyPointerToListArrayNforce,
			    &pyPointerToListArrayGrav,
			    &pyPointerToX,                     // Global arrays
			    &pyPointerToD,
			    &pyPointerToDeld,
			    &pyPointerToDcur,
			    &pyPointerToId,
			    &pyPointerToIwink,
			    &pyPointerToWink,
			    &pyPointerToListArrayNsysdat,
			    &pyPointerToListArrayIddmat,
			    &pyPointerToIbond,                 // Boundary condition arrays
			    &pyPointerToBond,
			    &pyPointerToDx,                    // Slippery node arrays
			    &pyPointerToDeldx,
			    &pyPointerToDxcur,
			    &pyPointerToDiforc,
			    &pyPointerToIdx,
			    &pyPointerToIwinkx,
			    &pyPointerToWinkx,
			    &pyPointerToIdslp,
			    &pyPointerToIpslp,
			    &pyPointerToIdhist,
			    &pyPointerToFault,                 // Split node arrays
			    &pyPointerToNfault,
			    &pyPointerToDfault,
			    &pyPointerToTfault,
			    &pyPointerToS,                     // Local stiffness matrix arrays
			    &pyPointerToStemp,
			    &pyPointerToState,                 // Element arrays
			    &pyPointerToDstate,
			    &pyPointerToState0,
			    &pyPointerToDmat,
			    &pyPointerToIens,
			    &pyPointerToLm,
			    &pyPointerToLmx,
			    &pyPointerToLmf,
			    &pyPointerToIvfamily,
			    &pyPointerToListArrayNpar,
			    &pyPointerToIelno,                 // Traction BC arrays
			    &pyPointerToIside,
			    &pyPointerToIhistry,
			    &pyPointerToPres,
			    &pyPointerToPdir,
			    &pyPointerToListArrayPropertyList, // Material property arrays
			    &pyPointerToMaterialModelInfo,
			    &pyPointerToGauss,                 // Element type arrays
			    &pyPointerToSh,
			    &pyPointerToShj,
			    &pyPointerToListArrayElementTypeInfo,
			    &pyPointerToHistry,                // Time information
			    &pyPointerToListArrayRtimdat,
			    &pyPointerToListArrayNtimdat,
			    &pyPointerToListArrayNvisdat,
			    &pyPointerToMaxstp,
			    &pyPointerToDelt,
			    &pyPointerToAlfa,
			    &pyPointerToMaxit,
			    &pyPointerToNtdinit,
			    &pyPointerToLgdef,
			    &pyPointerToUtol,
			    &pyPointerToFtol,
			    &pyPointerToEtol,
			    &pyPointerToItmax,
			    &pyPointerToListArrayRgiter,       // Iterative solution information
			    &pyPointerToSkew,                  // Skew rotation information
			    &pyPointerToListArrayNcodat,       // Input/output information
			    &pyPointerToListArrayNunits,
			    &pyPointerToListArrayNprint,
			    &pyPointerToIstatout,
			    &pyPointerToNstatout,
			    &asciiOutputFile,                  // Output file names
			    &plotOutputFile,
			    &ucdOutputRoot,
			    &autoprestrStage,                  // PETSc logging
                            &iterateEvent);


  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  Mat      A = (Mat) PyCObject_AsVoidPtr(pyA);
  double*  pointerToBextern = (double*) PyCObject_AsVoidPtr(pyPointerToBextern);
  double*  pointerToBtraction = (double*) PyCObject_AsVoidPtr(pyPointerToBtraction);
  double*  pointerToBgravity = (double*) PyCObject_AsVoidPtr(pyPointerToBgravity);
  double*  pointerToBconcForce = (double*) PyCObject_AsVoidPtr(pyPointerToBconcForce);
  double*  pointerToBintern = (double*) PyCObject_AsVoidPtr(pyPointerToBintern);
  double*  pointerToBresid = (double*) PyCObject_AsVoidPtr(pyPointerToBresid);
  double*  pointerToBwink = (double*) PyCObject_AsVoidPtr(pyPointerToBwink);
  double*  pointerToBwinkx = (double*) PyCObject_AsVoidPtr(pyPointerToBwinkx);
  double*  pointerToDispVec = (double*) PyCObject_AsVoidPtr(pyPointerToDispVec);
  double*  pointerToDprev = (double*) PyCObject_AsVoidPtr(pyPointerToDprev);
  int*  pointerToListArrayNforce = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNforce);
  double*  pointerToListArrayGrav = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayGrav);
  double*  pointerToX = (double*) PyCObject_AsVoidPtr(pyPointerToX);
  double*  pointerToD = (double*) PyCObject_AsVoidPtr(pyPointerToD);
  double*  pointerToDeld = (double*) PyCObject_AsVoidPtr(pyPointerToDeld);
  double*  pointerToDcur = (double*) PyCObject_AsVoidPtr(pyPointerToDcur);
  int*  pointerToId = (int*) PyCObject_AsVoidPtr(pyPointerToId);
  int*  pointerToIwink = (int*) PyCObject_AsVoidPtr(pyPointerToIwink);
  double*  pointerToWink = (double*) PyCObject_AsVoidPtr(pyPointerToWink);
  int*  pointerToListArrayNsysdat = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNsysdat);
  int*  pointerToListArrayIddmat = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayIddmat);
  int*  pointerToIbond = (int*) PyCObject_AsVoidPtr(pyPointerToIbond);
  double*  pointerToBond = (double*) PyCObject_AsVoidPtr(pyPointerToBond);
  double*  pointerToDx = (double*) PyCObject_AsVoidPtr(pyPointerToDx);
  double*  pointerToDeldx = (double*) PyCObject_AsVoidPtr(pyPointerToDeldx);
  double*  pointerToDxcur = (double*) PyCObject_AsVoidPtr(pyPointerToDxcur);
  double*  pointerToDiforc = (double*) PyCObject_AsVoidPtr(pyPointerToDiforc);
  int*  pointerToIdx = (int*) PyCObject_AsVoidPtr(pyPointerToIdx);
  int*  pointerToIwinkx = (int*) PyCObject_AsVoidPtr(pyPointerToIwinkx);
  double*  pointerToWinkx = (double*) PyCObject_AsVoidPtr(pyPointerToWinkx);
  int*  pointerToIdslp = (int*) PyCObject_AsVoidPtr(pyPointerToIdslp);
  int*  pointerToIpslp = (int*) PyCObject_AsVoidPtr(pyPointerToIpslp);
  int*  pointerToIdhist = (int*) PyCObject_AsVoidPtr(pyPointerToIdhist);
  double*  pointerToFault = (double*) PyCObject_AsVoidPtr(pyPointerToFault);
  int*  pointerToNfault = (int*) PyCObject_AsVoidPtr(pyPointerToNfault);
  double*  pointerToDfault = (double*) PyCObject_AsVoidPtr(pyPointerToDfault);
  double*  pointerToTfault = (double*) PyCObject_AsVoidPtr(pyPointerToTfault);
  double*  pointerToS = (double*) PyCObject_AsVoidPtr(pyPointerToS);
  double*  pointerToStemp = (double*) PyCObject_AsVoidPtr(pyPointerToStemp);
  double*  pointerToState = (double*) PyCObject_AsVoidPtr(pyPointerToState);
  double*  pointerToDstate = (double*) PyCObject_AsVoidPtr(pyPointerToDstate);
  double*  pointerToState0 = (double*) PyCObject_AsVoidPtr(pyPointerToState0);
  double*  pointerToDmat = (double*) PyCObject_AsVoidPtr(pyPointerToDmat);
  int*  pointerToIens = (int*) PyCObject_AsVoidPtr(pyPointerToIens);
  int*  pointerToLm = (int*) PyCObject_AsVoidPtr(pyPointerToLm);
  int*  pointerToLmx = (int*) PyCObject_AsVoidPtr(pyPointerToLmx);
  int*  pointerToLmf = (int*) PyCObject_AsVoidPtr(pyPointerToLmf);
  int*  pointerToIvfamily = (int*) PyCObject_AsVoidPtr(pyPointerToIvfamily);
  int*  pointerToListArrayNpar = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNpar);
  int*  pointerToIelno = (int*) PyCObject_AsVoidPtr(pyPointerToIelno);
  int*  pointerToIside = (int*) PyCObject_AsVoidPtr(pyPointerToIside);
  int*  pointerToIhistry = (int*) PyCObject_AsVoidPtr(pyPointerToIhistry);
  double*  pointerToPres = (double*) PyCObject_AsVoidPtr(pyPointerToPres);
  double*  pointerToPdir = (double*) PyCObject_AsVoidPtr(pyPointerToPdir);
  double*  pointerToListArrayPropertyList = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayPropertyList);
  int*  pointerToMaterialModelInfo = (int*) PyCObject_AsVoidPtr(pyPointerToMaterialModelInfo);
  double*  pointerToGauss = (double*) PyCObject_AsVoidPtr(pyPointerToGauss);
  double*  pointerToSh = (double*) PyCObject_AsVoidPtr(pyPointerToSh);
  double*  pointerToShj = (double*) PyCObject_AsVoidPtr(pyPointerToShj);
  int*  pointerToListArrayElementTypeInfo = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayElementTypeInfo);
  double*  pointerToHistry = (double*) PyCObject_AsVoidPtr(pyPointerToHistry);
  double*  pointerToListArrayRtimdat = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayRtimdat);
  int*  pointerToListArrayNtimdat = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNtimdat);
  int*  pointerToListArrayNvisdat = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNvisdat);
  int*  pointerToMaxstp = (int*) PyCObject_AsVoidPtr(pyPointerToMaxstp);
  double*  pointerToDelt = (double*) PyCObject_AsVoidPtr(pyPointerToDelt);
  double*  pointerToAlfa = (double*) PyCObject_AsVoidPtr(pyPointerToAlfa);
  int*  pointerToMaxit = (int*) PyCObject_AsVoidPtr(pyPointerToMaxit);
  int*  pointerToNtdinit = (int*) PyCObject_AsVoidPtr(pyPointerToNtdinit);
  int*  pointerToLgdef = (int*) PyCObject_AsVoidPtr(pyPointerToLgdef);
  double*  pointerToUtol = (double*) PyCObject_AsVoidPtr(pyPointerToUtol);
  double*  pointerToFtol = (double*) PyCObject_AsVoidPtr(pyPointerToFtol);
  double*  pointerToEtol = (double*) PyCObject_AsVoidPtr(pyPointerToEtol);
  int*  pointerToItmax = (int*) PyCObject_AsVoidPtr(pyPointerToItmax);
  double*  pointerToListArrayRgiter = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayRgiter);
  double*  pointerToSkew = (double*) PyCObject_AsVoidPtr(pyPointerToSkew);
  int*  pointerToListArrayNcodat = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNcodat);
  int*  pointerToListArrayNunits = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNunits);
  int*  pointerToListArrayNprint = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNprint);
  int*  pointerToIstatout = (int*) PyCObject_AsVoidPtr(pyPointerToIstatout);
  int*  pointerToNstatout = (int*) PyCObject_AsVoidPtr(pyPointerToNstatout);

  autoprestr_f(&A,                                // Sparse matrix arrays
	       pointerToBextern,                  // Force vectors
	       pointerToBtraction,
	       pointerToBgravity,
	       pointerToBconcForce,
	       pointerToBintern,
	       pointerToBresid,
	       pointerToBwink,
	       pointerToBwinkx,
	       pointerToDispVec,
	       pointerToDprev,
	       pointerToListArrayNforce,
	       pointerToListArrayGrav,
	       pointerToX,                        // Global arrays
	       pointerToD,
	       pointerToDeld,
	       pointerToDcur,
	       pointerToId,
	       pointerToIwink,
	       pointerToWink,
	       pointerToListArrayNsysdat,
	       pointerToListArrayIddmat,
	       pointerToIbond,                    // Boundary condition arrays
	       pointerToBond,
	       pointerToDx,                       // Slippery node arrays
	       pointerToDeldx,
	       pointerToDxcur,
	       pointerToDiforc,
	       pointerToIdx,
	       pointerToIwinkx,
	       pointerToWinkx,
	       pointerToIdslp,
	       pointerToIpslp,
	       pointerToIdhist,
	       pointerToFault,                    // Split node arrays
	       pointerToNfault,
	       pointerToDfault,
	       pointerToTfault,
	       pointerToS,                        // Local stiffness matrix arrays
	       pointerToStemp,
	       pointerToState,                    // Element arrays
	       pointerToDstate,
	       pointerToState0,
	       pointerToDmat,
	       pointerToIens,
	       pointerToLm,
	       pointerToLmx,
	       pointerToLmf,
	       pointerToIvfamily,
	       pointerToListArrayNpar,
	       pointerToIelno,                    // Traction BC arrays
	       pointerToIside,
	       pointerToIhistry,
	       pointerToPres,
	       pointerToPdir,
	       pointerToListArrayPropertyList,    // Material property arrays
	       pointerToMaterialModelInfo,
	       pointerToGauss,                    // Element type arrays
	       pointerToSh,
	       pointerToShj,
	       pointerToListArrayElementTypeInfo,
	       pointerToHistry,                   // Time information
	       pointerToListArrayRtimdat,
	       pointerToListArrayNtimdat,
	       pointerToListArrayNvisdat,
	       pointerToMaxstp,
	       pointerToDelt,
	       pointerToAlfa,
	       pointerToMaxit,
	       pointerToNtdinit,
	       pointerToLgdef,
	       pointerToUtol,
	       pointerToFtol,
	       pointerToEtol,
	       pointerToItmax,
	       pointerToListArrayRgiter,          // Iterative solution information
	       pointerToSkew,                     // Skew rotation information
	       pointerToListArrayNcodat,          // Input/output information
	       pointerToListArrayNunits,
	       pointerToListArrayNprint,
	       pointerToIstatout,
	       pointerToNstatout,
	       asciiOutputFile,                   // Output file names
	       plotOutputFile,
	       ucdOutputRoot,
               &autoprestrStage,
               &iterateEvent,
	       &errorcode,                        // Error codes
	       errorstring,
	       strlen(asciiOutputFile),           // String lengths
	       strlen(plotOutputFile),
	       strlen(ucdOutputRoot),
	       strlen(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }
  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "asciiOutputFile:" << asciiOutputFile
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// version
// $Id: autoprestr.cc,v 1.9 2005/05/03 18:48:19 willic3 Exp $

// End of file
