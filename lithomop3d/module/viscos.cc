// -*- C++ -*-
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//                               Charles A. Williams
//                        Rensselaer Polytechnic Institute
//                        (C) 2004 All Rights Reserved
// 
//  Copyright 2004 Rensselaer Polytechnic Institute.
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

#include "viscos.h"
#include "exceptionhandler.h"
#include "lithomop3d_externs.h"
#include <stdio.h>
#include <string.h>


// Perform the time-dependent solution.

char pylithomop3d_viscos__doc__[] = "";
char pylithomop3d_viscos__name__[] = "viscos";

PyObject * pylithomop3d_viscos(PyObject *, PyObject *args)
{
  PyObject* pyPointerToAlnz;                  // Sparse matrix arrays
  PyObject* pyPointerToPcg;
  PyObject* pyPointerToZcg;
  PyObject* pyPointerToDprev;
  PyObject* pyPointerToJa;
  PyObject* pyPointerToBextern;               // Force vectors
  PyObject* pyPointerToBtraction;
  PyObject* pyPointerToBgravity;
  PyObject* pyPointerToBconcForce;
  // PyObject* pyPointerToBprestress;
  PyObject* pyPointerToBintern;
  PyObject* pyPointerToBresid;
  PyObject* pyPointerToBwork;
  PyObject* pyPointerToDispVec;
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
  PyObject* pyPointerToIen;
  PyObject* pyPointerToLm;
  PyObject* pyPointerToLmx;
  PyObject* pyPointerToLmf;
  PyObject* pyPointerToInfiel;
  PyObject* pyPointerToListArrayIddmat;
  PyObject* pyPointerToListArrayNpar;
  PyObject* pyPointerToIelno;                 // Traction BC arrays
  PyObject* pyPointerToIside;
  PyObject* pyPointerToIhistry;
  PyObject* pyPointerToPres;
  PyObject* pyPointerToPdir;
  PyObject* pyPointerToListArrayPropertyList; // Material property arrays
  PyObject* pyPointerToMhist;
  PyObject* pyPointerToMaterialInfo;
  PyObject* pyPointerToMaterialModelInfo;
  PyObject* pyPointerToMaterialModelStates;
  PyObject* pyPointerGauss;                   // Element type arrays
  PyObject* pyPointerToSh;
  PyObject* pyPointerToShj;
  PyObject* pyPointerToElementTypeInfo;
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
  PyObject* pyPointerToListArrayGcurr;
  PyObject* pyPointerToListArrayGi;
  PyObject* pyPointerToListArrayGprev;
  PyObject* pyPointerToListArrayGtol;
  PyObject* pyPointerToListArrayRmin;
  PyObject* pyPointerToListArrayRmult;
  PyObject* pyPointerToListArrayNsiter;
  PyObject* pyPointerToSkew;                  // Skew rotation information
  PyObject* pyPointerToIprint;                // Input/output information
  PyObject* pyPointerToListArrayNcodat;
  PyObject* pyPointerToListArrayNunits;
  PyObject* pyPointerToListArrayNprint;
  PyObject* pyPointerToIstatout;
  char* asciiOutputFile;                      // Output file names
  char* plotOutputFile;
  char* ucdOutputRoot;

  int ok = PyArg_ParseTuple(args, "OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOsss:viscos",
			    &pyPointerToAlnz,                  // Sparse matrix arrays
			    &pyPointerToPcg,
			    &pyPointerToZcg,
			    &pyPointerToDprev,
			    &pyPointerToJa,
			    &pyPointerToBextern,               // Force vectors
			    &pyPointerToBtraction,
			    &pyPointerToBgravity,
			    &pyPointerToBconcForce,
			    // &pyPointerToBprestress,
			    &pyPointerToBintern,
			    &pyPointerToBresid,
			    &pyPointerToBwork,
			    &pyPointerToDispVec,
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
			    &pyPointerToIen,
			    &pyPointerToLm,
			    &pyPointerToLmx,
			    &pyPointerToLmf,
			    &pyPointerToInfiel,
			    &pyPointerToListArrayIddmat,
			    &pyPointerToListArrayNpar,
			    &pyPointerToIelno,                 // Traction BC arrays
			    &pyPointerToIside,
			    &pyPointerToIhistry,
			    &pyPointerToPres,
			    &pyPointerToPdir,
			    &pyPointerToListArrayPropertyList, // Material property arrays
			    &pyPointerToMhist,
			    &pyPointerToMaterialInfo,
			    &pyPointerToMaterialModelInfo,
			    &pyPointerToMaterialModelStates,
			    &pyPointerGauss,                   // Element type arrays
			    &pyPointerToSh,
			    &pyPointerToShj,
			    &pyPointerToElementTypeInfo,
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
			    &pyPointerToListArrayGcurr,
			    &pyPointerToListArrayGi,
			    &pyPointerToListArrayGprev,
			    &pyPointerToListArrayGtol,
			    &pyPointerToListArrayRmin,
			    &pyPointerToListArrayRmult,
			    &pyPointerToListArrayNsiter,
			    &pyPointerToSkew,                  // Skew rotation information
			    &pyPointerToIprint,                // Input/output information
			    &pyPointerToListArrayNcodat,
			    &pyPointerToListArrayNunits,
			    &pyPointerToListArrayNprint,
			    &pyPointerToIstatout,
			    &asciiOutputFile,                  // Output file names
			    &plotOutputFile,
			    &ucdOutputRoot);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 1024;
  char errorstring[maxsize];
  double*  pointerToAlnz = (double*) PyCObject_AsVoidPtr(pyPointerToAlnz);
  double*  pointerToPcg = (double*) PyCObject_AsVoidPtr(pyPointerToPcg);
  double*  pointerToZcg = (double*) PyCObject_AsVoidPtr(pyPointerToZcg);
  double*  pointerToDprev = (double*) PyCObject_AsVoidPtr(pyPointerToDprev);
  int*  pointerToJa = (int*) PyCObject_AsVoidPtr(pyPointerToJa);
  double*  pointerToBextern = (double*) PyCObject_AsVoidPtr(pyPointerToBextern);
  double*  pointerToBtraction = (double*) PyCObject_AsVoidPtr(pyPointerToBtraction);
  double*  pointerToBgravity = (double*) PyCObject_AsVoidPtr(pyPointerToBgravity);
  double*  pointerToBconcForce = (double*) PyCObject_AsVoidPtr(pyPointerToBconcForce);
  // double*  pointerToBprestress = (double*) PyCObject_AsVoidPtr(pyPointerToBprestress);
  double*  pointerToBintern = (double*) PyCObject_AsVoidPtr(pyPointerToBintern);
  double*  pointerToBresid = (double*) PyCObject_AsVoidPtr(pyPointerToBresid);
  double*  pointerToBwork = (double*) PyCObject_AsVoidPtr(pyPointerToBwork);
  double*  pointerToDispVec = (double*) PyCObject_AsVoidPtr(pyPointerToDispVec);
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
  int*  pointerToIen = (int*) PyCObject_AsVoidPtr(pyPointerToIen);
  int*  pointerToLm = (int*) PyCObject_AsVoidPtr(pyPointerToLm);
  int*  pointerToLmx = (int*) PyCObject_AsVoidPtr(pyPointerToLmx);
  int*  pointerToLmf = (int*) PyCObject_AsVoidPtr(pyPointerToLmf);
  int*  pointerToInfiel = (int*) PyCObject_AsVoidPtr(pyPointerToInfiel);
  int*  pointerToListArrayIddmat = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayIddmat);
  int*  pointerToListArrayNpar = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNpar);
  int*  pointerToIelno = (int*) PyCObject_AsVoidPtr(pyPointerToIelno);
  int*  pointerToIside = (int*) PyCObject_AsVoidPtr(pyPointerToIside);
  int*  pointerToIhistry = (int*) PyCObject_AsVoidPtr(pyPointerToIhistry);
  double*  pointerToPres = (double*) PyCObject_AsVoidPtr(pyPointerToPres);
  double*  pointerToPdir = (double*) PyCObject_AsVoidPtr(pyPointerToPdir);
  double*  pointerToListArrayPropertyList = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayPropertyList);
  int*  pointerToMhist = (int*) PyCObject_AsVoidPtr(pyPointerToMhist);
  int*  pointerToMaterialInfo = (int*) PyCObject_AsVoidPtr(pyPointerToMaterialInfo);
  int*  pointerToMaterialModelInfo = (int*) PyCObject_AsVoidPtr(pyPointerToMaterialModelInfo);
  int*  pointerToMaterialModelStates = (int*) PyCObject_AsVoidPtr(pyPointerToMaterialModelStates);
  double*  pointerGauss = (double*) PyCObject_AsVoidPtr(pyPointerGauss);
  double*  pointerToSh = (double*) PyCObject_AsVoidPtr(pyPointerToSh);
  double*  pointerToShj = (double*) PyCObject_AsVoidPtr(pyPointerToShj);
  int*  pointerToElementTypeInfo = (int*) PyCObject_AsVoidPtr(pyPointerToElementTypeInfo);
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
  double*  pointerToListArrayGcurr = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayGcurr);
  double*  pointerToListArrayGi = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayGi);
  double*  pointerToListArrayGprev = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayGprev);
  double*  pointerToListArrayGtol = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayGtol);
  double*  pointerToListArrayRmin = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayRmin);
  double*  pointerToListArrayRmult = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayRmult);
  int*  pointerToListArrayNsiter = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNsiter);
  double*  pointerToSkew = (double*) PyCObject_AsVoidPtr(pyPointerToSkew);
  int*  pointerToIprint = (int*) PyCObject_AsVoidPtr(pyPointerToIprint);
  int*  pointerToListArrayNcodat = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNcodat);
  int*  pointerToListArrayNunits = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNunits);
  int*  pointerToListArrayNprint = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNprint);
  int*  pointerToIstatout = (int*) PyCObject_AsVoidPtr(pyPointerToIstatout);


  viscos_f(pointerToAlnz,                     // Sparse matrix arrays
	   pointerToPcg,
	   pointerToZcg,
	   pointerToDprev,
	   pointerToJa,
	   pointerToBextern,                  // Force vectors
	   pointerToBtraction,
	   pointerToBgravity,
	   pointerToBconcForce,
	   // pointerToBprestress,
	   pointerToBintern,
	   pointerToBresid,
	   pointerToBwork,
	   pointerToDispVec,
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
	   pointerToIen,
	   pointerToLm,
	   pointerToLmx,
	   pointerToLmf,
	   pointerToInfiel,
	   pointerToListArrayIddmat,
	   pointerToListArrayNpar,
	   pointerToIelno,                    // Traction BC arrays
	   pointerToIside,
	   pointerToIhistry,
	   pointerToPres,
	   pointerToPdir,
	   pointerToListArrayPropertyList,    // Material property arrays
	   pointerToMhist,
	   pointerToMaterialInfo,
	   pointerToMaterialModelInfo,
	   pointerToMaterialModelStates,
	   pointerGauss,                      // Element type arrays
	   pointerToSh,
	   pointerToShj,
	   pointerToElementTypeInfo,
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
	   pointerToListArrayGcurr,
	   pointerToListArrayGi,
	   pointerToListArrayGprev,
	   pointerToListArrayGtol,
	   pointerToListArrayRmin,
	   pointerToListArrayRmult,
	   pointerToListArrayNsiter,
	   pointerToSkew,                     // Skew rotation information
	   pointerToIprint,                   // Input/output information
	   pointerToListArrayNcodat,
	   pointerToListArrayNunits,
	   pointerToListArrayNprint,
	   pointerToIstatout,
	   asciiOutputFile,                   // Output file names
	   plotOutputFile,
	   ucdOutputRoot,
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
// $Id: viscos.cc,v 1.7 2005/02/24 00:41:01 willic3 Exp $

// End of file
