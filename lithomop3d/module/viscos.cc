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
#include "lithomop3d_externs.h"


// Perform the time-dependent solution.

char pylithomop3d_viscos__doc__[] = "";
char pylithomop3d_viscos__name__[] = "viscos";

PyObject * pylithomop3d_viscos(PyObject *, PyObject *args)
{
  PyObject* pyPointerToAlnz;              // Sparse matrix arrays
  PyObject* pyPointerToPcg;
  PyObject* pyPointerToZcg;
  PyObject* pyPointerToJa;
  PyObject* pyPointerToB;                 // Force vectors
  PyObject* pyPointerToBtot;
  PyObject* pyPointerToBres;
  PyObject* pyPointerToPvec;
  PyObject* pyPointerToGvec1;
  PyObject* pyPointerToGvec2;
  PyObject* pyPointerToX;                 // Global arrays
  PyObject* pyPointerToD;
  PyObject* pyPointerToDx;
  PyObject* pyPointerToDeld;
  PyObject* pyPointerToDeldx;
  PyObject* pyPointerToDprev;
  PyObject* pyPointerToDcur;
  PyObject* pyPointerToDxcur;
  PyObject* pyPointerToId;
  PyObject* pyPointerToIdx;
  PyObject* pyPointerToSkew;
  PyObject* pyPointerToHistry;
  PyObject* pyPointerToIen;                // Element arrays
  PyObject* pyPointerToInfin;
  PyObject* pyPointerToMat;
  PyObject* pyPointerToLm;
  PyObject* pyPointerToLmx;
  PyObject* pyPointerToLmf;
  PyObject* pyPointerToProp;
  PyObject* pyPointerToListArrayGauss;
  PyObject* pyPointerToIbond;              // Boundary condition arrays
  PyObject* pyPointerToBond;
  PyObject* pyPointerToDmat;               // Gauss point arrays
  PyObject* pyPointerToStn;
  PyObject* pyPointerToScur;
  PyObject* pyPointerToSt0;
  PyObject* pyPointerToEps;
  PyObject* pyPointerToDeps;
  PyObject* pyPointerToBeta;
  PyObject* pyPointerToDbeta;
  PyObject* pyPointerToBetb;
  PyObject* pyPointerToDbetb;
  PyObject* pyPointerToListArrayIddmat;
  PyObject* pyPointerToIelno;               // Traction BC arrays
  PyObject* pyPointerToIside;
  PyObject* pyPointerToIhistry;
  PyObject* pyPointerToPres;
  PyObject* pyPointerToPdir;
  PyObject* pyPointerToMaxstp;              // Time step arrays
  PyObject* pyPointerToDelt;
  PyObject* pyPointerToAlfa;
  PyObject* pyPointerToMaxit;
  PyObject* pyPointerToMaxitc;
  PyObject* pyPointerToLgdef;
  PyObject* pyPointerToIbbar;
  PyObject* pyPointerToUtol;
  PyObject* pyPointerToFtol;
  PyObject* pyPointerToEtol;
  PyObject* pyPointerToItmax;
  PyObject* pyPointerToIprint;              // Time step printing info
  PyObject* pyPointerToFault;               // Split node arrays
  PyObject* pyPointerToNfault;
  PyObject* pyPointerToDfault;
  PyObject* pyPointerToTfault;
  PyObject* pyPointerToIdftn;
  PyObject* pyPointerToIdslp;               // Slippery node arrays
  PyObject* pyPointerToIpslp;
  PyObject* pyPointerToDiforc;
  PyObject* pyPointerToIdhist;
  PyObject* pyPointerToIwink;               // Winkler arrays
  PyObject* pyPointerToWink;
  PyObject* pyPointerToIwinkx;
  PyObject* pyPointerToWinkx;
  PyObject* pyPointerToS;                   // Local stiffness matrix arrays
  PyObject* pyPointerToStemp;
  PyObject* pyPointerToListArrayGcurr;      // Global variables stored as lists
  PyObject* pyPointerToListArrayGi;         // in python and as arrays in f77.
  PyObject* pyPointerToListArrayGprev;
  PyObject* pyPointerToListArrayGrav;
  PyObject* pyPointerToListArrayGtol;
  PyObject* pyPointerToListArrayNcodat;
  PyObject* pyPointerToListArrayNconsts;
  PyObject* pyPointerToListArrayNdimens;
  PyObject* pyPointerToListArrayNpar;
  PyObject* pyPointerToListArrayNprint;
  PyObject* pyPointerToListArrayNsiter;
  PyObject* pyPointerToListArrayNsysdat;
  PyObject* pyPointerToListArrayNtimdat;
  PyObject* pyPointerToListArrayNunits;
  PyObject* pyPointerToListArrayNvisdat;
  PyObject* pyPointerToListArrayRconsts;
  PyObject* pyPointerToListArrayRgiter;
  PyObject* pyPointerToListArrayRmin;
  PyObject* pyPointerToListArrayRmult;
  PyObject* pyPointerToListArrayRtimdat;
  char* asciiOutputFile;                     // Output file names
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOss:viscos",
			    &pyPointerToAlnz,              // Sparse matrix arrays
			    &pyPointerToPcg,
			    &pyPointerToZcg,
			    &pyPointerToJa,
			    &pyPointerToB,                 // Force vectors
			    &pyPointerToBtot,
			    &pyPointerToBres,
			    &pyPointerToPvec,
			    &pyPointerToGvec1,
			    &pyPointerToGvec2,
			    &pyPointerToX,                 // Global arrays
			    &pyPointerToD,
			    &pyPointerToDx,
			    &pyPointerToDeld,
			    &pyPointerToDeldx,
			    &pyPointerToDprev,
			    &pyPointerToDcur,
			    &pyPointerToDxcur,
			    &pyPointerToId,
			    &pyPointerToIdx,
			    &pyPointerToSkew,
			    &pyPointerToHistry,
			    &pyPointerToIen,                // Element arrays
			    &pyPointerToInfin,
			    &pyPointerToMat,
			    &pyPointerToLm,
			    &pyPointerToLmx,
			    &pyPointerToLmf,
			    &pyPointerToProp,
			    &pyPointerToListArrayGauss,
			    &pyPointerToIbond,              // Boundary condition arrays
			    &pyPointerToBond,
			    &pyPointerToDmat,               // Gauss point arrays
			    &pyPointerToStn,
			    &pyPointerToScur,
			    &pyPointerToSt0,
			    &pyPointerToEps,
			    &pyPointerToDeps,
			    &pyPointerToBeta,
			    &pyPointerToDbeta,
			    &pyPointerToBetb,
			    &pyPointerToDbetb,
			    &pyPointerToListArrayIddmat,
			    &pyPointerToIelno,               // Traction BC arrays
			    &pyPointerToIside,
			    &pyPointerToIhistry,
			    &pyPointerToPres,
			    &pyPointerToPdir,
			    &pyPointerToMaxstp,              // Time step arrays
			    &pyPointerToDelt,
			    &pyPointerToAlfa,
			    &pyPointerToMaxit,
			    &pyPointerToMaxitc,
			    &pyPointerToLgdef,
			    &pyPointerToIbbar,
			    &pyPointerToUtol,
			    &pyPointerToFtol,
			    &pyPointerToEtol,
			    &pyPointerToItmax,
			    &pyPointerToIprint,              // Time step printing info
			    &pyPointerToFault,               // Split node arrays
			    &pyPointerToNfault,
			    &pyPointerToDfault,
			    &pyPointerToTfault,
			    &pyPointerToIdftn,
			    &pyPointerToIdslp,               // Slippery node arrays
			    &pyPointerToIpslp,
			    &pyPointerToDiforc,
			    &pyPointerToIdhist,
			    &pyPointerToIwink,               // Winkler arrays
			    &pyPointerToWink,
			    &pyPointerToIwinkx,
			    &pyPointerToWinkx,
			    &pyPointerToS,                   // Local stiffness matrix arrays
			    &pyPointerToStemp,
			    &pyPointerToListArrayGcurr,      // Global variables stored as lists
			    &pyPointerToListArrayGi,         // in python and as arrays in f77.
			    &pyPointerToListArrayGprev,
			    &pyPointerToListArrayGrav,
			    &pyPointerToListArrayGtol,
			    &pyPointerToListArrayNcodat,
			    &pyPointerToListArrayNconsts,
			    &pyPointerToListArrayNdimens,
			    &pyPointerToListArrayNpar,
			    &pyPointerToListArrayNprint,
			    &pyPointerToListArrayNsiter,
			    &pyPointerToListArrayNsysdat,
			    &pyPointerToListArrayNtimdat,
			    &pyPointerToListArrayNunits,
			    &pyPointerToListArrayNvisdat,
			    &pyPointerToListArrayRconsts,
			    &pyPointerToListArrayRgiter,
			    &pyPointerToListArrayRmin,
			    &pyPointerToListArrayRmult,
			    &pyPointerToListArrayRtimdat,
			    &asciiOutputFile,                     // Output file names
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  double* pointerToAlnz = (double*) PyCObject_AsVoidPtr(pyPointerToAlnz);
  double* pointerToPcg = (double*) PyCObject_AsVoidPtr(pyPointerToPcg);
  double* pointerToZcg = (double*) PyCObject_AsVoidPtr(pyPointerToZcg);
  int* pointerToJa = (int*) PyCObject_AsVoidPtr(pyPointerToJa);
  double* pointerToB  = (double*) PyCObject_AsVoidPtr(pyPointerToB);
  double* pointerToBtot = (double*) PyCObject_AsVoidPtr(pyPointerToBtot);
  double* pointerToBres = (double*) PyCObject_AsVoidPtr(pyPointerToBres);
  double* pointerToPvec = (double*) PyCObject_AsVoidPtr(pyPointerToPvec);
  double* pointerToGvec1 = (double*) PyCObject_AsVoidPtr(pyPointerToGvec1);
  double* pointerToGvec2 = (double*) PyCObject_AsVoidPtr(pyPointerToGvec2);
  double* pointerToX = (double*) PyCObject_AsVoidPtr(pyPointerToX);
  double* pointerToD = (double*) PyCObject_AsVoidPtr(pyPointerToD);
  double* pointerToDx = (double*) PyCObject_AsVoidPtr(pyPointerToDx);
  double* pointerToDeld = (double*) PyCObject_AsVoidPtr(pyPointerToDeld);
  double* pointerToDeldx = (double*) PyCObject_AsVoidPtr(pyPointerToDeldx);
  double* pointerToDprev = (double*) PyCObject_AsVoidPtr(pyPointerToDprev);
  double* pointerToDcur = (double*) PyCObject_AsVoidPtr(pyPointerToDcur);
  double* pointerToDxcur = (double*) PyCObject_AsVoidPtr(pyPointerToDxcur);
  int* pointerToId = (int*) PyCObject_AsVoidPtr(pyPointerToId);
  int* pointerToIdx = (int*) PyCObject_AsVoidPtr(pyPointerToIdx);
  double* pointerToSkew = (double*) PyCObject_AsVoidPtr(pyPointerToSkew);
  double* pointerToHistry = (double*) PyCObject_AsVoidPtr(pyPointerToHistry);
  int* pointerToIen = (int*) PyCObject_AsVoidPtr(pyPointerToIen);
  int* pointerToInfin = (int*) PyCObject_AsVoidPtr(pyPointerToInfin);
  int* pointerToMat = (int*) PyCObject_AsVoidPtr(pyPointerToMat);
  int* pointerToLm = (int*) PyCObject_AsVoidPtr(pyPointerToLm);
  int* pointerToLmx = (int*) PyCObject_AsVoidPtr(pyPointerToLmx);
  int* pointerToLmf = (int*) PyCObject_AsVoidPtr(pyPointerToLmf);
  double* pointerToProp = (double*) PyCObject_AsVoidPtr(pyPointerToProp);
  double* pointerToListArrayGauss = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayGauss);
  int* pointerToIbond = (int*) PyCObject_AsVoidPtr(pyPointerToIbond);
  double* pointerToBond = (double*) PyCObject_AsVoidPtr(pyPointerToBond);
  double* pointerToDmat = (double*) PyCObject_AsVoidPtr(pyPointerToDmat);
  double* pointerToStn = (double*) PyCObject_AsVoidPtr(pyPointerToStn);
  double* pointerToScur = (double*) PyCObject_AsVoidPtr(pyPointerToScur);
  double* pointerToSt0 = (double*) PyCObject_AsVoidPtr(pyPointerToSt0);
  double* pointerToEps = (double*) PyCObject_AsVoidPtr(pyPointerToEps);
  double* pointerToDeps = (double*) PyCObject_AsVoidPtr(pyPointerToDeps);
  double* pointerToBeta = (double*) PyCObject_AsVoidPtr(pyPointerToBeta);
  double* pointerToDbeta = (double*) PyCObject_AsVoidPtr(pyPointerToDbeta);
  double* pointerToBetb = (double*) PyCObject_AsVoidPtr(pyPointerToBetb);
  double* pointerToDbetb = (double*) PyCObject_AsVoidPtr(pyPointerToDbetb);
  int* pointerToListArrayIddmat = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayIddmat);
  int* pointerToIelno = (int*) PyCObject_AsVoidPtr(pyPointerToIelno);
  int* pointerToIside = (int*) PyCObject_AsVoidPtr(pyPointerToIside);
  int* pointerToIhistry = (int*) PyCObject_AsVoidPtr(pyPointerToIhistry);
  double* pointerToPres = (double*) PyCObject_AsVoidPtr(pyPointerToPres);
  double* pointerToPdir = (double*) PyCObject_AsVoidPtr(pyPointerToPdir);
  int* pointerToMaxstp = (int*) PyCObject_AsVoidPtr(pyPointerToMaxstp);
  double* pointerToDelt = (double*) PyCObject_AsVoidPtr(pyPointerToDelt);
  double* pointerToAlfa = (double*) PyCObject_AsVoidPtr(pyPointerToAlfa);
  int* pointerToMaxit = (int*) PyCObject_AsVoidPtr(pyPointerToMaxit);
  int* pointerToMaxitc = (int*) PyCObject_AsVoidPtr(pyPointerToMaxitc);
  int* pointerToLgdef = (int*) PyCObject_AsVoidPtr(pyPointerToLgdef);
  int* pointerToIbbar = (int*) PyCObject_AsVoidPtr(pyPointerToIbbar);
  double* pointerToUtol = (double*) PyCObject_AsVoidPtr(pyPointerToUtol);
  double* pointerToFtol = (double*) PyCObject_AsVoidPtr(pyPointerToFtol);
  double* pointerToEtol = (double*) PyCObject_AsVoidPtr(pyPointerToEtol);
  int* pointerToItmax = (int*) PyCObject_AsVoidPtr(pyPointerToItmax);
  int* pointerToIprint = (int*) PyCObject_AsVoidPtr(pyPointerToIprint);
  double* pointerToFault = (double*) PyCObject_AsVoidPtr(pyPointerToFault);
  int* pointerToNfault = (int*) PyCObject_AsVoidPtr(pyPointerToNfault);
  double* pointerToDfault = (double*) PyCObject_AsVoidPtr(pyPointerToDfault);
  double* pointerToTfault = (double*) PyCObject_AsVoidPtr(pyPointerToTfault);
  int* pointerToIdftn = (int*) PyCObject_AsVoidPtr(pyPointerToIdftn);
  int* pointerToIdslp = (int*) PyCObject_AsVoidPtr(pyPointerToIdslp);
  int* pointerToIpslp = (int*) PyCObject_AsVoidPtr(pyPointerToIpslp);
  double* pointerToDiforc = (double*) PyCObject_AsVoidPtr(pyPointerToDiforc);
  int* pointerToIdhist = (int*) PyCObject_AsVoidPtr(pyPointerToIdhist);
  int* pointerToIwink = (int*) PyCObject_AsVoidPtr(pyPointerToIwink);
  double* pointerToWink = (double*) PyCObject_AsVoidPtr(pyPointerToWink);
  int* pointerToIwinkx = (int*) PyCObject_AsVoidPtr(pyPointerToIwinkx);
  double* pointerToWinkx = (double*) PyCObject_AsVoidPtr(pyPointerToWinkx);
  double* pointerToS = (double*) PyCObject_AsVoidPtr(pyPointerToS);
  double* pointerToStemp = (double*) PyCObject_AsVoidPtr(pyPointerToStemp);
  double* pointerToListArrayGcurr = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayGcurr);
  double* pointerToListArrayGi = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayGi);
  double* pointerToListArrayGprev = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayGprev);
  double* pointerToListArrayGrav = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayGrav);
  double* pointerToListArrayGtol = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayGtol);
  int* pointerToListArrayNcodat = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNcodat);
  int* pointerToListArrayNconsts = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNconsts);
  int* pointerToListArrayNdimens = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNdimens);
  int* pointerToListArrayNpar = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNpar);
  int* pointerToListArrayNprint = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNprint);
  int* pointerToListArrayNsiter = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNsiter);
  int* pointerToListArrayNsysdat = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNsysdat);
  int* pointerToListArrayNtimdat = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNtimdat);
  int* pointerToListArrayNunits = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNunits);
  int* pointerToListArrayNvisdat = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNvisdat);
  double* pointerToListArrayRconsts = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayRconsts);
  double* pointerToListArrayRgiter = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayRgiter);
  double* pointerToListArrayRmin = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayRmin);
  double* pointerToListArrayRmult = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayRmult);
  double* pointerToListArrayRtimdat = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayRtimdat);

  viscos_f(pointerToAlnz,              // Sparse matrix arrays
	   pointerToPcg,
	   pointerToZcg,
	   pointerToJa,
	   pointerToB,                 // Force vectors
	   pointerToBtot,
	   pointerToBres,
	   pointerToPvec,
	   pointerToGvec1,
	   pointerToGvec2,
	   pointerToX,                 // Global arrays
	   pointerToD,
	   pointerToDx,
	   pointerToDeld,
	   pointerToDeldx,
	   pointerToDprev,
	   pointerToDcur,
	   pointerToDxcur,
	   pointerToId,
	   pointerToIdx,
	   pointerToSkew,
	   pointerToHistry,
	   pointerToIen,                // Element arrays
	   pointerToInfin,
	   pointerToMat,
	   pointerToLm,
	   pointerToLmx,
	   pointerToLmf,
	   pointerToProp,
	   pointerToListArrayGauss,
	   pointerToIbond,              // Boundary condition arrays
	   pointerToBond,
	   pointerToDmat,               // Gauss point arrays
	   pointerToStn,
	   pointerToScur,
	   pointerToSt0,
	   pointerToEps,
	   pointerToDeps,
	   pointerToBeta,
	   pointerToDbeta,
	   pointerToBetb,
	   pointerToDbetb,
	   pointerToListArrayIddmat,
	   pointerToIelno,               // Traction BC arrays
	   pointerToIside,
	   pointerToIhistry,
	   pointerToPres,
	   pointerToPdir,
	   pointerToMaxstp,              // Time step arrays
	   pointerToDelt,
	   pointerToAlfa,
	   pointerToMaxit,
	   pointerToMaxitc,
	   pointerToLgdef,
	   pointerToIbbar,
	   pointerToUtol,
	   pointerToFtol,
	   pointerToEtol,
	   pointerToItmax,
	   pointerToIprint,              // Time step printing info
	   pointerToFault,               // Split node arrays
	   pointerToNfault,
	   pointerToDfault,
	   pointerToTfault,
	   pointerToIdftn,
	   pointerToIdslp,               // Slippery node arrays
	   pointerToIpslp,
	   pointerToDiforc,
	   pointerToIdhist,
	   pointerToIwink,               // Winkler arrays
	   pointerToWink,
	   pointerToIwinkx,
	   pointerToWinkx,
	   pointerToS,                   // Local stiffness matrix arrays
	   pointerToStemp,
	   pointerToListArrayGcurr,      // Global variables stored as lists
	   pointerToListArrayGi,         // in python and as arrays in f77.
	   pointerToListArrayGprev,
	   pointerToListArrayGrav,
	   pointerToListArrayGtol,
	   pointerToListArrayNcodat,
	   pointerToListArrayNconsts,
	   pointerToListArrayNdimens,
	   pointerToListArrayNpar,
	   pointerToListArrayNprint,
	   pointerToListArrayNsiter,
	   pointerToListArrayNsysdat,
	   pointerToListArrayNtimdat,
	   pointerToListArrayNunits,
	   pointerToListArrayNvisdat,
	   pointerToListArrayRconsts,
	   pointerToListArrayRgiter,
	   pointerToListArrayRmin,
	   pointerToListArrayRmult,
	   pointerToListArrayRtimdat,
	   asciiOutputFile,              // Output file names
	   plotOutputFile,
	   strlen(asciiOutputFile),
	   strlen(plotOutputFile));

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
// $Id: viscos.cc,v 1.1 2004/04/14 21:24:47 willic3 Exp $

// End of file
