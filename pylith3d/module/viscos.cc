// -*- C++ -*-
// 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  PyLith by Charles A. Williams
//  Copyright (c) 2003-2006 Rensselaer Polytechnic Institute
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

#include <petscmat.h>
#include <portinfo>
#include "journal/debug.h"

#include <Python.h>

#include "viscos.h"
#include "exceptionhandler.h"
#include "pylith3d_externs.h"
#include <stdio.h>
#include <string.h>


// Perform the time-dependent solution.

char pypylith3d_viscos__doc__[] = "";
char pypylith3d_viscos__name__[] = "viscos";

PyObject * pypylith3d_viscos(PyObject *, PyObject *args)
{
  PyObject* pyA;                              // Sparse matrix arrays
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
  PyObject* pyPointerToIen;
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
  PyObject* pyPointerToIprint;                // Input/output information
  PyObject* pyPointerToListArrayNcodat;
  PyObject* pyPointerToListArrayNunits;
  PyObject* pyPointerToListArrayNprint;
  PyObject* pyPointerToIstatout;
  PyObject* pyPointerToNstatout;
  char* asciiOutputFile;                      // Output file names
  char* plotOutputFile;
  char* ucdOutputRoot;
  int viscousStage;
  int iterateEvent;

  int ok = PyArg_ParseTuple(args, "OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOsssii:viscos",
			    &pyA,                              // Sparse matrix arrays
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
			    &pyPointerToIen,
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
			    &pyPointerToIprint,                // Input/output information
			    &pyPointerToListArrayNcodat,
			    &pyPointerToListArrayNunits,
			    &pyPointerToListArrayNprint,
			    &pyPointerToIstatout,
			    &pyPointerToNstatout,
			    &asciiOutputFile,                  // Output file names
			    &plotOutputFile,
			    &ucdOutputRoot,
                            &viscousStage,
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
  int*  pointerToIen = (int*) PyCObject_AsVoidPtr(pyPointerToIen);
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
  int*  pointerToIprint = (int*) PyCObject_AsVoidPtr(pyPointerToIprint);
  int*  pointerToListArrayNcodat = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNcodat);
  int*  pointerToListArrayNunits = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNunits);
  int*  pointerToListArrayNprint = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNprint);
  int*  pointerToIstatout = (int*) PyCObject_AsVoidPtr(pyPointerToIstatout);
  int*  pointerToNstatout = (int*) PyCObject_AsVoidPtr(pyPointerToNstatout);


  viscos_f(&A,                                // Sparse matrix arrays
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
	   pointerToIen,
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
	   pointerToIprint,                   // Input/output information
	   pointerToListArrayNcodat,
	   pointerToListArrayNunits,
	   pointerToListArrayNprint,
	   pointerToIstatout,
	   pointerToNstatout,
	   asciiOutputFile,                   // Output file names
	   plotOutputFile,
	   ucdOutputRoot,
           &viscousStage,
           &iterateEvent,
	   &errorcode,                        // Error codes
	   errorstring,
	   strlen(asciiOutputFile),           // String lengths
	   strlen(plotOutputFile),
	   strlen(ucdOutputRoot),
	   sizeof(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }


  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "asciiOutputFile:" << asciiOutputFile
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_None;
}


// version
// $Id: viscos.cc,v 1.15 2005/06/07 19:39:11 willic3 Exp $

// End of file
