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
  PyObject* pyPointerToIelindx;
  PyObject* pyPointerToTractionverts;         // Traction BC arrays
  PyObject* pyPointerToTractionvals;
  PyObject* pyPointerToGauss2d;
  PyObject* pyPointerToSh2d;
  PyObject* pyPointerToListArrayElementTypeInfo2d;
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

  int ok = PyArg_ParseTuple(args, "OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOsssii:autoprestr",
                            &pyA,                             // Sparse matrix arrays
                            &pyRhs,
                            &pySol,
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
                            &pyPointerToIelindx,
                            &pyPointerToTractionverts,         // Traction BC arrays
                            &pyPointerToTractionvals,
                            &pyPointerToGauss2d,
                            &pyPointerToSh2d,
                            &pyPointerToListArrayElementTypeInfo2d,
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
  Vec      rhs = (Vec) PyCObject_AsVoidPtr(pyRhs);
  Vec      sol = (Vec) PyCObject_AsVoidPtr(pySol);
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
  int*  pointerToIelindx = (int*) PyCObject_AsVoidPtr(pyPointerToIelindx);
  int*  pointerToTractionverts = (int*) PyCObject_AsVoidPtr(pyPointerToTractionverts);
  double*  pointerToTractionvals = (double*) PyCObject_AsVoidPtr(pyPointerToTractionvals);
  double*  pointerToGauss2d = (double*) PyCObject_AsVoidPtr(pyPointerToGauss2d);
  double*  pointerToSh2d = (double*) PyCObject_AsVoidPtr(pyPointerToSh2d);
  int*  pointerToListArrayElementTypeInfo2d = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayElementTypeInfo2d);
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
               &rhs,
               &sol,
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
               pointerToIelindx,
               pointerToTractionverts,            // Traction BC arrays
               pointerToTractionvals,
               pointerToGauss2d,
               pointerToSh2d,
               pointerToListArrayElementTypeInfo2d,
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



// Perform the elastic solution.

char pypylith3d_elastc__doc__[] = "";
char pypylith3d_elastc__name__[] = "elastc";

PyObject * pypylith3d_elastc(PyObject *, PyObject *args)
{
  PyObject* pyA;                              // Sparse matrix arrays
  PyObject* pyRhs;
  PyObject* pySol;
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
  PyObject* pyPointerToIelindx;
  PyObject* pyPointerToTractionverts;         // Traction BC arrays
  PyObject* pyPointerToTractionvals;
  PyObject* pyPointerToGauss2d;
  PyObject* pyPointerToSh2d;
  PyObject* pyPointerToListArrayElementTypeInfo2d;
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
  int elasticStage, iterateEvent;             // PETSc logging

  int ok = PyArg_ParseTuple(args, "OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOsssii:elastc",
			    &pyA,                              // Sparse matrix arrays
			    &pyRhs,
			    &pySol,
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
                            &pyPointerToIelindx,
                            &pyPointerToTractionverts,         // Traction BC arrays
                            &pyPointerToTractionvals,
                            &pyPointerToGauss2d,
                            &pyPointerToSh2d,
                            &pyPointerToListArrayElementTypeInfo2d,
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
  int*  pointerToIelindx = (int*) PyCObject_AsVoidPtr(pyPointerToIelindx);
  int*  pointerToTractionverts = (int*) PyCObject_AsVoidPtr(pyPointerToTractionverts);
  double*  pointerToTractionvals = (double*) PyCObject_AsVoidPtr(pyPointerToTractionvals);
  double*  pointerToGauss2d = (double*) PyCObject_AsVoidPtr(pyPointerToGauss2d);
  double*  pointerToSh2d = (double*) PyCObject_AsVoidPtr(pyPointerToSh2d);
  int*  pointerToListArrayElementTypeInfo2d = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayElementTypeInfo2d);
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

  elastc_f(&A,                                // Sparse matrix arrays
           &rhs,
	   &sol,
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
	   pointerToIelindx,
	   pointerToTractionverts,            // Traction BC arrays
	   pointerToTractionvals,
	   pointerToGauss2d,
	   pointerToSh2d,
	   pointerToListArrayElementTypeInfo2d,
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
           &elasticStage,
           &iterateEvent,
	   &errorcode,                        // Error codes
	   errorstring,
	   strlen(asciiOutputFile),           // String lengths
	   strlen(plotOutputFile),
	   strlen(ucdOutputRoot),
	   sizeof(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    //int *fail = NULL;
    //int i = fail[0];
    return NULL;
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




// Assign equation numbers to Iwink and Iwinkx arrays, and compact
// Wink and Winkx arrays to correspond to assigned BC.

char pypylith3d_assign_wink__doc__[] = "";
char pypylith3d_assign_wink__name__[] = "assign_wink";

PyObject * pypylith3d_assign_wink(PyObject *, PyObject *args)
{
  PyObject* pyPointerToWinkdef;
  PyObject* pyPointerToWink;
  PyObject* pyPointerToIwinkdef;
  PyObject* pyPointerToIwinkid;
  PyObject* pyPointerToIwink;
  PyObject* pyPointerToId;
  int numberNodes;
  int numberWinklerForces;
  int numberWinklerEntries;

  int ok = PyArg_ParseTuple(args, "OOOOOOiii:assign_wink",
			    &pyPointerToWinkdef,
			    &pyPointerToWink,
			    &pyPointerToIwinkdef,
			    &pyPointerToIwinkid,
			    &pyPointerToIwink,
			    &pyPointerToId,
			    &numberNodes,
			    &numberWinklerForces,
			    &numberWinklerEntries);

  if (!ok) {
    return 0;
  }

  double* pointerToWinkdef = (double*) PyCObject_AsVoidPtr(pyPointerToWinkdef);
  double* pointerToWink = (double*) PyCObject_AsVoidPtr(pyPointerToWink);
  int* pointerToIwinkdef = (int*) PyCObject_AsVoidPtr(pyPointerToIwinkdef);
  int* pointerToIwinkid = (int*) PyCObject_AsVoidPtr(pyPointerToIwinkid);
  int* pointerToIwink = (int*) PyCObject_AsVoidPtr(pyPointerToIwink);
  int* pointerToId = (int*) PyCObject_AsVoidPtr(pyPointerToId);

  assign_wink_f(pointerToWinkdef,
		pointerToWink,
		pointerToIwinkdef,
		pointerToIwinkid,
		pointerToIwink,
		pointerToId,
		&numberNodes,
		&numberWinklerForces,
		&numberWinklerEntries);

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberWinklerForces:" << numberWinklerForces
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
  PyObject* pyPointerToId;
  PyObject* pyPointerToIdx;
  PyObject* pyPointerToIbond;
  PyObject* pyPointerToNslip;
  PyObject* pyPointerToIdslp;
  int numberSlipperyNodeEntries;
  int numberNodes;
  int totalNumberSlipperyNodes;

  int ok = PyArg_ParseTuple(args, "OOOOOiii:create_id",
			    &pyPointerToId,
			    &pyPointerToIdx,
			    &pyPointerToIbond,
			    &pyPointerToNslip,
			    &pyPointerToIdslp,
			    &numberSlipperyNodeEntries,
			    &numberNodes,
			    &totalNumberSlipperyNodes);

  if (!ok) {
    return 0;
  }

  int* pointerToId = (int*) PyCObject_AsVoidPtr(pyPointerToId);
  int* pointerToIdx = (int*) PyCObject_AsVoidPtr(pyPointerToIdx);
  int* pointerToIbond = (int*) PyCObject_AsVoidPtr(pyPointerToIbond);
  int* pointerToNslip = (int*) PyCObject_AsVoidPtr(pyPointerToNslip);
  int* pointerToIdslp = (int*) PyCObject_AsVoidPtr(pyPointerToIdslp);
  int numberGlobalEquations = 0;

  create_id_f(pointerToId,
	  pointerToIdx,
	  pointerToIbond,
	  pointerToNslip,
	  pointerToIdslp,
	  &numberSlipperyNodeEntries,
	  &numberNodes,
	  &totalNumberSlipperyNodes,
	  &numberGlobalEquations);

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberGlobalEquations:" << numberGlobalEquations
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i", numberGlobalEquations);
}


// Form id array for split nodes

char pypylith3d_id_split__doc__[] = "";
char pypylith3d_id_split__name__[] = "id_split";

PyObject * pypylith3d_id_split(PyObject *, PyObject *args)
{
  PyObject* pyPointerToNfault;
  PyObject* pyPointerToIdftn;
  int numberNodes;
  int numberSplitNodeEntries;
  int totalNumberSplitNodes;

  int ok = PyArg_ParseTuple(args, "OOiii:id_split",
			    &pyPointerToNfault,
			    &pyPointerToIdftn,
			    &numberNodes,
			    &numberSplitNodeEntries,
			    &totalNumberSplitNodes);

  if (!ok) {
    return 0;
  }

  int* pointerToNfault = (int*) PyCObject_AsVoidPtr(pyPointerToNfault);
  int* pointerToIdftn = (int*) PyCObject_AsVoidPtr(pyPointerToIdftn);

  id_split_f(pointerToNfault,
	     pointerToIdftn,
	     &numberNodes,
	     &numberSplitNodeEntries,
	     &totalNumberSplitNodes);

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberSplitNodeEntries:" << numberSplitNodeEntries
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
  PyObject* pyPointerToX;
  PyObject* pyPointerToXtmp;
  PyObject* pyPointerToIdslp;
  PyObject* pyPointerToIpslp;
  PyObject* pyPointerToItmp;
  PyObject* pyPointerToItmp1;
  PyObject* pyPointerToItmp2;
  PyObject* pyPointerToNslip;
  int numberSlipperyNodeEntries;
  int totalNumberSlipperyNodes;
  int numberNodes;

  int ok = PyArg_ParseTuple(args, "OOOOOOOOiii:nfind",
			    &pyPointerToX,
			    &pyPointerToXtmp,
			    &pyPointerToIdslp,
			    &pyPointerToIpslp,
			    &pyPointerToItmp,
			    &pyPointerToItmp1,
			    &pyPointerToItmp2,
			    &pyPointerToNslip,
			    &numberSlipperyNodeEntries,
			    &totalNumberSlipperyNodes,
			    &numberNodes);

  if (!ok) {
    return 0;
  }

  double* pointerToX = (double*) PyCObject_AsVoidPtr(pyPointerToX);
  double* pointerToXtmp = (double*) PyCObject_AsVoidPtr(pyPointerToXtmp);
  int* pointerToIdslp = (int*) PyCObject_AsVoidPtr(pyPointerToIdslp);
  int* pointerToIpslp = (int*) PyCObject_AsVoidPtr(pyPointerToIpslp);
  int* pointerToItmp = (int*) PyCObject_AsVoidPtr(pyPointerToItmp);
  int* pointerToItmp1 = (int*) PyCObject_AsVoidPtr(pyPointerToItmp1);
  int* pointerToItmp2 = (int*) PyCObject_AsVoidPtr(pyPointerToItmp2);
  int* pointerToNslip = (int*) PyCObject_AsVoidPtr(pyPointerToNslip);

  nfind_f(pointerToX,
	  pointerToXtmp,
	  pointerToIdslp,
	  pointerToIpslp,
	  pointerToItmp,
	  pointerToItmp1,
	  pointerToItmp2,
	  pointerToNslip,
	  &numberSlipperyNodeEntries,
	  &totalNumberSlipperyNodes,
	  &numberNodes);

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "pointerToX:" << pointerToX
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
  PyObject* pyPointerToBond;
  double displacementScaleFactor;
  double velocityScaleFactor;
  double forceScaleFactor;
  PyObject* pyPointerToIbond;
  int numberNodes;
  int numberBcEntries;
  int f77FileInput;
  char* bcInputFile;

  int ok = PyArg_ParseTuple(args, "OdddOiiis:read_bc",
			    &pyPointerToBond,
			    &displacementScaleFactor,
			    &velocityScaleFactor,
			    &forceScaleFactor,
			    &pyPointerToIbond,
			    &numberNodes,
			    &numberBcEntries,
			    &f77FileInput,
			    &bcInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToBond = (double*) PyCObject_AsVoidPtr(pyPointerToBond);
  int* pointerToIbond = (int*) PyCObject_AsVoidPtr(pyPointerToIbond);
  int numberConcForces = 0;

  read_bc_f(pointerToBond,
	    &displacementScaleFactor,
	    &velocityScaleFactor,
	    &forceScaleFactor,
	    pointerToIbond,
	    &numberNodes,
	    &numberBcEntries,
	    &numberConcForces,
	    &f77FileInput,
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
    << "numberBcEntries:" << numberBcEntries
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
  PyObject* pyPointerToIen;
  PyObject* pyPointerToMat;
  int numberVolumeElementNodes;
  int numberVolumeElements;
  int numberNodes;
  int numberVolumeElementFamilies;
  int f77FileInput;
  char* connectivityInputFile;

  int ok = PyArg_ParseTuple(args, "OOiiiiis:read_connect",
                            &pyPointerToIen,
                            &pyPointerToMat,
                            &numberVolumeElementNodes,
                            &numberVolumeElements,
                            &numberNodes,
                            &numberVolumeElementFamilies,
                            &f77FileInput,
                            &connectivityInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* pointerToIen = (int*) PyCObject_AsVoidPtr(pyPointerToIen);
  int* pointerToMat = (int*) PyCObject_AsVoidPtr(pyPointerToMat);

  read_connect_f(pointerToIen,
		 pointerToMat,
		 &numberVolumeElementNodes,
                 &numberVolumeElements,
                 &numberNodes,
                 &numberVolumeElementFamilies,
		 &f77FileInput,
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
    << "numberVolumeElementFamilies:" << numberVolumeElementFamilies
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
  PyObject* pyPointerToX;
  double coordinateScaleFactor;
  int numberNodes;
  int f77FileInput;
  char* coordinateInputFile;

  int ok = PyArg_ParseTuple(args, "Odiis:read_coords",
			    &pyPointerToX,
			    &coordinateScaleFactor,
			    &numberNodes,
			    &f77FileInput,
			    &coordinateInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToX = (double*) PyCObject_AsVoidPtr(pyPointerToX);

  read_coords_f(pointerToX,
		&coordinateScaleFactor,
		&numberNodes,
		&f77FileInput,
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
    << "coordinateScaleFactor:" << coordinateScaleFactor
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
  PyObject* pyPointerToDiforc;
  PyObject* pyPointerToNslip;
  PyObject* pyPointerToIdhist;
  int numberSlipperyNodeEntries;
  int numberDifferentialForceEntries;
  int numberNodes;
  int f77FileInput;
  char* differentialForceInputFile;

  int ok = PyArg_ParseTuple(args, "OOOiiiis:read_diff",
			    &pyPointerToDiforc,
			    &pyPointerToNslip,
			    &pyPointerToIdhist,
			    &numberSlipperyNodeEntries,
			    &numberDifferentialForceEntries,
			    &numberNodes,
			    &f77FileInput,
			    &differentialForceInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToDiforc = (double*) PyCObject_AsVoidPtr(pyPointerToDiforc);
  int* pointerToNslip = (int*) PyCObject_AsVoidPtr(pyPointerToNslip);
  int* pointerToIdhist = (int*) PyCObject_AsVoidPtr(pyPointerToIdhist);

  read_diff_f(pointerToDiforc,
	      pointerToNslip,
	      pointerToIdhist,
	      &numberSlipperyNodeEntries,
	      &numberDifferentialForceEntries,
	      &numberNodes,
	      &f77FileInput,
	      differentialForceInputFile,
	      &errorcode,
	      errorstring,
	      strlen(differentialForceInputFile),
	      sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberDifferentialForceEntries:" << numberDifferentialForceEntries
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
  PyObject* pyPointerToIprint;
  int numberFullOutputs;
  int analysisTypeInt;
  int numberCycles;
  int totalNumberTimeSteps;
  int f77FileInput;
  char* fullOutputInputFile;

  int ok = PyArg_ParseTuple(args, "Oiiiiis:read_fuldat",
			    &pyPointerToIprint,
			    &numberFullOutputs,
			    &analysisTypeInt,
			    &numberCycles,
			    &totalNumberTimeSteps,
			    &f77FileInput,
			    &fullOutputInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* pointerToIprint = (int*) PyCObject_AsVoidPtr(pyPointerToIprint);

  read_fuldat_f(pointerToIprint,
		&numberFullOutputs,
		&analysisTypeInt,
		&numberCycles,
		&totalNumberTimeSteps,
		&f77FileInput,
		fullOutputInputFile,
		&errorcode,
		errorstring,
		strlen(fullOutputInputFile),
		sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberFullOutputs:" << numberFullOutputs
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
  PyObject* pyPointerToHistry;
  PyObject* pyPointerToTimes;
  int numberLoadHistories;
  int totalNumberTimeSteps;
  int f77FileInput;
  char* loadHistoryInputFile;

  int ok = PyArg_ParseTuple(args, "OOiiis:read_hist",
			    &pyPointerToHistry,
			    &pyPointerToTimes,
			    &numberLoadHistories,
			    &totalNumberTimeSteps,
			    &f77FileInput,
			    &loadHistoryInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToHistry = (double*) PyCObject_AsVoidPtr(pyPointerToHistry);
  double* pointerToTimes = (double*) PyCObject_AsVoidPtr(pyPointerToTimes);

  read_hist_f(pointerToHistry,
	      pointerToTimes,
	      &numberLoadHistories,
	      &totalNumberTimeSteps,
	      &f77FileInput,
	      loadHistoryInputFile,
	      &errorcode,
	      errorstring,
	      strlen(loadHistoryInputFile),
	      sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberLoadHistories:" << numberLoadHistories
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
  // PyObject* pyPointerToStn;
  // PyObject* pyPointerToSt0;
  // PyObject* pyPointerToListArrayPrscal;
  // int numberStressComponents;
  // int numberGaussPoints;
  // int numberPrestressGaussPoints;
  // int numberElements;
  // int numberPrestressEntries;
  // int prestressAutoComputeInt;
  // int asciiOutputInt;
  // int f77FileInput;
  // int f77AsciiOutput;
  // char* prestressInputFile;
  // char* asciiOutputFile;

  // int ok = PyArg_ParseTuple(args, "OOOiiiiiiiiiss:read_prestr",
			    // &pyPointerToStn,
			    // &pyPointerToSt0,
			    // &pyPointerToListArrayPrscal,
			    // &numberStressComponents,
			    // &numberGaussPoints,
			    // &numberPrestressGaussPoints,
			    // &numberElements,
			    // &numberPrestressEntries,
			    // &prestressAutoComputeInt,
			    // &asciiOutputInt,
			    // &f77FileInput,
			    // &f77AsciiOutput,
			    // &prestressInputFile,
			    // &asciiOutputFile);

  // if (!ok) {
    // return 0;
  // }

  // int errorcode = 0;
  // const int maxsize = 4096;
  // char errorstring[maxsize];
  // double* pointerToStn = (double*) PyCObject_AsVoidPtr(pyPointerToStn);
  // double* pointerToSt0 = (double*) PyCObject_AsVoidPtr(pyPointerToSt0);
  // double* pointerToListArrayPrscal = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayPrscal);

  // read_prestr_f(pointerToStn,
		// pointerToSt0,
		// pointerToListArrayPrscal,
		// &numberStressComponents,
		// &numberGaussPoints,
		// &numberPrestressGaussPoints,
		// &numberElements,
		// &numberPrestressEntries,
		// &prestressAutoComputeInt,
		// &asciiOutputInt,
		// &f77FileInput,
		// &f77AsciiOutput,
		// &errorcode,
		// prestressInputFile,
		// asciiOutputFile,
		// strlen(prestressInputFile),
		// strlen(asciiOutputFile));
    
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
  PyObject* pyPointerToSkew;
  double rotationScaleFactor;
  int numberRotationEntries;
  int numberNodes;
  int autoRotateSlipperyNodesInt;
  int f77FileInput;
  char* rotationInputFile;

  int ok = PyArg_ParseTuple(args, "Odiiiis:read_skew",
			    &pyPointerToSkew,
			    &rotationScaleFactor,
			    &numberRotationEntries,
			    &numberNodes,
			    &autoRotateSlipperyNodesInt,
			    &f77FileInput,
			    &rotationInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToSkew = (double*) PyCObject_AsVoidPtr(pyPointerToSkew);

  read_skew_f(pointerToSkew,
	      &rotationScaleFactor,
	      &numberRotationEntries,
	      &numberNodes,
	      &autoRotateSlipperyNodesInt,
	      &f77FileInput,
	      rotationInputFile,
	      &errorcode,
	      errorstring,
	      strlen(rotationInputFile),
	      sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberRotationEntries:" << numberRotationEntries
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
  PyObject* pyPointerToNslip;
  int numberSlipperyNodeEntries;
  int numberNodes;
  int autoRotateSlipperyNodesInt;
  int f77FileInput;
  char* slipperyNodeInputFile;

  int ok = PyArg_ParseTuple(args, "Oiiiis:read_slip",
			    &pyPointerToNslip,
			    &numberSlipperyNodeEntries,
			    &numberNodes,
			    &autoRotateSlipperyNodesInt,
			    &f77FileInput,
			    &slipperyNodeInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* pointerToNslip = (int*) PyCObject_AsVoidPtr(pyPointerToNslip);
  int totalNumberSlipperyNodes = 0;

  read_slip_f(pointerToNslip,
	      &numberSlipperyNodeEntries,
	      &totalNumberSlipperyNodes,
	      &numberNodes,
	      &autoRotateSlipperyNodesInt,
	      &f77FileInput,
	      slipperyNodeInputFile,
	      &errorcode,
	      errorstring,
	      strlen(slipperyNodeInputFile),
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
  PyObject* pyPointerToFault;
  PyObject* pyPointerToNfault;
  int numberSplitNodeEntries;
  int numberNodes;
  int numberVolumeElements;
  int f77FileInput;
  char* splitNodeInputFile;

  int ok = PyArg_ParseTuple(args, "OOiiiis:read_split",
			    &pyPointerToFault,
			    &pyPointerToNfault,
			    &numberSplitNodeEntries,
			    &numberNodes,
			    &numberVolumeElements,
			    &f77FileInput,
			    &splitNodeInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToFault = (double*) PyCObject_AsVoidPtr(pyPointerToFault);
  int* pointerToNfault = (int*) PyCObject_AsVoidPtr(pyPointerToNfault);
  int totalNumberSplitNodes = 0;

  read_split_f(pointerToFault,
	       pointerToNfault,
	       &numberSplitNodeEntries,
	       &totalNumberSplitNodes,
	       &numberNodes,
	       &numberVolumeElements,
	       &f77FileInput,
	       splitNodeInputFile,
	       &errorcode,
	       errorstring,
	       strlen(splitNodeInputFile),
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
  PyObject* pyPointerToIstatout;
  PyObject* pyPointerToNstatout;
  int f77FileInput;
  char* stateVariableInputFile;

  int ok = PyArg_ParseTuple(args, "OOis:read_stateout",
			    &pyPointerToIstatout,
			    &pyPointerToNstatout,
			    &f77FileInput,
			    &stateVariableInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* pointerToIstatout = (int*) PyCObject_AsVoidPtr(pyPointerToIstatout);
  int* pointerToNstatout = (int*) PyCObject_AsVoidPtr(pyPointerToNstatout);

  read_stateout_f(pointerToIstatout,
	          pointerToNstatout,
	          &f77FileInput,
	          stateVariableInputFile,
	          &errorcode,
	          errorstring,
	          strlen(stateVariableInputFile),
	          sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "pointerToIstatout:" << pointerToIstatout
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
  PyObject* pyPointerToDelt;
  PyObject* pyPointerToAlfa;
  PyObject* pyPointerToUtol;
  PyObject* pyPointerToFtol;
  PyObject* pyPointerToEtol;
  PyObject* pyPointerToTimes;
  double timeScaleFactor;
  PyObject* pyPointerToMaxstp;
  PyObject* pyPointerToMaxit;
  PyObject* pyPointerToNtdinit;
  PyObject* pyPointerToLgdef;
  PyObject* pyPointerToItmax;
  int numberTimeStepGroups;
  int totalNumberTimeSteps;
  int f77FileInput;
  char* timeStepInputFile;

  int ok = PyArg_ParseTuple(args, "OOOOOOdOOOOOiiis:read_timdat",
			    &pyPointerToDelt,
			    &pyPointerToAlfa,
			    &pyPointerToUtol,
			    &pyPointerToFtol,
			    &pyPointerToEtol,
			    &pyPointerToTimes,
			    &timeScaleFactor,
			    &pyPointerToMaxstp,
			    &pyPointerToMaxit,
			    &pyPointerToNtdinit,
			    &pyPointerToLgdef,
			    &pyPointerToItmax,
			    &numberTimeStepGroups,
			    &totalNumberTimeSteps,
			    &f77FileInput,
			    &timeStepInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToDelt = (double*) PyCObject_AsVoidPtr(pyPointerToDelt);
  double* pointerToAlfa = (double*) PyCObject_AsVoidPtr(pyPointerToAlfa);
  double* pointerToUtol = (double*) PyCObject_AsVoidPtr(pyPointerToUtol);
  double* pointerToFtol = (double*) PyCObject_AsVoidPtr(pyPointerToFtol);
  double* pointerToEtol = (double*) PyCObject_AsVoidPtr(pyPointerToEtol);
  double* pointerToTimes = (double*) PyCObject_AsVoidPtr(pyPointerToTimes);
  int* pointerToMaxstp = (int*) PyCObject_AsVoidPtr(pyPointerToMaxstp);
  int* pointerToMaxit = (int*) PyCObject_AsVoidPtr(pyPointerToMaxit);
  int* pointerToNtdinit = (int*) PyCObject_AsVoidPtr(pyPointerToNtdinit);
  int* pointerToLgdef = (int*) PyCObject_AsVoidPtr(pyPointerToLgdef);
  int* pointerToItmax = (int*) PyCObject_AsVoidPtr(pyPointerToItmax);

  read_timdat_f(pointerToDelt,
		pointerToAlfa,
		pointerToUtol,
		pointerToFtol,
		pointerToEtol,
		pointerToTimes,
		&timeScaleFactor,
		pointerToMaxstp,
		pointerToMaxit,
		pointerToNtdinit,
		pointerToLgdef,
		pointerToItmax,
		&numberTimeStepGroups,
		&totalNumberTimeSteps,
		&f77FileInput,
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
    << "totalNumberTimeSteps:" << totalNumberTimeSteps
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
  PyObject* pyPointerToTractionverts;
  PyObject* pyPointerToTractionvals;
  double tractionScaleFactor;
  int numberTractionBc;
  int numberSurfaceElementNodes;
  int f77FileInput;
  char* tractionInputFile;

  int ok = PyArg_ParseTuple(args, "OOdiiis:read_tractions",
			    &pyPointerToTractionverts,
			    &pyPointerToTractionvals,
			    &tractionScaleFactor,
			    &numberTractionBc,
			    &numberSurfaceElementNodes,
			    &f77FileInput,
			    &tractionInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* pointerToTractionverts = (int*) PyCObject_AsVoidPtr(pyPointerToTractionverts);
  double* pointerToTractionvals = (double*) PyCObject_AsVoidPtr(pyPointerToTractionvals);

  read_tractions_f(pointerToTractionverts,
		   pointerToTractionvals,
		   &tractionScaleFactor,
		   &numberTractionBc,
		   &numberSurfaceElementNodes,
		   &f77FileInput,
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
    << "tractionScaleFactor:" << tractionScaleFactor
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
  PyObject* pyPointerToWinkdef;
  PyObject* pyPointerToListArrayWscal;
  PyObject* pyPointerToIwinkdef;
  PyObject* pyPointerToIwinkid;
  int numberWinklerForces;
  int numberWinklerEntries;
  int f77FileInput;
  char* winklerInputFile;

  int ok = PyArg_ParseTuple(args, "OOOOiiis:read_wink",
			    &pyPointerToWinkdef,
			    &pyPointerToListArrayWscal,
			    &pyPointerToIwinkdef,
			    &pyPointerToIwinkid,
			    &numberWinklerForces,
			    &numberWinklerEntries,
			    &f77FileInput,
			    &winklerInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToWinkdef = (double*) PyCObject_AsVoidPtr(pyPointerToWinkdef);
  double* pointerToListArrayWscal = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayWscal);
  int* pointerToIwinkdef = (int*) PyCObject_AsVoidPtr(pyPointerToIwinkdef);
  int* pointerToIwinkid = (int*) PyCObject_AsVoidPtr(pyPointerToIwinkid);

  read_wink_f(pointerToWinkdef,
	      pointerToListArrayWscal,
	      pointerToIwinkdef,
	      pointerToIwinkid,
	      &numberWinklerForces,
	      &numberWinklerEntries,
	      &f77FileInput,
	      winklerInputFile,
	      &errorcode,
	      errorstring,
	      strlen(winklerInputFile),
	      sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberWinklerForces:" << numberWinklerForces
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
  PyObject* pyPointerToIen;
  PyObject* pyPointerToMat;
  PyObject* pyPointerToMaterialModelInfo;
  PyObject* pyPointerToVolumeElementFamilyList;
  PyObject* pyPointerToIvfamily;
  PyObject* pyPointerToIens;
  PyObject* pyPointerToIvftmp;
  PyObject* pyPointerToIndxiel;
  PyObject* pyPointerToIelindx;
  int numberVolumeElementNodes;
  int numberVolumeElementGaussPoints;
  int maxNumberVolumeElementFamilies;
  int numberVolumeElementFamilies;
  int prestressFlag;
  int numberVolumeElements;
  int numberNodes;

  int ok = PyArg_ParseTuple(args, "OOOOOOOOOiiiiiii:sort_elements",
			    &pyPointerToIen,
			    &pyPointerToMat,
			    &pyPointerToMaterialModelInfo,
			    &pyPointerToVolumeElementFamilyList,
			    &pyPointerToIvfamily,
			    &pyPointerToIens,
			    &pyPointerToIvftmp,
			    &pyPointerToIndxiel,
			    &pyPointerToIelindx,
			    &numberVolumeElementNodes,
			    &numberVolumeElementGaussPoints,
			    &maxNumberVolumeElementFamilies,
			    &numberVolumeElementFamilies,
			    &prestressFlag,
			    &numberVolumeElements,
			    &numberNodes);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* pointerToIen = (int*) PyCObject_AsVoidPtr(pyPointerToIen);
  int* pointerToMat = (int*) PyCObject_AsVoidPtr(pyPointerToMat);
  int* pointerToMaterialModelInfo = (int*) PyCObject_AsVoidPtr(pyPointerToMaterialModelInfo);
  int* pointerToVolumeElementFamilyList = (int*) PyCObject_AsVoidPtr(pyPointerToVolumeElementFamilyList);
  int* pointerToIvfamily = (int*) PyCObject_AsVoidPtr(pyPointerToIvfamily);
  int* pointerToIens = (int*) PyCObject_AsVoidPtr(pyPointerToIens);
  int* pointerToIvftmp = (int*) PyCObject_AsVoidPtr(pyPointerToIvftmp);
  int* pointerToIndxiel = (int*) PyCObject_AsVoidPtr(pyPointerToIndxiel);
  int* pointerToIelindx = (int*) PyCObject_AsVoidPtr(pyPointerToIelindx);
  int stateSize = 0;
  int state0Size = 0;
  int propertySize = 0;

  sort_elements_f(pointerToIen,
		  pointerToMat,
		  pointerToMaterialModelInfo,
		  pointerToVolumeElementFamilyList,
		  pointerToIvfamily,
		  pointerToIens,
		  pointerToIvftmp,
		  pointerToIndxiel,
		  pointerToIelindx,
		  &numberVolumeElementNodes,
		  &numberVolumeElementGaussPoints,
		  &maxNumberVolumeElementFamilies,
		  &numberVolumeElementFamilies,
		  &prestressFlag,
		  &numberVolumeElements,
		  &numberNodes,
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
  PyObject* pyPointerToNslip;
  PyObject* pyPointerToIndxiel;
  int numberSlipperyNodeEntries;
  int numberVolumeElements;

  int ok = PyArg_ParseTuple(args, "OOii:sort_slip_nodes",
			    &pyPointerToNslip,
			    &pyPointerToIndxiel,
			    &numberSlipperyNodeEntries,
			    &numberVolumeElements);

  if (!ok) {
    return 0;
  }

  int* pointerToNslip = (int*) PyCObject_AsVoidPtr(pyPointerToNslip);
  int* pointerToIndxiel = (int*) PyCObject_AsVoidPtr(pyPointerToIndxiel);

  sort_slip_nodes_f(pointerToNslip,
		    pointerToIndxiel,
		    &numberSlipperyNodeEntries,
		    &numberVolumeElements);

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberSlipperyNodeEntries:" << numberSlipperyNodeEntries
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
  PyObject* pyPointerToNfault;
  PyObject* pyPointerToIndxiel;
  int numberSplitNodeEntries;
  int numberVolumeElements;

  int ok = PyArg_ParseTuple(args, "OOii:sort_split_nodes",
			    &pyPointerToNfault,
			    &pyPointerToIndxiel,
			    &numberSplitNodeEntries,
			    &numberVolumeElements);

  if (!ok) {
    return 0;
  }

  int* pointerToNfault = (int*) PyCObject_AsVoidPtr(pyPointerToNfault);
  int* pointerToIndxiel = (int*) PyCObject_AsVoidPtr(pyPointerToIndxiel);

  sort_split_nodes_f(pointerToNfault,
		    pointerToIndxiel,
		    &numberSplitNodeEntries,
		    &numberVolumeElements);

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberSplitNodeEntries:" << numberSplitNodeEntries
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
  int numberGlobalEquations;
  PyObject* pyPointerToLm;
  PyObject* pyPointerToLmx;
  int numberVolumeElements;
  int totalNumberSlipperyNodes;
  int numberVolumeElementNodes;

  int ok = PyArg_ParseTuple(args, "iOOiii:cmp_stiffsz",
			    &numberGlobalEquations,
			    &pyPointerToLm,
			    &pyPointerToLmx,
			    &numberVolumeElements,
			    &totalNumberSlipperyNodes,
			    &numberVolumeElementNodes);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* pointerToLm = (int*) PyCObject_AsVoidPtr(pyPointerToLm);
  int* pointerToLmx = (int*) PyCObject_AsVoidPtr(pyPointerToLmx);
  int workingArraySize = 0;

  cmp_stiffsz_f(&numberGlobalEquations,
	   	pointerToLm,
	   	pointerToLmx,
	   	&numberVolumeElements,
	   	&workingArraySize,
	   	&totalNumberSlipperyNodes,
	   	&numberVolumeElementNodes,
	   	&errorcode,
	   	errorstring,
	   	sizeof(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }
		  
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "workingArraySize:" << workingArraySize
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i", workingArraySize);
}



// Create linked list of nonzero row and column entries in the stiffness matrix

char pypylith3d_lnklst__doc__[] = "";
char pypylith3d_lnklst__name__[] = "lnklst";

PyObject * pypylith3d_lnklst(PyObject *, PyObject *args)
{
  int numberGlobalEquations;
  PyObject* pyPointerToLm;
  PyObject* pyPointerToLmx;
  int numberVolumeElements;
  int numberVolumeElementNodes;
  int numberVolumeElementEquations;
  PyObject* pyPointerToIndx;
  PyObject* pyPointerToLink;
  PyObject* pyPointerToNbrs;
  int workingArraySize;
  int totalNumberSlipperyNodes;

  int ok = PyArg_ParseTuple(args, "iOOiiiOOOii:lnklst",
			    &numberGlobalEquations,
			    &pyPointerToLm,
			    &pyPointerToLmx,
			    &numberVolumeElements,
			    &numberVolumeElementNodes,
			    &numberVolumeElementEquations,
			    &pyPointerToIndx,
			    &pyPointerToLink,
			    &pyPointerToNbrs,
			    &workingArraySize,
			    &totalNumberSlipperyNodes);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* pointerToLm = (int*) PyCObject_AsVoidPtr(pyPointerToLm);
  int* pointerToLmx = (int*) PyCObject_AsVoidPtr(pyPointerToLmx);
  int* pointerToIndx = (int*) PyCObject_AsVoidPtr(pyPointerToIndx);
  int* pointerToLink = (int*) PyCObject_AsVoidPtr(pyPointerToLink);
  int* pointerToNbrs = (int*) PyCObject_AsVoidPtr(pyPointerToNbrs);
  int stiffnessMatrixSize =0;
  int stiffnessOffDiagonalSize =0;

  lnklst_f(&numberGlobalEquations,
	   pointerToLm,
	   pointerToLmx,
	   &numberVolumeElements,
	   &numberVolumeElementNodes,
	   &numberVolumeElementEquations,
	   pointerToIndx,
	   pointerToLink,
	   pointerToNbrs,
	   &workingArraySize,
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
  PyObject* pyPointerToId;
  int numberNodes;
  PyObject* pyPointerToIens;
  PyObject* pyPointerToLm;
  int numberVolumeElements;
  int numberVolumeElementNodes;

  int ok = PyArg_ParseTuple(args, "OiOOii:local",
			    &pyPointerToId,
			    &numberNodes,
			    &pyPointerToIens,
			    &pyPointerToLm,
			    &numberVolumeElements,
			    &numberVolumeElementNodes);

  if (!ok) {
    return 0;
  }

  int* pointerToId = (int*) PyCObject_AsVoidPtr(pyPointerToId);
  int* pointerToIens = (int*) PyCObject_AsVoidPtr(pyPointerToIens);
  int* pointerToLm = (int*) PyCObject_AsVoidPtr(pyPointerToLm);

  local_f(pointerToId,
	  &numberNodes,
	  pointerToIens,
	  pointerToLm,
	  &numberVolumeElements,
	  &numberVolumeElementNodes);

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberVolumeElements:" << numberVolumeElements
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
  PyObject* pyPointerToIens;
  PyObject* pyPointerToLmf;
  int numberVolumeElements;
  PyObject* pyPointerToNfault;
  int numberSplitNodeEntries;
  int numberVolumeElementNodes;

  int ok = PyArg_ParseTuple(args, "OOiOii:localf",
			    &pyPointerToIens,
			    &pyPointerToLmf,
			    &numberVolumeElements,
			    &pyPointerToNfault,
			    &numberSplitNodeEntries,
			    &numberVolumeElementNodes);

  if (!ok) {
    return 0;
  }

  int* pointerToIens = (int*) PyCObject_AsVoidPtr(pyPointerToIens);
  int* pointerToLmf = (int*) PyCObject_AsVoidPtr(pyPointerToLmf);
  int* pointerToNfault = (int*) PyCObject_AsVoidPtr(pyPointerToNfault);

  localf_f(pointerToIens,
	   pointerToLmf,
	   &numberVolumeElements,
	   pointerToNfault,
	   &numberSplitNodeEntries,
	   &numberVolumeElementNodes);
		  
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberSplitNodeEntries:" << numberSplitNodeEntries
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
  PyObject* pyPointerToIdx;
  int numberNodes;
  PyObject* pyPointerToIens;
  PyObject* pyPointerToLmx;
  int numberVolumeElements;
  PyObject* pyPointerToNslip;
  int numberSlipperyNodeEntries;
  int numberVolumeElementNodes;

  int ok = PyArg_ParseTuple(args, "OiOOiOii:localx",
			    &pyPointerToIdx,
  			    &numberNodes,
  			    &pyPointerToIens,
  			    &pyPointerToLmx,
  			    &numberVolumeElements,
  			    &pyPointerToNslip,
			    &numberSlipperyNodeEntries,
  			    &numberVolumeElementNodes);

  if (!ok) {
    return 0;
  }

  int* pointerToIdx = (int*) PyCObject_AsVoidPtr(pyPointerToIdx);
  int* pointerToIens = (int*) PyCObject_AsVoidPtr(pyPointerToIens);
  int* pointerToLmx = (int*) PyCObject_AsVoidPtr(pyPointerToLmx);
  int* pointerToNslip = (int*) PyCObject_AsVoidPtr(pyPointerToNslip);

  localx_f(pointerToIdx,
	   &numberNodes,
	   pointerToIens,
	   pointerToLmx,
	   &numberVolumeElements,
	   pointerToNslip,
	   &numberSlipperyNodeEntries,
	   &numberVolumeElementNodes);

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberVolumeElements:" << numberVolumeElements
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
  PyObject* pyPointerToIndx;
  PyObject* pyPointerToLink;
  PyObject* pyPointerToNbrs;
  int numberGlobalEquations;
  int stiffnessMatrixSize;
  int workingArraySize;

  int ok = PyArg_ParseTuple(args, "OOOOiii:makemsr",
                            &pyA,
			    &pyPointerToIndx,
			    &pyPointerToLink,
			    &pyPointerToNbrs,
			    &numberGlobalEquations,
			    &stiffnessMatrixSize,
			    &workingArraySize);

  if (!ok) {
    return 0;
  }
  Mat A = (Mat) PyCObject_AsVoidPtr(pyA);
  int* pointerToIndx = (int*) PyCObject_AsVoidPtr(pyPointerToIndx);
  int* pointerToLink = (int*) PyCObject_AsVoidPtr(pyPointerToLink);
  int* pointerToNbrs = (int*) PyCObject_AsVoidPtr(pyPointerToNbrs);
  int minimumNonzeroTermsPerRow = 0;
  int maximumNonzeroTermsPerRow = 0;
  double averageNonzeroTermsPerRow = 0.0;

  makemsr_f(&A,
	    pointerToIndx,
	    pointerToLink,
	    pointerToNbrs,
	    &numberGlobalEquations,
	    &stiffnessMatrixSize,
	    &workingArraySize,
	    &minimumNonzeroTermsPerRow,
	    &maximumNonzeroTermsPerRow,
	    &averageNonzeroTermsPerRow);

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "workingArraySize:" << workingArraySize
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  // return Py_None;
  return Py_BuildValue("iid",
		       minimumNonzeroTermsPerRow,
		       maximumNonzeroTermsPerRow,
		       averageNonzeroTermsPerRow);
}



// Perform the time-dependent solution.

char pypylith3d_viscos__doc__[] = "";
char pypylith3d_viscos__name__[] = "viscos";

PyObject * pypylith3d_viscos(PyObject *, PyObject *args)
{
  PyObject* pyA;                              // Sparse matrix arrays
  PyObject* pyRhs;
  PyObject* pySol;
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
  PyObject* pyPointerToIelindx;
  PyObject* pyPointerToTractionverts;         // Traction BC arrays
  PyObject* pyPointerToTractionvals;
  PyObject* pyPointerToGauss2d;
  PyObject* pyPointerToSh2d;
  PyObject* pyPointerToListArrayElementTypeInfo2d;
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

  int ok = PyArg_ParseTuple(args, "OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOsssii:viscos",
			    &pyA,                              // Sparse matrix arrays
			    &pyRhs,
			    &pySol,
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
			    &pyPointerToIelindx,
                            &pyPointerToTractionverts,         // Traction BC arrays
                            &pyPointerToTractionvals,
                            &pyPointerToGauss2d,
                            &pyPointerToSh2d,
                            &pyPointerToListArrayElementTypeInfo2d,
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
  Vec      rhs = (Vec) PyCObject_AsVoidPtr(pyRhs);
  Vec      sol = (Vec) PyCObject_AsVoidPtr(pySol);
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
  int*  pointerToIelindx = (int*) PyCObject_AsVoidPtr(pyPointerToIelindx);
  int*  pointerToTractionverts = (int*) PyCObject_AsVoidPtr(pyPointerToTractionverts);
  double*  pointerToTractionvals = (double*) PyCObject_AsVoidPtr(pyPointerToTractionvals);
  double*  pointerToGauss2d = (double*) PyCObject_AsVoidPtr(pyPointerToGauss2d);
  double*  pointerToSh2d = (double*) PyCObject_AsVoidPtr(pyPointerToSh2d);
  int*  pointerToListArrayElementTypeInfo2d = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayElementTypeInfo2d);
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


  viscos_f(&A,                            // 1: Sparse matrix arrays
	   &rhs,
	   &sol,
	   pointerToBextern,                  // Force vectors
	   pointerToBtraction,
	   pointerToBgravity,
	   pointerToBconcForce,
	   pointerToBintern,
	   pointerToBresid,
	   pointerToBwink,                    // 10
	   pointerToBwinkx,
	   pointerToDispVec,
	   pointerToDprev,
	   pointerToListArrayNforce,
	   pointerToListArrayGrav,
	   pointerToX,                        // Global arrays
	   pointerToD,
	   pointerToDeld,
	   pointerToDcur,
	   pointerToId,                       // 20
	   pointerToIwink,
	   pointerToWink,
	   pointerToListArrayNsysdat,
	   pointerToListArrayIddmat,
	   pointerToIbond,                    // Boundary condition arrays
	   pointerToBond,
	   pointerToDx,                       // Slippery node arrays
	   pointerToDeldx,
	   pointerToDxcur,
	   pointerToDiforc,                   // 30
	   pointerToIdx,
	   pointerToIwinkx,
	   pointerToWinkx,
	   pointerToIdslp,
	   pointerToIpslp,
	   pointerToIdhist,
	   pointerToFault,                    // Split node arrays
	   pointerToNfault,
	   pointerToDfault,
	   pointerToTfault,                   // 40
	   pointerToS,                        // Local stiffness matrix arrays
	   pointerToStemp,
	   pointerToState,                    // Element arrays
	   pointerToDstate,
	   pointerToState0,
	   pointerToDmat,
	   pointerToIen,
	   pointerToLm,
	   pointerToLmx,
	   pointerToLmf,                      // 50
	   pointerToIvfamily,
	   pointerToListArrayNpar,
	   pointerToIelindx,
	   pointerToTractionverts,            // Traction BC arrays
	   pointerToTractionvals,
	   pointerToGauss2d,
	   pointerToSh2d,
	   pointerToListArrayElementTypeInfo2d,
	   pointerToListArrayPropertyList,    // Material property arrays
	   pointerToMaterialModelInfo,        // 60
	   pointerToGauss,                    // Element type arrays
	   pointerToSh,
	   pointerToShj,
	   pointerToListArrayElementTypeInfo,
	   pointerToHistry,                   // Time information
	   pointerToListArrayRtimdat,
	   pointerToListArrayNtimdat,
	   pointerToListArrayNvisdat,
	   pointerToMaxstp,
	   pointerToDelt,                     // 70
	   pointerToAlfa,
	   pointerToMaxit,
	   pointerToNtdinit,
	   pointerToLgdef,
	   pointerToUtol,
	   pointerToFtol,
	   pointerToEtol,
	   pointerToItmax,
	   pointerToListArrayRgiter,          // Iterative solution information
	   pointerToSkew,                     // 80: Skew rotation information
	   pointerToIprint,                   // Input/output information
	   pointerToListArrayNcodat,
	   pointerToListArrayNunits,
	   pointerToListArrayNprint,
	   pointerToIstatout,
	   pointerToNstatout,
	   asciiOutputFile,                   // Output file names
	   plotOutputFile,
	   ucdOutputRoot,
	   &viscousStage,                     // 90: PETSC logging
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

// Perform the time-dependent solution.

char pypylith3d_viscos_step__doc__[] = "";
char pypylith3d_viscos_step__name__[] = "viscos_step";

PyObject * pypylith3d_viscos_step(PyObject *, PyObject *args)
{
  PyObject* pyA;                              // Sparse matrix arrays
  PyObject* pyRhs;
  PyObject* pySol;
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
  PyObject* pyPointerToIelindx;
  PyObject* pyPointerToTractionverts;         // Traction BC arrays
  PyObject* pyPointerToTractionvals;
  PyObject* pyPointerToGauss2d;
  PyObject* pyPointerToSh2d;
  PyObject* pyPointerToListArrayElementTypeInfo2d;
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
  int ntot, ltim, indexx, jcyc, i, j, skc, nfirst, nstep, lgdefp;
  double time, deltp;
  PyObject *pyGtol;

  int ok = PyArg_ParseTuple(args, "OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOsssiiiiiiiiiiiddiO!:viscos_step",
			    &pyA,                              // Sparse matrix arrays
			    &pyRhs,
			    &pySol,
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
			    &pyPointerToIelindx,
                            &pyPointerToTractionverts,         // Traction BC arrays
                            &pyPointerToTractionvals,
                            &pyPointerToGauss2d,
                            &pyPointerToSh2d,
                            &pyPointerToListArrayElementTypeInfo2d,
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
  int*  pointerToIelindx = (int*) PyCObject_AsVoidPtr(pyPointerToIelindx);
  int*  pointerToTractionverts = (int*) PyCObject_AsVoidPtr(pyPointerToTractionverts);
  double*  pointerToTractionvals = (double*) PyCObject_AsVoidPtr(pyPointerToTractionvals);
  double*  pointerToGauss2d = (double*) PyCObject_AsVoidPtr(pyPointerToGauss2d);
  double*  pointerToSh2d = (double*) PyCObject_AsVoidPtr(pyPointerToSh2d);
  int*  pointerToListArrayElementTypeInfo2d = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayElementTypeInfo2d);
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
	   pointerToBextern,                  // Force vectors
	   pointerToBtraction,
	   pointerToBgravity,
	   pointerToBconcForce,
	   pointerToBintern,
	   pointerToBresid,
	   pointerToBwink,                    // 10
	   pointerToBwinkx,
	   pointerToDispVec,
	   pointerToDprev,
	   pointerToListArrayNforce,
	   pointerToListArrayGrav,
	   pointerToX,                        // Global arrays
	   pointerToD,
	   pointerToDeld,
	   pointerToDcur,
	   pointerToId,                       // 20
	   pointerToIwink,
	   pointerToWink,
	   pointerToListArrayNsysdat,
	   pointerToListArrayIddmat,
	   pointerToIbond,                    // Boundary condition arrays
	   pointerToBond,
	   pointerToDx,                       // Slippery node arrays
	   pointerToDeldx,
	   pointerToDxcur,
	   pointerToDiforc,                   // 30
	   pointerToIdx,
	   pointerToIwinkx,
	   pointerToWinkx,
	   pointerToIdslp,
	   pointerToIpslp,
	   pointerToIdhist,
	   pointerToFault,                    // Split node arrays
	   pointerToNfault,
	   pointerToDfault,
	   pointerToTfault,                   // 40
	   pointerToS,                        // Local stiffness matrix arrays
	   pointerToStemp,
	   pointerToState,                    // Element arrays
	   pointerToDstate,
	   pointerToState0,
	   pointerToDmat,
	   pointerToIen,
	   pointerToLm,
	   pointerToLmx,
	   pointerToLmf,                      // 50
	   pointerToIvfamily,
	   pointerToListArrayNpar,
		pointerToIelindx,
		pointerToTractionverts,            // Traction BC arrays
	   pointerToTractionvals,
	   pointerToGauss2d,
	   pointerToSh2d,
	   pointerToListArrayElementTypeInfo2d,
	   pointerToListArrayPropertyList,    // Material property arrays
	   pointerToMaterialModelInfo,        // 60
	   pointerToGauss,                    // Element type arrays
	   pointerToSh,
	   pointerToShj,
	   pointerToListArrayElementTypeInfo,
	   pointerToHistry,                   // Time information
	   pointerToListArrayRtimdat,
	   pointerToListArrayNtimdat,
	   pointerToListArrayNvisdat,
	   pointerToMaxstp,
	   pointerToDelt,                     // 70
	   pointerToAlfa,
	   pointerToMaxit,
	   pointerToNtdinit,
	   pointerToLgdef,
	   pointerToUtol,
	   pointerToFtol,
	   pointerToEtol,
	   pointerToItmax,
	   pointerToListArrayRgiter,          // Iterative solution information
	   pointerToSkew,                     // 80: Skew rotation information
	   pointerToIprint,                   // Input/output information
	   pointerToListArrayNcodat,
	   pointerToListArrayNunits,
	   pointerToListArrayNprint,
	   pointerToIstatout,
	   pointerToNstatout,
	   asciiOutputFile,                   // Output file names
	   plotOutputFile,
	   ucdOutputRoot,
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
	   strlen(asciiOutputFile),           // String lengths
	   strlen(plotOutputFile),
	   strlen(ucdOutputRoot),
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
  PyObject* pyPointerToListArrayNprint;
  PyObject* pyPointerToListArrayNunits;
  char* asciiOutputFile;                      // Output file names
  char* plotOutputFile;
  int viscousStage;

  int ok = PyArg_ParseTuple(args, "OOssi:viscos_setup",
			    &pyPointerToListArrayNprint,
			    &pyPointerToListArrayNunits,
			    &asciiOutputFile,                  // Output file names
                &plotOutputFile,
                &viscousStage);

  if (!ok) {
    return 0;
  }

  int  errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int *pointerToListArrayNprint = (int *) PyCObject_AsVoidPtr(pyPointerToListArrayNprint);
  int *pointerToListArrayNunits = (int *) PyCObject_AsVoidPtr(pyPointerToListArrayNunits);

  viscos_setup_f(pointerToListArrayNprint,
                 pointerToListArrayNunits,
                 asciiOutputFile,                   // Output file names
                 plotOutputFile,
                 &viscousStage,                     // PETSC logging
                 &errorcode,                        // Error codes
                 errorstring,
                 strlen(asciiOutputFile),           // String lengths
                 strlen(plotOutputFile),
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

// Cleanup the time-dependent solution.

char pypylith3d_viscos_cleanup__doc__[] = "";
char pypylith3d_viscos_cleanup__name__[] = "viscos_cleanup";

PyObject * pypylith3d_viscos_cleanup(PyObject *, PyObject *args)
{
  PyObject* pyPointerToListArrayNtimdat;
  PyObject* pyPointerToListArrayNprint;
  PyObject* pyPointerToListArrayNunits;

  int ok = PyArg_ParseTuple(args, "OOO:viscos_cleanup",&pyPointerToListArrayNtimdat, &pyPointerToListArrayNprint, &pyPointerToListArrayNunits);

  if (!ok) {
    return 0;
  }

  int  errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int *pointerToListArrayNtimdat = (int *) PyCObject_AsVoidPtr(pyPointerToListArrayNtimdat);
  int *pointerToListArrayNprint = (int *) PyCObject_AsVoidPtr(pyPointerToListArrayNprint);
  int *pointerToListArrayNunits  = (int *) PyCObject_AsVoidPtr(pyPointerToListArrayNunits);

  viscos_cleanup_f(pointerToListArrayNtimdat,
                   pointerToListArrayNprint,
                   pointerToListArrayNunits,
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
  PyObject* pyPointerToBond;
  PyObject* pyPointerToIbond;
  int numberNodes;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOiiis:write_bc",
			    &pyPointerToBond,
			    &pyPointerToIbond,
			    &numberNodes,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToBond = (double*) PyCObject_AsVoidPtr(pyPointerToBond);
  int* pointerToIbond = (int*) PyCObject_AsVoidPtr(pyPointerToIbond);

  write_bc_f(pointerToBond,
	     pointerToIbond,
	     &numberNodes,
	     &f77AsciiOutput,
	     &asciiOutputInt,
	     asciiOutputFile,
	     &errorcode,
	     errorstring,
	     strlen(asciiOutputFile),
	     sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberNodes:" << numberNodes
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
  PyObject* pyPointerToIens;
  PyObject* pyPointerToIvfamily;
  PyObject* pyPointerToIndxiel;
  int numberVolumeElementNodes;
  int numberVolumeElementGaussPoints;
  int numberVolumeElements;
  int volumeElementType;
  int numberVolumeElementFamilies;
  int f77AsciiOutput;
  int f77PlotOutput;
  int asciiOutputInt;
  int plotOutputInt;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "OOOiiiiiiiiiss:write_connect",
                            &pyPointerToIens,
                            &pyPointerToIvfamily,
                            &pyPointerToIndxiel,
                            &numberVolumeElementNodes,
                            &numberVolumeElementGaussPoints,
                            &numberVolumeElements,
                            &volumeElementType,
                            &numberVolumeElementFamilies,
                            &f77AsciiOutput,
                            &f77PlotOutput,
                            &asciiOutputInt,
                            &plotOutputInt,
                            &asciiOutputFile,
                            &plotOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* pointerToIens = (int*) PyCObject_AsVoidPtr(pyPointerToIens);
  int* pointerToIvfamily = (int*) PyCObject_AsVoidPtr(pyPointerToIvfamily);
  int* pointerToIndxiel = (int*) PyCObject_AsVoidPtr(pyPointerToIndxiel);

  write_connect_f(pointerToIens,
		  pointerToIvfamily,
		  pointerToIndxiel,
		  &numberVolumeElementNodes,
		  &numberVolumeElementGaussPoints,
		  &numberVolumeElements,
		  &volumeElementType,
		  &numberVolumeElementFamilies,
		  &f77AsciiOutput,
		  &f77PlotOutput,
		  &asciiOutputInt,
		  &plotOutputInt,
		  asciiOutputFile,
		  plotOutputFile,
		  &errorcode,
		  errorstring,
		  strlen(asciiOutputFile),
		  strlen(plotOutputFile),
		  sizeof(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberVolumeElements:" << numberVolumeElements
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
  PyObject* pyPointerToX;
  int numberNodes;
  int f77AsciiOutput;
  int f77PlotOutput;
  int asciiOutputInt;
  int plotOutputInt;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "Oiiiiiss:write_coords",
			    &pyPointerToX,
			    &numberNodes,
			    &f77AsciiOutput,
			    &f77PlotOutput,
			    &asciiOutputInt,
			    &plotOutputInt,
			    &asciiOutputFile,
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToX = (double*) PyCObject_AsVoidPtr(pyPointerToX);

  write_coords_f(pointerToX,
		 &numberNodes,
		 &f77AsciiOutput,
		 &f77PlotOutput,
		 &asciiOutputInt,
		 &plotOutputInt,
		 asciiOutputFile,
		 plotOutputFile,
		 &errorcode,
		 errorstring,
		 strlen(asciiOutputFile),
		 strlen(plotOutputFile),
		 sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberNodes:" << numberNodes
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
  PyObject* pyPointerToDiforc;
  PyObject* pyPointerToNslip;
  PyObject* pyPointerToIdhist;
  int numberSlipperyNodeEntries;
  int numberDifferentialForceEntries;
  int numberNodes;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOOiiiiis:write_diff",
			    &pyPointerToDiforc,
			    &pyPointerToNslip,
			    &pyPointerToIdhist,
			    &numberSlipperyNodeEntries,
			    &numberDifferentialForceEntries,
			    &numberNodes,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToDiforc = (double*) PyCObject_AsVoidPtr(pyPointerToDiforc);
  int* pointerToNslip = (int*) PyCObject_AsVoidPtr(pyPointerToNslip);
  int* pointerToIdhist = (int*) PyCObject_AsVoidPtr(pyPointerToIdhist);

  write_diff_f(pointerToDiforc,
	       pointerToNslip,
	       pointerToIdhist,
	       &numberSlipperyNodeEntries,
	       &numberDifferentialForceEntries,
	       &numberNodes,
	       &f77AsciiOutput,
	       &asciiOutputInt,
	       asciiOutputFile,
	       &errorcode,
	       errorstring,
	       strlen(asciiOutputFile),
	       sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberDifferentialForceEntries:" << numberDifferentialForceEntries
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
  int numberVolumeElements;
  int numberVolumeElementNodes;
  int numberVolumeElementGaussPoints;
  int volumeElementType;
  int quadratureOrderInt;
  int prestressAutoComputeInt;
  int prestressAutoChangeElasticPropsInt;
  double prestressAutoComputePoisson;
  double prestressAutoComputeYoungs;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "iiiiiiiddiis:write_element_info",
			    &numberVolumeElements,
			    &numberVolumeElementNodes,
			    &numberVolumeElementGaussPoints,
			    &volumeElementType,
			    &quadratureOrderInt,
			    &prestressAutoComputeInt,
			    &prestressAutoChangeElasticPropsInt,
			    &prestressAutoComputePoisson,
			    &prestressAutoComputeYoungs,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  write_element_info_f(&numberVolumeElements,
		       &numberVolumeElementNodes,
		       &numberVolumeElementGaussPoints,
		       &volumeElementType,
		       &quadratureOrderInt,
		       &prestressAutoComputeInt,
		       &prestressAutoChangeElasticPropsInt,
		       &prestressAutoComputePoisson,
		       &prestressAutoComputeYoungs,
		       &f77AsciiOutput,
		       &asciiOutputInt,
		       asciiOutputFile,strlen(asciiOutputFile));

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberVolumeElements:" << numberVolumeElements
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
  PyObject* pyPointerToIprint;
  int numberFullOutputs;
  int analysisTypeInt;
  int numberCycles;
  int totalNumberTimeSteps;
  int f77AsciiOutput;
  int f77PlotOutput;
  int asciiOutputInt;
  int plotOutputInt;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "Oiiiiiiiiss:write_fuldat",
			    &pyPointerToIprint,
			    &numberFullOutputs,
			    &analysisTypeInt,
			    &numberCycles,
			    &totalNumberTimeSteps,
			    &f77AsciiOutput,
			    &f77PlotOutput,
			    &asciiOutputInt,
			    &plotOutputInt,
			    &asciiOutputFile,
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* pointerToIprint = (int*) PyCObject_AsVoidPtr(pyPointerToIprint);

  write_fuldat_f(pointerToIprint,
		 &numberFullOutputs,
		 &analysisTypeInt,
		 &numberCycles,
		 &totalNumberTimeSteps,
		 &f77AsciiOutput,
		 &f77PlotOutput,
		 &asciiOutputInt,
		 &plotOutputInt,
		 asciiOutputFile,
		 plotOutputFile,
		 &errorcode,
		 errorstring,
		 strlen(asciiOutputFile),
		 strlen(plotOutputFile),
		 sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberFullOutputs:" << numberFullOutputs
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
  int asciiOutputInt;
  int plotOutputInt;
  int numberNodes;
  int analysisTypeInt;
  int debuggingOutputInt;
  int f77AsciiOutput;
  int f77PlotOutput;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "siiiiiiiss:write_global_info",
			    &title,
			    &asciiOutputInt,
			    &plotOutputInt,
			    &numberNodes,
			    &analysisTypeInt,
			    &debuggingOutputInt,
			    &f77AsciiOutput,
			    &f77PlotOutput,
			    &asciiOutputFile,
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  write_global_info_f(title,
		      &asciiOutputInt,
		      &plotOutputInt,
		      &numberNodes,
		      &analysisTypeInt,
		      &debuggingOutputInt,
		      &f77AsciiOutput,
		      &f77PlotOutput,
		      asciiOutputFile,
		      plotOutputFile,
		      strlen(title),
		      strlen(asciiOutputFile),
		      strlen(plotOutputFile));

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "debuggingOutputInt:" << debuggingOutputInt
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
  PyObject* pyPointerToHistry;
  PyObject* pyPointerToTimes;
  int numberLoadHistories;
  int totalNumberTimeSteps;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOiiiis:write_hist",
			    &pyPointerToHistry,
			    &pyPointerToTimes,
			    &numberLoadHistories,
			    &totalNumberTimeSteps,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToHistry = (double*) PyCObject_AsVoidPtr(pyPointerToHistry);
  double* pointerToTimes = (double*) PyCObject_AsVoidPtr(pyPointerToTimes);

  write_hist_f(pointerToHistry,
	       pointerToTimes,
	       &numberLoadHistories,
	       &totalNumberTimeSteps,
	       &f77AsciiOutput,
	       &asciiOutputInt,
	       asciiOutputFile,
	       &errorcode,
	       errorstring,
	       strlen(asciiOutputFile),
	       sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberLoadHistories:" << numberLoadHistories
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
  PyObject* pyPointerToListArrayPropertyList;
  PyObject* pyPointerToListArrayGrav;
  PyObject* pyPointerToIvfamily;
  PyObject* pyPointerToMaterialModelInfo;
  int numberVolumeElementFamilies;
  int propertySize;
  int asciiOutputInt;
  int plotOutputInt;
  int f77AsciiOutput;
  int f77PlotOutput;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "OOOOiiiiiiss:write_props",
			    &pyPointerToListArrayPropertyList,
			    &pyPointerToListArrayGrav,
			    &pyPointerToIvfamily,
			    &pyPointerToMaterialModelInfo,
			    &numberVolumeElementFamilies,
			    &propertySize,
			    &asciiOutputInt,
			    &plotOutputInt,
			    &f77AsciiOutput,
			    &f77PlotOutput,
			    &asciiOutputFile,
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToListArrayPropertyList = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayPropertyList);
  double* pointerToListArrayGrav = (double*) PyCObject_AsVoidPtr(pyPointerToListArrayGrav);
  int* pointerToIvfamily = (int*) PyCObject_AsVoidPtr(pyPointerToIvfamily);
  int* pointerToMaterialModelInfo = (int*) PyCObject_AsVoidPtr(pyPointerToMaterialModelInfo);

  write_props_f(pointerToListArrayPropertyList,
		pointerToListArrayGrav,
		pointerToIvfamily,
		pointerToMaterialModelInfo,
		&numberVolumeElementFamilies,
		&propertySize,
		&asciiOutputInt,
		&plotOutputInt,
		&f77AsciiOutput,
		&f77PlotOutput,
		asciiOutputFile,
		plotOutputFile,
		&errorcode,
		errorstring,
		strlen(asciiOutputFile),
		strlen(plotOutputFile),
		sizeof(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberVolumeElementFamilies:" << numberVolumeElementFamilies
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
  PyObject* pyPointerToSkew;
  int numberRotationEntries;
  int autoRotateSlipperyNodesInt;
  int numberNodes;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "Oiiiiis:write_skew",
			    &pyPointerToSkew,
			    &numberRotationEntries,
			    &autoRotateSlipperyNodesInt,
			    &numberNodes,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToSkew = (double*) PyCObject_AsVoidPtr(pyPointerToSkew);

  write_skew_f(pointerToSkew,
	       &numberRotationEntries,
	       &autoRotateSlipperyNodesInt,
	       &numberNodes,
	       &f77AsciiOutput,
	       &asciiOutputInt,
	       asciiOutputFile,
	       &errorcode,
	       errorstring,
	       strlen(asciiOutputFile),
	       sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberRotationEntries:" << numberRotationEntries
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
  PyObject* pyPointerToNslip;
  int numberSlipperyNodeEntries;
  int totalNumberSlipperyNodes;
  int f77AsciiOutput;
  int f77PlotOutput;
  int asciiOutputInt;
  int plotOutputInt;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "Oiiiiiiss:write_slip",
			    &pyPointerToNslip,
			    &numberSlipperyNodeEntries,
			    &totalNumberSlipperyNodes,
			    &f77AsciiOutput,
			    &f77PlotOutput,
			    &asciiOutputInt,
			    &plotOutputInt,
			    &asciiOutputFile,
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* pointerToNslip = (int*) PyCObject_AsVoidPtr(pyPointerToNslip);

  write_slip_f(pointerToNslip,
	       &numberSlipperyNodeEntries,
	       &totalNumberSlipperyNodes,
	       &f77AsciiOutput,
	       &f77PlotOutput,
	       &asciiOutputInt,
	       &plotOutputInt,
	       asciiOutputFile,
	       plotOutputFile,
	       &errorcode,
	       errorstring,
	       strlen(asciiOutputFile),
	       strlen(plotOutputFile),
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
  int numberGlobalEquations;
  int stiffnessMatrixSize;
  int minimumNonzeroTermsPerRow;
  int maximumNonzeroTermsPerRow;
  double averageNonzeroTermsPerRow;
  int asciiOutputInt;
  int f77AsciiOutput;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "iiiidiis:write_sparse_info",
			    &numberGlobalEquations,
			    &stiffnessMatrixSize,
			    &minimumNonzeroTermsPerRow,
			    &maximumNonzeroTermsPerRow,
			    &averageNonzeroTermsPerRow,
			    &asciiOutputInt,
			    &f77AsciiOutput,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  write_sparse_info_f(&numberGlobalEquations,
		      &stiffnessMatrixSize,
		      &minimumNonzeroTermsPerRow,
		      &maximumNonzeroTermsPerRow,
		      &averageNonzeroTermsPerRow,
		      &asciiOutputInt,
		      &f77AsciiOutput,
		      asciiOutputFile,strlen(asciiOutputFile));
		  
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "minimumNonzeroTermsPerRow:" << minimumNonzeroTermsPerRow
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
  PyObject* pyPointerToFault;
  PyObject* pyPointerToNfault;
  int numberSplitNodeEntries;
  int f77AsciiOutput;
  int f77PlotOutput;
  int asciiOutputInt;
  int plotOutputInt;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "OOiiiiiss:write_split",
			    &pyPointerToFault,
			    &pyPointerToNfault,
			    &numberSplitNodeEntries,
			    &f77AsciiOutput,
			    &f77PlotOutput,
			    &asciiOutputInt,
			    &plotOutputInt,
			    &asciiOutputFile,
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToFault = (double*) PyCObject_AsVoidPtr(pyPointerToFault);
  int* pointerToNfault = (int*) PyCObject_AsVoidPtr(pyPointerToNfault);

  write_split_f(pointerToFault,
		pointerToNfault,
		&numberSplitNodeEntries,
		&f77AsciiOutput,
		&f77PlotOutput,
		&asciiOutputInt,
		&plotOutputInt,
		asciiOutputFile,
		plotOutputFile,
		&errorcode,
		errorstring,
		strlen(asciiOutputFile),
		strlen(plotOutputFile),
		sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberSplitNodeEntries:" << numberSplitNodeEntries
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
  PyObject* pyPointerToIdftn;
  int totalNumberSplitNodes;
  int f77PlotOutput;
  int plotOutputInt;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "Oiiis:write_split_plot",
			    &pyPointerToIdftn,
			    &totalNumberSplitNodes,
			    &f77PlotOutput,
			    &plotOutputInt,
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  int* pointerToIdftn = (int*) PyCObject_AsVoidPtr(pyPointerToIdftn);

  write_split_plot_f(pointerToIdftn,
		     &totalNumberSplitNodes,
		     &f77PlotOutput,
		     &plotOutputInt,
		     plotOutputFile,
		     strlen(plotOutputFile));
    

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
  PyObject* pyPointerToIstatout;
  PyObject* pyPointerToNstatout;
  int f77AsciiOutput;
  int f77PlotOutput;
  int asciiOutputInt;
  int plotOutputInt;
  char* asciiOutputFile;
  char* plotOutputFile;

  int ok = PyArg_ParseTuple(args, "OOiiiiss:write_stateout",
			    &pyPointerToIstatout,
			    &pyPointerToNstatout,
			    &f77AsciiOutput,
			    &f77PlotOutput,
			    &asciiOutputInt,
			    &plotOutputInt,
			    &asciiOutputFile,
			    &plotOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* pointerToIstatout = (int*) PyCObject_AsVoidPtr(pyPointerToIstatout);
  int* pointerToNstatout = (int*) PyCObject_AsVoidPtr(pyPointerToNstatout);

  write_stateout_f(pointerToIstatout,
		   pointerToNstatout,
		   &f77AsciiOutput,
		   &f77PlotOutput,
		   &asciiOutputInt,
		   &plotOutputInt,
		   asciiOutputFile,
		   plotOutputFile,
		   &errorcode,
		   errorstring,
		   strlen(asciiOutputFile),
		   strlen(plotOutputFile),
		   sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "pointerToIstatout:" << pointerToIstatout
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
  double stressTolerance;
  double minimumStrainPerturbation;
  double initialStrainPerturbation;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "dddiis:write_strscomp",
			    &stressTolerance,
			    &minimumStrainPerturbation,
			    &initialStrainPerturbation,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  write_strscomp_f(&stressTolerance,
		   &minimumStrainPerturbation,
		   &initialStrainPerturbation,
		   &f77AsciiOutput,
		   &asciiOutputInt,
		   asciiOutputFile,strlen(asciiOutputFile));

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "initialStrainPerturbation:" << initialStrainPerturbation
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
  int usePreviousDisplacementFlag;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "iiis:write_subiter",
			    &usePreviousDisplacementFlag,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }


  write_subiter_f(&usePreviousDisplacementFlag,
		  &f77AsciiOutput,
		  &asciiOutputInt,
		  asciiOutputFile,strlen(asciiOutputFile));
		  
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "usePreviousDisplacementFlag:" << usePreviousDisplacementFlag
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
  PyObject* pyPointerToDelt;
  PyObject* pyPointerToAlfa;
  PyObject* pyPointerToUtol;
  PyObject* pyPointerToFtol;
  PyObject* pyPointerToEtol;
  PyObject* pyPointerToTimes;
  PyObject* pyPointerToMaxstp;
  PyObject* pyPointerToMaxit;
  PyObject* pyPointerToNtdinit;
  PyObject* pyPointerToLgdef;
  PyObject* pyPointerToItmax;
  int numberTimeStepGroups;
  int totalNumberTimeSteps;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOOOOOOOOOOiiiis:write_timdat",
			    &pyPointerToDelt,
			    &pyPointerToAlfa,
			    &pyPointerToUtol,
			    &pyPointerToFtol,
			    &pyPointerToEtol,
			    &pyPointerToTimes,
			    &pyPointerToMaxstp,
			    &pyPointerToMaxit,
			    &pyPointerToNtdinit,
			    &pyPointerToLgdef,
			    &pyPointerToItmax,
			    &numberTimeStepGroups,
			    &totalNumberTimeSteps,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToDelt = (double*) PyCObject_AsVoidPtr(pyPointerToDelt);
  double* pointerToAlfa = (double*) PyCObject_AsVoidPtr(pyPointerToAlfa);
  double* pointerToUtol = (double*) PyCObject_AsVoidPtr(pyPointerToUtol);
  double* pointerToFtol = (double*) PyCObject_AsVoidPtr(pyPointerToFtol);
  double* pointerToEtol = (double*) PyCObject_AsVoidPtr(pyPointerToEtol);
  double* pointerToTimes = (double*) PyCObject_AsVoidPtr(pyPointerToTimes);
  int* pointerToMaxstp = (int*) PyCObject_AsVoidPtr(pyPointerToMaxstp);
  int* pointerToMaxit = (int*) PyCObject_AsVoidPtr(pyPointerToMaxit);
  int* pointerToNtdinit = (int*) PyCObject_AsVoidPtr(pyPointerToNtdinit);
  int* pointerToLgdef = (int*) PyCObject_AsVoidPtr(pyPointerToLgdef);
  int* pointerToItmax = (int*) PyCObject_AsVoidPtr(pyPointerToItmax);

  write_timdat_f(pointerToDelt,
		 pointerToAlfa,
		 pointerToUtol,
		 pointerToFtol,
		 pointerToEtol,
		 pointerToTimes,
		 pointerToMaxstp,
		 pointerToMaxit,
		 pointerToNtdinit,
		 pointerToLgdef,
		 pointerToItmax,
		 &numberTimeStepGroups,
		 &totalNumberTimeSteps,
		 &f77AsciiOutput,
		 &asciiOutputInt,
		 asciiOutputFile,
		 &errorcode,
		 errorstring,
		 strlen(asciiOutputFile),
		 sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "totalNumberTimeSteps:" << totalNumberTimeSteps
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
  PyObject* pyPointerToTractionverts;
  PyObject* pyPointerToTractionvals;
  int numberTractionBc;
  int numberSurfaceElementNodes;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOiiiis:write_tractions",
			    &pyPointerToTractionverts,
			    &pyPointerToTractionvals,
			    &numberTractionBc,
			    &numberSurfaceElementNodes,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int* pointerToTractionverts = (int*) PyCObject_AsVoidPtr(pyPointerToTractionverts);
  double* pointerToTractionvals = (double*) PyCObject_AsVoidPtr(pyPointerToTractionvals);

  write_tractions_f(pointerToTractionverts,
		    pointerToTractionvals,
		    &numberTractionBc,
		    &numberSurfaceElementNodes,
		    &f77AsciiOutput,
		    &asciiOutputInt,
		    asciiOutputFile,
		    &errorcode,
		    errorstring,
		    strlen(asciiOutputFile),
		    sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberTractionBc:" << numberTractionBc
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
  PyObject* pyPointerToX;
  int numberNodes;
  PyObject* pyPointerToIens;
  PyObject* pyPointerToIvfamily;
  int numberVolumeElements;
  int numberVolumeElementFamilies;
  PyObject* pyPointerToSh;
  int numberVolumeElementNodes;
  int numberVolumeElementGaussPoints;
  int volumeElementType;
  PyObject* pyPointerToIstatout;
  PyObject* pyPointerToNstatout;
  int f77UcdOutput;
  int ucdOutputInt;
  char* ucdOutputRoot;

  int ok = PyArg_ParseTuple(args, "OiOOiiOiiiOOiis:write_ucd_mesh",
			    &pyPointerToX,
			    &numberNodes,
			    &pyPointerToIens,
			    &pyPointerToIvfamily,
			    &numberVolumeElements,
			    &numberVolumeElementFamilies,
			    &pyPointerToSh,
			    &numberVolumeElementNodes,
			    &numberVolumeElementGaussPoints,
			    &volumeElementType,
			    &pyPointerToIstatout,
			    &pyPointerToNstatout,
			    &f77UcdOutput,
			    &ucdOutputInt,
			    &ucdOutputRoot);

  if (!ok) {
    return 0;
  }

  double* pointerToX = (double*) PyCObject_AsVoidPtr(pyPointerToX);
  int* pointerToIens = (int*) PyCObject_AsVoidPtr(pyPointerToIens);
  int* pointerToIvfamily = (int*) PyCObject_AsVoidPtr(pyPointerToIvfamily);
  double* pointerToSh = (double*) PyCObject_AsVoidPtr(pyPointerToSh);
  int* pointerToIstatout = (int*) PyCObject_AsVoidPtr(pyPointerToIstatout);
  int* pointerToNstatout = (int*) PyCObject_AsVoidPtr(pyPointerToNstatout);

  write_ucd_mesh_f(pointerToX,
		   &numberNodes,
		   pointerToIens,
		   pointerToIvfamily,
		   &numberVolumeElements,
		   &numberVolumeElementFamilies,
		   pointerToSh,
		   &numberVolumeElementNodes,
		   &numberVolumeElementGaussPoints,
		   &volumeElementType,
		   pointerToIstatout,
		   pointerToNstatout,
		   &f77UcdOutput,
		   &ucdOutputInt,
		   ucdOutputRoot,strlen(ucdOutputRoot));
		  
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberVolumeElements:" << numberVolumeElements
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
  PyObject* pyPointerToWinkdef;
  PyObject* pyPointerToIwinkdef;
  PyObject* pyPointerToIwinkid;
  int numberWinklerEntries;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOOiiis:write_wink",
			    &pyPointerToWinkdef,
			    &pyPointerToIwinkdef,
			    &pyPointerToIwinkid,
			    &numberWinklerEntries,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToWinkdef = (double*) PyCObject_AsVoidPtr(pyPointerToWinkdef);
  int* pointerToIwinkdef = (int*) PyCObject_AsVoidPtr(pyPointerToIwinkdef);
  int* pointerToIwinkid = (int*) PyCObject_AsVoidPtr(pyPointerToIwinkid);

  write_wink_f(pointerToWinkdef,
	       pointerToIwinkdef,
	       pointerToIwinkid,
	       &numberWinklerEntries,
	       &f77AsciiOutput,
	       &asciiOutputInt,
	       asciiOutputFile,
	       &errorcode,
	       errorstring,
	       strlen(asciiOutputFile),
	       sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberWinklerEntries:" << numberWinklerEntries
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
  PyObject* pyPointerToWinkxdef;
  PyObject* pyPointerToIwinkxdef;
  PyObject* pyPointerToIwinkxid;
  int numberSlipperyWinklerEntries;
  int f77AsciiOutput;
  int asciiOutputInt;
  char* asciiOutputFile;

  int ok = PyArg_ParseTuple(args, "OOOiiis:write_winkx",
			    &pyPointerToWinkxdef,
			    &pyPointerToIwinkxdef,
			    &pyPointerToIwinkxid,
			    &numberSlipperyWinklerEntries,
			    &f77AsciiOutput,
			    &asciiOutputInt,
			    &asciiOutputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  double* pointerToWinkxdef = (double*) PyCObject_AsVoidPtr(pyPointerToWinkxdef);
  int* pointerToIwinkxdef = (int*) PyCObject_AsVoidPtr(pyPointerToIwinkxdef);
  int* pointerToIwinkxid = (int*) PyCObject_AsVoidPtr(pyPointerToIwinkxid);

  write_winkx_f(pointerToWinkxdef,
		pointerToIwinkxdef,
		pointerToIwinkxid,
		&numberSlipperyWinklerEntries,
		&f77AsciiOutput,
		&asciiOutputInt,
		asciiOutputFile,
		&errorcode,
		errorstring,
		strlen(asciiOutputFile),
		sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberSlipperyWinklerEntries:" << numberSlipperyWinklerEntries
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
  int f77FileInput;
  char* displacementUnits;
  char* velocityUnits;
  char* forceUnits;
  char* bcInputFile;

  int ok = PyArg_ParseTuple(args, (char *) "issss:scan_bc",
			    &f77FileInput,
			    &displacementUnits,
			    &velocityUnits,
			    &forceUnits,
			    &bcInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int numberBcEntries = 0;

  scan_bc_f(&numberBcEntries,
	    &f77FileInput,
	    displacementUnits,
	    velocityUnits,
	    forceUnits,
	    bcInputFile,
	    &errorcode,
	    errorstring,
	    strlen(displacementUnits),
	    strlen(velocityUnits),
	    strlen(forceUnits),
	    strlen(bcInputFile),
	    sizeof(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberBcEntries:" << numberBcEntries
    << journal::endl;

  // return
  return Py_BuildValue((char *) "i", numberBcEntries);
}


// Scan connectivities

char pypylith3d_scan_connect__doc__[] = "";
char pypylith3d_scan_connect__name__[] = "scan_connect";

PyObject * pypylith3d_scan_connect(PyObject *, PyObject *args)
{
  PyObject* pyPointerToListArrayNumberElementNodesBase;
  PyObject* pyPointerToMaterialModelInfo;
  PyObject* pyPointerToListArrayMaterialModel;
  PyObject* pyPointerToVolumeElementFamilyList;
  int maxNumberVolumeElementFamilies;
  int numberMaterials;
  int f77FileInput;
  char* connectivityInputFile;

  int ok = PyArg_ParseTuple(args, (char *) "OOOOiiis:scan_connect",
			    &pyPointerToListArrayNumberElementNodesBase,
			    &pyPointerToMaterialModelInfo,
                            &pyPointerToListArrayMaterialModel,
			    &pyPointerToVolumeElementFamilyList,
			    &maxNumberVolumeElementFamilies,
			    &numberMaterials,
			    &f77FileInput,
			    &connectivityInputFile);

  if (!ok) {
    return 0;
  }

  int* pointerToListArrayNumberElementNodesBase = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayNumberElementNodesBase);
  int* pointerToMaterialModelInfo = (int*) PyCObject_AsVoidPtr(pyPointerToMaterialModelInfo);
  int* pointerToListArrayMaterialModel = (int*) PyCObject_AsVoidPtr(pyPointerToListArrayMaterialModel);
  int* pointerToVolumeElementFamilyList = (int*) PyCObject_AsVoidPtr(pyPointerToVolumeElementFamilyList);
  int numberVolumeElements = 0;
  int numberVolumeElementFamilies = 0;
  int volumeElementType = 0;
  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];

  scan_connect_f(pointerToListArrayNumberElementNodesBase,
		 pointerToMaterialModelInfo,
		 pointerToListArrayMaterialModel,
		 pointerToVolumeElementFamilyList,
		 &maxNumberVolumeElementFamilies,
		 &numberMaterials,
		 &numberVolumeElements,
		 &numberVolumeElementFamilies,
		 &volumeElementType,
		 &f77FileInput,
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
    << "numberVolumeElementFamilies:" << numberVolumeElementFamilies
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue((char *) "iii", numberVolumeElements,
		       numberVolumeElementFamilies,
		       volumeElementType);
}


// Scan coordinates

char pypylith3d_scan_coords__doc__[] = "";
char pypylith3d_scan_coords__name__[] = "scan_coords";

PyObject * pypylith3d_scan_coords(PyObject *, PyObject *args)
{
  int f77FileInput;
  char *coordinateUnits;
  char *coordinateInputFile;

  int ok = PyArg_ParseTuple(args, (char *) "iss:scan_coords",
			    &f77FileInput,
			    &coordinateUnits,
			    &coordinateInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int numberNodes = 0;

  scan_coords_f(&numberNodes,
		&f77FileInput,
		coordinateUnits,
		coordinateInputFile,
		&errorcode,
		errorstring,
		strlen(coordinateUnits),
		strlen(coordinateInputFile),
		sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberNodes:" << numberNodes
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue((char *) "i", numberNodes);
}


// Scan differential forces

char pypylith3d_scan_diff__doc__[] = "";
char pypylith3d_scan_diff__name__[] = "scan_diff";

PyObject * pypylith3d_scan_diff(PyObject *, PyObject *args)
{
  int numberSlipperyNodeEntries;
  int f77FileInput;
  char* differentialForceInputFile;

  int ok = PyArg_ParseTuple(args, (char *) "iis:scan_diff",
			    &numberSlipperyNodeEntries,
			    &f77FileInput,
			    &differentialForceInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int numberDifferentialForceEntries = 0;

  scan_diff_f(&numberSlipperyNodeEntries,
	      &numberDifferentialForceEntries,
	      &f77FileInput,
	      differentialForceInputFile,
	      &errorcode,
	      errorstring,
	      strlen(differentialForceInputFile),
	      sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberDifferentialForceEntries:" << numberDifferentialForceEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue((char *) "i", numberDifferentialForceEntries);
}


// Scan time steps at which full output is desired

char pypylith3d_scan_fuldat__doc__[] = "";
char pypylith3d_scan_fuldat__name__[] = "scan_fuldat";

PyObject * pypylith3d_scan_fuldat(PyObject *, PyObject *args)
{
  int analysisTypeInt;
  int totalNumberTimeSteps;
  int f77FileInput;
  char* fullOutputInputFile;

  int ok = PyArg_ParseTuple(args, (char *) "iiis:scan_fuldat",
			    &analysisTypeInt,
			    &totalNumberTimeSteps,
			    &f77FileInput,
			    &fullOutputInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int numberFullOutputs = 0;

  scan_fuldat_f(&analysisTypeInt,
		&totalNumberTimeSteps,
		&numberFullOutputs,
		&f77FileInput,
		fullOutputInputFile,
		&errorcode,
		errorstring,
		strlen(fullOutputInputFile),
		sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberFullOutputs:" << numberFullOutputs
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue((char *) "i", numberFullOutputs);
}


// Scan load histories

char pypylith3d_scan_hist__doc__[] = "";
char pypylith3d_scan_hist__name__[] = "scan_hist";

PyObject * pypylith3d_scan_hist(PyObject *, PyObject *args)
{
  int f77FileInput;
  char* loadHistoryInputFile;

  int ok = PyArg_ParseTuple(args, (char *) "is:scan_hist",
			    &f77FileInput,
			    &loadHistoryInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int numberLoadHistories = 0;

  scan_hist_f(&numberLoadHistories,
	      &f77FileInput,
	      loadHistoryInputFile,
	      &errorcode,
	      errorstring,
	      strlen(loadHistoryInputFile),
	      sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberLoadHistories:" << numberLoadHistories
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue((char *) "i", numberLoadHistories);
}


// Scan element prestresses

  // char pypylith3d_scan_prestr__doc__[] = "";
  // char pypylith3d_scan_prestr__name__[] = "scan_prestr";

  // PyObject * pypylith3d_scan_prestr(PyObject *, PyObject *args)
  // {
  //   int numberStressComponents;
  //   int numberPrestressGaussPoints;
  //   int numberElements;
  //   int prestressAutoComputeInt;
  //   int f77FileInput;
  //   char* prestressInputFile;

  //   int ok = PyArg_ParseTuple(args, "iiiiis:scan_prestr",
  // 			    &numberStressComponents,
  // 			    &numberPrestressGaussPoints,
  // 			    &numberElements,
  // 			    &prestressAutoComputeInt,
  // 			    &f77FileInput,
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
  // 		&prestressAutoComputeInt,
  // 		&f77FileInput,
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
  int f77FileInput;
  char* rotationUnits;
  char* rotationInputFile;

  int ok = PyArg_ParseTuple(args, (char *) "iss:scan_skew",
			    &f77FileInput,
			    &rotationUnits,
			    &rotationInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int numberRotationEntries = 0;

  scan_skew_f(&numberRotationEntries,
	      &f77FileInput,
	      rotationUnits,
	      rotationInputFile,
	      &errorcode,
	      errorstring,
	      strlen(rotationUnits),
	      strlen(rotationInputFile),
	      sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberRotationEntries:" << numberRotationEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue((char *) "i", numberRotationEntries);
}


// Scan slippery node entries

char pypylith3d_scan_slip__doc__[] = "";
char pypylith3d_scan_slip__name__[] = "scan_slip";

PyObject * pypylith3d_scan_slip(PyObject *, PyObject *args)
{
  int f77FileInput;
  char* slipperyNodeInputFile;

  int ok = PyArg_ParseTuple(args, (char *) "is:scan_slip",
			    &f77FileInput,
			    &slipperyNodeInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int numberSlipperyNodeEntries = 0;

  scan_slip_f(&numberSlipperyNodeEntries,
	      &f77FileInput,
	      slipperyNodeInputFile,
	      &errorcode,
	      errorstring,
	      strlen(slipperyNodeInputFile),
	      sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberSlipperyNodeEntries:" << numberSlipperyNodeEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue((char *) "i", numberSlipperyNodeEntries);
}


// Scan split node entries

char pypylith3d_scan_split__doc__[] = "";
char pypylith3d_scan_split__name__[] = "scan_split";

PyObject * pypylith3d_scan_split(PyObject *, PyObject *args)
{
  int f77FileInput;
  char* splitNodeInputFile;

  int ok = PyArg_ParseTuple(args, (char *) "is:scan_split",
			    &f77FileInput,
			    &splitNodeInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int numberSplitNodeEntries = 0;

  scan_split_f(&numberSplitNodeEntries,
	       &f77FileInput,
	       splitNodeInputFile,
	       &errorcode,
	       errorstring,
	       strlen(splitNodeInputFile),
	       sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberSplitNodeEntries:" << numberSplitNodeEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue((char *) "i", numberSplitNodeEntries);
}


// Scan time step data

char pypylith3d_scan_timdat__doc__[] = "";
char pypylith3d_scan_timdat__name__[] = "scan_timdat";

PyObject * pypylith3d_scan_timdat(PyObject *, PyObject *args)
{
  int f77FileInput;
  char* timeUnits;
  char* timeStepInputFile;

  int ok = PyArg_ParseTuple(args, (char *) "iss:scan_timdat",
			    &f77FileInput,
			    &timeUnits,
			    &timeStepInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int numberTimeStepGroups = 0;
  int totalNumberTimeSteps = 0;

  scan_timdat_f(&totalNumberTimeSteps,
		&numberTimeStepGroups,
		&f77FileInput,
		timeUnits,
		timeStepInputFile,
		&errorcode,
		errorstring,
		strlen(timeUnits),
		strlen(timeStepInputFile),
		sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "totalNumberTimeSteps:" << totalNumberTimeSteps
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue((char *) "ii", numberTimeStepGroups,
		       totalNumberTimeSteps);
}



// Scan traction BC

char pypylith3d_scan_tractions__doc__[] = "";
char pypylith3d_scan_tractions__name__[] = "scan_tractions";

PyObject * pypylith3d_scan_tractions(PyObject *, PyObject *args)
{
  int maxElementNodes2d;
  int f77FileInput;
  char *tractionUnits;
  char *tractionInputFile;

  int ok = PyArg_ParseTuple(args, (char *) "iiss:scan_tractions",
                            &maxElementNodes2d,
			    &f77FileInput,
			    &tractionUnits,
			    &tractionInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int numberTractionBc = 0;

  scan_tractions_f(&numberTractionBc,
		   &maxElementNodes2d,
		   &f77FileInput,
		   tractionUnits,
		   tractionInputFile,
		   &errorcode,
		   errorstring,
		   strlen(tractionUnits),
		   strlen(tractionInputFile),
		   sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberTractionBc:" << numberTractionBc
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue((char *) "i", numberTractionBc);
}

// Read winkler BC

char pypylith3d_scan_wink__doc__[] = "";
char pypylith3d_scan_wink__name__[] = "scan_wink";

PyObject * pypylith3d_scan_wink(PyObject *, PyObject *args)
{
  int f77FileInput;
  char* winklerInputFile;

  int ok = PyArg_ParseTuple(args, (char *) "is:scan_wink",
			    &f77FileInput,
			    &winklerInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int numberWinklerEntries = 0;
  int numberWinklerForces = 0;

  scan_wink_f(&numberWinklerEntries,
	      &numberWinklerForces,
	      &f77FileInput,
	      winklerInputFile,
	      &errorcode,
	      errorstring,
	      strlen(winklerInputFile),
	      sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberWinklerForces:" << numberWinklerForces
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue((char *) "ii", numberWinklerEntries,
	               numberWinklerForces);
}


// Read winkler BC for slippery nodes

char pypylith3d_scan_winkx__doc__[] = "";
char pypylith3d_scan_winkx__name__[] = "scan_winkx";

PyObject * pypylith3d_scan_winkx(PyObject *, PyObject *args)
{
  int numberSlipperyNodeEntries;
  int f77FileInput;
  char* slipperyWinklerInputFile;

  int ok = PyArg_ParseTuple(args, (char *) "iis:scan_winkx",
			    &numberSlipperyNodeEntries,
			    &f77FileInput,
			    &slipperyWinklerInputFile);

  if (!ok) {
    return 0;
  }

  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];
  int numberSlipperyWinklerEntries = 0;
  int numberSlipperyWinklerForces = 0;

  scan_winkx_f(&numberSlipperyNodeEntries,
	       &numberSlipperyWinklerEntries,
	       &numberSlipperyWinklerForces,
	       &f77FileInput,
	       slipperyWinklerInputFile,
	       &errorcode,
	       errorstring,
	       strlen(slipperyWinklerInputFile),
	       sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberSlipperyWinklerForces:" << numberSlipperyWinklerForces
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue((char *) "ii", numberSlipperyWinklerEntries,
		       numberSlipperyWinklerForces);
}
    

// Initialize material model info

char pypylith3d_matmod_def__doc__[] = "";
char pypylith3d_matmod_def__name__[] = "matmod_def";

PyObject * pypylith3d_matmod_def(PyObject *, PyObject *args)
{
  PyObject* pyPointerToMaterialModelInfo;

  int ok = PyArg_ParseTuple(args, "O:matmod_def",
			    &pyPointerToMaterialModelInfo);

  if (!ok) {
    return 0;
  }

  int* pointerToMaterialModelInfo = (int*) PyCObject_AsVoidPtr(pyPointerToMaterialModelInfo);

  matmod_def_f(pointerToMaterialModelInfo);
    
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "pointerToMaterialModelInfo:" << pointerToMaterialModelInfo
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
  PyObject* pyPointerToSh;
  PyObject* pyPointerToShj;
  PyObject* pyPointerToGauss;
  int quadratureOrderInt;
  int elementType;
  int numberVolumeElementNodes;
  int numberVolumeElementGaussPoints;

  int ok = PyArg_ParseTuple(args, "OOOiiii:preshape",
			    &pyPointerToSh,
			    &pyPointerToShj,
			    &pyPointerToGauss,
			    &quadratureOrderInt,
			    &elementType,
			    &numberVolumeElementNodes,
			    &numberVolumeElementGaussPoints);

  if (!ok) {
    return 0;
  }

  double* pointerToSh = (double*) PyCObject_AsVoidPtr(pyPointerToSh);
  double* pointerToShj = (double*) PyCObject_AsVoidPtr(pyPointerToShj);
  double* pointerToGauss = (double*) PyCObject_AsVoidPtr(pyPointerToGauss);
  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];

  preshape_f(pointerToSh,
	     pointerToShj,
	     pointerToGauss,
	     &quadratureOrderInt,
	     &elementType,
	     &numberVolumeElementNodes,
	     &numberVolumeElementGaussPoints,
	     &errorcode,
	     errorstring,
	     sizeof(errorstring));

  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }
    
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "quadratureOrderInt:" << quadratureOrderInt
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
  PyObject* pyPointerToSh2d;
  PyObject* pyPointerToGauss2d;
  int quadratureOrderInt;
  int elementType;
  int numberSurfaceElementNodes;
  int numberSurfaceElementGaussPoints;

  int ok = PyArg_ParseTuple(args, "OOiiii:preshape2d",
			    &pyPointerToSh2d,
			    &pyPointerToGauss2d,
			    &quadratureOrderInt,
			    &elementType,
			    &numberSurfaceElementNodes,
			    &numberSurfaceElementGaussPoints);

  if (!ok) {
    return 0;
  }

  double* pointerToSh2d = (double*) PyCObject_AsVoidPtr(pyPointerToSh2d);
  double* pointerToGauss2d = (double*) PyCObject_AsVoidPtr(pyPointerToGauss2d);
  int errorcode = 0;
  const int maxsize = 4096;
  char errorstring[maxsize];

  preshape2d_f(pointerToSh2d,
	       pointerToGauss2d,
	       &quadratureOrderInt,
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
    << "quadratureOrderInt:" << quadratureOrderInt
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
