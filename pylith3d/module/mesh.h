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

#if !defined(pypylith3d_mesh_h)
#define pypylith3d_mesh_h

// process mesh
extern char pypylith3d_processMesh__name__[];
extern char pypylith3d_processMesh__doc__[];
extern "C"
PyObject * pypylith3d_processMesh(PyObject *, PyObject *);

// create a PETSc Mat
extern char pypylith3d_createPETScMat__name__[];
extern char pypylith3d_createPETScMat__doc__[];
extern "C"
PyObject * pypylith3d_createPETScMat(PyObject *, PyObject *);

// destroy a PETSc Mat
extern char pypylith3d_destroyPETScMat__name__[];
extern char pypylith3d_destroyPETScMat__doc__[];
extern "C"
PyObject * pypylith3d_destroyPETScMat(PyObject *, PyObject *);

// output a PETSc Mesh and Fields
extern char pypylith3d_outputMesh__name__[];
extern char pypylith3d_outputMesh__doc__[];
extern "C"
PyObject * pypylith3d_outputMesh(PyObject *, PyObject *);

#endif

// End of file
