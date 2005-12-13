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

#include <petscmat.h>
#include <petscmesh.h>
#include <portinfo>
#include "journal/debug.h"

#include <Python.h>

#include "scanner.h"
#include "exceptionhandler.h"
#include "lithomop3d_externs.h"
#include <stdio.h>
#include <string.h>

#undef __FUNCT__
#define __FUNCT__ "IgnoreComments_PyLith"
PetscErrorCode IgnoreComments_PyLith(char *buf, PetscInt bufSize, FILE *f)
{
  PetscFunctionBegin;
  while((fgets(buf, bufSize, f) != NULL) && (buf[0] == '#')) {}
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ReadBoundary_PyLith"
PetscErrorCode ReadBoundary_PyLith(const char *baseFilename, PetscTruth useZeroBase, PetscInt *numBoundaryVertices, PetscInt *numBoundaryComponents, PetscInt **boundaryVertices, PetscScalar **boundaryValues)
{
  FILE          *f;
  PetscInt       maxVerts= 1024, vertexCount = 0;
  PetscInt       numComp = 3;
  PetscInt      *verts;
  PetscScalar   *values;
  char           bcFilename[2048];
  char           buf[2048];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscStrcpy(bcFilename, baseFilename);
  ierr = PetscStrcat(bcFilename, ".bc");
  f = fopen(bcFilename, "r");CHKERRQ(ierr);
  IgnoreComments_PyLith(buf, 2048, f);
  /* Ignore displacement units */
  fgets(buf, 2048, f);
  /* Ignore velocity units */
  fgets(buf, 2048, f);
  /* Ignore force units */
  fgets(buf, 2048, f);
  IgnoreComments_PyLith(buf, 2048, f);
  ierr = PetscMalloc(maxVerts*(numComp+1) * sizeof(PetscInt),    &verts);CHKERRQ(ierr);
  ierr = PetscMalloc(maxVerts*numComp * sizeof(PetscScalar), &values);CHKERRQ(ierr);
  do {
    const char *v = strtok(buf, " ");
    int vertex = atoi(v);
        
    if (!useZeroBase) vertex -= 1;
    if (vertexCount == maxVerts) {
      PetscInt *vtmp;
      PetscScalar *ctmp;

      vtmp = verts;
      ierr = PetscMalloc(maxVerts*2*(numComp+1) * sizeof(PetscInt), &verts);CHKERRQ(ierr);
      ierr = PetscMemcpy(verts, vtmp, maxVerts*(numComp+1) * sizeof(PetscInt));CHKERRQ(ierr);
      ierr = PetscFree(vtmp);CHKERRQ(ierr);
      ctmp = values;
      ierr = PetscMalloc(maxVerts*2*numComp * sizeof(PetscScalar), &values);CHKERRQ(ierr);
      ierr = PetscMemcpy(values, ctmp, maxVerts*numComp * sizeof(PetscScalar));CHKERRQ(ierr);
      ierr = PetscFree(ctmp);CHKERRQ(ierr);
      maxVerts *= 2;
    }
    verts[vertexCount*(numComp+1)+0] = vertex;
    /* X boundary condition*/
    v = strtok(NULL, " ");
    verts[vertexCount*(numComp+1)+1] = atoi(v);
    /* Y boundary condition*/
    v = strtok(NULL, " ");
    verts[vertexCount*(numComp+1)+2] = atoi(v);
    /* Z boundary condition*/
    v = strtok(NULL, " ");
    verts[vertexCount*(numComp+1)+3] = atoi(v);
    /* X boundary value */
    v = strtok(NULL, " ");
    values[vertexCount*numComp+0] = atof(v);
    /* Y boundary value */
    v = strtok(NULL, " ");
    values[vertexCount*numComp+1] = atof(v);
    /* Z boundary value */
    v = strtok(NULL, " ");
    values[vertexCount*numComp+2] = atof(v);
    vertexCount++;
  } while(fgets(buf, 2048, f) != NULL);
  fclose(f);
  *numBoundaryVertices = vertexCount;
  *numBoundaryComponents = numComp;
  *boundaryVertices = verts;
  *boundaryValues   = values;
  PetscFunctionReturn(0);
}

#include <src/dm/mesh/sieve/ALE_exception.hh>
#include <src/dm/mesh/sieve/ALE_mem.hh>
#include <src/dm/mesh/sieve/ALE_containers.hh>
#include <src/dm/mesh/sieve/ALE.hh>
#include <src/dm/mesh/sieve/PreSieve.hh>
#include <src/dm/mesh/sieve/Sieve.hh>
#include <src/dm/mesh/sieve/Stack.hh>
#include <src/dm/mesh/sieve/IndexBundle.hh>

#undef __FUNCT__
#define __FUNCT__ "WriteBoundary_PyLith"
PetscErrorCode WriteBoundary_PyLith(const char *baseFilename, ALE::Sieve *topology, ALE::Sieve *boundary, ALE::IndexBundle *boundaryBundle, Vec boundaryVec)
{
  FILE          *f;
  char           bcFilename[2048];
  PetscScalar   *boundaryValues;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ALE::IndexBundle vertexBundle(topology);
  vertexBundle.setFiberDimensionByDepth(0, 1);
  vertexBundle.computeOverlapIndices();
  vertexBundle.computeGlobalIndices();

  ierr = PetscStrcpy(bcFilename, baseFilename);
  ierr = PetscStrcat(bcFilename, ".bc");
  f = fopen(bcFilename, "w");CHKERRQ(ierr);
  fprintf(f, "displacement_units = m\n");
  fprintf(f, "velocity_units = m/s\n");
  fprintf(f, "force_units = newton\n");
  fprintf(f, "#\n");
  fprintf(f, "# The last row for each node applies\n");
  fprintf(f, "#\n");
  fprintf(f, "#  Node X BC Y BC Z BC   X Value          Y Value          Z Value\n");
  fprintf(f, "#\n");
  ALE::Obj<ALE::Point_set> vertices = boundary->cap();

  ierr = VecGetArray(boundaryVec, &boundaryValues);
  for(ALE::Point_set::iterator v_itor = vertices->begin(); v_itor != vertices->end(); v_itor++) {
    ALE::Point                 vertex = *v_itor;
    ALE::Obj<ALE::Point_array> intervals = vertexBundle.getLocalOrderedClosureIndices(ALE::Point_set(vertex));
    ALE::Obj<ALE::Point_set>   support = boundary->support(vertex);
    ALE::Point                 interval = boundaryBundle->getFiberInterval(vertex);
    int                        constraints[3] = {0, 0, 0};

    for(ALE::Point_set::iterator s_itor = support->begin(); s_itor != support->end(); s_itor++) {
      ALE::Point boundaryPoint = *s_itor;

      constraints[boundaryPoint.prefix] = boundaryPoint.index;
    }
    fprintf(f, "%7d %4d %4d %4d % 16.8E % 16.8E % 16.8E\n", intervals->begin()->prefix+1,
            constraints[0], constraints[1], constraints[2],
            boundaryValues[interval.prefix+0], boundaryValues[interval.prefix+1], boundaryValues[interval.prefix+2]);
  }
  ierr = VecRestoreArray(boundaryVec, &boundaryValues);
  fclose(f);
  PetscFunctionReturn(0);
}

// Process mesh

char pylithomop3d_processMesh__doc__[] = "";
char pylithomop3d_processMesh__name__[] = "processMesh";

PyObject * pylithomop3d_processMesh(PyObject *, PyObject *args)
{
  char *meshInputFile;
  char  meshOutputFile[2048];

  int ok = PyArg_ParseTuple(args, "s:processMesh", &meshInputFile);

  if (!ok) {
    return 0;
  }

  MPI_Comm          comm = PETSC_COMM_WORLD;
  PetscMPIInt       rank;
  Mesh              mesh;
  ALE::Sieve       *topology;
  ALE::PreSieve    *orientation;
  ALE::Sieve       *boundary;
  ALE::IndexBundle *boundaryBundle;
  Vec               boundaryVec;
  PetscViewer       viewer;
  PetscInt         *boundaryVertices;
  PetscScalar      *boundaryValues;
  PetscInt          numBoundaryVertices, numBoundaryComponents;
  PetscErrorCode    ierr;

  ierr = MeshCreatePyLith(comm, meshInputFile, &mesh);
  ierr = MeshDistribute(mesh);
  ierr = ReadBoundary_PyLith(meshInputFile, PETSC_FALSE, &numBoundaryVertices, &numBoundaryComponents, &boundaryVertices, &boundaryValues);
  ierr = MeshCreateBoundary(mesh, numBoundaryVertices, numBoundaryComponents, boundaryVertices, boundaryValues, (void **) &boundaryBundle, &boundaryVec);
  ierr = PetscViewerCreate(comm, &viewer);
  ierr = PetscViewerSetType(viewer, PETSC_VIEWER_ASCII);
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_PYLITH_LOCAL);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ);
  ierr = PetscExceptionTry1(PetscViewerFileSetName(viewer, meshInputFile), PETSC_ERR_FILE_OPEN);
  if (PetscExceptionValue(ierr)) {
    /* this means that a caller above me has also tryed this exception so I don't handle it here, pass it up */
  } else if (PetscExceptionCaught(ierr, PETSC_ERR_FILE_OPEN)) {
    ierr = 0;
  } 
  ierr = MeshView(mesh, viewer);
  ierr = PetscViewerDestroy(viewer);

  ierr = MPI_Comm_rank(comm, &rank);
  sprintf(meshOutputFile, "%s.%d", meshInputFile, rank);
  ierr = MeshGetTopology(mesh, (void **) &topology);
  ierr = MeshGetOrientation(mesh, (void **) &orientation);
  ierr = MeshGetBoundary(mesh, (void **) &boundary);
  ierr = WriteBoundary_PyLith(meshOutputFile, topology, boundary, boundaryBundle, boundaryVec);

  ALE::IndexBundle        *fieldBundle = new ALE::IndexBundle(topology);
  ALE::Obj<ALE::Point_set> cap = boundary->cap();

  fieldBundle->setFiberDimensionByDepth(0, 3);
  for(ALE::Point_set::iterator b_itor = cap->begin(); b_itor != cap->end(); b_itor++) {
    ALE::Point vertex = *b_itor;
#if 0
    int numConstraints = boundaryBundle->getFiberInterval(vertex).index;
#else
    ALE::Obj<ALE::Point_set> support = boundary->support(vertex);
    int numConstraints = 0;

    for(ALE::Point_set::iterator s_itor = support->begin(); s_itor != support->end(); s_itor++) {
      numConstraints += (*s_itor).index;
    }
#endif
    printf("Setting dimension of (%d, %d) to %d\n", vertex.prefix, vertex.index, 3 - numConstraints);
    fieldBundle->setFiberDimension(vertex, 3 - numConstraints);
  }
  fieldBundle->computeOverlapIndices();
  fieldBundle->computeGlobalIndices();
  ierr = MeshSetBundle(mesh, (void *) fieldBundle);

  ierr = VecDestroy(boundaryVec);
  delete boundaryBundle;

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "Output new mesh into: " << meshOutputFile
    << journal::endl;

  // return
  PyObject *pyMesh = PyCObject_FromVoidPtr(mesh, NULL);
  return Py_BuildValue("sN", meshOutputFile, pyMesh);
}

// Scan boundary conditions

char pylithomop3d_scan_bc__doc__[] = "";
char pylithomop3d_scan_bc__name__[] = "scan_bc";

PyObject * pylithomop3d_scan_bc(PyObject *, PyObject *args)
{
  int f77FileInput;
  char* displacementUnits;
  char* velocityUnits;
  char* forceUnits;
  char* bcInputFile;

  int ok = PyArg_ParseTuple(args, "issss:scan_bc",
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
	    strlen(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberBcEntries:" << numberBcEntries
    << journal::endl;

  // return
  return Py_BuildValue("i", numberBcEntries);
}


// Scan connectivities

char pylithomop3d_scan_connect__doc__[] = "";
char pylithomop3d_scan_connect__name__[] = "scan_connect";

PyObject * pylithomop3d_scan_connect(PyObject *, PyObject *args)
{
  PyObject* pyPointerToListArrayNumberElementNodesBase;
  PyObject* pyPointerToMaterialModelInfo;
  PyObject* pyPointerToListArrayMaterialModel;
  PyObject* pyPointerToVolumeElementFamilyList;
  int maxNumberVolumeElementFamilies;
  int numberMaterials;
  int f77FileInput;
  char* connectivityInputFile;

  int ok = PyArg_ParseTuple(args, "OOOOiiis:scan_connect",
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
		 strlen(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberVolumeElementFamilies:" << numberVolumeElementFamilies
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("iii", numberVolumeElements,
		       numberVolumeElementFamilies,
		       volumeElementType);
}


// Read coordinates

char pylithomop3d_scan_coords__doc__[] = "";
char pylithomop3d_scan_coords__name__[] = "scan_coords";

PyObject * pylithomop3d_scan_coords(PyObject *, PyObject *args)
{
  int f77FileInput;
  char *coordinateUnits;
  char *coordinateInputFile;

  int ok = PyArg_ParseTuple(args, "iss:scan_coords",
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
		strlen(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberSpaceDimensions:" << numberNodes
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i", numberNodes);
}


// Read differential forces

char pylithomop3d_scan_diff__doc__[] = "";
char pylithomop3d_scan_diff__name__[] = "scan_diff";

PyObject * pylithomop3d_scan_diff(PyObject *, PyObject *args)
{
  int numberSlipperyNodeEntries;
  int f77FileInput;
  char* differentialForceInputFile;

  int ok = PyArg_ParseTuple(args, "iis:scan_diff",
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
	      strlen(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberDifferentialForceEntries:" << numberDifferentialForceEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i",numberDifferentialForceEntries);
}


// Read time steps at which full output is desired

char pylithomop3d_scan_fuldat__doc__[] = "";
char pylithomop3d_scan_fuldat__name__[] = "scan_fuldat";

PyObject * pylithomop3d_scan_fuldat(PyObject *, PyObject *args)
{
  int analysisTypeInt;
  int totalNumberTimeSteps;
  int f77FileInput;
  char* fullOutputInputFile;

  int ok = PyArg_ParseTuple(args, "iiis:scan_fuldat",
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
		strlen(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberFullOutputs:" << numberFullOutputs
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i",numberFullOutputs);
}


// Read load histories

char pylithomop3d_scan_hist__doc__[] = "";
char pylithomop3d_scan_hist__name__[] = "scan_hist";

PyObject * pylithomop3d_scan_hist(PyObject *, PyObject *args)
{
  int f77FileInput;
  char* loadHistoryInputFile;

  int ok = PyArg_ParseTuple(args, "is:scan_hist",
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
	      strlen(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberLoadHistories:" << numberLoadHistories
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i",numberLoadHistories);
}


// Read element prestresses

  // char pylithomop3d_scan_prestr__doc__[] = "";
  // char pylithomop3d_scan_prestr__name__[] = "scan_prestr";

  // PyObject * pylithomop3d_scan_prestr(PyObject *, PyObject *args)
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

  //   journal::debug_t debug("lithomop3d");
  //   debug
  //     << journal::at(__HERE__)
  //     << "numberPrestressEntries:" << numberPrestressEntries
  //     << journal::endl;

  // return
  //   Py_INCREF(Py_None);
  //   return Py_BuildValue("i",numberPrestressEntries);
  // }


// Read local coordinate rotations

char pylithomop3d_scan_skew__doc__[] = "";
char pylithomop3d_scan_skew__name__[] = "scan_skew";

PyObject * pylithomop3d_scan_skew(PyObject *, PyObject *args)
{
  int f77FileInput;
  char* rotationUnits;
  char* rotationInputFile;

  int ok = PyArg_ParseTuple(args, "iss:scan_skew",
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
	      strlen(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberRotationEntries:" << numberRotationEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i",numberRotationEntries);
}


// Read slippery node entries

char pylithomop3d_scan_slip__doc__[] = "";
char pylithomop3d_scan_slip__name__[] = "scan_slip";

PyObject * pylithomop3d_scan_slip(PyObject *, PyObject *args)
{
  int f77FileInput;
  char* slipperyNodeInputFile;

  int ok = PyArg_ParseTuple(args, "is:scan_slip",
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
	      strlen(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberSlipperyNodeEntries:" << numberSlipperyNodeEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i",numberSlipperyNodeEntries);
}


// Read split node entries

char pylithomop3d_scan_split__doc__[] = "";
char pylithomop3d_scan_split__name__[] = "scan_split";

PyObject * pylithomop3d_scan_split(PyObject *, PyObject *args)
{
  int f77FileInput;
  char* splitNodeInputFile;

  int ok = PyArg_ParseTuple(args, "is:scan_split",
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
	       strlen(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberSplitNodeEntries:" << numberSplitNodeEntries
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i",numberSplitNodeEntries);
}


// Read time step data

char pylithomop3d_scan_timdat__doc__[] = "";
char pylithomop3d_scan_timdat__name__[] = "scan_timdat";

PyObject * pylithomop3d_scan_timdat(PyObject *, PyObject *args)
{
  int f77FileInput;
  char* timeUnits;
  char* timeStepInputFile;

  int ok = PyArg_ParseTuple(args, "iss:scan_timdat",
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
		strlen(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "totalNumberTimeSteps:" << totalNumberTimeSteps
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("ii",numberTimeStepGroups,
		       totalNumberTimeSteps);
}


// Read traction BC

// char pylithomop3d_scan_traction__doc__[] = "";
// char pylithomop3d_scan_traction__name__[] = "scan_traction";

// PyObject * pylithomop3d_scan_traction(PyObject *, PyObject *args)
// {
  // char* tractionBcUnits;
  // int f77FileInput;
  // char* tractionInputFile;

  // int ok = PyArg_ParseTuple(args, "sis:scan_traction",
			    // &tractionBcUnits,
			    // &f77FileInput,
			    // &tractionInputFile);

  // if (!ok) {
    // return 0;
  // }

  // int errorcode = 0;
  // const int maxsize = 4096;
  // char errorstring[maxsize];
  // int numberTractionBc = 0;

  // scan_traction_f(&numberTractionBc,
		  // &f77FileInput,
		  // tractionBcUnits,
		  // tractionInputFile,
		  // &errorcode,
		  // errorstring,
		  // strlen(tractionBcUnits),
		  // strlen(tractionInputFile),
		  // strlen(errorstring));
    
  // if(0 != exceptionhandler(errorcode, errorstring)) {
    // return 0;
  // }

  // journal::debug_t debug("lithomop3d");
  // debug
    // << journal::at(__HERE__)
    // << "numberTractionBc:" << numberTractionBc
    // << journal::endl;

  // return
  // Py_INCREF(Py_None);
  // return Py_BuildValue("i",numberTractionBc);
// }


// Read winkler BC

char pylithomop3d_scan_wink__doc__[] = "";
char pylithomop3d_scan_wink__name__[] = "scan_wink";

PyObject * pylithomop3d_scan_wink(PyObject *, PyObject *args)
{
  int f77FileInput;
  char* winklerInputFile;

  int ok = PyArg_ParseTuple(args, "is:scan_wink",
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
	      strlen(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberWinklerForces:" << numberWinklerForces
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("ii",numberWinklerEntries,
	               numberWinklerForces);
}


// Read winkler BC for slippery nodes

char pylithomop3d_scan_winkx__doc__[] = "";
char pylithomop3d_scan_winkx__name__[] = "scan_winkx";

PyObject * pylithomop3d_scan_winkx(PyObject *, PyObject *args)
{
  int numberSlipperyNodeEntries;
  int f77FileInput;
  char* slipperyWinklerInputFile;

  int ok = PyArg_ParseTuple(args, "iis:scan_winkx",
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
	       strlen(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("lithomop3d");
  debug
    << journal::at(__HERE__)
    << "numberSlipperyWinklerForces:" << numberSlipperyWinklerForces
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("ii",numberSlipperyWinklerEntries,
		       numberSlipperyWinklerForces);
}
    
// version
// $Id: scanner.cc,v 1.7 2005/06/07 19:39:11 willic3 Exp $

// End of file
