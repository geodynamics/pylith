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

#include <petscmesh.h>
#include <petscmat.h>
#include <portinfo>
#include "journal/debug.h"

#include <Python.h>

#include "scanner.h"
#include "exceptionhandler.h"
#include "pylith3d_externs.h"
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

#undef __FUNCT__
#define __FUNCT__ "WriteBoundary_PyLith"
PetscErrorCode WriteBoundary_PyLith(const char *baseFilename, ALE::Obj<ALE::Mesh> mesh)
{
  FILE                       *f;
  char                        bcFilename[2048];
  typedef std::pair<ALE::Mesh::field_type::patch_type,int> patch_type;
  ALE::Obj<ALE::Mesh::foliation_type> boundaries = mesh->getBoundaries();
  ALE::Obj<ALE::Mesh::bundle_type>    vertexBundle = mesh->getBundle(0);
  ALE::Mesh::field_type::patch_type patch;
  PetscErrorCode              ierr;

  PetscFunctionBegin;
  if (mesh->debug) {
    boundaries->view("PyLith boundaries");
  }
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
  ALE::Obj<ALE::Mesh::sieve_type::traits::depthSequence> vertices = boundaries->getTopology()->depthStratum(0);

  for(ALE::Mesh::sieve_type::traits::depthSequence::iterator v_itor = vertices->begin(); v_itor != vertices->end(); v_itor++) {
    int    constraints[3];
    double values[3] = {0.0, 0.0, 0.0};

    for(int c = 0; c < 3; c++) {
      ALE::Mesh::foliation_type::patch_type p(patch, c+1);

      constraints[c] = boundaries->getFiberDimension(p, *v_itor);
      if (constraints[c]) {
        values[c] = boundaries->restrict(p, *v_itor)[0];
      }
    }

    if (constraints[0] || constraints[1] || constraints[2]) {
      fprintf(f, "%7d %4d %4d %4d % 16.8E % 16.8E % 16.8E\n", vertexBundle->getIndex(patch, *v_itor).prefix+1,
              constraints[0], constraints[1], constraints[2], values[0], values[1], values[2]);
    }
  }
  fclose(f);
  PetscFunctionReturn(0);
}

// Process mesh

PetscErrorCode MeshView_Sieve_Newer(ALE::Obj<ALE::Mesh> mesh, PetscViewer viewer);

char pypylith3d_processMesh__doc__[] = "";
char pypylith3d_processMesh__name__[] = "processMesh";

PyObject * pypylith3d_processMesh(PyObject *, PyObject *args)
{
  char *meshInputFile;
  char  meshOutputFile[2048];
  int   interpolateMesh;

  int ok = PyArg_ParseTuple(args, "si:processMesh", &meshInputFile, &interpolateMesh);

  if (!ok) {
    return 0;
  }

  journal::debug_t  debug("pylith3d");
  MPI_Comm          comm = PETSC_COMM_WORLD;
  PetscMPIInt       rank;
  ALE::Obj<ALE::Mesh> mesh;
  PetscViewer       viewer;
  PetscInt         *boundaryVertices;
  PetscScalar      *boundaryValues;
  PetscInt          numBoundaryVertices, numBoundaryComponents;
  PetscErrorCode    ierr;

  ierr = MPI_Comm_rank(comm, &rank);
  sprintf(meshOutputFile, "%s.%d", meshInputFile, rank);
  mesh = ALE::PyLithBuilder::createNew(comm, meshInputFile, (bool) interpolateMesh);
  debug << journal::at(__HERE__) << "[" << rank << "]Created new PETSc Mesh for " << meshInputFile << journal::endl;
  mesh = mesh->distribute();
  debug << journal::at(__HERE__) << "[" << rank << "]Distributed PETSc Mesh"  << journal::endl;
  ierr = ReadBoundary_PyLith(meshInputFile, PETSC_FALSE, &numBoundaryVertices, &numBoundaryComponents, &boundaryVertices, &boundaryValues);

  typedef std::pair<ALE::Mesh::field_type::patch_type,int> patch_type;
  ALE::Obj<ALE::Mesh::foliation_type> boundaries = mesh->getBoundaries();
  ALE::Mesh::field_type::patch_type patch;
  std::set<int> seen;
  int numElements = mesh->getBundle(mesh->getTopology()->depth())->getGlobalOffsets()[mesh->commSize()];

  boundaries->setTopology(mesh->getTopology());
  for(int c = 0; c < numBoundaryComponents; c++) {
    boundaries->setPatch(mesh->getTopology()->leaves(), ALE::Mesh::foliation_type::patch_type(patch, c+1));
  }
  // Reverse order allows newer conditions to override older, as required by PyLith
  for(int v = numBoundaryVertices-1; v >= 0; v--) {
    ALE::Mesh::point_type vertex(0, boundaryVertices[v*(numBoundaryComponents+1)] + numElements);

    if (seen.find(vertex.index) == seen.end()) {
      for(int c = 0; c < numBoundaryComponents; c++) {
        if (boundaryVertices[v*(numBoundaryComponents+1)+c+1]) {
          boundaries->setFiberDimension(ALE::Mesh::foliation_type::patch_type(patch, c+1), vertex, 1);
        }
      }
      seen.insert(vertex.index);
    }
  }
  boundaries->orderPatches();
  for(int v = 0; v < numBoundaryVertices; v++) {
    ALE::Mesh::point_type vertex(0, boundaryVertices[v*(numBoundaryComponents+1)] + numElements);

    for(int c = 0; c < numBoundaryComponents; c++) {
      boundaries->update(ALE::Mesh::foliation_type::patch_type(patch, c+1), vertex, &boundaryValues[v*numBoundaryComponents+c]);
    }
  }
  debug << journal::at(__HERE__) << "[" << rank << "]Created boundary conditions"  << journal::endl;

  bool refineMesh = false;
  if (refineMesh) {
    double refinementLimit = 2.4e4*2.4e4*2.4e4*0.01*0.5;
    bool   interpolate     = true;
    mesh = ALE::Generator::refine(mesh, refinementLimit, interpolate);
  }

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
  ierr = MeshView_Sieve_Newer(mesh, viewer);
  ierr = PetscViewerDestroy(viewer);
  debug << journal::at(__HERE__) << "[" << rank << "]Output new PyLith mesh into: " << meshOutputFile << journal::endl;

  ierr = WriteBoundary_PyLith(meshOutputFile, mesh);
  debug << journal::at(__HERE__) << "[" << rank << "]Wrote PyLith boundary conditions"  << journal::endl;

  ALE::Obj<ALE::Mesh::field_type> field = mesh->getField("displacement");
  ALE::Obj<ALE::Mesh::sieve_type::traits::depthSequence> vertices = mesh->getTopology()->depthStratum(0);

  field->setPatch(mesh->getTopology()->leaves(), patch);
  field->setFiberDimensionByDepth(patch, 0, 3);
  for(ALE::Mesh::sieve_type::traits::depthSequence::iterator v_itor = vertices->begin(); v_itor != vertices->end(); v_itor++) {
    int numConstraints = 0;

    for(int c = 0; c < numBoundaryComponents; c++) {
      numConstraints += boundaries->getFiberDimension(ALE::Mesh::foliation_type::patch_type(patch, c+1), *v_itor);
    }

    if (numConstraints > 0) {
      if (mesh->debug) {
        std::cout << "[" << rank << "]Setting dimension of " << *v_itor << " to " << 3 - numConstraints << std::endl;
      }
      field->setFiberDimension(patch, *v_itor, 3 - numConstraints);
    }
  }
  field->orderPatches();
  field->createGlobalOrder();
  if (mesh->debug) {
    field->view("Displacement field");
  }
  ALE::Obj<ALE::Mesh::sieve_type::traits::heightSequence> elements = mesh->getTopology()->heightStratum(0);
  ALE::Obj<ALE::Mesh::bundle_type> vertexBundle = mesh->getBundle(0);
  std::string orderName("element");

  for(ALE::Mesh::sieve_type::traits::heightSequence::iterator e_iter = elements->begin(); e_iter != elements->end(); e_iter++) {
    // setFiberDimensionByDepth() does not work here since we only want it to apply to the patch cone
    //   What we really need is the depthStratum relative to the patch
    ALE::Obj<ALE::Mesh::bundle_type::order_type::coneSequence> cone = vertexBundle->getPatch(orderName, *e_iter);

    field->setPatch(orderName, cone, *e_iter);
    for(ALE::Mesh::bundle_type::order_type::coneSequence::iterator c_iter = cone->begin(); c_iter != cone->end(); ++c_iter) {
      field->setFiberDimension(orderName, *e_iter, *c_iter, field->getFiberDimension(patch, *c_iter));
    }
  }
  field->orderPatches(orderName);
  debug << journal::at(__HERE__) << "[" << rank << "]Created displacement Field"  << journal::endl;

  // return
  PyObject *pyMesh = PyCObject_FromVoidPtr(mesh.ptr(), NULL);
  mesh.int_allocator->del(mesh.refCnt);
  mesh.refCnt = NULL;
  return Py_BuildValue("sN", meshOutputFile, pyMesh);
}

// Create a PETSc Mat
char pypylith3d_createPETScMat__doc__[] = "";
char pypylith3d_createPETScMat__name__[] = "createPETScMat";

PyObject * pypylith3d_createPETScMat(PyObject *, PyObject *args)
{
  PyObject *pyMesh, *pyA, *pyRhs, *pySol;
  MPI_Comm comm = PETSC_COMM_WORLD;
  Mat      A;
  Vec      rhs, sol;

  int ok = PyArg_ParseTuple(args, "O:createPETScMat", &pyMesh);
  if (!ok) {
    return 0;
  }

  ALE::Mesh *mesh = (ALE::Mesh *) PyCObject_AsVoidPtr(pyMesh);
  const int *offsets = mesh->getField("displacement")->getGlobalOffsets();
  int size = offsets[mesh->commRank()+1] - offsets[mesh->commRank()];

  if (MatCreate(comm, &A)) {
    PyErr_SetString(PyExc_RuntimeError, "Could not create PETSc Mat");
    return 0;
  }
  if (MatSetSizes(A, size, size, PETSC_DETERMINE, PETSC_DETERMINE)) {
    PyErr_SetString(PyExc_RuntimeError, "Could not set sizes for PETSc Mat");
    return 0;
  }
  if (VecCreate(comm, &rhs)) {
    PyErr_SetString(PyExc_RuntimeError, "Could not create PETSc Rhs");
    return 0;
  }
  if (VecSetSizes(rhs, size, PETSC_DETERMINE)) {
    PyErr_SetString(PyExc_RuntimeError, "Could not set sizes for PETSc Rhs");
    return 0;
  }
  if (VecSetFromOptions(rhs)) {
    PyErr_SetString(PyExc_RuntimeError, "Could not set options for PETSc Rhs");
    return 0;
  }
  if (VecDuplicate(rhs, &sol)) {
    PyErr_SetString(PyExc_RuntimeError, "Could not create PETSc Sol");
    return 0;
  }

  PetscObjectContainer c;
  PetscErrorCode       ierr;

  ierr = PetscObjectContainerCreate(comm, &c);
  ierr = PetscObjectContainerSetPointer(c, mesh);
  ierr = PetscObjectCompose((PetscObject) A, "mesh", (PetscObject) c);
  ierr = PetscObjectContainerDestroy(c);

  VecScatter injection;
  ierr = MeshGetGlobalScatter(mesh, "displacement", rhs, &injection);
  ierr = PetscObjectContainerCreate(comm, &c);
  ierr = PetscObjectContainerSetPointer(c, mesh);
  ierr = PetscObjectCompose((PetscObject) rhs, "mesh", (PetscObject) c);
  ierr = PetscObjectContainerDestroy(c);
  ierr = PetscObjectCompose((PetscObject) rhs, "injection", (PetscObject) injection);
  ierr = PetscObjectContainerCreate(comm, &c);
  ierr = PetscObjectContainerSetPointer(c, mesh);
  ierr = PetscObjectCompose((PetscObject) sol, "mesh", (PetscObject) c);
  ierr = PetscObjectContainerDestroy(c);
  ierr = PetscObjectCompose((PetscObject) sol, "injection", (PetscObject) injection);

  ierr = MatSetFromOptions(A);
  ierr = preallocateMatrix(A, mesh, mesh->getField("displacement"));

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "Created PETSc Mat:" << size
    << journal::endl;

  // return Py_None;
  pyA = PyCObject_FromVoidPtr(A, NULL);
  pyRhs = PyCObject_FromVoidPtr(rhs, NULL);
  pySol = PyCObject_FromVoidPtr(sol, NULL);
  return Py_BuildValue("NNN", pyA, pyRhs, pySol);
}

// Destroy a PETSc Mat

char pypylith3d_destroyPETScMat__doc__[] = "";
char pypylith3d_destroyPETScMat__name__[] = "destroyPETScMat";

PyObject * pypylith3d_destroyPETScMat(PyObject *, PyObject *args)
{
  PyObject *pyA,*pyRhs, *pySol;
  Mat A;
  Vec rhs, sol;

  int ok = PyArg_ParseTuple(args, "OOO:destroyPETScMat", &pyA, &pyRhs, &pySol);
  if (!ok) {
    return 0;
  }

  A = (Mat) PyCObject_AsVoidPtr(pyA);
  if (MatDestroy(A)) {
    PyErr_SetString(PyExc_RuntimeError, "Could not destroy PETSc Mat");
    return 0;
  }
  rhs = (Vec) PyCObject_AsVoidPtr(pyRhs);
  if (VecDestroy(rhs)) {
    PyErr_SetString(PyExc_RuntimeError, "Could not destroy PETSc Rhs");
    return 0;
  }
  sol = (Vec) PyCObject_AsVoidPtr(pySol);
  if (VecDestroy(sol)) {
    PyErr_SetString(PyExc_RuntimeError, "Could not destroy PETSc Sol");
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "Destroyed PETSc Mat"
    << journal::endl;

  // return Py_None;
  Py_INCREF(Py_None);
  return Py_None;
}

PetscErrorCode FieldView_Sieve(ALE::Obj<ALE::Mesh>, const std::string&, PetscViewer);

char pypylith3d_outputMesh__doc__[] = "";
char pypylith3d_outputMesh__name__[] = "outputMesh";

PyObject * pypylith3d_outputMesh(PyObject *, PyObject *args)
{
  PyObject *pyMesh, *pySol;
  char     *meshBaseFile;

  int ok = PyArg_ParseTuple(args, "sOO:outputMesh", &meshBaseFile, &pyMesh, &pySol);
  if (!ok) {
    return 0;
  }

  ALE::Mesh  *mesh = (ALE::Mesh *) PyCObject_AsVoidPtr(pyMesh);
  Vec         sol  = (Vec)         PyCObject_AsVoidPtr(pySol);
  std::string filename(meshBaseFile);
  PetscViewer viewer;
  Vec         partition;

  ALE::Obj<ALE::Mesh> m(mesh);

  // Injection Vec in to Field
  ALE::Obj<ALE::Mesh::field_type> displacement = m->getField("displacement");
  ALE::Mesh::field_type::patch_type patch;
  Vec        l;
  VecScatter injection;

  VecCreateSeqWithArray(PETSC_COMM_SELF, displacement->getSize(patch), displacement->restrict(patch), &l);
  PetscObjectQuery((PetscObject) sol, "injection", (PetscObject *) &injection);
  VecScatterBegin(sol, l, INSERT_VALUES, SCATTER_REVERSE, injection);
  VecScatterEnd(sol, l, INSERT_VALUES, SCATTER_REVERSE, injection);
  VecDestroy(l);

  // Create complete field by adding BC
  ALE::Obj<ALE::Mesh::field_type> full_displacement = m->getField("full_displacement");
  ALE::Obj<ALE::Mesh::field_type::order_type::baseSequence> patches = displacement->getPatches();
  ALE::Obj<ALE::Mesh::foliation_type> boundaries = m->getBoundaries();

  // This is wrong if the domain changes
  if (!full_displacement->getGlobalOrder()) {
    for(ALE::Mesh::field_type::order_type::baseSequence::iterator p_iter = patches->begin(); p_iter != patches->end(); ++p_iter) {
      full_displacement->setPatch(displacement->getPatch(*p_iter), *p_iter);
      full_displacement->setFiberDimensionByDepth(*p_iter, 0, 3);
    }
    full_displacement->orderPatches();
    full_displacement->createGlobalOrder();
  }
  for(ALE::Mesh::field_type::order_type::baseSequence::iterator p_iter = patches->begin(); p_iter != patches->end(); ++p_iter) {
    ALE::Obj<ALE::Mesh::field_type::order_type::coneSequence> elements = full_displacement->getPatch(*p_iter);

    for(ALE::Mesh::field_type::order_type::coneSequence::iterator e_iter = elements->begin(); e_iter != elements->end(); ++e_iter) {
      const int     dim   = displacement->getIndex(*p_iter, *e_iter).index;
      const double *array = displacement->restrict(*p_iter, *e_iter);
      int           v     = 0;
      double        values[3];

      for(int c = 0; c < 3; c++) {
        ALE::Mesh::foliation_type::patch_type        bPatch(*p_iter, c+1);
        const ALE::Mesh::foliation_type::index_type& idx = boundaries->getIndex(bPatch, *e_iter);

        if (idx.index > 0) {
          values[c] = 0.0;
        } else if (dim > 0) {
          values[c] = array[v++];
        }
      }
      if (v != dim) {
        std::cout << "ERROR: Invalid size " << v << " used for " << *e_iter << " with index " << displacement->getIndex(*p_iter, *e_iter) << std::endl;
      }
      full_displacement->updateAdd(*p_iter, *e_iter, values);
    }
  }
  ALE::Obj<ALE::Mesh::foliation_type::order_type::baseSequence> bdPatches = boundaries->getPatches();

  for(ALE::Mesh::foliation_type::order_type::baseSequence::iterator p_iter = bdPatches->begin(); p_iter != bdPatches->end(); ++p_iter) {
    ALE::Obj<ALE::Mesh::foliation_type::order_type::coneSequence> elements = boundaries->getPatch(*p_iter);

    for(ALE::Mesh::foliation_type::order_type::coneSequence::iterator e_iter = elements->begin(); e_iter != elements->end(); ++e_iter) {
      const double *array = full_displacement->restrict((*p_iter).first, *e_iter);
      double        values[3];

      if (boundaries->getIndex(*p_iter, *e_iter).index > 0) {
        for(int c = 0; c < 3; c++) {
          values[c] = array[c];
        }
        values[(*p_iter).second-1] = boundaries->restrict(*p_iter, *e_iter)[0];
        full_displacement->update((*p_iter).first, *e_iter, values);
      }
    }
  }

  PetscViewerCreate(mesh->comm(), &viewer);
  PetscViewerSetType(viewer, PETSC_VIEWER_ASCII);
  PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_VTK);
  filename += ".vtk";
  PetscViewerFileSetName(viewer, filename.c_str());
  MeshView_Sieve_Newer(m, viewer);
  FieldView_Sieve(m, "full_displacement", viewer);
  PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_VTK_CELL);
  FieldView_Sieve(m, "material", viewer);
  //FieldView_Sieve(partition, viewer);
  PetscViewerPopFormat(viewer);
  PetscViewerDestroy(viewer);

  m.int_allocator->del(m.refCnt);
  m.refCnt = NULL;

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "Output PETSc Mesh and Solution"
    << journal::endl;

  // return Py_None;
  Py_INCREF(Py_None);
  return Py_None;
}

char pypylith3d_interpolatePoints__doc__[] = "";
char pypylith3d_interpolatePoints__name__[] = "interpolatePoints";

PyObject * pypylith3d_interpolatePoints(PyObject *, PyObject *args)
{
  PyObject *pyMesh, *pySol, *pyPoints,*pyValues;

  int ok = PyArg_ParseTuple(args, "OOO:outputMesh", &pyMesh, &pySol, &pyPoints);
  if (!ok) {
    return 0;
  }

  ALE::Mesh  *mesh = (ALE::Mesh *) PyCObject_AsVoidPtr(pyMesh);
  Vec         sol  = (Vec)         PyCObject_AsVoidPtr(pySol);
  double     *points;

  // Convert Numeric matrix to C array

  ALE::Obj<ALE::Mesh> m(mesh);

  // Injection Vec in to Field
  ALE::Obj<ALE::Mesh::field_type> displacement = m->getField("displacement");
  ALE::Mesh::field_type::patch_type patch;
  Vec        l;
  VecScatter injection;

  VecCreateSeqWithArray(PETSC_COMM_SELF, displacement->getSize(patch), displacement->restrict(patch), &l);
  PetscObjectQuery((PetscObject) sol, "injection", (PetscObject *) &injection);
  VecScatterBegin(sol, l, INSERT_VALUES, SCATTER_REVERSE, injection);
  VecScatterEnd(sol, l, INSERT_VALUES, SCATTER_REVERSE, injection);
  VecDestroy(l);

  // Create complete field by adding BC
  ALE::Obj<ALE::Mesh::field_type> full_displacement = m->getField("full_displacement");
  ALE::Obj<ALE::Mesh::field_type::order_type::baseSequence> patches = displacement->getPatches();
  ALE::Obj<ALE::Mesh::foliation_type> boundaries = m->getBoundaries();

  // This is wrong if the domain changes
  if (!full_displacement->getGlobalOrder()) {
    for(ALE::Mesh::field_type::order_type::baseSequence::iterator p_iter = patches->begin(); p_iter != patches->end(); ++p_iter) {
      full_displacement->setPatch(displacement->getPatch(*p_iter), *p_iter);
      full_displacement->setFiberDimensionByDepth(*p_iter, 0, 3);
    }
    full_displacement->orderPatches();
    full_displacement->createGlobalOrder();
  }
  for(ALE::Mesh::field_type::order_type::baseSequence::iterator p_iter = patches->begin(); p_iter != patches->end(); ++p_iter) {
    ALE::Obj<ALE::Mesh::field_type::order_type::coneSequence> elements = full_displacement->getPatch(*p_iter);

    for(ALE::Mesh::field_type::order_type::coneSequence::iterator e_iter = elements->begin(); e_iter != elements->end(); ++e_iter) {
      const int     dim   = displacement->getIndex(*p_iter, *e_iter).index;
      const double *array = displacement->restrict(*p_iter, *e_iter);
      int           v     = 0;
      double        values[3];

      for(int c = 0; c < 3; c++) {
        ALE::Mesh::foliation_type::patch_type        bPatch(*p_iter, c+1);
        const ALE::Mesh::foliation_type::index_type& idx = boundaries->getIndex(bPatch, *e_iter);

        if (idx.index > 0) {
          values[c] = 0.0;
        } else if (dim > 0) {
          values[c] = array[v++];
        }
      }
      if (v != dim) {
        std::cout << "ERROR: Invalid size " << v << " used for " << *e_iter << " with index " << displacement->getIndex(*p_iter, *e_iter) << std::endl;
      }
      full_displacement->updateAdd(*p_iter, *e_iter, values);
    }
  }
  ALE::Obj<ALE::Mesh::foliation_type::order_type::baseSequence> bdPatches = boundaries->getPatches();

  for(ALE::Mesh::foliation_type::order_type::baseSequence::iterator p_iter = bdPatches->begin(); p_iter != bdPatches->end(); ++p_iter) {
    ALE::Obj<ALE::Mesh::foliation_type::order_type::coneSequence> elements = boundaries->getPatch(*p_iter);

    for(ALE::Mesh::foliation_type::order_type::coneSequence::iterator e_iter = elements->begin(); e_iter != elements->end(); ++e_iter) {
      const double *array = full_displacement->restrict((*p_iter).first, *e_iter);
      double        values[3];

      if (boundaries->getIndex(*p_iter, *e_iter).index > 0) {
        for(int c = 0; c < 3; c++) {
          values[c] = array[c];
        }
        values[(*p_iter).second-1] = boundaries->restrict(*p_iter, *e_iter)[0];
        full_displacement->update((*p_iter).first, *e_iter, values);
      }
    }
  }

  //interpolatePoints(m, full_displacement, points, values);

  // Convert C array to Numeric matrix

  m.int_allocator->del(m.refCnt);
  m.refCnt = NULL;

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "Interpolated points"
    << journal::endl;

  return Py_BuildValue("N", pyValues);
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
  return Py_BuildValue("i", numberBcEntries);
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
  return Py_BuildValue("iii", numberVolumeElements,
		       numberVolumeElementFamilies,
		       volumeElementType);
}


// Read coordinates

char pypylith3d_scan_coords__doc__[] = "";
char pypylith3d_scan_coords__name__[] = "scan_coords";

PyObject * pypylith3d_scan_coords(PyObject *, PyObject *args)
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
		sizeof(errorstring));
    
  if(0 != exceptionhandler(errorcode, errorstring)) {
    return 0;
  }

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "numberSpaceDimensions:" << numberNodes
    << journal::endl;

  // return
  Py_INCREF(Py_None);
  return Py_BuildValue("i", numberNodes);
}


// Read differential forces

char pypylith3d_scan_diff__doc__[] = "";
char pypylith3d_scan_diff__name__[] = "scan_diff";

PyObject * pypylith3d_scan_diff(PyObject *, PyObject *args)
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
  return Py_BuildValue("i",numberDifferentialForceEntries);
}


// Read time steps at which full output is desired

char pypylith3d_scan_fuldat__doc__[] = "";
char pypylith3d_scan_fuldat__name__[] = "scan_fuldat";

PyObject * pypylith3d_scan_fuldat(PyObject *, PyObject *args)
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
  return Py_BuildValue("i",numberFullOutputs);
}


// Read load histories

char pypylith3d_scan_hist__doc__[] = "";
char pypylith3d_scan_hist__name__[] = "scan_hist";

PyObject * pypylith3d_scan_hist(PyObject *, PyObject *args)
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
  return Py_BuildValue("i",numberLoadHistories);
}


// Read element prestresses

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


// Read local coordinate rotations

char pypylith3d_scan_skew__doc__[] = "";
char pypylith3d_scan_skew__name__[] = "scan_skew";

PyObject * pypylith3d_scan_skew(PyObject *, PyObject *args)
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
  return Py_BuildValue("i",numberRotationEntries);
}


// Read slippery node entries

char pypylith3d_scan_slip__doc__[] = "";
char pypylith3d_scan_slip__name__[] = "scan_slip";

PyObject * pypylith3d_scan_slip(PyObject *, PyObject *args)
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
  return Py_BuildValue("i",numberSlipperyNodeEntries);
}


// Read split node entries

char pypylith3d_scan_split__doc__[] = "";
char pypylith3d_scan_split__name__[] = "scan_split";

PyObject * pypylith3d_scan_split(PyObject *, PyObject *args)
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
  return Py_BuildValue("i",numberSplitNodeEntries);
}


// Read time step data

char pypylith3d_scan_timdat__doc__[] = "";
char pypylith3d_scan_timdat__name__[] = "scan_timdat";

PyObject * pypylith3d_scan_timdat(PyObject *, PyObject *args)
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
  return Py_BuildValue("ii",numberTimeStepGroups,
		       totalNumberTimeSteps);
}


// Read traction BC

// char pypylith3d_scan_traction__doc__[] = "";
// char pypylith3d_scan_traction__name__[] = "scan_traction";

// PyObject * pypylith3d_scan_traction(PyObject *, PyObject *args)
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
		  // sizeof(errorstring));
    
  // if(0 != exceptionhandler(errorcode, errorstring)) {
    // return 0;
  // }

  // journal::debug_t debug("pylith3d");
  // debug
    // << journal::at(__HERE__)
    // << "numberTractionBc:" << numberTractionBc
    // << journal::endl;

  // return
  // Py_INCREF(Py_None);
  // return Py_BuildValue("i",numberTractionBc);
// }


// Read winkler BC

char pypylith3d_scan_wink__doc__[] = "";
char pypylith3d_scan_wink__name__[] = "scan_wink";

PyObject * pypylith3d_scan_wink(PyObject *, PyObject *args)
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
  return Py_BuildValue("ii",numberWinklerEntries,
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
  return Py_BuildValue("ii",numberSlipperyWinklerEntries,
		       numberSlipperyWinklerForces);
}
    
// version
// $Id: scanner.cc,v 1.7 2005/06/07 19:39:11 willic3 Exp $

// End of file
