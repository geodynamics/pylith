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

#include <Distribution.hh>
#include <petscmesh.h>
#include <src/dm/mesh/meshpylith.h>
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
  while((fgets(buf, bufSize, f) != NULL) && ((buf[0] == '#') || (buf[0] == '\0'))) {}
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
  /* Don't need this one since IgnoreComments reads until "failure".
   * fgets(buf, 2048, f);
   */
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
PetscErrorCode WriteBoundary_PyLith(const char *baseFilename, const ALE::Obj<ALE::Mesh>& mesh)
{
  ALE::Mesh::foliated_section_type::patch_type      patch      = 0;
  const ALE::Obj<ALE::Mesh::numbering_type>&        vNumbering = mesh->getFactory()->getLocalNumbering(mesh->getTopology(), patch, 0);
  const ALE::Obj<ALE::Mesh::foliated_section_type>& boundaries = mesh->getBoundariesNew();
  FILE          *f;
  char           bcFilename[2048];
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (mesh->debug()) {
    boundaries->view("PyLith boundaries");
  }
  ierr = PetscStrcpy(bcFilename, baseFilename);
  ierr = PetscStrcat(bcFilename, ".bc");

  // Determine if we have bc stuff
  const ALE::Obj<ALE::Mesh::topology_type::label_sequence>& vertices = boundaries->getTopology()->depthStratum(patch, 0);
  bool haveBC = false;
  for(ALE::Mesh::topology_type::label_sequence::iterator v_iter = vertices->begin(); v_iter != vertices->end(); ++v_iter)
    if (boundaries->getFiberDimension(patch, *v_iter) > 0) {
      haveBC = true;
      break;
    } // if

  f = fopen(bcFilename, "w");CHKERRQ(ierr);
  if (haveBC) {
    // Only write header if bc file contains information
    // If header is written and rest of file is empty,
    // then we have a problem and reading will fail with
    // error message.
    fprintf(f, "displacement_units = m\n");
    fprintf(f, "velocity_units = m/s\n");
    fprintf(f, "force_units = newton\n");
    fprintf(f, "#\n");
    fprintf(f, "# The last row for each node applies\n");
    fprintf(f, "#\n");
    fprintf(f, "#  Node X BC Y BC Z BC   X Value          Y Value          Z Value\n");
    fprintf(f, "#\n");
  } // if

  for(ALE::Mesh::topology_type::label_sequence::iterator v_iter = vertices->begin(); v_iter != vertices->end(); ++v_iter) {
    int    constraints[3] = {0, 0, 0};
    double values[3] = {0.0, 0.0, 0.0};
    int    size = boundaries->getFiberDimension(patch, *v_iter);
    const ALE::Mesh::foliated_section_type::value_type *array = boundaries->restrict(patch, *v_iter);

    for(int c = 0; c < size; c++) {
      constraints[array[c].first] = 1;
      values[array[c].first]      = array[c].second;
    }

    if (constraints[0] || constraints[1] || constraints[2]) {
      fprintf(f, "%7d %4d %4d %4d % 16.8E % 16.8E % 16.8E\n", vNumbering->getIndex(*v_iter)+1,
              constraints[0], constraints[1], constraints[2], values[0], values[1], values[2]);
    }
  }
  fclose(f);
  PetscFunctionReturn(0);
}

// Process mesh

PetscErrorCode MeshView_Sieve(const ALE::Obj<ALE::Mesh>& mesh, PetscViewer viewer);

char pypylith3d_processMesh__doc__[] = "";
char pypylith3d_processMesh__name__[] = "processMesh";

PyObject * pypylith3d_processMesh(PyObject *, PyObject *args)
{
  char *meshInputFile;
  char  *meshOutputFile;
  int   interpolateMesh;
  char* partitioner;

  int ok = PyArg_ParseTuple(args, (char *) "ssis:processMesh", &meshOutputFile, &meshInputFile, &interpolateMesh, &partitioner);

  if (!ok) {
    return 0;
  }

  using ALE::Obj;
  journal::debug_t  debug("pylith3d");
  MPI_Comm          comm = PETSC_COMM_WORLD;
  PetscMPIInt       rank;
  Mesh              mesh;
  Obj<ALE::Mesh>    m;
  PetscViewer       viewer;
  PetscInt         *boundaryVertices;
  PetscScalar      *boundaryValues;
  PetscInt          numBoundaryVertices, numBoundaryComponents;
  int               debugFlag = 0;
  PetscErrorCode    ierr;

  ierr = MPI_Comm_rank(comm, &rank);
  ierr = MeshCreatePyLith(comm, 3, meshInputFile, PETSC_FALSE, (PetscTruth) interpolateMesh, &mesh);
  ierr = MeshGetMesh(mesh, m);
  m->setDebug(debugFlag);
  int numElements = m->getTopology()->heightStratum(0, 0)->size();
  debug << journal::at(__HERE__) << "[" << rank << "]Created new PETSc Mesh for " << meshInputFile << journal::endl;
  m = ALE::New::Distribution<ALE::Mesh::topology_type>::distributeMesh(m, partitioner);
  ierr = MeshSetMesh(mesh, m);
  debug << journal::at(__HERE__) << "[" << rank << "]Distributed PETSc Mesh"  << journal::endl;
  ierr = ReadBoundary_PyLith(meshInputFile, PETSC_FALSE, &numBoundaryVertices, &numBoundaryComponents, &boundaryVertices, &boundaryValues);

  const Obj<ALE::Mesh::foliated_section_type>& boundaries = m->getBoundariesNew();
  ALE::Mesh::foliated_section_type::patch_type patch      = 0;
  std::set<int> seen;

  ierr = MPI_Bcast(&numElements, 1, MPI_INT, 0, comm);
  boundaries->setTopology(m->getTopology());
  // Reverse order allows newer conditions to override older, as required by PyLith
  for(int v = numBoundaryVertices-1; v >= 0; v--) {
    ALE::Mesh::point_type vertex(boundaryVertices[v*(numBoundaryComponents+1)] + numElements);
    int size = 0;

    if (seen.find(vertex) == seen.end()) {
      for(int c = 0; c < numBoundaryComponents; c++) {
        size += boundaryVertices[v*(numBoundaryComponents+1)+c+1];
      }
      boundaries->setFiberDimension(patch, vertex, size);
      seen.insert(vertex);
    }
  }
  boundaries->allocate();
  for(int v = 0; v < numBoundaryVertices; v++) {
    ALE::Mesh::point_type vertex(boundaryVertices[v*(numBoundaryComponents+1)] + numElements);
    ALE::Mesh::foliated_section_type::value_type values[3];

    for(int c = 0, i = 0; c < numBoundaryComponents; c++) {
      if (boundaryVertices[v*(numBoundaryComponents+1)+c+1]) {
        values[i].first  = c;
        values[i].second = boundaryValues[v*numBoundaryComponents+c];
        i++;
      }
    }
    boundaries->update(patch, vertex, values);
  }
  debug << journal::at(__HERE__) << "[" << rank << "]Created boundary conditions"  << journal::endl;

#if 0
  bool refineMesh = false;
  if (refineMesh) {
    double refinementLimit = 2.4e4*2.4e4*2.4e4*0.01*0.5;
    bool   interpolate     = true;
    mesh = ALE::Generator::refine(mesh, refinementLimit, interpolate);
  }
#endif

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
  debug << journal::at(__HERE__) << "[" << rank << "]Output new PyLith mesh into: " << meshOutputFile << journal::endl;

  ierr = WriteBoundary_PyLith(meshOutputFile, m);
  debug << journal::at(__HERE__) << "[" << rank << "]Wrote PyLith boundary conditions"  << journal::endl;

  const Obj<ALE::Mesh::topology_type::label_sequence>& vertices = m->getTopology()->depthStratum(0, 0);
  SectionReal                       section;
  Obj<ALE::Mesh::real_section_type> s;

  ierr = MeshGetSectionReal(mesh, "default", &section);
  ierr = PetscObjectSetName((PetscObject) section, "displacement");
  ierr = MeshSetSectionReal(mesh, section);
  ierr = SectionRealGetSection(section, s);
  s->setFiberDimensionByDepth(0, 0, 3);
  for(ALE::Mesh::topology_type::label_sequence::iterator v_iter = vertices->begin(); v_iter != vertices->end(); ++v_iter) {
    int numConstraints = boundaries->getFiberDimension(patch, *v_iter);

    if (numConstraints > 0) {
      if (m->debug()) {
        std::cout << "[" << rank << "]Setting dimension of " << *v_iter << " to " << 3 - numConstraints << std::endl;
      }
      s->setFiberDimension(0, *v_iter, 3 - numConstraints);
    }
  }
  s->allocate();
  if (m->debug()) {
    s->view("Displacement field");
  }
  ierr = SectionRealDestroy(section);
  debug << journal::at(__HERE__) << "[" << rank << "]Created displacement Field"  << journal::endl;

  m->getFactory()->constructInverseOrder(m->getFactory()->getLocalNumbering(m->getTopology(), 0, m->getTopology()->depth()));

  // return
  PyObject *pyMesh = PyCObject_FromVoidPtr(mesh, NULL);
  //PyObject *pyMesh = PyCObject_FromVoidPtr(mesh.ptr(), NULL);
  //mesh.int_allocator->del(mesh.refCnt);
  //mesh.refCnt = NULL;
  return pyMesh;
}

// Create a PETSc Mat
char pypylith3d_createPETScMat__doc__[] = "";
char pypylith3d_createPETScMat__name__[] = "createPETScMat";

PyObject * pypylith3d_createPETScMat(PyObject *, PyObject *args)
{
  using ALE::Obj;
  PyObject *pyMesh, *pyA, *pyRhs, *pySol;
  MPI_Comm comm = PETSC_COMM_WORLD;
  Mat      A;
  Vec      rhs, sol;
  PetscErrorCode ierr;

  int ok = PyArg_ParseTuple(args, (char *) "O:createPETScMat", &pyMesh);
  if (!ok) {
    return 0;
  }

  Mesh mesh = (Mesh) PyCObject_AsVoidPtr(pyMesh);
  Obj<ALE::Mesh> m;

  ierr = MeshGetMesh(mesh, m);
  const ALE::Obj<ALE::Mesh::order_type>& offsets = m->getFactory()->getGlobalOrder(m->getTopology(), 0, "displacement", m->getRealSection("displacement")->getAtlas());
  int localSize = offsets->getLocalSize();
  int globalSize = offsets->getGlobalSize();

  if (MatCreate(comm, &A)) {
    PyErr_SetString(PyExc_RuntimeError, "Could not create PETSc Mat");
    return 0;
  }
  if (MatSetSizes(A, localSize, localSize, globalSize, globalSize)) {
    PyErr_SetString(PyExc_RuntimeError, "Could not set sizes for PETSc Mat");
    return 0;
  }
  if (VecCreate(comm, &rhs)) {
    PyErr_SetString(PyExc_RuntimeError, "Could not create PETSc Rhs");
    return 0;
  }
  if (VecSetSizes(rhs, localSize, globalSize)) {
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

  ierr = PetscObjectCompose((PetscObject) A, "mesh", (PetscObject) mesh);

  VecScatter injection = NULL;
  ierr = MeshGetGlobalScatter(mesh, &injection);
  ierr = PetscObjectCompose((PetscObject) rhs, "mesh",      (PetscObject) mesh);
  ierr = PetscObjectCompose((PetscObject) rhs, "injection", (PetscObject) injection);
  ierr = PetscObjectCompose((PetscObject) sol, "mesh",      (PetscObject) mesh);
  ierr = PetscObjectCompose((PetscObject) sol, "injection", (PetscObject) injection);

  ierr = MatSetFromOptions(A);
  ierr = preallocateMatrix(m->getTopology(), m->getRealSection("displacement")->getAtlas(), m->getFactory()->getGlobalOrder(m->getTopology(), 0, "displacement", m->getRealSection("displacement")->getAtlas()), A);

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "Created PETSc Mat: " << localSize << " " << globalSize
    << journal::endl;

  // return Py_None;
  pyA = PyCObject_FromVoidPtr(A, NULL);
  pyRhs = PyCObject_FromVoidPtr(rhs, NULL);
  pySol = PyCObject_FromVoidPtr(sol, NULL);
  return Py_BuildValue((char *) "NNN", pyA, pyRhs, pySol);
}

// Destroy a PETSc Mat

char pypylith3d_destroyPETScMat__doc__[] = "";
char pypylith3d_destroyPETScMat__name__[] = "destroyPETScMat";

PyObject * pypylith3d_destroyPETScMat(PyObject *, PyObject *args)
{
  PyObject *pyA,*pyRhs, *pySol;
  Mat A;
  Vec rhs, sol;

  int ok = PyArg_ParseTuple(args, (char *) "OOO:destroyPETScMat", &pyA, &pyRhs, &pySol);
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

PetscErrorCode updateDisplacement(SectionReal displacement, Vec sol) {
  VecScatter     injection;
  Vec            lv;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = SectionRealGetLocalVector(displacement, &lv);CHKERRQ(ierr);
  ierr = PetscObjectQuery((PetscObject) sol, "injection", (PetscObject *) &injection);CHKERRQ(ierr);
  ierr = VecScatterBegin(sol, lv, INSERT_VALUES, SCATTER_REVERSE, injection);CHKERRQ(ierr);
  ierr = VecScatterEnd(sol, lv, INSERT_VALUES, SCATTER_REVERSE, injection);CHKERRQ(ierr);
  ierr = VecDestroy(lv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// Create complete displacement field by adding BC
PetscErrorCode createFullDisplacement(Mesh mesh, SectionReal *fullDisplacement) {
  using ALE::Obj;
  Obj<ALE::Mesh> m;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MeshGetMesh(mesh, m);CHKERRQ(ierr);
  if (!m->hasRealSection("full_displacement")) {
    const Obj<ALE::Mesh::topology_type::sheaf_type>& patches    = m->getTopology()->getPatches();
    const Obj<ALE::Mesh::real_section_type>&         full       = m->getRealSection("full_displacement");

    for(ALE::Mesh::topology_type::sheaf_type::iterator p_iter = patches->begin(); p_iter != patches->end(); ++p_iter) {
      full->setFiberDimensionByDepth(p_iter->first, 0, 3);
    }
    full->allocate();
  }
  const Obj<ALE::Mesh::foliated_section_type>&     boundaries = m->getBoundariesNew();
  const Obj<ALE::Mesh::topology_type::sheaf_type>& patches    = m->getTopology()->getPatches();
  const Obj<ALE::Mesh::real_section_type>&         disp       = m->getRealSection("displacement");
  const Obj<ALE::Mesh::real_section_type>&         full       = m->getRealSection("full_displacement");

  for(ALE::Mesh::topology_type::sheaf_type::iterator p_iter = patches->begin(); p_iter != patches->end(); ++p_iter) {
    const ALE::Obj<ALE::Mesh::topology_type::label_sequence>& vertices = m->getTopology()->depthStratum(p_iter->first, 0);

    for(ALE::Mesh::topology_type::label_sequence::iterator v_iter = vertices->begin(); v_iter != vertices->end(); ++v_iter) {
      const ALE::Mesh::foliated_section_type::value_type *constVal = boundaries->restrict(p_iter->first, *v_iter);
      const ALE::Mesh::real_section_type::value_type     *array    = disp->restrict(p_iter->first, *v_iter);
      const int numConst = boundaries->size(p_iter->first, *v_iter);
      const int dim      = disp->getFiberDimension(p_iter->first, *v_iter);
      double    values[3];
      int       v = 0;

      for(int c = 0; c < 3; c++) {
        int i;

        for(i = 0; i < numConst; i++) {
          if (constVal[i].first == c) {
            values[c] = constVal[i].second;
            break;
          }
        }
        if (i == numConst) {
          values[c] = array[v++];
        }
      }
      if (v != dim) {
        std::cout << "ERROR: Invalid size " << v << " used for " << *v_iter << " with index " << disp->getIndex(p_iter->first, *v_iter) << std::endl;
      }
      full->updateAdd(p_iter->first, *v_iter, values);
    }
  }
  ierr = MeshGetSectionReal(mesh, "full_displacement", fullDisplacement);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

char pypylith3d_outputMesh__doc__[] = "";
char pypylith3d_outputMesh__name__[] = "outputMesh";

PyObject * pypylith3d_outputMesh(PyObject *, PyObject *args)
{
  using ALE::Obj;
  PyObject *pyMesh, *pySol;
  char     *meshBaseFile;

  int ok = PyArg_ParseTuple(args, (char *) "sOO:outputMesh", &meshBaseFile, &pyMesh, &pySol);
  if (!ok) {
    return 0;
  }

  Mesh           mesh = (Mesh) PyCObject_AsVoidPtr(pyMesh);
  Vec            sol  = (Vec)  PyCObject_AsVoidPtr(pySol);
  MPI_Comm       comm;
  SectionReal    displacement, fullDisplacement;
  std::string    filename(meshBaseFile);
  PetscViewer    viewer;
  //Vec       partition;
  PetscTruth     hasMaterial;
  PetscErrorCode ierr;

  ierr = PetscObjectGetComm((PetscObject) mesh, &comm);
  ierr = MeshGetSectionReal(mesh, "displacement", &displacement);
  ierr = updateDisplacement(displacement, sol);
  ierr = createFullDisplacement(mesh, &fullDisplacement);
  ierr = PetscViewerCreate(comm, &viewer);
  ierr = PetscViewerSetType(viewer, PETSC_VIEWER_ASCII);
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_VTK);
  filename += ".vtk";
  ierr = PetscViewerFileSetName(viewer, filename.c_str());
  ierr = MeshView(mesh, viewer);
  ierr = SectionRealView(fullDisplacement, viewer);
  ierr = MeshHasSectionInt(mesh, "material", &hasMaterial);
  if (hasMaterial) {
    SectionInt material;

    ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_VTK_CELL);
    ierr = MeshGetSectionInt(mesh, "material", &material);
    ierr = SectionIntView(material, viewer);
    ierr = PetscViewerPopFormat(viewer);
    ierr = SectionIntDestroy(material);
  }
  //SectionIntView(partition, viewer);
  ierr = PetscViewerDestroy(viewer);
  ierr = SectionRealDestroy(displacement);
  ierr = SectionRealDestroy(fullDisplacement);

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "Output PETSc Mesh and Solution"
    << journal::endl;

  // return Py_None;
  Py_INCREF(Py_None);
  return Py_None;
}

#if 0
#include "Numeric/arrayobject.h"

char pypylith3d_interpolatePoints__doc__[] = "";
char pypylith3d_interpolatePoints__name__[] = "interpolatePoints";

PyObject * pypylith3d_interpolatePoints(PyObject *, PyObject *args)
{
  using ALE::Obj;
  PyObject      *pyMesh, *pySol;
  PyArrayObject *pyPoints;

  int ok = PyArg_ParseTuple(args, (char *) "OOO!:interpolatePoints", &pyMesh, &pySol, &PyArray_Type, &pyPoints);
  if (!ok) {
    return 0;
  }
  if ((pyPoints->nd != 2) || (pyPoints->descr->type_num != PyArray_DOUBLE)) {
    PyErr_SetString(PyExc_ValueError, "points must be a 2d array with double values");
    return 0;
  }
  if (pyPoints->dimensions[1] != 3) {
    PyErr_SetString(PyExc_ValueError, "points must be a 3d");
    return 0;
  }
  if ((pyPoints->strides[0] != 3 * sizeof(double)) || (pyPoints->strides[1] != sizeof(double))) {
    PyErr_SetString(PyExc_ValueError, "points must be a contiguous array");
    return 0;
  }

  Mesh           mesh = (Mesh) PyCObject_AsVoidPtr(pyMesh);
  Vec            sol  = (Vec)  PyCObject_AsVoidPtr(pySol);
  SectionReal    displacement, fullDisplacement;
  const int      numPoints = pyPoints->dimensions[0];
  double        *values;
  PetscErrorCode ierr;

  ierr = MeshGetSectionReal(mesh, "displacement", &displacement);
  ierr = updateDisplacement(displacement, sol);
  ierr = createFullDisplacement(mesh, &fullDisplacement);
  ierr = MeshInterpolatePoints(mesh, fullDisplacement, numPoints, (double *) pyPoints->data, &values);
  ierr = SectionRealDestroy(displacement);
  ierr = PetscFree(values);

  int            dims[2]  = {numPoints, 3};
  PyArrayObject *pyValues = (PyArrayObject *) PyArray_FromDims(2, dims, PyArray_DOUBLE);
  double        *data     = (double *) pyValues->data;

  for(int p = 0; p < numPoints; ++p) {
    for(int d = 0; d < 3; d++) {
      data[p*3+d] = values[p*3+d];
    }
  }

  ierr = PetscFree(values);
  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "Interpolated points"
    << journal::endl;

  return Py_BuildValue((char *) "N", pyValues);
}
#endif

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
    
// version
// $Id: scanner.cc,v 1.7 2005/06/07 19:39:11 willic3 Exp $

// End of file
