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

#include "mesh.h"

#include "config.h"

#include <Distribution.hh>
#include <petscmesh.h>
#include <src/dm/mesh/meshpylith.h>
#include <petscmat.h>
#include "journal/debug.h"

#include <Python.h>

#include <stdio.h>
#include <string.h>

#include "PyLithMeshLib.h"

#undef __FUNCT__
#define __FUNCT__ "IgnoreComments_PyLith"
static PetscErrorCode IgnoreComments_PyLith(char *buf, PetscInt bufSize, FILE *f)
{
  PetscFunctionBegin;
  while((fgets(buf, bufSize, f) != NULL) && ((buf[0] == '#') || (buf[0] == '\0'))) {}
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ReadBoundary_PyLith"
static PetscErrorCode ReadBoundary_PyLith(const char *baseFilename, PetscTruth useZeroBase, PetscInt *numBoundaryVertices, PetscInt *numBoundaryComponents, PetscInt **boundaryVertices, PetscScalar **boundaryValues)
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
  // ierr = PetscStrcat(bcFilename, ".bc");
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
static PetscErrorCode WriteBoundary_PyLith(const char *baseFilename, const ALE::Obj<ALECompat::Mesh>& mesh)
{
  ALECompat::Mesh::foliated_section_type::patch_type      patch      = 0;
  const ALE::Obj<ALECompat::Mesh::numbering_type>&        vNumbering = mesh->getFactory()->getLocalNumbering(mesh->getTopology(), patch, 0);
  const ALE::Obj<ALECompat::Mesh::foliated_section_type>& boundaries = mesh->getBoundariesNew();
  FILE          *f;
  char           bcFilename[2048];
  PetscErrorCode ierr;
  MPI_Comm          comm = PETSC_COMM_WORLD;
  PetscMPIInt       rank;
  char           suff[9];

  ierr = MPI_Comm_rank(comm, &rank);
  sprintf(suff, "%s%d%s", ".", rank, ".bc");

  PetscFunctionBegin;
  if (mesh->debug()) {
    boundaries->view("PyLith boundaries");
  }
  int slen = std::strlen(baseFilename);
  std::strncpy(bcFilename, baseFilename, slen-3);
  bcFilename[slen-3] = '\0';
  // ierr = PetscStrcpy(bcFilename, baseFilename.substr(0,slen-4));
  ierr = PetscStrcat(bcFilename, suff);

  // Determine if we have bc stuff
  const ALE::Obj<ALECompat::Mesh::topology_type::label_sequence>& vertices = boundaries->getTopology()->depthStratum(patch, 0);
  bool haveBC = false;
  for(ALECompat::Mesh::topology_type::label_sequence::iterator v_iter = vertices->begin(); v_iter != vertices->end(); ++v_iter)
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

  for(ALECompat::Mesh::topology_type::label_sequence::iterator v_iter = vertices->begin(); v_iter != vertices->end(); ++v_iter) {
    int    constraints[3] = {0, 0, 0};
    double values[3] = {0.0, 0.0, 0.0};
    int    size = boundaries->getFiberDimension(patch, *v_iter);
    const ALECompat::Mesh::foliated_section_type::value_type *array = boundaries->restrict(patch, *v_iter);

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

PyObject *PyLithMeshLib::Mesh::_processMesh(PyMeshObject *self)
{
  using ALE::Obj;
  journal::debug_t  debug("pylith3d");
  MPI_Comm          comm = PETSC_COMM_WORLD;
  PetscMPIInt       rank;
  Obj<ALECompat::Mesh>    m;
  PetscViewer       viewer;
  PetscInt         *boundaryVertices;
  PetscScalar      *boundaryValues;
  PetscInt          numBoundaryVertices, numBoundaryComponents;
  int               debugFlag = 0;
  PetscErrorCode    ierr;

  ierr = MPI_Comm_rank(comm, &rank);
  ierr = MeshCompatCreatePyLith(comm, 3, self->meshInputFile, PETSC_FALSE, (PetscTruth) self->interpolateMesh, &self->mesh);
  ierr = MeshCompatGetMesh(self->mesh, m);
  m->setDebug(debugFlag);
  int numElements = m->getTopology()->heightStratum(0, 0)->size();
  debug << journal::at(__HERE__) << "[" << rank << "]Created new PETSc Mesh for " << self->meshInputFile << journal::endl;
  m = ALECompat::New::Distribution<ALECompat::Mesh::topology_type>::distributeMesh(m, self->partitioner);
  ierr = MeshCompatSetMesh(self->mesh, m);
  debug << journal::at(__HERE__) << "[" << rank << "]Distributed PETSc Mesh"  << journal::endl;
  ierr = ReadBoundary_PyLith(self->meshBcFile, PETSC_FALSE, &numBoundaryVertices, &numBoundaryComponents, &boundaryVertices, &boundaryValues);

  const Obj<ALECompat::Mesh::foliated_section_type>& boundaries = m->getBoundariesNew();
  ALECompat::Mesh::foliated_section_type::patch_type patch      = 0;
  std::set<int> seen;

  ierr = MPI_Bcast(&numElements, 1, MPI_INT, 0, comm);
  boundaries->setTopology(m->getTopology());
  // Reverse order allows newer conditions to override older, as required by PyLith
  for(int v = numBoundaryVertices-1; v >= 0; v--) {
    ALECompat::Mesh::point_type vertex(boundaryVertices[v*(numBoundaryComponents+1)] + numElements);
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
    ALECompat::Mesh::point_type vertex(boundaryVertices[v*(numBoundaryComponents+1)] + numElements);
    ALECompat::Mesh::foliated_section_type::value_type values[3];

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

  ierr = PetscViewerCreate(comm, &viewer);
  ierr = PetscViewerSetType(viewer, PETSC_VIEWER_ASCII);
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_PYLITH_LOCAL);
  ierr = PetscViewerFileSetMode(viewer, FILE_MODE_READ);
  ierr = PetscExceptionTry1(PetscViewerFileSetName(viewer, self->meshInputFile), PETSC_ERR_FILE_OPEN);
  if (PetscExceptionValue(ierr)) {
    /* this means that a caller above me has also tryed this exception so I don't handle it here, pass it up */
  } else if (PetscExceptionCaught(ierr, PETSC_ERR_FILE_OPEN)) {
    ierr = 0;
  } 
  ierr = MeshView(self->mesh, viewer);
  ierr = PetscViewerDestroy(viewer);

  char           bcFilename[2048];
  char           suff[9];

  sprintf(suff, "%s%d%s", ".", rank, ".bc");

  int slen = std::strlen(self->meshBcFile);
  std::strncpy(bcFilename, self->meshBcFile, slen-3);
  // ierr = PetscStrcpy(bcFilename, self->meshBcFile.substr(0,slen-4));
  ierr = PetscStrcat(bcFilename, suff);

  debug << journal::at(__HERE__) << "[" << rank << "]Output new PyLith mesh into: " << bcFilename << journal::endl;

  ierr = WriteBoundary_PyLith(self->meshBcFile, m);
  debug << journal::at(__HERE__) << "[" << rank << "]Wrote PyLith boundary conditions"  << journal::endl;

  const Obj<ALECompat::Mesh::topology_type::label_sequence>& vertices = m->getTopology()->depthStratum(0, 0);
  Obj<ALECompat::Mesh::real_section_type> s = m->getRealSection("default");

  m->setRealSection("displacement", s);
  s->setFiberDimensionByDepth(0, 0, 3);
  for(ALECompat::Mesh::topology_type::label_sequence::iterator v_iter = vertices->begin(); v_iter != vertices->end(); ++v_iter) {
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
  debug << journal::at(__HERE__) << "[" << rank << "]Created displacement Field"  << journal::endl;

  m->getFactory()->constructInverseOrder(m->getFactory()->getLocalNumbering(m->getTopology(), 0, m->getTopology()->depth()));

  Py_INCREF(Py_None);
  return Py_None;
}

// Create a PETSc Mat

PyObject *PyLithMeshLib::Mesh::_createMat(PyMeshObject *self)
{
  using ALE::Obj;
  MPI_Comm comm = PETSC_COMM_WORLD;
  PetscErrorCode ierr;

  Obj<ALECompat::Mesh> m;

  ierr = MeshCompatGetMesh(self->mesh, m);
  const ALE::Obj<ALECompat::Mesh::order_type>& offsets = m->getFactory()->getGlobalOrder(m->getTopology(), 0, "displacement", m->getRealSection("displacement")->getAtlas());
  int localSize = offsets->getLocalSize();
  int globalSize = offsets->getGlobalSize();

  if (MatCreate(comm, &self->A)) {
    PyErr_SetString(PyExc_RuntimeError, "Could not create PETSc Mat");
    return 0;
  }
  if (MatSetSizes(self->A, localSize, localSize, globalSize, globalSize)) {
    PyErr_SetString(PyExc_RuntimeError, "Could not set sizes for PETSc Mat");
    return 0;
  }
  if (VecCreate(comm, &self->rhs)) {
    PyErr_SetString(PyExc_RuntimeError, "Could not create PETSc Rhs");
    return 0;
  }
  if (VecSetSizes(self->rhs, localSize, globalSize)) {
    PyErr_SetString(PyExc_RuntimeError, "Could not set sizes for PETSc Rhs");
    return 0;
  }
  if (VecSetFromOptions(self->rhs)) {
    PyErr_SetString(PyExc_RuntimeError, "Could not set options for PETSc Rhs");
    return 0;
  }
  if (VecDuplicate(self->rhs, &self->sol)) {
    PyErr_SetString(PyExc_RuntimeError, "Could not create PETSc Sol");
    return 0;
  }

  ierr = PetscObjectCompose((PetscObject) self->A, "mesh", (PetscObject) self->mesh);

  VecScatter injection = NULL;
  ierr = MeshCompatGetGlobalScatter(self->mesh, &injection);
  ierr = PetscObjectCompose((PetscObject) self->rhs, "mesh",      (PetscObject) self->mesh);
  ierr = PetscObjectCompose((PetscObject) self->rhs, "injection", (PetscObject) injection);
  ierr = PetscObjectCompose((PetscObject) self->sol, "mesh",      (PetscObject) self->mesh);
  ierr = PetscObjectCompose((PetscObject) self->sol, "injection", (PetscObject) injection);

  ierr = MatSetFromOptions(self->A);
  ierr = preallocateMatrixCompat(m->getTopology(), m->getRealSection("displacement")->getAtlas(), m->getFactory()->getGlobalOrder(m->getTopology(), 0, "displacement", m->getRealSection("displacement")->getAtlas()), self->A);

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "Created PETSc Mat: " << localSize << " " << globalSize
    << journal::endl;

  Py_INCREF(Py_None);
  return Py_None;
}

// Destroy a PETSc Mat

PyObject *PyLithMeshLib::Mesh::_destroyMat(PyMeshObject *self)
{
  if (MatDestroy(self->A)) {
    PyErr_SetString(PyExc_RuntimeError, "Could not destroy PETSc Mat");
    return 0;
  }
  self->A = 0;
  
  if (VecDestroy(self->rhs)) {
    PyErr_SetString(PyExc_RuntimeError, "Could not destroy PETSc Rhs");
    return 0;
  }
  self->rhs = 0;
  
  if (VecDestroy(self->sol)) {
    PyErr_SetString(PyExc_RuntimeError, "Could not destroy PETSc Sol");
    return 0;
  }
  self->sol = 0;

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "Destroyed PETSc Mat"
    << journal::endl;

  Py_INCREF(Py_None);
  return Py_None;
}

PyObject *PyLithMeshLib::Mesh::_destroyMesh(PyMeshObject *self)
    
{
    if (MeshDestroy(self->mesh)) {
        PyErr_SetString(PyExc_RuntimeError, "Could not destroy PETSc Mesh");
        return 0;
    }
    self->mesh = 0;
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
  ierr = VecScatterBegin(injection, sol, lv, INSERT_VALUES, SCATTER_REVERSE);CHKERRQ(ierr);
  ierr = VecScatterEnd(injection, sol, lv, INSERT_VALUES, SCATTER_REVERSE);CHKERRQ(ierr);
  ierr = VecDestroy(lv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// Create complete displacement field by adding BC
PetscErrorCode createFullDisplacement(Mesh mesh, SectionReal *fullDisplacement) {
  using ALE::Obj;
  Obj<ALECompat::Mesh> m;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MeshCompatGetMesh(mesh, m);CHKERRQ(ierr);
  if (!m->hasRealSection("full_displacement")) {
    const Obj<ALECompat::Mesh::topology_type::sheaf_type>& patches    = m->getTopology()->getPatches();
    const Obj<ALECompat::Mesh::real_section_type>&         full       = m->getRealSection("full_displacement");

    for(ALECompat::Mesh::topology_type::sheaf_type::iterator p_iter = patches->begin(); p_iter != patches->end(); ++p_iter) {
      full->setFiberDimensionByDepth(p_iter->first, 0, 3);
    }
    full->allocate();
  }
  const Obj<ALECompat::Mesh::foliated_section_type>&     boundaries = m->getBoundariesNew();
  const Obj<ALECompat::Mesh::topology_type::sheaf_type>& patches    = m->getTopology()->getPatches();
  const Obj<ALECompat::Mesh::real_section_type>&         disp       = m->getRealSection("displacement");
  const Obj<ALECompat::Mesh::real_section_type>&         full       = m->getRealSection("full_displacement");

  for(ALECompat::Mesh::topology_type::sheaf_type::iterator p_iter = patches->begin(); p_iter != patches->end(); ++p_iter) {
    const ALE::Obj<ALECompat::Mesh::topology_type::label_sequence>& vertices = m->getTopology()->depthStratum(p_iter->first, 0);

    for(ALECompat::Mesh::topology_type::label_sequence::iterator v_iter = vertices->begin(); v_iter != vertices->end(); ++v_iter) {
      const ALECompat::Mesh::foliated_section_type::value_type *constVal = boundaries->restrict(p_iter->first, *v_iter);
      const ALECompat::Mesh::real_section_type::value_type     *array    = disp->restrict(p_iter->first, *v_iter);
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

PyObject *PyLithMeshLib::Mesh::_outputMesh(PyMeshObject *self, char *fileRoot)
{
  using ALE::Obj;

  MPI_Comm       comm;
  SectionReal    displacement, fullDisplacement;
  std::string    filename(fileRoot);
  PetscViewer    viewer;
  //Vec       partition;
  PetscTruth     hasMaterial;
  PetscErrorCode ierr;

#if 0
  ierr = PetscObjectGetComm((PetscObject) self->mesh, &comm);
  ierr = MeshGetSectionReal(self->mesh, "displacement", &displacement);
  ierr = updateDisplacement(displacement, self->sol);
  ierr = createFullDisplacement(self->mesh, &fullDisplacement);
  ierr = PetscViewerCreate(comm, &viewer);
  ierr = PetscViewerSetType(viewer, PETSC_VIEWER_ASCII);
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_VTK);
  filename += ".vtk";
  ierr = PetscViewerFileSetName(viewer, filename.c_str());
  ierr = MeshView(self->mesh, viewer);
  ierr = SectionRealView(fullDisplacement, viewer);
  ierr = MeshHasSectionInt(self->mesh, "material", &hasMaterial);
  if (hasMaterial) {
    SectionInt material;

    ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_VTK_CELL);
    ierr = MeshGetSectionInt(self->mesh, "material", &material);
    ierr = SectionIntView(material, viewer);
    ierr = PetscViewerPopFormat(viewer);
    ierr = SectionIntDestroy(material);
  }
  //SectionIntView(partition, viewer);
  ierr = PetscViewerDestroy(viewer);
  ierr = SectionRealDestroy(displacement);
  ierr = SectionRealDestroy(fullDisplacement);
#endif

  journal::debug_t debug("pylith3d");
  debug
    << journal::at(__HERE__)
    << "DISABLED: Output PETSc Mesh and Solution"
    << journal::endl;

  Py_INCREF(Py_None);
  return Py_None;
}


// End of file
