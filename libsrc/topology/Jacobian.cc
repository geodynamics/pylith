// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include <portinfo>

#include "Jacobian.hh" // implementation of class methods

#include "Mesh.hh" // USES Mesh
#include "SolutionFields.hh" // USES SolutionFields
#include "Field.hh" // USES Field

#include "pylith/utils/petscerror.h" // USES CHECK_PETSC_ERROR

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::Jacobian::Jacobian(const SolutionFields& fields) :
  _fields(fields),
  _matrix(0)
{ // constructor
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = fields.mesh().sieveMesh();
  const ALE::Obj<Mesh::RealSection>& solnSection = fields.solution().section();

  PetscErrorCode err = MeshCreateMatrix(sieveMesh, solnSection, 
					MATAIJ, &_matrix);
  CHECK_PETSC_ERROR_MSG(err, "Could not create PETSc sparse matrix "
			"associated with system Jacobian.");
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::topology::Jacobian::~Jacobian(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::topology::Jacobian::deallocate(void)
{ // deallocate
  if (0 != _matrix) {
    PetscErrorCode err = MatDestroy(_matrix); _matrix = 0;
    CHECK_PETSC_ERROR(err);
  } // if
} // deallocate

// ----------------------------------------------------------------------
// Get PETSc matrix.
const PetscMat
pylith::topology::Jacobian::matrix(void) const
{ // matrix
  return _matrix;
} // matrix

// ----------------------------------------------------------------------
// Get PETSc matrix.
PetscMat
pylith::topology::Jacobian::matrix(void)
{ // matrix
  return _matrix;
} // matrix

// ----------------------------------------------------------------------
// Assemble matrix.
void
pylith::topology::Jacobian::assemble(const char* mode)
{ // assemble
  PetscErrorCode err = 0;
  if (0 == strcmp(mode, "final_assembly")) {
    err = MatAssemblyBegin(_matrix, MAT_FINAL_ASSEMBLY); 
    CHECK_PETSC_ERROR(err);
    err = MatAssemblyEnd(_matrix, MAT_FINAL_ASSEMBLY);
    CHECK_PETSC_ERROR(err);
  } else if (0 == strcmp(mode, "flush_assembly")) {
    err = MatAssemblyBegin(_matrix, MAT_FLUSH_ASSEMBLY);
    CHECK_PETSC_ERROR(err);
    err = MatAssemblyEnd(_matrix, MAT_FLUSH_ASSEMBLY);
    CHECK_PETSC_ERROR(err);
  } else
    throw std::runtime_error("Unknown mode for assembly of sparse matrix "
			     "associated with system Jacobian.");
} // assemble

// ----------------------------------------------------------------------
// Set entries in matrix to zero (retain structure).
void
pylith::topology::Jacobian::zero(void)
{ // zero
  PetscErrorCode err = MatZeroEntries(_matrix);
  CHECK_PETSC_ERROR(err);
} // zero

// ----------------------------------------------------------------------
// View matrix to stdout.
void
pylith::topology::Jacobian::view(void)
{ // view
  PetscErrorCode err = MatView(_matrix, PETSC_VIEWER_STDOUT_WORLD);
  CHECK_PETSC_ERROR(err);
} // view

// ----------------------------------------------------------------------
// Write matrix to binary file.
void
pylith::topology::Jacobian::write(const char* filename)
{ // write
  PetscViewer viewer;

  const MPI_Comm comm = _fields.mesh().comm();

  PetscErrorCode err = 
    PetscViewerBinaryOpen(comm, filename, FILE_MODE_WRITE, &viewer);
  CHECK_PETSC_ERROR(err);

  err = MatView(_matrix, viewer); CHECK_PETSC_ERROR(err);
  err = PetscViewerDestroy(viewer); CHECK_PETSC_ERROR(err);
} // write


// End of file 
