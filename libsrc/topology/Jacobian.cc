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

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::Jacobian::Jacobian(const SolutionFields& fields) :
  _fields(fields),
  _matrix(0)
{ // constructor
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = fields.mesh().sieveMesh();
  const ALE::Obj<Mesh::RealSection>& solnSection = fields.solution().section();

  _matrix = new PetscMat;
  assert(0 != _matrix);
  PetscErrorCode err = MeshCreateMatrix(sieveMesh, solnSection, 
					MATAIJ, _matrix);
  if (err) {
    PetscError(__LINE__,__FUNCT__,__FILE__,__SDIR__,err,0," ");
    throw std::runtime_error("Could not create PETSc sparse matrix "
			     "associated with system Jacobian.");
  } // if
  

} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::topology::Jacobian::~Jacobian(void)
{ // destructor
  MatDestroy(*_matrix);
  delete _matrix; _matrix = 0;
} // destructor

// ----------------------------------------------------------------------
// Get PETSc matrix.
const PetscMat*
pylith::topology::Jacobian::matrix(void) const
{ // matrix
  return _matrix;
} // matrix

// ----------------------------------------------------------------------
// Get PETSc matrix.
PetscMat*
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
    err = MatAssemblyBegin(*_matrix, MAT_FINAL_ASSEMBLY);
    if (err) {
      PetscError(__LINE__,__FUNCT__,__FILE__,__SDIR__,err,0," ");
      throw std::runtime_error("Error beginning final assembly of sparse "
			       "matrix associated with system Jacobian.");
    } // if
    err = MatAssemblyEnd(*_matrix, MAT_FINAL_ASSEMBLY);
    if (err) {
      PetscError(__LINE__,__FUNCT__,__FILE__,__SDIR__,err,0," ");
      throw std::runtime_error("Error ending final assembly of sparse "
			       "matrix associated with system Jacobian.");
    } // if
  } else if (0 == strcmp(mode, "flush_assembly")) {
    err = MatAssemblyBegin(*_matrix, MAT_FLUSH_ASSEMBLY);
    if (err) {
      PetscError(__LINE__,__FUNCT__,__FILE__,__SDIR__,err,0," ");
      throw std::runtime_error("Error beginning flush assembly of sparse "
			       "matrix associated with system Jacobian.");
    } // if
    err = MatAssemblyEnd(*_matrix, MAT_FLUSH_ASSEMBLY);
    if (err) {
      PetscError(__LINE__,__FUNCT__,__FILE__,__SDIR__,err,0," ");
      throw std::runtime_error("Error ending flush assembly of sparse "
			       "matrix associated with system Jacobian.");
    } // if
  } else
    throw std::runtime_error("Unknown mode for assembly of sparse matrix "
			     "associated with system Jacobian.");
} // assemble

// ----------------------------------------------------------------------
// Set entries in matrix to zero (retain structure).
void
pylith::topology::Jacobian::zero(void)
{ // zero
  PetscErrorCode err = MatZeroEntries(*_matrix);
  if (err) {
    PetscError(__LINE__,__FUNCT__,__FILE__,__SDIR__,err,0," ");
    throw std::runtime_error("Error zeroing entries of sparse matrix "
			     "associated with system Jacobian.");
  } // if
} // zero

// ----------------------------------------------------------------------
// View matrix to stdout.
void
pylith::topology::Jacobian::view(void)
{ // view
  PetscErrorCode err = MatView(*_matrix, PETSC_VIEWER_STDOUT_WORLD);
  if (err) {
    PetscError(__LINE__,__FUNCT__,__FILE__,__SDIR__,err,0," ");
    throw std::runtime_error("Error viewing sparse matrix associatd "
			     "with system Jacobian.");
  } // if
} // view

// ----------------------------------------------------------------------
// Write matrix to binary file.
void
pylith::topology::Jacobian::write(const char* filename)
{ // write
  PetscViewer viewer;

  const MPI_Comm comm = _fields.mesh().comm();

  PetscErrorCode err = 
    PetscViewerBinaryOpen(comm, filename,
			  FILE_MODE_WRITE, &viewer);
  if (err) {
    PetscError(__LINE__,__FUNCT__,__FILE__,__SDIR__,err,0," ");
    throw std::runtime_error("Could not create PETSc binary viewer for "
			     "sparse matrix associated with system Jacobian.");
  } // if

  err = MatView(*_matrix, viewer);
  if (err) {
    PetscError(__LINE__,__FUNCT__,__FILE__,__SDIR__,err,0," ");
    throw std::runtime_error("Could not view PETSc sparse matrix associated "
			     "with system Jacobian.");
  } // if

  err = PetscViewerDestroy(viewer);
  if (err) {
    PetscError(__LINE__,__FUNCT__,__FILE__,__SDIR__,err,0," ");
    throw std::runtime_error("Could not destroy PETSc binary viewer for "
			     "sparse matrix associated with system Jacobian.");
  } // if

} // write


// End of file 
