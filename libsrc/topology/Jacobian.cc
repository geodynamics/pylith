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
pylith::topology::Jacobian::Jacobian(const SolutionFields& fields,
				     const char* matrixType,
				     const bool blockOkay) :
  _fields(fields),
  _matrix(0)
{ // constructor
  const ALE::Obj<Mesh::SieveMesh>& sieveMesh = fields.mesh().sieveMesh();
  const ALE::Obj<Mesh::RealSection>& solnSection = fields.solution().section();
  ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
  logger.stagePush("Jacobian");

  // Set blockFlag to -1 if okay to set block size equal to fiber
  // dimension, otherwise use a block size of 1.
  const int blockFlag = (blockOkay) ? -1 : 1;

  PetscErrorCode err = MeshCreateMatrix(sieveMesh, solnSection, 
					matrixType, &_matrix, blockFlag);
  CHECK_PETSC_ERROR_MSG(err, "Could not create PETSc sparse matrix "
			"associated with system Jacobian.");
  logger.stagePop();
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

// ----------------------------------------------------------------------
// Verify symmetry of matrix.
void
pylith::topology::Jacobian::verifySymmetry(void) const
{ // verifySymmetry
  const PetscMat matSparse = _matrix;

  int nrows = 0;
  int ncols = 0;
  MatGetSize(matSparse, &nrows, &ncols);

  PetscMat matDense;
  PetscMat matSparseAIJ;
  MatConvert(matSparse, MATSEQAIJ, MAT_INITIAL_MATRIX, &matSparseAIJ);
  MatConvert(matSparseAIJ, MATSEQDENSE, MAT_INITIAL_MATRIX, &matDense);

  double_array vals(nrows*ncols);
  int_array rows(nrows);
  int_array cols(ncols);
  for (int iRow=0; iRow < nrows; ++iRow)
    rows[iRow] = iRow;
  for (int iCol=0; iCol < ncols; ++iCol)
    cols[iCol] = iCol;
  MatGetValues(matDense, nrows, &rows[0], ncols, &cols[0], &vals[0]);
  const double tolerance = 1.0e-06;
  bool isSymmetric = true;
  for (int iRow=0; iRow < nrows; ++iRow)
    for (int iCol=0; iCol < ncols; ++iCol) {
      const int indexIJ = ncols*iRow+iCol;
      const int indexJI = nrows*iCol+iRow;
      const double valIJ = vals[indexIJ];
      const double valJI = vals[indexJI];
      if (fabs(valIJ) > 1.0)
	if (fabs(1.0 - valJI/valIJ) > tolerance) {
	  std::cerr << "Mismatch: " 
		    << "(" << iRow << ", " << iCol << ") = " << valIJ
		    << ", (" << iCol << ", " << iRow << ") = " << valJI
		    << std::endl;
	  isSymmetric = false;
	} // if
      else
	if (fabs(valJI - valIJ) > tolerance) {
	  std::cerr << "Mismatch: " 
		    << "(" << iRow << ", " << iCol << ") = " << valIJ
		    << ", (" << iCol << ", " << iRow << ") = " << valJI
		    << std::endl;
	  isSymmetric = false;
	} // if
    } // for
  MatDestroy(matDense);
  MatDestroy(matSparseAIJ);
  if (!isSymmetric)
    throw std::runtime_error("Jacobian matrix is not symmetric.");
} // verifySymmetry


// End of file 
