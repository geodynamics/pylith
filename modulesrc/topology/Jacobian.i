// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

/**
 * @file modulesrc/topology/Jacobian.i
 *
 * @brief Python interface to C++ Jacobian object.
 */

namespace pylith {
  namespace topology {

    class Jacobian
    { // Jacobian

      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :

      /** Default constructor.
       *
       * @param field Field associated with mesh and solution of the problem.
       * @param matrixType Type of PETSc sparse matrix.
       * @param blockOkay True if okay to use block size equal to fiberDim
       * (all or none of the DOF at each point are constrained).
       */
      Jacobian(const Field& field,
	       const char* matrixType ="aij",
	       const bool blockOkay =false);

      /// Destructor.
      ~Jacobian(void);
      
      /// Deallocate PETSc and local data structures.
      void deallocate(void);

      /** Get PETSc matrix.
       *
       * @returns PETSc sparse matrix.
       */
      const PetscMat* matrix(void) const;
      
      /** Get PETSc matrix.
       *
       * @returns PETSc sparse matrix.
       */
      PetscMat* matrix(void);
      
      /** Get matrix type.
       *
       * @returns Matrix type.
       */
      const char* matrixType(void) const;

      /** Assemble matrix.
       *
       * @param mode Assembly mode.
       */
      void assemble(const char* mode);
      
      /// Set entries in matrix to zero (retain structure).
      void zero(void);
      
      /// View matrix to stdout.
      void view(void) const;
      
      /** Write matrix to binary file.
       *
       * @param filename Name of file.
       * @param comm MPI communicator.
       */
      void write(const char* filename,
		 const MPI_Comm comm);

      /// Verify symmetry of matrix. For debugger purposes only.
      void verifySymmetry(void) const;

    }; // Jacobian

  } // topology
} // pylith

// End of file 
