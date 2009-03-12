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
       * @param fields Fields associated with mesh and solution of the problem.
       */
      Jacobian(const SolutionFields& fields);
      
      /// Destructor.
      ~Jacobian(void);
      
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
      
      /** Assemble matrix.
       *
       * @param mode Assembly mode.
       */
      void assemble(const char* mode);
      
      /// Set entries in matrix to zero (retain structure).
      void zero(void);
      
      /// View matrix to stdout.
      void view(void);
      
      /** Write matrix to binary file.
       *
       * @param filename Name of file.
       */
      void write(const char* filename);

    }; // Jacobian

  } // topology
} // pylith

// End of file 
