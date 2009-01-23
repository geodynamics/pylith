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
 * @file modulesrc/topology/Mesh.hh
 *
 * @brief Python interface to C++ PyLith Mesh object.
 */

namespace pylith {
  namespace topology {

    // Mesh -------------------------------------------------------------
    class Mesh
    { // Mesh

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /** Default constructor.
       *
       * @param comm MPI communicator for mesh.
       * @param dim Dimension associated with mesh cells.
       */
      Mesh(const MPI_Comm& comm =PETSC_COMM_WORLD,
	   const int dim =3); 

      /// Default destructor
      ~Mesh(void);

      /** Set coordinate system.
       *
       * @param cs Coordinate system.
       */
      void coordsys(const spatialdata::geocoords::CoordSys* cs);
      
      /** Get coordinate system.
       *
       * @returns Coordinate system.
       */
      const spatialdata::geocoords::CoordSys* coordsys(void) const;
      
      /** Set debug flag.
       *
       * @param value Turn on debugging if true.
       */
      void debug(const bool value);
      
      /** Get debug flag.
       *
       * @param Get debugging flag.
       */
      bool debug(void) const;

      /** Get dimension of mesh.
       *
       * @returns Dimension of mesh.
       */
      int dimension(void) const;
      
      /** Get MPI communicator associated with mesh.
       *
       * @returns MPI communicator.
       */
      const MPI_Comm comm(void) const;
    
      /// Initialize the finite-element mesh.
      void initialize(void);
      
      /** Print mesh to stdout.
       *
       * @param label Label for mesh.
       */
      void view(const char* label);

    }; // Mesh

  } // topology
} // pylith


// End of file
