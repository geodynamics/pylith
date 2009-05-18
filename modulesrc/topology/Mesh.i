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

      /// Default constructor.
      Mesh(void);

      /** Constructor with dimension and communicator.
       *
       * @param dim Dimension associated with mesh cells.
       * @param comm MPI communicator for mesh.
       */
      Mesh(const int dim,
	   const MPI_Comm& comm =PETSC_COMM_WORLD); 

      /// Default destructor
      ~Mesh(void);
      
      /** Create Sieve mesh.
       *
       * @param dim Dimension associated with mesh cells.
       */
      void createSieveMesh(const int dim=3); 
      
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

      /** Get representative cone size of mesh.
       *
       * @returns Representative cone size of mesh.
       */
      int coneSize(void) const;
      
      /** Get number of vertices in mesh.
       *
       * @returns Number of vertices in mesh.
       */
      int numVertices(void) const;
      
      /** Get number of cells in mesh.
       *
       * @returns Number of cells in mesh.
       */
      int numCells(void) const;
      
      /** Set MPI communicator associated with mesh.
       *
       * @param value MPI communicator.
       */
      void comm(const MPI_Comm value);
    
      /** Get MPI communicator associated with mesh.
       *
       * @returns MPI communicator.
       */
      const MPI_Comm comm(void) const;
    
      /** Initialize the finite-element mesh.
       *
       * @param normalizer Nondimensionalizer.
       */
      void nondimensionalize(const spatialdata::units::Nondimensional& normalizer);

      /** Print mesh to stdout.
       *
       * @param label Label for mesh.
       */
      void view(const char* label);

      /** Return the names of all vertex groups.
       *
       * @returns an array of all vertex group names
       */
      void groups(int *numNames, char ***outNames);

      /** Return the size of a group.
       *
       * @returns the number of vertices in the group
       */
      int groupSize(const char *name);

    }; // Mesh

  } // topology
} // pylith


// End of file
