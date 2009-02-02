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
 * @file modulesrc/topology/SubMesh.i
 *
 * @brief Python interface to C++ Mesh object.
 */

namespace pylith {
  namespace topology {

    class SubMesh
    { // SubMesh

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      SubMesh(void);
      
      /** Constructor with mesh and label for vertices marking boundary.
       *
       * @param mesh Finite-element mesh over domain.
       * @param label Label for vertices marking boundary.
       */
      SubMesh(const Mesh& mesh,
	      const char* label);

      /// Default destructor
      ~SubMesh(void);
      
      /** Create Sieve mesh.
       *
       * @param mesh Finite-element mesh over domain.
       * @param label Label for vertices marking boundary.
       */
      void createSubMesh(const Mesh& mesh,
			 const char* label); 

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
      
    }; // SubMesh

  } // topology
} // pylith


// End of file
