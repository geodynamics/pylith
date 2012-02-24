// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
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
      SubMesh(const pylith::topology::Mesh& mesh,
	      const char* label);

      /// Default destructor
      ~SubMesh(void);
      
      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);

      /** Create Sieve mesh.
       *
       * @param mesh Finite-element mesh over domain.
       * @param label Label for vertices marking boundary.
       */
      void createSubMesh(const pylith::topology::Mesh& mesh,
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

      /** Get representative cone size for mesh.
       *
       * @returns Representative cone size for mesh.
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
