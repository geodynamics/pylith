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
 * @file modulesrc/topology/Distributor.i
 *
 * @brief Python interface to C++ Distributor object.
 */

namespace pylith {
  namespace topology {

    class Distributor
    { // Distributor

      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :

      /// Constructor
      Distributor(void);
      
      /// Destructor
      ~Distributor(void);
      
      /** Distribute mesh among processors.
       *
       * @param newMesh Distributed mesh (result).
       * @param origMesh Mesh to distribute.
   * @param partitionerName Name of PETSc partitioner to use in distributing mesh.
       */
      static
      void distribute(pylith::topology::Mesh* const newMesh,
		      const pylith::topology::Mesh& origMesh,
		      const char* partitionerName);

      /** Write partitioning info for distributed mesh.
       *
       * @param writer Data writer for partition information.
       * @param mesh Distributed mesh.
       */
      static
      void write(pylith::meshio::DataWriter* const writer,
		 const pylith::topology::Mesh& mesh);

    }; // Distributor

  } // topology
} // pylith


// End of file 
