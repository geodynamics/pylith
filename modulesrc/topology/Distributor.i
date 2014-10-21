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
// Copyright (c) 2010-2014 University of California, Davis
//
// See COPYING for license information.
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
       */
      static
      void distribute(pylith::topology::Mesh* const newMesh,
		      const pylith::topology::Mesh& origMesh);

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
