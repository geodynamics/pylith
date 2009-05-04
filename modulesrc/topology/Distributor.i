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
       * @param partitioner Name of partitioner to use in distributing mesh.
       */
      static
      void distribute(pylith::topology::Mesh* const newMesh,
		      const pylith::topology::Mesh& origMesh,
		      const char* partitioner);

      /** Write partitioning info for distributed mesh.
       *
       * @param writer Data writer for partition information.
       * @param mesh Distributed mesh.
       * @param cs Coordinate system for mesh.
       */
      static
      void write(pylith::meshio::DataWriter<pylith::topology::Mesh, pylith::topology::Field<pylith::topology::Mesh> >* const writer,
		 const pylith::topology::Mesh& mesh);

    }; // Distributor

  } // topology
} // pylith


// End of file 
