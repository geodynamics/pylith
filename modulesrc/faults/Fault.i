// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

/** @file modulesrc/faults/Fault.i
 *
 * @brief Python interface to C++ Fault object.
 */

namespace pylith {
  namespace faults {

    class Fault
    { // class Fault

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      Fault(void);

      /// Destructor.
      virtual
      ~Fault(void);
      
      /** Set material identifier of fault.
       *
       * @param value Fault identifier
       */
      void id(const int value);
      
      /** Get material identifier of fault.
       *
       * @returns Fault identifier
       */
      int id(void) const;
      
      /** Set label of group of vertices associated with fault.
       *
       * @param value Label of fault
       */
      void label(const char* value);
      
      /** Get label of group of vertices associated with fault.
       *
       * @returns Label of fault
       */
      const char* label(void) const;
      
      /** Adjust mesh topology for fault implementation.
       *
       * @param mesh PETSc mesh
       */
      virtual
      void adjustTopology(pylith::topology::Mesh* mesh,
			  const bool flipFault =false) = 0;
      
      /** Initialize fault. Determine orientation and setup boundary
       * condition parameters.
       *
       * @param mesh PETSc mesh
       * @param cs Coordinate system for mesh
       * @param upDir Direction perpendicular to along-strike direction that is 
       *   not collinear with fault normal (usually "up" direction but could 
       *   be up-dip direction; only applies to fault surfaces in a 3-D domain).
       * @param normalDir General preferred direction for fault normal
       *   (used to pick which of two possible normal directions for
       *   interface; only applies to fault surfaces in a 3-D domain).
       * @param matDB Database of bulk elastic properties for fault region
       *   (used to improve conditioning of Jacobian matrix)
       */
      virtual
      void initialize(const pylith::topology::Mesh& mesh,
		      const double upDir[3],
		      const double normalDir[3],
		      spatialdata::spatialdb::SpatialDB* matDB) = 0;
      
      /** Get mesh associated with fault fields.
       *
       * @returns PETSc mesh object
       */
      const pylith::topology::SubMesh& faultMesh(void) const;
      
      /** Get vertex field associated with integrator.
       *
       * @param name Name of vertex field.
       * @param fields Solution fields.
       * @returns Vertex field.
       */
      virtual
      const pylith::topology::Field<pylith::topology::SubMesh>&
      vertexField(const char* name,
		  const pylith::topology::SolutionFields& fields) = 0;
      
      /** Get cell field associated with integrator.
       *
       * @param name Name of cell field.
       * @param fields Solution fields.
       * @returns Cell field.
       */
      virtual
      const pylith::topology::Field<pylith::topology::SubMesh>&
      cellField(const char* name,
		const pylith::topology::SolutionFields& fields) = 0;
      
    }; // class Fault
    
  } // faults
} // pylith


// End of file 
