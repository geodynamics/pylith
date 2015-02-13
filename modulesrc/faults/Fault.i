// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
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
      
      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
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

      /** Set label of group of vertices defining buried edge of fault.
       *
       * @param value Label of fault
       */
      void edge(const char* value);
      
      /** Get label of group of vertices defining buried edge of fault.
       *
       * @returns Label of fault
       */
      const char* edge(void) const;

      /** Get dimension of mesh.
       *
       * @returns Dimension of mesh.
       */
      int dimension(void) const;

      /** Get number of vertices per cell for mesh.
       *
       * @returns Number of vertices per cell for mesh.
       */
      int numCorners(void) const;
  
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

      /** Get the number of vertices associated with the fault (before
       * fault mesh exists).
       *
       * @param mesh PETSc mesh
       * @return Number of vertices on the fault.
       */
      virtual
      int numVerticesNoMesh(const pylith::topology::Mesh& mesh) const = 0;

      /** Adjust mesh topology for fault implementation.
       *
       * @param mesh PETSc mesh
       */
      %apply int *INOUT {int *firstFaultVertex, int *firstFaultCell};
      virtual
      void adjustTopology(pylith::topology::Mesh* const mesh,
                          int *firstFaultVertex,
                          int *firstLagrangeVertex,
                          int *firstFaultCell) = 0;
      %clear int *firstFaultVertex, int *firstFaultCell;
      
      /** Initialize fault. Determine orientation and setup boundary
       * condition parameters.
       *
       * @param mesh PETSc mesh
       * @param cs Coordinate system for mesh
       * @param upDir Direction perpendicular to along-strike direction that is 
       *   not collinear with fault normal (usually "up" direction but could 
       *   be up-dip direction; applies to fault surfaces in 2-D and 3-D).
       */
      virtual
      void initialize(const pylith::topology::Mesh& mesh,
		      const PylithScalar upDir[3]) = 0;
      
      /** Get mesh associated with fault fields.
       *
       * @returns PETSc mesh object
       */
      const pylith::topology::Mesh& faultMesh(void) const;
      
      /** Get vertex field associated with integrator.
       *
       * @param name Name of vertex field.
       * @param fields Solution fields.
       * @returns Vertex field.
       */
      virtual
      const pylith::topology::Field& vertexField(const char* name,
						 const pylith::topology::SolutionFields* fields =0) = 0;
      
      /** Get cell field associated with integrator.
       *
       * @param name Name of cell field.
       * @param fields Solution fields.
       * @returns Cell field.
       */
      virtual
      const pylith::topology::Field& cellField(const char* name,
					       const pylith::topology::SolutionFields* fields =0) = 0;
      
    }; // class Fault
    
  } // faults
} // pylith


// End of file 
