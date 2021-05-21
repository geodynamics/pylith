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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file modulesrc/topology/Fields.i
 *
 * @brief Python interface to C++ Fields object.
 */

namespace pylith {
  namespace topology {

    class Fields
    { // Fields

      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :

      /** Default constructor.
       *
       * @param mesh Finite-element mesh.
       */
      Fields(const pylith::topology::Mesh& mesh);

      /// Destructor.
      virtual
      ~Fields(void);

      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);

      /** Check if fields contains a given field.
       *
       * @param name Name of field.
       * @return True if fields contains field, false otherwise.
       */
      bool hasField(const char* name) const;

      /** Add field.
       *
       * @param name Name of field.
       * @param label Label for field.
       */
      void add(const char* name,
	       const char* label);

      /** Get field.
       *
       * @param name Name of field.
       */
      const pylith::topology::Field& get(const char* name) const;
	   
      /** Get field.
       *
       * @param name Name of field.
       */
      pylith::topology::Field& get(const char* name);
	   
      /** Get mesh associated with fields.
       *
       * @returns Finite-element mesh.
       */
      const pylith::topology::Mesh& mesh(void) const;

    }; // Fields

  } // topology
} // pylith


// End of file 
