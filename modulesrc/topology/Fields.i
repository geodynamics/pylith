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
// Copyright (c) 2010-2015 University of California, Davis
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

      /** Add field.
       *
       * @param name Name of field.
       * @param label Label for field.
       * @param domain Type of points over which to define field.
       * @param fiberDim Fiber dimension for field.
       */
      void add(const char* name,
	       const char* label,
	       const pylith::topology::FieldBase::DomainEnum domain,
	       const int fiberDim);

      /** Delete field.
       *
       * @param name Name of field.
       */
      void delField(const char* name);

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
	   
      /** Copy layout to other fields.
       *
       * @param name Name of field to use as template for layout.
       */
      void copyLayout(const char* name);
      
      /** Get mesh associated with fields.
       *
       * @returns Finite-element mesh.
       */
      const pylith::topology::Mesh& mesh(void) const;

      /** Return the names of all fields.
       *
       * @param numValues Number of field values [output].
       * @param values Values of field values [output].
       */
      void fieldNames(int* numValues, 
		      char*** values) const;

    }; // Fields

  } // topology
} // pylith


// End of file 
