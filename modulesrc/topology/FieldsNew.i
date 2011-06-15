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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file modulesrc/topology/FieldsNew.i
 *
 * @brief Python interface to C++ FieldsNew object.
 */

namespace pylith {
  namespace topology {

    template<typename mesh_type>
    class FieldsNew
    { // FieldsNew

      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :

      /** Default constructor.
       *
       * @param mesh Finite-element mesh.
       */
      FieldsNew(const mesh_type& mesh);

      /// Destructor.
      virtual
      ~FieldsNew(void);

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
       * @param fiberDim Fiber dimension for field.
       */
      void add(const char* name,
	       const char* label,
	       const int fiberDim,
	       FieldBase::VectorFieldEnum vectorFieldType =FieldBase::OTHER,
	       const double scale =1.0,
	       const bool dimsOkay =false);

      /** Create and allocate Sieve section.
       *
       * @param points Points over which to define section.
       */
      void allocate(const ALE::Obj<typename mesh_type::SieveMesh::label_sequence>& points);

      /** Create and allocate Sieve section.
       *
       * @param points Points over which to define section.
       */
      void allocate(const pylith::int_array& points);

      /** Create and allocate Sieve section.
       *
       * @param domain Type of points over which to define section.
       * @param stratum Stratum depth (for vertices) and height (for cells).
       */
      void allocate(const FieldBase::DomainEnum domain,
		    const int stratum =0);

      /** Get field.
       *
       * @param name Name of field.
       * @returns Field.
       */
      Field<mesh_type>& get(const char* name);
	   
      /** Get mesh associated with fields.
       *
       * @returns Finite-element mesh.
       */
      const mesh_type& mesh(void) const;

      /** Get section containing fields.
       *
       * @returns Sieve section
       */
      const ALE::Obj<typename mesh_type::RealUniformSection>& section(void) const;

      /** Compute total fiber dimension for section.
       *
       * @returns Fiber dimension.
       */
      int fiberDim(void) const;

      /** Get index of first value of field in section.
       *
       * @param name Name of field.
       * @returns Index of first value of field in section.
       */
      int sectionIndex(const char* name) const;

      /** Get fiber dimension of field in section.
       *
       * @param name Name of field.
       * @returns Fiber dimension of field in section.
       */
      int sectionFiberDim(const char* name) const;

      /// Complete section by assembling across processors.
      void complete(void);

      /** Return the names of all fields.
       *
       * @param numValues Number of fields,
       * @param values Names of fields.
       */
      void fieldNames(int* numValues,
		      char*** values) const;

      /** View fields and section.
       *
       * @param label Label for fields.
       */
      void view(const char* label);

    }; // FieldsNew

  } // topology
} // pylith


// End of file 
