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
 * @file libsrc/topology/Fields.hh
 *
 * @brief Container for managing multiple fields over a finite-element
 * mesh.
 */

#if !defined(pylith_topology_fields_hh)
#define pylith_topology_fields_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

#include "pylith/topology/FieldBase.hh" // USES FieldBase::DomainEnum

#include <string> // USES std::string
#include <map> // USES std::map

// Fields ---------------------------------------------------------------
/// Container for managing multiple fields over a finite-element mesh.
class pylith::topology::Fields { // Fields
    friend class TestFieldsMesh; // unit testing
    friend class TestFieldsSubmesh; // unit testing

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /** Default constructor.
     *
     * @param mesh Finite-element mesh.
     */
    Fields(const Mesh& mesh);

    /// Destructor.
    virtual ~Fields(void);

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
     * @param label Label for field (used in output).
     */
    void add(const char* name,
             const char* label);

    /** Delete field.
     *
     * @param name Name of field.
     */
    void del(const char* name);

    /** Get field.
     *
     * @param name Name of field.
     */
    const Field& get(const char* name) const;

    /** Get field.
     *
     * @param name Name of field.
     */
    Field& get(const char* name);

    /** Get mesh associated with fields.
     *
     * @returns Finite-element mesh.
     */
    const Mesh& mesh(void) const; // :TODO: @brad Remove obsolete.

    // PROTECTED TYPEDEFS ///////////////////////////////////////////////////
protected:

    typedef std::map< std::string, Field* > map_type;

    // PROTECTED MEMBERS ////////////////////////////////////////////////////
protected:

    map_type _fields;
    const Mesh& _mesh;

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    Fields(const Fields&); ///< Not implemented
    const Fields& operator=(const Fields&); ///< Not implemented

}; // Fields

#endif // pylith_topology_fields_hh

// End of file
