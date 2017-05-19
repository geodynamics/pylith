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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file unittests/libtests/meshio/FieldFactory.hh
 *
 * @brief C++ helper class for setting up solution fields for testing.
 */

#if !defined(pylith_problems_fieldfactory_hh)
#define pylith_problems_fieldfactory_hh

/// Namespace for pylith package
namespace pylith {
    namespace meshio {
        class FieldFactory;
    } // meshio
} // pylith

// Include directives ---------------------------------------------------
#include "pylith/topology/topologyfwd.hh" // HOLDSA Fields
#include "pylith/topology/FieldBase.hh" // USES FieldBase::Descretization

// FieldFactory-----------------------------------------------
/// @brief C++ helper class for setting up auxiliary fields for materials.
class pylith::meshio::FieldFactory {
    friend class TestFieldFactory;   // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /** Default constructor.
     *
     * @param[inout] fields Container with fields.
     */
    FieldFactory(pylith::topology::Fields& fields);

    /// Destructor.
    ~FieldFactory(void);

    /** Add scalar field and set field values.
     *
     * @param[in] discretization Discretization for scalar field.
     * @param[in] values Array of values for field (numPoints, numComponents).
     * @param[in] numPoints Number of points associated with field.
     * @param[in] numComponents Number of values per point.
     */
    void scalar(const pylith::topology::FieldBase::Discretization& discretization,
                const PylithScalar* values,
                const PylithInt numPoints,
                const PylithInt numComponents);

    /** Add vector field.
     *
     * @param[in] discretization Discretization for vector field.
     * @param[in] values Array of values for field (numPoints, numComponents).
     * @param[in] numPoints Number of points associated with field.
     * @param[in] numComponents Number of values per point.
     */
    void vector(const pylith::topology::FieldBase::Discretization& discretization,
                const PylithScalar* values,
                const PylithInt numPoints,
                const PylithInt numComponents);

    /** Add tensor field.
     *
     * @param[in] discretization Discretization for tensor field.
     * @param[in] values Array of values for field (numPoints, numComponents).
     * @param[in] numPoints Number of points associated with field.
     * @param[in] numComponents Number of values per point.
     */
    void tensor(const pylith::topology::FieldBase::Discretization& discretization,
                const PylithScalar* values,
                const PylithInt numPoints,
                const PylithInt numComponents);

    /** Add other field.
     *
     * @param[in] discretization Discretization for other field.
     * @param[in] values Array of values for field (numPoints, numComponents).
     * @param[in] numPoints Number of points associated with field.
     * @param[in] numComponents Number of values per point.
     */
    void other(const pylith::topology::FieldBase::Discretization& discretization,
               const PylithScalar* values,
               const PylithInt numPoints,
               const PylithInt numComponents);

    /** Set values in field.
     *
     * @param field Field to populate.
     * @param values Array of values.
     * @param numPoints Number of points associated with values.
     * @param numComponents Number of values per point.
     */
    void setField(pylith::topology::Field* field,
                  const PylithScalar* values,
                  const PylithInt numPoints,
                  const PylithInt numComponents);


    // PRIVATE MEMBERS ////////////////////////////////////////////////////
private:

    pylith::topology::Fields& _fields; ///< Contrainer with fields.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    FieldFactory(const FieldFactory &);   ///< Not implemented.
    const FieldFactory& operator=(const FieldFactory&);   ///< Not implemented

}; // class FieldFactory

#endif // pylith_problems_fieldfactory_hh


// End of file
