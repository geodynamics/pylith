// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

/** @file tests/libtests/meshio/FieldFactory.hh
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
#include "pylith/topology/topologyfwd.hh" // HOLDSA Field
#include "pylith/topology/FieldBase.hh" // USES FieldBase::Descretization

// FieldFactory-----------------------------------------------
/// @brief C++ helper class for setting up auxiliary fields for materials.
class pylith::meshio::FieldFactory {
    friend class TestFieldFactory; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /** Default constructor.
     *
     * @param[inout] fields Container with fields.
     */
    FieldFactory(pylith::topology::Field& field);

    /// Destructor.
    ~FieldFactory(void);

    /** Add scalar field and set field values.
     *
     * @param[in] discretization Discretization for scalar field.
     */
    void addScalar(const pylith::topology::FieldBase::Discretization& discretization);

    /** Add vector field.
     *
     * @param[in] discretization Discretization for vector field.
     */
    void addVector(const pylith::topology::FieldBase::Discretization& discretization);

    /** Add tensor field.
     *
     * @param[in] discretization Discretization for tensor field.
     */
    void addTensor(const pylith::topology::FieldBase::Discretization& discretization);

    /** Add other field.
     *
     * @param[in] discretization Discretization for other field.
     */
    void addOther(const pylith::topology::FieldBase::Discretization& discretization);

    /** Set values in field.
     *
     * @param[in] values Array of values for all subfields.
     * @param[in] numPoints Number of points associated with values.
     * @param[in] numDOF Total number of values per point.
     */
    void setValues(const PylithScalar* values,
                   const PylithInt numPoints,
                   const PylithInt numDOF);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////
private:

    pylith::topology::Field& _field; ///< Field witn subfields.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    FieldFactory(const FieldFactory &); ///< Not implemented.
    const FieldFactory& operator=(const FieldFactory&); ///< Not implemented

}; // class FieldFactory

#endif // pylith_problems_fieldfactory_hh

// End of file
