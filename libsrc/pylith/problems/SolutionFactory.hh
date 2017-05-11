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

/** @file libsrc/problems/SolutionFactory.hh
 *
 * @brief C++ helper class for setting up solution fields for testing.
 */

#if !defined(pylith_problems_solutionfactory_hh)
#define pylith_materials_solutionfactory_hh

// Include directives ---------------------------------------------------
#include "problemsfwd.hh" // forward declarations
#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/topology/FieldBase.hh" // USES FieldBase::Descretization

#include "spatialdata/units/unitsfwd.hh" // HASA Nondimensional
#include "spatialdata/spatialdb/spatialdbfwd.hh" // HASA Nondimensional

// SolutionFactory-----------------------------------------------
/// @brief C++ helper class for setting up auxiliary fields for materials.
class pylith::problems::SolutionFactory : public pylith::utils::GenericComponent {
    friend class TestSolutionFactory;   // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /** Default constructor.
     *
     * @param[inout] solution Solution field.
     * @param[in] normalizer Nondimensionalizer for problem.
     * @param[in] spaceDim Spatial dimension of problem.
     */
    SolutionFactory(pylith::topology::Field& solution,
                    const spatialdata::units::Nondimensional& normalizer);

    /// Destructor.
    ~SolutionFactory(void);

    /** Add displacement field to solution field.
     *
     * @param[in] discretization Discretization for displacement field.
     */
    void displacement(const pylith::topology::FieldBase::Discretization& discretization);

    /** Add velocity field to solution field.
     *
     * @param[in] discretization Discretization for velocity field.
     */
    void velocity(const pylith::topology::FieldBase::Discretization& discretization);

    /** Add pressure field to solution field.
     *
     * @param[in] discretization Discretization for pressure field.
     */
    void pressure(const pylith::topology::FieldBase::Discretization& discretization);

    /** Add fluid pressure field to solution field.
     *
     * @param[in] discretization Discretization for fluid pressure field.
     */
    void fluidPressure(const pylith::topology::FieldBase::Discretization& discretization);

    /** Add temperature field to solution field.
     *
     * @param[in] discretization Discretization for temperature field.
     */
    void temperature(const pylith::topology::FieldBase::Discretization& discretization);

    /** Add time derivative of displacement field to solution field.
     *
     * @param[in] discretization Discretization for time derivative of displacement field.
     */
    void displacementDot(const pylith::topology::FieldBase::Discretization& discretization);

    /** Add time derivative of velocity field to solution field.
     *
     * @param[in] discretization Discretization for time derivative of velocity field.
     */
    void velocityDot(const pylith::topology::FieldBase::Discretization& discretization);

    /** Add time derivative of pressure field to solution field.
     *
     * @param[in] discretization Discretization for time derivative of pressure field.
     */
    void pressureDot(const pylith::topology::FieldBase::Discretization& discretization);

    /** Add time derivative of fluid pressure field to solution field.
     *
     * @param[in] discretization Discretization for time derivative of fluid pressure field.
     */
    void fluidPressureDot(const pylith::topology::FieldBase::Discretization& discretization);

    /** Add time derivative of temperature field to solution field.
     *
     * @param[in] discretization Discretization for time derivative of temperature field.
     */
    void temperatureDot(const pylith::topology::FieldBase::Discretization& discretization);

    /** Allocate and populate field with values using spatial database.
     *
     * @param[in] db Spatial database.
     */
    void setValues(spatialdata::spatialdb::SpatialDB* db);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////
private:

    pylith::topology::Field& _solution; ///< Solution field.
    const spatialdata::units::Nondimensional& _normalizer; ///< Nondimensionalizer.
    const int _spaceDim; ///< Spatal dimension of problem.

    static const char* _genericComponent; ///< Name of Generic component.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    SolutionFactory(const SolutionFactory &);   ///< Not implemented.
    const SolutionFactory& operator=(const SolutionFactory&);   ///< Not implemented

}; // class SolutionFactory

#endif // pylith_problems_solutionfactory_hh


// End of file
