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

/** @file libsrc/problems/SolutionFactory.hh
 *
 * @brief C++ helper class for setting up solution fields for testing.
 */

#if !defined(pylith_problems_solutionfactory_hh)
#define pylith_problems_solutionfactory_hh

// Include directives ---------------------------------------------------
#include "problemsfwd.hh" // forward declarations
#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/topology/FieldBase.hh" // USES FieldBase::Descretization

#include "spatialdata/units/unitsfwd.hh" // HASA Nondimensional
#include "spatialdata/spatialdb/spatialdbfwd.hh" // HASA Nondimensional

// SolutionFactory-----------------------------------------------
/// @brief C++ helper class for setting up solution subfields for unit tests.
class pylith::problems::SolutionFactory : public pylith::utils::GenericComponent {
    friend class TestSolutionFactory; // unit testing

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

    /** Add displacement subfield to solution field.
     *
     * @param[in] discretization Discretization for displacement subfield.
     */
    void addDisplacement(const pylith::topology::FieldBase::Discretization& discretization);

    /** Add velocity subfield to solution field.
     *
     * @param[in] discretization Discretization for velocity subfield.
     */
    void addVelocity(const pylith::topology::FieldBase::Discretization& discretization);

    /** Add pressure subfield to solution field.
     *
     * @param[in] discretization Discretization for pressure subfield.
     */
    void addPressure(const pylith::topology::FieldBase::Discretization& discretization);

    /** Add trace strain subfield to solution field.
     *
     * @param[in] discretization Discretization for trace strain subfield.
     */
    void addTraceStrain(const pylith::topology::FieldBase::Discretization& discretization);

    /** Add fault Lagrange multiplier subfield to solution field.
     *
     * @param[in] discretization Discretization for fault Lagrange multiplier subfield.
     */
    void addLagrangeMultiplierFault(const pylith::topology::FieldBase::Discretization& discretization);

    /** Add temperature subfield to solution field.
     *
     * @param[in] discretization Discretization for temperature subfield.
     */
    void addTemperature(const pylith::topology::FieldBase::Discretization& discretization);

    /** Allocate and populate subfield with values using spatial database.
     *
     * @param[in] db Spatial database.
     */
    void setValues(spatialdata::spatialdb::SpatialDB* db);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////
private:

    pylith::topology::Field& _solution; ///< Solution field.
    const spatialdata::units::Nondimensional& _normalizer; ///< Nondimensionalizer.
    const int _spaceDim; ///< Spatal dimension of problem.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    SolutionFactory(const SolutionFactory &); ///< Not implemented.
    const SolutionFactory& operator=(const SolutionFactory&); ///< Not implemented

}; // class SolutionFactory

#endif // pylith_problems_solutionfactory_hh

// End of file
