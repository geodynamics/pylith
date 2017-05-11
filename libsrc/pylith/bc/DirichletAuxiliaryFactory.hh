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

/** @file libsrc/bc/DirichletAuxiliaryFactory.hh
 *
 * @brief C++ helper class for setting up auxiliary fields for boundary conditions.
 */

#if !defined(pylith_bc_dirichletauxiliaryfactory_hh)
#define pylith_bc_dirichletauxiliaryfactory_hh

// Include directives ---------------------------------------------------
#include "bcfwd.hh" // forward declarations
#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/topology/FieldBase.hh" // USES FieldBase::Discretization

// DirichletAuxiliaryFactory-----------------------------------------------
/// @brief C++ helper class for setting up auxiliary fields for bnoundary conditions.
class pylith::bc::DirichletAuxiliaryFactory : public pylith::utils::GenericComponent {
    friend class TestDirichletAuxiliaryFactory;   // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /** Default constructor.
     *
     * @param[inout] bc Dirichlet boundary condition.
     * @param[in] solution Solution field.
     * @param[in] timeScale Time scale associated with nondimensionalization.
     */
    DirichletAuxiliaryFactory(const DirichletNew& bc,
			      const pylith::topology::Field& solution,
	const PylithReal timeScale);

    /// Destructor.
    ~DirichletAuxiliaryFactory(void);

    /// Add initial amplitude field to auxiliary fields.
    void initialAmplitude(void) const;

    /// Add rate amplitude field to auxiliary fields.
    void rateAmplitude(void) const;

    /// Add rate start time amplitude field to auxiliary fields.
    void rateStartTime(void) const;

    /// Add time history amplitude field to auxiliary fields.
    void timeHistoryAmplitude(void) const;

    /// Add time history start time field to auxiliary fields.
    void timeHistoryStartTime(void) const;

    /// Add time history value field to auxiliary fields.
    void timeHistoryValue(void) const;

    // PRIVATE MEMBERS ////////////////////////////////////////////////////
private:

    const DirichletNew& _bc; ///< Boundary condition with auxiliary fields.
    const pylith::topology::FieldBase::Description& _description; ///< Description for constrained field.
    const int _spaceDim; ///< Spatial dimension.
    const PylithReal _timeScale; ///< Time scale associated with nondimensionalization.

    static const char* _genericComponent; ///< Name of generic component.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    DirichletAuxiliaryFactory(const DirichletAuxiliaryFactory &);   ///< Not implemented.
    const DirichletAuxiliaryFactory& operator=(const DirichletAuxiliaryFactory&);   ///< Not implemented

}; // class DirichletAuxiliaryFactory

#endif // pylith_bc_dirichletauxiliaryfactory_hh


// End of file
