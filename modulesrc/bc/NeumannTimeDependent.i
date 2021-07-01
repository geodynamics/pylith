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

/** @file modulesrc/bc/NeumannTimeDependent.i
 *
 * @brief Python interface to C++ NeumannTimeDependent object.
 */

namespace pylith {
    namespace bc {
        class NeumannTimeDependent : public pylith::bc::BoundaryCondition {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            NeumannTimeDependent(void);

            /// Destructor.
            ~NeumannTimeDependent(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Set time history database.
             *
             * @param[in] db Time history database.
             */
            void setTimeHistoryDB(spatialdata::spatialdb::TimeHistory* th);

            /** Get time history database.
             *
             * @preturns Time history database.
             */
            const spatialdata::spatialdb::TimeHistory* getTimeHistoryDB(void);

            /** Use initial value term in time history expression.
             *
             * @param[in] value True if using initial value term in expression.
             */
            void useInitial(const bool value);

            /** Get flag associated with using initial value term in time history expression.
             *
             * @returns True if using initial value term in expression, false otherwise.
             */
            bool useInitial(void) const;

            /** Use rate value term in time history expression.
             *
             * @param[in] value True if using rate value term in expression.
             */
            void useRate(const bool value);

            /** Get flag associated with using rate value term in time history expression.
             *
             * @returns True if using rate value term in expression, false otherwise.
             */
            bool useRate(void) const;

            /** Use time history term in time history expression.
             *
             * @param[in] value True if using time history term in expression.
             */
            void useTimeHistory(const bool value);

            /** Get flag associated with using time history term in time history expression.
             *
             * @returns True if using time history term in expression, false otherwise.
             */
            bool useTimeHistory(void) const;

            /** Set name of scale associated with Neumann boundary
             * condition (e.g., 'pressure' for elasticity).
             *
             * A Neumann boundary condition constrains the gradient in
             * a solution subfield. In some cases the constraint is
             * actually on a scaled version of the gradient as is the
             * case of a Neumann boundary condition for elasticity
             * that constrains boundary tractions.
             *
             * @param value Name of scale for nondimensionalizing Neumann boundary condition.
             */
            void setScaleName(const char* value);

            /** Create integrator and set kernels.
             *
             * @param[in] solution Solution field.
             * @returns Integrator if applicable, otherwise NULL.
             */
            pylith::feassemble::Integrator* createIntegrator(const pylith::topology::Field& solution);

            /** Create constraint and set kernels.
             *
             * @param[in] solution Solution field.
             * @returns Constraint if applicable, otherwise NULL.
             */
            pylith::feassemble::Constraint* createConstraint(const pylith::topology::Field& solution);

            /** Create auxiliary field.
             *
             * @param[in] solution Solution field.
             * @param[in\ domainMesh Finite-element mesh associated with integration domain.
             *
             * @returns Auxiliary field if applicable, otherwise NULL.
             */
            pylith::topology::Field* createAuxiliaryField(const pylith::topology::Field& solution,
                                                          const pylith::topology::Mesh& domainMesh);

            /** Create derived field.
             *
             * @param[in] solution Solution field.
             * @param[in\ domainMesh Finite-element mesh associated with integration domain.
             *
             * @returns Derived field if applicable, otherwise NULL.
             */
            pylith::topology::Field* createDerivedField(const pylith::topology::Field& solution,
                                                        const pylith::topology::Mesh& domainMesh);

            /** Update auxiliary subfields at beginning of time step.
             *
             * @param[out] auxiliaryField Auxiliary field.
             * @param[in] t Current time.
             */
            void updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                                      const double t);

            // PROTECTED METHODS
            // ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

            /** Get auxiliary factory associated with physics.
             *
             * @return Auxiliary factory for physics object.
             */
            pylith::feassemble::AuxiliaryFactory* _getAuxiliaryFactory(void);

            /** Update kernel constants.
             *
             * @param[in] dt Current time step.
             */
            void _updateKernelConstants(const PylithReal dt);

        };

        // NeumannTimeDependent

    } // bc
} // pylith

// End of file
