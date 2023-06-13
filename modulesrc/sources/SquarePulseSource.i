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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file modulesrc/sources/SquarePulseSource.i
 *
 * Python interface to C++ SquarePulseSource.
 */

namespace pylith {
    namespace sources {
        class SquarePulseSource : public pylith::sources::Source {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            SquarePulseSource(void);

            /// Destructor.
            ~SquarePulseSource(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Set time history database.
             *
             * @param[in] db Time history database.
             */
            void setTimeHistoryDB(spatialdata::spatialdb::TimeHistory* th);

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

            /** Get time history database.
             *
             * @preturns Time history database.
             */
            const spatialdata::spatialdb::TimeHistory* getTimeHistoryDB(void);

            /** Verify configuration is acceptable.
             *
             * @param[in] solution Solution field.
             */
            void verifyConfiguration(const pylith::topology::Field& solution) const;

            /** Create integrator and set kernels.
             *
             * @param[in] solution Solution field.
             *
             *  @returns Integrator if applicable, otherwise NULL.
             */
            pylith::feassemble::Integrator* createIntegrator(const pylith::topology::Field& solution);

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

            // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////
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

        }; // class SquarePulseSource

    } // sources
} // pylith

// End of file