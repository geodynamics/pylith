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

/** @file modulesrc/sources/TimeHistoryWavelet.i
 *
 * Python interface to C++ TimeHistoryWavelet.
 */

namespace pylith {
    namespace sources {
        class TimeHistoryWavelet: public pylith::sources::SourceTimeFunctionMomentTensorForce {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            TimeHistoryWavelet(void);

            /// Destructor.
            ~TimeHistoryWavelet(void);

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

            /** Get auxiliary factory associated with physics.
             *
             * @return Auxiliary factory for physics object.
             */
            pylith::sources::AuxiliaryFactoryMomentTensorForce* getAuxiliaryFactory(void);

            /** Add source time function subfields to auxiliary field.
             *
             * @param[inout] auxiliaryField Auxiliary field.
             */
            void addAuxiliarySubfields(void);

            /** Get g1v kernel for residual, G(t,s).
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return residual kernel for g1v
             * .
             */
            PetscPointFunc getKernelg1v_explicit(const spatialdata::geocoords::CoordSys* coordsys) const;

            /** Update auxiliary field for current time.
             *
             * @param[inout] auxiliaryField Auxiliary field to update.
             * @param[in] t Current time.
             * @param[in] timeScale Time scale for nondimensionalization.
             */
            void updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                                      const PylithReal t,
                                      const PylithReal timeScale);

        };

        // class TimeHistoryWavelet

    } // sources
} // pylith

// End of file
