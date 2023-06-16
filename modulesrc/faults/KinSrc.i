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
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/** @file modulesrc/faults/KinSrc.i
 *
 * @brief Python interface to C++ KinSrc object.
 */

namespace pylith {
    namespace faults {
        class KinSrc : public pylith::utils::PyreComponent {
            // PUBLIC METHODS /////////////////////////////////////////////////
public:

            /// Default constructor.
            KinSrc(void);

            /// Destructor.
            ~KinSrc(void);

            /// Deallocate PETSc and local data structures.
            virtual
            void deallocate(void);

            /** Set origin time for earthquake source.
             *
             * @param value Origin time for earthquake source.
             */
            void setOriginTime(const PylithReal value);

            /** Get origin time for earthquake source.
             *
             * @returns Origin time for earthquake source.
             */
            PylithReal getOriginTime(void) const;

            /** Get auxiliary field associated with the kinematic source.
             *
             * @return field Auxiliary field for the kinematic source.
             */
            const pylith::topology::Field& auxField(void) const;

            /** Set the spatial database for filling auxiliary subfields.
             *
             * @param[in] value Pointer to database.
             */
            void auxFieldDB(spatialdata::spatialdb::SpatialDB* value);

            /** Initialize kinematic (prescribed slip) earthquake source.
             *
             * @param[in] auxField Auxiliary field associated with fault finite-element integration.
             * @param[in] normalizer Normalizer for nondimensionalizing values.
             * @param[in] cs Coordinate system for problem.
             */
            void initialize(const pylith::topology::Field& auxField,
                            const spatialdata::units::Nondimensional& normalizer,
                            const spatialdata::geocoords::CoordSys* cs);

            /** Get requested slip subfields at time t.
             *
             * @param[inout] slipLocalVec Local PETSc vector for slip, slip rate, or slip accelerationvalues.
             * @param[in] faultAuxiliaryField Auxiliary field for fault.
             * @param[in] t Time t.
             * @param[in] timeScale Time scale for nondimensionalization.
             * @param[in] bitSlipSubfields Slip subfields to compute.
             */
            virtual
            void getSlipSubfields(PetscVec slipLocalVec,
                                  pylith::topology::Field* faultAuxiliaryField,
                                  const PylithScalar t,
                                  const PylithScalar timeScale,
                                  const int bitSlipSubfields);

            // PROTECTED METHODS //////////////////////////////////////////////////
protected:

            /** Setup auxiliary subfields (discretization and query fns).
             *
             * Create subfields in auxiliary fields (includes name of the field,
             * vector field type, discretization, and scale for
             * nondimensionalization) and set query functions for filling them
             * from a spatial database.
             *
             * @attention The order of the calls to subfieldAdd() must match the
             * order of the auxiliary fields in the FE kernels.
             *
             * @param[in] normalizer Normalizer for nondimensionalizing values.
             * @param[in] cs Coordinate system for problem.
             */
            virtual
            void _auxiliaryFieldSetup(const spatialdata::units::Nondimensional& normalizer,
                                      const spatialdata::geocoords::CoordSys* cs) = 0;

        }; // class KinSrc

    } // faults
} // pylith

// End of file
