// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** @file modulesrc/materials/IsotropicLinearElasticity.i
 *
 * Python interface to C++ IsotropicLinearElasticity.
 */

namespace pylith {
    namespace materials {
        class IsotropicLinearElasticity: public pylith::materials::RheologyElasticity {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            IsotropicLinearElasticity(void);

            /// Destructor.
            ~IsotropicLinearElasticity(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Use reference stress and strain in computation of stress and
             * strain?
             *
             * @param[in] value Flag indicating to include reference stress and strain.
             */
            void useReferenceState(const bool value);

            /** Use reference stress and strain in computation of stress and
             * strain?
             *
             * @returns True if using reference stress and strain, false otherwise.
             */
            bool useReferenceState(void) const;

            /** Get auxiliary factory associated with physics.
             *
             * @return Auxiliary factory for physics object.
             */
            pylith::materials::AuxiliaryFactoryElasticity* getAuxiliaryFactory(void);

            /** Add rheology subfields to auxiliary field.
             *
             * @param[inout] auxiliaryField Auxiliary field.
             */
            void addAuxiliarySubfields(void);

            /** Get stress kernel for RHS residual, G(t,s).
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return RHS residual kernel for stress.
             */
            PetscPointFnWrapper getKernelf1v(const spatialdata::geocoords::CoordSys* coordsys) const;

            /** Get elastic constants kernel for RHS Jacobian G(t,s).
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return RHS Jacobian kernel for elastic constants.
             */
            PetscPointJacFnWrapper getKernelJf3vu(const spatialdata::geocoords::CoordSys* coordsys) const;

            /** Get f0 kernel for LHS interface residual, F(t,s,dot{s}), for negative fault face.
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return LHS residual f0 kernel.
             */
            PetscBdPointFnWrapper getKernelf0Neg(const spatialdata::geocoords::CoordSys* coordsys) const;

            /** Get f0 kernel for LHS interface residual, F(t,s,dot{s}), for positive fault face.
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return LHS residual f0 kernel.
             */
            PetscBdPointFnWrapper getKernelf0Pos(const spatialdata::geocoords::CoordSys* coordsys) const;

            /** Get stress kernel for derived field.
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return Project kernel for computing stress subfield in derived field.
             */
            PetscPointFnWrapper getKernelCauchyStressVector(const spatialdata::geocoords::CoordSys* coordsys) const;

        };

        // class IsotropicLinearElasticity

    } // materials
} // pylith

// End of file
