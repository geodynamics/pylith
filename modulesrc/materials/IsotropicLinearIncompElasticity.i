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

/** @file modulesrc/materials/IsotropicLinearIncompElasticity.i
 *
 * Python interface to C++ IsotropicLinearIncompElasticity.
 */

namespace pylith {
    namespace materials {
        class IsotropicLinearIncompElasticity : public pylith::materials::RheologyIncompressibleElasticity {
            // PUBLIC METHODS /////////////////////////////////////////////////////
public:

            /// Default constructor.
            IsotropicLinearIncompElasticity(void);

            /// Destructor.
            ~IsotropicLinearIncompElasticity(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Include reference stress/strain?
             *
             * @param value Flag indicating to include reference stress/strain.
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

            /** Get f0p kernel for LHS residual, F(t,s,\dot{s}).
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return RHS residual kernel for pressure.
             */
            PetscPointFunc getKernelf0p(const spatialdata::geocoords::CoordSys* coordsys) const;

            /** Get f1u kernel for LHS residual, F(t,s,\dot{s}).
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return RHS residual kernel for stress.
             */
            PetscPointFunc getKernelf1u(const spatialdata::geocoords::CoordSys* coordsys) const;

            /** Get Jf0pp kernel for LHS Jacobian, F(t,s,\dot{s}).
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return LHS Jf0pp kernel.
             */
            PetscPointJac getKernelJf0pp(const spatialdata::geocoords::CoordSys* coordsys) const;

            /** Get Jf3uu kernel for LHS Jacobian F(t,s,\dot{s}).
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return LHS Jf3uu kernel.
             */
            PetscPointJac getKernelJf3uu(const spatialdata::geocoords::CoordSys* coordsys) const;

            /** Get stress kernel for derived field.
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return Project kernel for computing stress subfield in derived field.
             */
            PetscPointFunc getKernelCauchyStressVector(const spatialdata::geocoords::CoordSys* coordsys) const;

        };

        // class IsotropicLinearIncompElasticity

    } // materials
} // pylith

// End of file
