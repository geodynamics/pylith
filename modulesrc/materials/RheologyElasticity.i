// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** @file modulesrc/materials/RheologyElasticity.i
 *
 * Python interface to C++ abstract base class RheologyElasticity.
 */

namespace pylith {
    namespace materials {
        class RheologyElasticity: public pylith::utils::PyreComponent {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            RheologyElasticity(void);

            /// Destructor.
            virtual ~RheologyElasticity(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Get auxiliary factory associated with physics.
             *
             * @return Auxiliary factory for physics object.
             */
            virtual
            pylith::materials::AuxiliaryFactoryElasticity* getAuxiliaryFactory(void) = 0;

            /// Add rheology subfields to auxiliary field.
            virtual
            void addAuxiliarySubfields(void) = 0;

            /** Get stress kernel for RHS residual, G(t,s).
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return RHS residual kernel for stress.
             */
            virtual
            PetscPointFunc getKernelf1v(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

            /** Get elastic constants kernel for RHS Jacobian G(t,s).
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return RHS Jacobian kernel for elastic constants.
             */
            virtual
            PetscPointJac getKernelJf3vu(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

            /** Get f0 kernel for LHS interface residual, F(t,s,dot{s}), for negative fault face.
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return LHS residual f0 kernel.
             */
            virtual
            PetscBdPointFunc getKernelf0Neg(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

            /** Get f0 kernel for LHS interface residual, F(t,s,dot{s}), for positive fault face.
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return LHS residual f0 kernel.
             */
            virtual
            PetscBdPointFunc getKernelf0Pos(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

            /** Get stress kernel for derived field.
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return Project kernel for computing stress subfield in derived field.
             */
            virtual
            PetscPointFunc getKernelCauchyStressVector(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

            /** Add kernels for updating state variables.
             *
             * @param[inout] kernels Array of kernels for updating state variables.
             * @param[in] coordsys Coordinate system.
             */
            virtual
            void addKernelsUpdateStateVars(std::vector < pylith::feassemble::IntegratorDomain::ProjectKernels > * kernels,
                                           const spatialdata::geocoords::CoordSys* coordsys) const;

            /** Update kernel constants.
             *
             * @param[inout] kernelConstants Array of constants used in integration kernels.
             * @param[in] dt Current time step.
             */
            virtual
            void updateKernelConstants(pylith::real_array* kernelConstants,
                                       const PylithReal dt) const;

        };

        // class RheologyElasticity

    } // materials
} // pylith

// End of file
