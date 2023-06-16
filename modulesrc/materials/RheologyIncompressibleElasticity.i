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

/** @file modulesrc/materials/RheologyIncompressibleElasticity.i
 *
 * Python interface to C++ abstract base class RheologyIncompressibleElasticity.
 */

namespace pylith {
    namespace materials {
        class RheologyIncompressibleElasticity : public pylith::utils::PyreComponent {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            RheologyIncompressibleElasticity(void);

            /// Destructor.
            virtual ~RheologyIncompressibleElasticity(void);

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

            /** Get f0p kernel for LHS residual, F(t,s,\dot{s}).
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return RHS residual kernel for pressure.
             */
            virtual
            PetscPointFunc getKernelf0p(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

            /** Get f1u kernel for LHS residual, F(t,s,\dot{s}).
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return LHS residual kernel for stress.
             */
            virtual
            PetscPointFunc getKernelf1u(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

            /** Get Jf0pp kernel for LHS Jacobian F(t,s,\dot{s}).
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return LHS Jf0pp kernel.
             */
            virtual
            PetscPointJac getKernelJf0pp(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

            /** Get Jf3uu kernel for LHS Jacobian F(t,s,\dot{s}).
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return LHS Jacobian kernel for elastic constants.
             */
            virtual
            PetscPointJac getKernelJf3uu(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

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
            void addKernelsUpdateStateVars(std::vector<pylith::feassemble::IntegratorDomain::ProjectKernels>* kernels,
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

        // class RheologyIncompressibleElasticity

    } // materials
} // pylith

// End of file
