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

/** @file modulesrc/materials/IncompressibleElasticity.i
 *
 * Python interface to C++ IncompressibleElasticity.
 */

namespace pylith {
    namespace materials {
        class IncompressibleElasticity : public pylith::materials::Material {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            IncompressibleElasticity(void);

            /// Destructor.
            ~IncompressibleElasticity(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Include body force?
             *
             * @param[in] value Flag indicating to include body force term.
             */
            void useBodyForce(const bool value);

            /** Include body force?
             *
             * @returns True if including body force term, false otherwise.
             */
            bool useBodyForce(void) const;

            /** Set bulk rheology.
             *
             * @param[in] rheology Bulk rheology for elasticity.
             */
            void setBulkRheology(pylith::materials::RheologyIncompressibleElasticity* const rheology);

            /** Get bulk rheology.
             *
             * @returns Bulk rheology for elasticity.
             */
            pylith::materials::RheologyIncompressibleElasticity* getBulkRheology(void) const;

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

            /** Get derived factory associated with physics.
             *
             * @return Derived factory for physics object.
             */
            pylith::topology::FieldFactory* _getDerivedFactory(void);

        }; // class Elasticity

    } // materials
} // pylith

// End of file
