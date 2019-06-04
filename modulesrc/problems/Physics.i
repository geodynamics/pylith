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

/** @file modulesrc/problems/Physics.i
 *
 * Python interface to C++ abstract base class Physics.
 */

namespace pylith {
    namespace problems {
        class Physics : public pylith::utils::PyreComponent {
            // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Constructor
            Physics(void);

            /// Destructor
            virtual ~Physics(void);

            /// Deallocate PETSc and local data structures.
            virtual
            void deallocate(void);

            /** Set manager of scales used to nondimensionalize problem.
             *
             * @param dim Nondimensionalizer.
             */
            void setNormalizer(const spatialdata::units::Nondimensional& dim);

            /** Set spatial database for populating auxiliary field.
             *
             * @param[in] value Spatial database with iniital values for auxiliary field.
             */
            void setAuxiliaryFieldDB(spatialdata::spatialdb::SpatialDB* const value);

            /** Set discretization information for auxiliary subfield.
             *
             * @param[in] subfieldName Name of auxiliary subfield.
             * @param[in] basisOrder Polynomial order for basis.
             * @param[in] quadOrder Order of quadrature rule.
             * @param[in] dimension Dimension of points for discretization.
             * @param[in] isBasisContinuous True if basis is continuous.
             * @param[in] feSpace Finite-element space.
             */
            void setAuxiliarySubfieldDiscretization(const char* subfieldName,
                                                    const int basisOrder,
                                                    const int quadOrder,
                                                    const int dimension,
                                                    const bool isBasisContinuous,
                                                    const pylith::topology::FieldBase::SpaceEnum feSpace);

            /** Set discretization information for derived subfield.
             *
             * @param[in] subfieldName Name of auxiliary subfield.
             * @param[in] basisOrder Polynomial order for basis.
             * @param[in] quadOrder Order of quadrature rule.
             * @param[in] dimension Dimension of points for discretization.
             * @param[in] isBasisContinuous True if basis is continuous.
             * @param[in] feSpace Finite-element space.
             */
            void setDerivedSubfieldDiscretization(const char* subfieldName,
                                                  const int basisOrder,
                                                  const int quadOrder,
                                                  const int dimension,
                                                  const bool isBasisContinuous,
                                                  const pylith::topology::FieldBase::SpaceEnum feSpace);

            /** Register observer to receive notifications.
             *
             * Observers are used for output.
             *
             * @param[in] observer Observer to receive notifications.
             */
            void registerObserver(pylith::problems::ObserverPhysics* observer);

            /** Remove observer from receiving notifications.
             *
             * @param[in] observer Observer to remove.
             */
            void removeObserver(pylith::problems::ObserverPhysics* observer);

            /** Get observers receiving notifications of physics updates.
             *
             * @returns Observers receiving notifications.
             */
            pylith::problems::ObserversPhysics* getObservers(void);

            /** Get constants used in kernels (point-wise functions).
             *
             * @param[in] solution Solution field.
             * @param[in] dt Current time step.
             *
             * @return Array of constants.
             */
            const pylith::real_array& getKernelConstants(const PylithReal dt);

            /** Verify configuration is acceptable.
             *
             * @param[in] solution Solution field.
             */
            virtual
            void verifyConfiguration(const pylith::topology::Field& solution) const = 0;

            /** Create integrator and set kernels.
             *
             * @param[in] solution Solution field.
             * @returns Integrator if applicable, otherwise NULL.
             */
            virtual
            pylith::feassemble::Integrator* createIntegrator(const pylith::topology::Field& solution) = 0;

            /** Create constraint and set kernels.
             *
             * @param[in] solution Solution field.
             * @returns Constraint if applicable, otherwise NULL.
             */
            virtual
            pylith::feassemble::Constraint* createConstraint(const pylith::topology::Field& solution) = 0;

            /** Create auxiliary field.
             *
             * @param[in] solution Solution field.
             * @param[in\ physicsMesh Finite-element mesh associated with physics.
             *
             * @returns Auxiliary field if applicable, otherwise NULL.
             */
            virtual
            pylith::topology::Field* createAuxiliaryField(const pylith::topology::Field& solution,
                                                          const pylith::topology::Mesh& physicsMesh) = 0;

            /** Create derived field.
             *
             * @param[in] solution Solution field.
             * @param[in\ physicsMesh Finite-element mesh associated with physics.
             *
             * @returns Derived field if applicable, otherwise NULL.
             */
            virtual
            pylith::topology::Field* createDerivedField(const pylith::topology::Field& solution,
                                                        const pylith::topology::Mesh& physicsMesh) = 0;

            /** Update time-dependent auxiliary field.
             *
             * @param[inout] auxiliaryField Auxiliary field.
             * @param[in] t Current time.
             */
            virtual
            void updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                                      const double t);

            // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////
protected:

            /** Get auxiliary factory associated with physics.
             *
             * @return Auxiliary factory for physics object.
             */
            virtual
            pylith::feassemble::AuxiliaryFactory* _getAuxiliaryFactory(void) = 0;

            /** Get derived factory associated with physics.
             *
             * @return Derived factory for physics object.
             */
            virtual
            pylith::topology::FieldFactory* _getDerivedFactory(void);

            /** Update kernel constants.
             *
             * @param[in] dt Current time step.
             */
            virtual
            void _updateKernelConstants(const PylithReal dt);

        };

        // class Physics

    } // problems
} // pylith

// End of file
