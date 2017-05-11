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

/** @file modulesrc/feassemble/Constraint.i
 *
 * @brief Python interface to C++ abstract base Constraint.
 */

namespace pylith {
    namespace feassemble {

        class ConstraintPointwise : public pylith::utils::PyreComponent
        { // class ConstraintPointwise

        // PUBLIC METHODS /////////////////////////////////////////////////
public:

        /// Default constructor.
        ConstraintPointwise(void);

        /// Destructor.
        virtual
        ~ConstraintPointwise(void);

        /// Deallocate PETSc and local data structures.
        virtual
        void deallocate(void);

        /** Set name of field in solution to constrain.
         *
         * @param[in] value Name of field in solution to constrain.
         */
        void field(const char* value);

        /** Get name of field in solution to constrain.
         *
         * @returns Name of field in solution to constrain.
         */
        const char* field(void) const;

        /** Set indices of constrained degrees of freedom at each location.
         *
         * Example: [0, 1] to apply forces to x and y degrees of freedom in
         * a Cartesian coordinate system.
         *
         * @param[in] dof Array of indices for constrained degrees of freedom.
         * @param[in] size Size of array
         */
	  %apply(int* IN_ARRAY1, int DIM1) {
	   (const int* flags,
	    const int size)
	     };
	 void constrainedDOF(const int* flags,
			     const int size);
	 %clear(const int* flags, const int size);

        /** Get indices of constrained degrees of freedom.
         *
         * @returns Array of indices for constrained degrees of freedom.
         */
	const pylith::int_array& constrainedDOF(void) const;

        /** Get auxiliary fields.
         *
         * @returns field Field over boundary.
         */
        const pylith::topology::Field& auxFields(void) const;

        /** Check whether constraint has a given auxiliary field.
         *
         * @param[in] name Name of field.
         *
         * @returns True if constraint has auxiliary field, false otherwise.
         */
        bool hasAuxField(const char* name);

        /** Get auxiliary field.
         *
         * @param[out] field Pointer to field.
         * @param[in] name Name of field to retrieve.
         */
        void getAuxField(pylith::topology::Field *field,
                         const char* name) const;

        /** Set spatial database for auxiliary fields.
         *
         * @param[in] value Pointer to database.
         */
        void auxFieldsDB(spatialdata::spatialdb::SpatialDB* value);

        /** Set discretization information for auxiliary subfield.
         *
         * @param[in] name Name of auxiliary subfield.
         * @param[in] basisOrder Polynomial order for basis.
         * @param[in] quadOrder Order of quadrature rule.
         * @param[in] isBasisContinuous True if basis is continuous.
         * @param[in] feSpace Finite-element space.
         */
        void auxFieldDiscretization(const char* name,
                                    const int basisOrder,
                                    const int quadOrder,
                                    const bool isBasisContinuous,
                                    const pylith::topology::FieldBase::SpaceEnum feSpace);

        /** Get discretization information for auxiliary subfield.
         *
         * @param[in] name Name of subfield.
         * @return Discretization information for auxiliary subfield. If
         * discretization information was not set, then use "default".
         */
        const pylith::topology::FieldBase::Discretization& auxFieldDiscretization(const char* name) const;

        /** Set manager of scales used to nondimensionalize problem.
         *
         * @param dim Nondimensionalizer.
         */
        void normalizer(const spatialdata::units::Nondimensional& dim);

        /** Verify configuration is acceptable.
         *
         * @param[in] solution Solution field.
         */
        virtual
        void verifyConfiguration(const pylith::topology::Field& solution) const;

        /** Initialize integrator.
         *
         * @param[in] solution Solution field (layout).
         */
        virtual
        void initialize(const pylith::topology::Field& solution) = 0;

        /** Update auxiliary fields at beginning of time step.
         *
         * @param[in] t Current time.
         * @param[in] dt Current time step.
         */
        virtual
        void prestep(const double t,
                     const double dt);

        /** Set constrained values in solution field.
         *
         * @param[out] solution Solution field.
         * @param[in] t Current time.
         */
        virtual
        void setSolution(pylith::topology::Field* solution,
			 const double t) = 0;

        }; // class ConstraintPointwise

    } // feassemble
} // pylith


// End of file
