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

/** @file modulesrc/bc/DirichletBoundary.i
 *
 * @brief Python interface to C++ DirichletBoundary object.
 */

namespace pylith {
    namespace bc {

        class DirichletNew :
  	public BoundaryConditionNew,
	public pylith::feassemble::ConstraintPointwise
        { // class DirichletNew

        // PUBLIC METHODS /////////////////////////////////////////////////
public:

        /// Default constructor.
        DirichletNew(void);

        /// Destructor.
        ~DirichletNew(void);

        /// Deallocate PETSc and local data structures.
        virtual
        void deallocate(void);

        /** Initialize boundary condition.
         *
         * @param[in] solution Solution field.
         */
        void initialize(const pylith::topology::Field& solution);

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
         */
        virtual
        void _auxFieldsSetup(void) = 0;

        /** Set kernels for RHS residual G(t,s).
         *
         * Potentially, there are g0 and g1 kernels for each equation. If no
         * kernel is needed, then set the kernel function to NULL.
         *
         * @param solution Solution field.
         */
        virtual
        void _setFEKernelsConstraint(const topology::Field& solution) = 0;
        
        }; // class DirichletNew

    } // bc
} // pylith


// End of file
