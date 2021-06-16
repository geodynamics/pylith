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

/** @file modulesrc/sources/Source.i
 *
 * Python interface to C++ abstract base class Source.
 */

namespace pylith {
    namespace sources {
        class Source : public pylith::problems::Physics {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /** Default constructor.
             *
             * @param dimension Spatial dimension associated with source.
             */
            Source(const int dimension);

            /// Destructor.
            virtual ~Source(void);

            /// Deallocate PETSc and local data structures.
            virtual
            void deallocate(void);

            /** Set value of label source-id used to identify source cells.
             *
             * @param value Source identifier
             */
            void setSourceId(const int value);

            /** Get value of label source-id used to identify source cells.
             *
             * @returns Source identifier
             */
            int getSourceId(void) const;

            /** Set descriptive label for source.
             *
             * @param value Label of source.
             */
            void setDescriptiveLabel(const char* value);

            /** Get descruptive label of source.
             *
             * @returns Label of source
             */
            const char* getDescriptiveLabel(void) const;

	    /** Create constraint and set kernels.
	     *
	     * @param[in] solution Solution field.
	     * @returns Constraint if applicable, otherwise NULL.
	     */
	    virtual
	    pylith::feassemble::Constraint* createConstraint(const pylith::topology::Field& solution);

            /** Set coordinates and names of points.
             *
             * @param[in] points Array of coordinates [numPoints * spaceDim].
             * @param[in] numPoints Number of points.
             * @param[in] spaceDim Spatial dimension for coordinates.
             * @param[in] pointNames Array with point names.
             * @param[in] numPointNames Number of point names.
             */
            %apply(double* IN_ARRAY2, int DIM1, int DIM2) {
	            (const PylithReal* pointCoords,
	            const int numPoints,
	            const int spaceDim)
	        };
            %apply(const char* const* string_list, const int list_len){
	            (const char* const* pointNames, const int numPointNames)
	        };
            void setPoints(const PylithReal* pointCoords,
                           const int numPoints,
                           const int spaceDim,
                           const char* const* pointNames,
                           const int numPointNames);
            %clear(const PylithReal* pointCoords, const int numPoints, const int spaceDim);
            %clear(const char* const* pointNames, const int numPointNames);
        }; // class Source

    } // sources
} // pylith

// End of file
