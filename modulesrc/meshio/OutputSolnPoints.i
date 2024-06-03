// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information. 
// =================================================================================================

/**
 * @file modulesrc/meshio/OutputSolnPoints.i
 *
 * @brief Python interface to C++ OutputSolnPoints object.
 */

namespace pylith {
    namespace meshio {
        class OutputSolnPoints : public pylith::meshio::OutputSoln {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Constructor.
            OutputSolnPoints(void);

            /// Destructor
            virtual ~OutputSolnPoints(void);

            /** Set coordinates and names of points.
             *
             * @param[in] points Array of coordinates [numPoints * spaceDim].
             * @param[in] numPoints Number of points.
             * @param[in] spaceDim Spatial dimension for coordinates.
             * @param[in] pointNames Array with point names.
             * @param[in] numPointNames Number of point banes.
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

            // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////
protected:

            /** Write solution at time step.
             *
             * @param[in] t Current time.
             * @param[in] tindex Current time step.
             * @param[in] solution Solution at time t.
             */
            void _writeSolnStep(const PylithReal t,
                                const PylithInt tindex,
                                const pylith::topology::Field& solution);

        }; // OutputSolnPoints

    } // meshio
} // pylith

// End of file
