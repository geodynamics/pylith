// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

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
            void setPoints(const PylithReal* pointCoords,
                           const PylithInt numPoints,
                           const PylithInt spaceDim,
                           const char* const* pointNames,
                           const PylithInt numPointNames);

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
