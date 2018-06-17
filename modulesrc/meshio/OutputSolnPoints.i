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
// Copyright (c) 2010-2017 University of California, Davis
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

        class pylith::meshio::OutputSolnPoints : public OutputSoln { // OutputSolnPoints

            // PUBLIC METHODS ///////////////////////////////////////////////////
public:

            /** Constructor
             *
             * @param[in] problem Problem to observe.
             */
            OutputSolnPoints(pylith::problems::Problem* const problem);

            /// Destructor
            ~OutputSolnPoints(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Set station coordinates and names.
             *
             * @param[in] points Array of station coordinates [numStations * spaceDim].
             * @param[in] numStations Number of stations.
             * @param[in] spaceDim Spatial dimension for coordinates.
             * @param[in] stationNames Array with station names.
             * @param[in] numStationNames Number of station banes.
             */
            %apply(PylithReal* IN_ARRAY2, int DIM1, int DIM2) {
                (const PylithReal* stationCoords,
                 const PylithInt numStations,
                 const PylithInt spaceDim)
            };
            %apply(const char* const* string_list, const int list_len){
                (const char* const* stationNames, const int numStationNames)
            };
            void stations(const PylithReal* stationCooords,
                          const PylithInt numStations,
                          const PylithInt spaceDim,
                          const char* const* stationNames,
                          const PylithInt numStationNames);
            %clear(const PylithReal* stationCoords, const PylthInt numStations, const PylithInt spaceDim);
            %clear(const char* const* stationNames, const PylithInt numStationNames);

            // PROTECTED MEMBERS ////////////////////////////////////////////////////
protected:

            /** Write solution at time step.
             *
             * @param[in] t Current time.
             * @param[in] tindex Current time step.
             * @param[in] solution Solution at time t.
             */
            void _writeDataStep(const PylithReal t,
                                const PylithInt tindex,
                                const pylith::topology::Field& solution);

        }; // OutputSolnPoints

    } // meshio
} // pylith


// End of file
