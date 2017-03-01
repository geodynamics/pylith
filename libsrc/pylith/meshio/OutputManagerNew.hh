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
 * @file libsrc/meshio/OutputManager.hh
 *
 * @brief Manager for output of finite-element data.
 */

#if !defined(pylith_meshio_outputmanagernew_hh)
#define pylith_meshio_outputmanagernew_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "spatialdata/geocoords/geocoordsfwd.hh" // USES CoordSys

// OutputManager --------------------------------------------------------
/// Manager for output of finite-element data.
class pylith::meshio::OutputManagerNew : public pylith::utils::PyreComponent
{ // OutputManager
friend class TestOutputManagerNew;   // unit testing

// PUBLIC ENUMS /////////////////////////////////////////////////////////
public:

enum TriggerEnum {
    SKIP_TIMESTEPS=0, ///< Skip X time steps between writes.
    ELAPSED_TIME=1, ///< Skip x time between writes.
}; // TriggerEnum

// PUBLIC METHODS ///////////////////////////////////////////////////////
public:

/// Constructor
OutputManagerNew(void);

/// Destructor
virtual
~OutputManagerNew(void);

/// Deallocate PETSc and local data structures.
virtual
void deallocate(void);

/** Set trigger for how often to write output.
 *
 * @param[in] flag Flag indicating which method to use for determining how often to write output.
 */
void trigger(TriggerEnum flag);

/** Get trigger for how often to write otuput.
 *
 * @returns Flag indicating which method to use for determining how often to write output.
 */
TriggerEnum trigger(void) const;

/** Set number of time steps to skip between writes.
 *
 * @param[in] Number of time steps to skip between writes.
 */
void numTimeStepsSkip(const int value);

/** Get number of time steps to skip between writes.
 *
 * @returns Number of time steps to skip between writes.
 */
int numTimeStepsSkip(void) const;

/** Set elapsed time between writes.
 *
 * @param[in] Elapsed time between writes.
 */
void timeSkip(const double value);

/** Get elapsed time between writes.
 *
 * @returns Elapsed time between writes.
 */
double timeSkip(void) const;

/** Set coordinate system in output. The vertex fields in the output
 * are not affected by any change in coordinates.
 *
 * @param[in] cs Coordinate system in output.
 */
void coordsys(const spatialdata::geocoords::CoordSys* cs);

/** Set writer to write data to file.
 *
 * @param[in] datawriter Writer for data.
 */
void writer(DataWriter* const datawriter);

/** Set filter for vertex data.
 *
 * @param[in] filter Filter to apply to vertex data before writing.
 */
void vertexFilter(VertexFilter* const filter);

/** Set filter for cell data.
 *
 * @param[in] filter Filter to apply to cell data before writing.
 */
void cellFilter(CellFilter* const filter);

/** Prepare for output.
 *
 * @param[in] mesh Finite-element mesh object.
 * @param[in] isInfo True if only writing info values.
 * @param[in] label Name of label defining cells to include in output (=0 means use all cells in mesh).
 * @param[in] labelId Value of label defining which cells to include.
 */
virtual
void open(const topology::Mesh& mesh,
          const bool isInfo,
          const char* label =NULL,
          const int labelId =0);

/// Close output files.
virtual
void close(void);

/** Setup file for writing fields at time step.
 *
 * @param[in] t Time of time step.
 * @param[in] mesh Finite-element mesh object.
 * @param[in] label Name of label defining cells to include in output (=0 means use all cells in mesh).
 * @param[in] labelId Value of label defining which cells to include.
 */
virtual
void openTimeStep(const PylithReal t,
                  const topology::Mesh& mesh,
                  const char* label =NULL,
                  const int labelId =0);

/// End writing fields at time step.
virtual
void closeTimeStep(void);

/** Append finite-element vertex field to file.
 *
 * @param[in] t Time associated with field.
 * @param[in] field Vertex field.
 * @param[in] mesh Mesh for output.
 */
virtual
void appendVertexField(const PylithReal t,
                       topology::Field& field,
                       const topology::Mesh& mesh);

/** Append finite-element cell field to file.
 *
 * @param[in] t Time associated with field.
 * @param[in] field Cell field.
 * @param[in] label Name of label defining cells to include in output (=0 means use all cells in mesh).
 * @param[in] labelId Value of label defining which cells to include.
 */
virtual
void appendCellField(const PylithReal t,
                     topology::Field& field,
                     const char* label =NULL,
                     const int labelId =0);

/** Check whether we want to write output at time t.
 *
 * @param[in] t Time of proposed write.
 * @param[in] timeStep Inxex of current time step.
 * @returns True if output should be written at time t, false otherwise.
 */
bool shouldWrite(const PylithReal t,
                 const PylithInt timeStep);

/** Get buffer for field.
 *
 * Find the most appropriate buffer that matches field, reusing and reallocating as necessary.
 *
 * @param[in] fieldIn Input field.
 * @param[in] name Name of subfield (optional).
 * @returns Field to use as buffer for outputting field.
 */
pylith::topology::Field& getBuffer(const pylith::topology::Field& fieldIn,
                                   const char* name =NULL);

// PROTECTED METHODS ////////////////////////////////////////////////////
protected:

/** Dimension field.
 *
 * @param[in] fieldIn Field to dimensionalize.
 */
topology::Field& _dimension(topology::Field& fieldIn);

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected:

topology::Fields* _fields;   ///< Buffer fields.
spatialdata::geocoords::CoordSys* _coordsys; ///< Coordinate system for output.
DataWriter* _writer;   ///< Writer for data.
VertexFilter* _vertexFilter;   ///< Filter applied to vertex data.
CellFilter* _cellFilter;   ///< Filter applied to cell data.

PylithReal _timeSkip; ///< Elapsed time between writes.
PylithReal _timeWrote; ///< Time when data was previously writtern.
PylithInt _numTimeStepsSkip; ///< Number of time steps to skip between writes.
PylithInt _timeStepWrote; ///< Time step when data was previously written.
TriggerEnum _trigger; ///< Flag indicating whether to use elapsed time or number of time steps when deciding when to write.

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

OutputManagerNew(const OutputManagerNew&);   ///< Not implemented.
const OutputManagerNew& operator=(const OutputManagerNew&);   ///< Not implemented

}; // OutputManagerNew

#endif // pylith_meshio_outputmanagernew_hh


// End of file
