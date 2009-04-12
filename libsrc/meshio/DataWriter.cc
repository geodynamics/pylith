// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include <portinfo>

// ----------------------------------------------------------------------
// Constructor
template<typename mesh_type, typename field_type>
pylith::meshio::DataWriter<mesh_type, field_type>::DataWriter(void) :
  _numTimeSteps(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
template<typename mesh_type, typename field_type>
pylith::meshio::DataWriter<mesh_type, field_type>::~DataWriter(void)
{ // destructor
} // destructor  

// ----------------------------------------------------------------------
// Prepare for writing files.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriter<mesh_type, field_type>::open(const mesh_type& mesh,
					    const int numTimeSteps,
					    const char* label,
					    const int labelId)
{ // open
  _numTimeSteps = numTimeSteps;
} // open

// ----------------------------------------------------------------------
// Close output files.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriter<mesh_type, field_type>::close(void)
{ // close
} // close

// ----------------------------------------------------------------------
// Prepare file for data at a new time step.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriter<mesh_type, field_type>::openTimeStep(const double t,
						    const mesh_type& mesh,
						    const char* label,
						    const int labelId)
{ // openTimeStep
} // openTimeStep

// ----------------------------------------------------------------------
// Cleanup after writing data for a time step.
template<typename mesh_type, typename field_type>
void
pylith::meshio::DataWriter<mesh_type, field_type>::closeTimeStep(void)
{ // closeTimeStep
} // closeTimeStep

// ----------------------------------------------------------------------
// Copy constructor.
template<typename mesh_type, typename field_type>
pylith::meshio::DataWriter<mesh_type, field_type>::DataWriter(const DataWriter& w)
{ // copy constructor
} // copy constructor


// End of file 
