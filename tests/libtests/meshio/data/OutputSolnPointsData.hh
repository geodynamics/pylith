// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

#if !defined(pylith_meshio_outputsolnpointsdata_hh)
#define pylith_meshio_outputsolnpointsdata_hh

#include "pylith/utils/types.hh" // HASA PylithScalar

namespace pylith {
    namespace meshio {
        class OutputSolnPointsData;
    } // pylith
} // meshio

class pylith::meshio::OutputSolnPointsData
{

// PUBLIC METHODS ///////////////////////////////////////////////////////
public:

/// Constructor
OutputSolnPointsData(void);

/// Destructor
~OutputSolnPointsData(void);

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

char* meshFilename;   ///< Filename for input mesh

/// @name Point information.
//@{
int spaceDim;   ///< Number of dimensions in vertex coordinates
int numPoints;   ///< Number of points.
PylithScalar* points;   ///< Coordinates of points [numPoints*spaceDim].
const char** names;   ///< Names of points (e.g., station names).
//@}

/// @name Field data.
//@{
int fiberDim;
PylithScalar* coefs;   ///< Polynomial coefficients [fiberDim*spaceDim] for computing field.
//@}

};

#endif // pylith_meshio_outputsolnpointsdata_hh

// End of file
