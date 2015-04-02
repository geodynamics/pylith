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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file modulesrc/meshio/Xdmf.i
 *
 * @brief Python interface to C++ Xdmf object.
 */

namespace pylith {
  namespace meshio {

    class pylith::meshio::Xdmf
    { // HDF5

      // PUBLIC METHODS -------------------------------------------------
    public :
  
      /// Default constructor.
      Xdmf(void);
      
      /// Destructor
      ~Xdmf(void);
      
      /** Write Xdmf file associated with HDF5 file.
       *
       * @param filenameXdmf Name of Xdmf file.
       * @param filenameHDF5 Name of HDF5 file.
       */
      void write(const char* filenameXdmf,
		 const char* filenameHDF5);
      
      
    }; // Xdmf
    
  } // meshio
} // pylith


// End of file 
