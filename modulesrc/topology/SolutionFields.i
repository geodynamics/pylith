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
 * @file modulesrc/topology/SolutionFields.i
 *
 * @brief Python interface to C++ SolutionFields object.
 */

namespace pylith {
  namespace topology {

    class SolutionFields : public Fields
    { // SolutionFields

      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :

      /** Default constructor.
       *
       * @param mesh Finite-element mesh.
       */
      SolutionFields(const Mesh& mesh);
      
      /// Destructor.
      ~SolutionFields(void);
      
      /** Set name of solution field.
       *
       * @param name Name of field that is the solution.
       */
      void solutionName(const char* name);
      
      /** Get solution field.
       *
       * @returns Solution field.
       */
      const Field& solution(void) const;
      
      /** Get solution field.
       *
       * @returns Solution field.
       */
      Field& solution(void);
      
    }; // SolutionFields

  } // topology
} // pylith

// End of file 
