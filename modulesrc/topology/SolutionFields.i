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

/**
 * @file modulesrc/topology/SolutionFields.i
 *
 * @brief Python interface to C++ SolutionFields object.
 */

%template(SolutionFieldsBase) pylith::topology::Fields<pylith::topology::Field<pylith::topology::Mesh> >;

namespace pylith {
  namespace topology {

    class SolutionFields : public Fields <Field <Mesh> >
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
      const Field<Mesh>& solution(void) const;
      
      /** Get solution field.
       *
       * @returns Solution field.
       */
      Field<Mesh>& solution(void);
      
    }; // SolutionFields

  } // topology
} // pylith

// End of file 
