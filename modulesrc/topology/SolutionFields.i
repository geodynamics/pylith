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
      
      /** Set name of field that will be used in the solve.
       *
       * @param name Name of field used in the solve.
       */
      void solveSolnName(const char* name);
      
      /** Get field used in the solve.
       *
       * @returns Field used in the solve.
       */
      const Field<Mesh>& solveSoln(void) const;
      
      /** Get field used in the solve.
       *
       * @returns Field used in the solve.
       */
      Field<Mesh>& solveSoln(void);

      /** Create history manager for a subset of the managed fields.
       *
       * @param fields Fields in history (first is most recent).
       * @param size Number of fields in history.
       */
      %apply(const char* const* string_list, const int list_len){
	(const char* const* fields,
	 const int size)
	  };
      void createHistory(const char* const* fields,
			 const int size);
      %clear(const char* const* fields, const int size);
      
      /** Shift fields in history. Handles to fields are shifted so that
       *  the most recent values become associated with the second most
       *  recent item in the history, etc.
       */
      void shiftHistory(void);
      
    }; // SolutionFields

  } // topology
} // pylith

// End of file 
