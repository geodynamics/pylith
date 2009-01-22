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
 * @file modulesrc/utils/EventLogger.i
 *
 * @brief Python interface to C++ EventLogger.
 */


namespace pylith {
  namespace utils {

    class EventLogger
    { // EventLogger

      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :

      /// Constructor
      EventLogger(void);

      /// Destructor
      ~EventLogger(void);

      /** Set name of logging class.
       *
       * @param name Name of logging class.
       */
      void className(const char* name);

      /** Get name of logging class.
       *
       * @returns Name of logging class.
       */
      const char* className(void) const;

      /// Setup logging class.
      void initialize(void);

      /** Register event.
       *
       * @prerequisite Must call initialize() before registerEvent().
       * 
       * @param name Name of event.
       * @returns Event identifier.
       */
      int registerEvent(const char* name);

      /** Get event identifier.
       *
       * @param name Name of event.
       * @returns Event identifier.
       */
      int eventId(const char* name);

      /** Log event begin.
       *
       * @param id Event identifier.
       */
      void eventBegin(const int id);

      /** Log event end.
       *
       * @param id Event identifier.
       */
      void eventEnd(const int id);

      /** Register stage.
       *
       * @prerequisite Must call initialize() before registerStage().
       * 
       * @param name Name of stage.
       * @returns Stage identifier.
       */
      int registerStage(const char* name);

      /** Get stage identifier.
       *
       * @param name Name of stage.
       * @returns Stage identifier.
       */
      int stageId(const char* name);

      /** Log stage begin.
       *
       * @param id Stage identifier.
       */
      void stagePush(const int id);
      
      /// Log stage end.
      void stagePop(void);

    }; // EventLogger

  } // utils
} // pylith


// End of file 
