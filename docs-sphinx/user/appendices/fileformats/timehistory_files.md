(sec:format:TimeHistoryIO)=
# *TimeHistory* Database Files

*TimeHistory* database files contain a header describing the number of points in the time history and the units for the time stamps followed by a list with pairs of time stamps and amplitude values.
The amplitude at an arbitrary point in time is computed via interpolation of the values in the database.
This means that the time history database must span the range of time values of interest.
The points in the time history must also be ordered in time.

```{code-block} cfg
// This time history database specifies temporal variation in
// amplitude. In this case we prescribe a triangular slip time
// history.
//
// Comments can appear almost anywhere in these files and are
// delimited with two slashes (//) just like in C++. All text and
// whitespace after the delimiter on a given line is ignored.
//
// The next line is the magic header for spatial database files
// in ASCII format.
#TIME HISTORY ascii
TimeHistory { // start specifying the database parameters
  num-points = 5 // number of points in time history

  // Specify the units used in the time stamps.
  time-units = year
} // end of TimeHistory header
// The time history values are listed after the parameters.
// Columns time and amplitude where the amplitude values are unitless.
 0.0 0.00
 2.0 1.00
 6.0 4.00
10.0 2.00
11.0 0.00
```
