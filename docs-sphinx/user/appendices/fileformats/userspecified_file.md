(sec:format:TimeStepUser)=
# User-Specified Time-Step File

This file lists the time-step sizes for nonuniform, user-specified time steps associated with the *TimeStepUser* object.
The file's format is an ASCII file that includes the units for the time-step sizes and then a list of the time steps.

```{code-block} cfg
// This time step file specifies five time steps with the units in years.
// Comments can appear almost anywhere in these files and are
// delimited with two slashes (//) just like in C++. All text and
// whitespace after the delimiter on a given line is ignored.
//
// Units for the time steps
units = year
1.0 // Comment
2.0
3.0
2.5
3.0
```
