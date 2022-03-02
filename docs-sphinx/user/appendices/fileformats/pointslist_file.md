(sec:format:PointsList)=
# *PointsList* File

This file lists the coordinates of the locations where output is requested for the `OutputSolnPoints` component.
The coordinate system is specified in the `OutputSolnPoints` component.

```{code-block} cfg
# Comments are limited to complete lines. The default delimiter for comments
# is '#', which can be changed via parameters. Additionally, the delimiter
# separating values can also be customized (default is whitespace).
#
# The first column is the station name. The coordinates of the points are given
# in the subsequent columns.
P0 1.0 -2.0 0.0
P1 2.0 -4.0 -0.1
P2 0.0 +2.0 0.0
P3 2.5 -0.2 -0.2
P4 0.0 2.0 +0.2
```
