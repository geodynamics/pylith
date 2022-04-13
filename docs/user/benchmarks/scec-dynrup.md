# SCEC Dynamic Rupture Benchmarks

The SCEC website <https://strike.scec.org/cvws/cgi-bin/cvws.cgi> includes a graphical user interface for examining the benchmark results.
Benchmark results for PyLith are available for TPV205-2D (horizontal slice through a vertical strike-slip fault), TPV205 (vertical strike-slip fault with high and low stress asperities), TPV210-2D (vertical slice through a 60-degree dipping normal fault), TPV210 (60-degree dipping normal fault), TPV11, TPV12, TPV13, TPV14-2D and TPV15-2D (horizontal slice through a verticel strike-slip fault with a branch), TPV14, TPV15, TPV 24, TPV25 (vertical strike-slip fault with a branch), TPV 16 and 17 (vertical strike-slip fault with spatially heterogeneous initial tractions), TPV 22 and 23 (vertical strike-slip fault with a stepover), TPV102 (vertical strike-slip fault with rate-state friction).

The benchmark results indicate that triangular and tetrahedral cells generate less numerical noise than quadrilateral or hexahedral cells.
The input files in the repository are updated for PyLith v2.0.0, so you will need to modify them if you use another version of PyLith.
