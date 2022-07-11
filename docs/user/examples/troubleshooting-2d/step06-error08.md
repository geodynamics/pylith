# Step 6: Error 8

## Error Message

```{code-block} console
---
caption: Output when running Step 6.
linenos: True
---
$ pylith step06_twofaults.cfg

 -- Solving problem.
0 TS dt 0.01 time 0.
    0 SNES Function norm 1.975141284264e-02 
    Linear solve converged due to CONVERGED_ATOL iterations 408
    1 SNES Function norm 7.771687291628e-13 
  Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 1
1 TS dt 0.01 time 0.01
 >> /software/baagaard/py38-venv/pylith-debug/lib/python3.8/site-packages/pylith/problems/Problem.py:201:finalize
 -- timedependent(info)
 -- Finalizing problem.
```

## Troubleshooting Strategy

The simulation ran without errors.
In visualizing the output we notice the slip distribution contains a sharp transition from 0 m to 2.0 m; we intended to prescribe slip that is uniform above y=-20 km and tapers linearly to 0 at y=-30 km.
We load the JSON parameter file into the PyLith Parameter Viewer and find that we are using the default `query_type` of `nearest` for the earthquake rupture parameters.
For our intended piecewise linear variation in slip, we need to use `linear` for the `query_type`.

## Resolution

```{code-block} cfg
---
caption: Correct error in `step06_twofaults.cfg`.
---
[pylithapp.problem.interfaces.fault]
...
db_auxiliary_field.query_type = linear
```
