# Step 6: Fault Slip Inversion

In this example we do a simple static slip inversion.
We treat the displacements at the fake GPS stations in Step 4 as the "observations" and use the Green's functions from Step 5 to invert for the fault slip that we prescribed in Step 4.

We use simple generalized inversion method with penalties to minimize the seismic moment.
The Python script `invert_slip.py` will load the observations from Step 4 and Green's functions and respones from Step 5 and invert for the slip.

```{code-block} console
---
caption: Run the fault slip inversion code.
---
$ ./invert_slip.py

# Show command line options for the inversion code
$ ./invert_slip.py --help
```

By default, the inversion code will write the results of the inversion to `output/step06_greensfns-inversion_results.txt`.

## Plotting the results

If you are using the PyLith binary, which includes the `matplotlib` Python module, or have it installed, then you can plot the results of the simulation using the `viz/plot_inversion_results.py` Python script.

```{code-block} console
---
caption: Plot the inversion results using the `matplotlib` Python module.
---
$ ./viz/plot_inversion_results.py
```

:::{figure-md} fig:example:strikeslip:2d:step06:solution
<img src="figs/step06-solution.*" alt="Results of slip inversion in Step 6." width="100%"/>

Results of slip inversion in Step 6.
The thick black line shows the prescribed slip in Step 4.
The thin colored lines show the slip from the inversion with different penalty factors.
:::
