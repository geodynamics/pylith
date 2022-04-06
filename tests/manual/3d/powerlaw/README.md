# 3D Power Law Manual Tests

These tests check the accuracy of the power-law solution for different time step sizes.
You can either choose to run all simulations (default) or a single simulation (axialtraction_powerlaw, axialtraction_powerlaw_n1, sheartraction_powerlaw).
You can also choose to test the Jacobian; however, this will take much longer to run.
The results are shown in a series of plots.

## Run all tests and check Jacobians

```bash
./run_tests.py --test_jacobian
```

## Run a single test (no Jacobian test)

```bash
./run_tests.py --sim axialtraction_powerlaw
```

## Computing Jacobians

The `jacobian_powerlaw.py` script is used to compute the Jacobian for a power-law material.

```
./jacobian_powerlaw.py
```
