(sec-user-governingeqns-elasticity-rheologies)=
# Bulk Rheologies for Elasticity

In this section we describe the mathematic formulations of the bulk rheologies for elasticity.
The bulk rheologies include

* Isotropic linear elasticity,
* Isotropic linear (Maxwell) viscoelasticity,
* Isotropic linear generalized Maxwell viscoelasticity, and
* Isotropic power-law viscoelasticity.

For the viscoelastic rheologies, we assume the viscous deformation is incompressible; consequently, we separate the stress and strain fields into deviatoric and volumetric components.

```{table} Mathematical notation for viscoelastic formulations.
:name: tab:viscoelasticity:notation

| Variable                   |   Symbol      | Definition                                        |
|:-------------------------------|:-----------------:|:-------------------------------------------------------|
| Mean stress                    |    $\mathit{P}$   | $\frac{\mathop{\mathrm{Tr}}(\boldsymbol{\sigma})}{3}$ |
| Mean strain                    | $\mathit{\theta}$ | $\frac{\mathop{\mathrm{Tr}}(\boldsymbol{\epsilon})}{3}$ |
| Deviatoric stress              | $\boldsymbol{\sigma}^{\mathit{dev}}$  | $\boldsymbol{\sigma} - \mathit{P}\mathit{\boldsymbol{I}}$ |
| Deviatoric strain              | $\boldsymbol{\epsilon}^{\mathit{dev}}$  | $\boldsymbol{\epsilon} - \mathit{\theta}\mathit{\boldsymbol{I}}$ |
| 2nd deviatoric stress invariant|    $\mathit{J}_{2}^{\prime}$    | $\frac{1}{2}\boldsymbol{\sigma}^{dev}:\boldsymbol{\sigma}^{dev}$              |
| 2nd deviatoric strain invariant|    $\mathit{L}_{2}^{\prime}$    | $\frac{1}{2}\boldsymbol{\epsilon}^{dev}:\boldsymbol{\epsilon}^{dev}$              |
```

:::{toctree}
linear-elastic.md
linear-maxwell.md
linear-genmax.md
effective-stress.md
powerlaw.md
alt-formulations
:::
