(sec-user-governingeqns-elasticity-linear-maxwell)=
# Linear Viscoelastic Models

{numref}`fig:viscoelasticity:models` shows schematic representations of the viscoelastic models. 
The linear Maxwell model can be represented by a spring in series with a linear dashpot.

:::{figure-md} fig:viscoelasticity:models
<img src="figs/viscoelastic-schematic.*" alt="Schematic representations of viscoelastic rheologies" width="100%"/>

Schematic representations of viscoelastic bulk rheologies available in PyLith.
Note that the models with linear (Newtonian) and power-law viscous behavior are treated as separate rheologies.
:::

For a one-dimensional model, the response is given by

```{math}
:label: eqn:maxwell:1d
\frac{\mathit{d\epsilon}_{Total}}{\mathit{dt}} = \frac{\mathit{d\epsilon}_{D}}{\mathit{dt}} + \frac{\mathit{d\epsilon}_{S}}{\mathit{dt}} = \frac{\mathit{\sigma}}{\mathit{\eta}} + \frac{1}{\mathit{E}}\frac{\mathit{d\sigma}}{\mathit{dt}},
```

where $\mathit{\epsilon}_{Total}$ is the total strain, $\mathit{\epsilon}_{D}$ is the strain in the dashpot, $\mathit{\epsilon}_{S}$ is the strain in the spring, $\mathit{\sigma}$ is the stress, $\mathit{\eta}$ is the viscosity of the dashpot, and $\mathit{E}$ is the spring constant.
When a Maxwell material is subjected to constant strain, the stresses relax exponentially with time.
When a Maxwell material is subjected to a constant stress, there is an immediate elastic strain, corresponding to the response of the spring, and a viscous strain that increases linearly with time.
Because the strain response is unbounded, the Maxwell model actually represents a fluid.

Another simple model is the Kelvin-Voigt model, which consists of a spring in parallel with a dashpot. In this case, the one-dimensional response is given by

```{math}
:label: eqn:kelvin:1d
\mathit{\sigma}(\mathit{t}) = \mathit{E\epsilon}(\mathit{t}) + \mathit{\eta}\frac{\mathit{d\epsilon}(\mathit{t})}{\mathit{dt}}.
```

As opposed to the Maxwell model, which represents a fluid, the Kelvin-Voigt model represents a solid undergoing reversible, viscoelastic strain.
If the material is subjected to a constant stress, it deforms at a decreasing rate, gradually approaching the strain that would occur for a purely elastic material. When the stress is released, the material gradually relaxes back to its undeformed state.
