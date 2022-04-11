(sec-user-governing-eqns-alternative-formulations)=
# Alternative Material Model Formulations

In some cases there are alternative formulations that may be used in future versions of PyLith, and we describe those here.

## Effective Stress Formulation for a Linear Maxwell Viscoelastic Material

A linear Maxwell viscoelastic material may be characterized by the same elastic parameters as an isotropic elastic material ($E$ and $\nu$), as well as the viscosity, $\eta$.
The creep strain increment is
```{math}
:label: eq:D1
\begin{gathered}
\underline{\Delta e}^{C}=\frac{\Delta t\phantom{}^{\tau}\underline{S}}{2\eta}\,\,.
\end{gathered}
```
Therefore,
```{math}
:label: eq:D2
\begin{gathered}
\Delta\overline{e}^{C}=\frac{\Delta t\sqrt{^{\tau}J_{2}^{\prime}}}{\sqrt{3\eta}}=\frac{\Delta t\phantom{}^{\tau}\overline{\sigma}}{3\eta}\,,\,\mathrm{and}\,^{\tau}\gamma=\frac{1}{2\eta}\,\,.
\end{gathered}
```
Substituting equations {math:numref}`eq:46`, {math:numref}`eq:D1`, and {math:numref}`eq:D2` into {math:numref}`eq:43`, we obtain
```{math}
:label: eq:D3
$$\begin{gathered}
^{t+\Delta t}\underline{S}=\frac{1}{a_{E}}\left\{ ^{t+\Delta t}\underline{e}^{\prime}-\frac{\Delta t}{2\eta}\left[(1-\alpha)^{t}\underline{S}+\alpha\phantom{}^{t+\Delta t}\underline{S}\right]\right\} +\underline{S}^{I}\,\,.
\end{gathered}
```
Solving for $^{t+\Delta t}\underline{S}$,
```{math}
:label: eq:D4
\begin{gathered}
^{t+\Delta t}\underline{S}=\frac{1}{a_{E}+\frac{\alpha\Delta t}{2\eta}}\left[^{t+\Delta t}\underline{e}^{\prime}-\frac{\Delta t}{2\eta}(1-\alpha)^{t}\underline{S}+\frac{1+\mathrm{\nu}}{E}\underline{S}^{I}\right]\,\,.
\end{gathered}
```
In this case it is possible to solve directly for the deviatoric stresses, and the effective stress function approach is not needed.
To obtain the total stress, we simply use
```{math}
:label: eq:D5
\begin{gathered}
^{t+\Delta t}\sigma_{ij}=\phantom{}^{t+\Delta t}S_{ij}+\frac{\mathit{1}}{a_{m}}\left(\,^{t+\Delta t}\theta-\theta^{I}\right)\delta_{ij}+P^{I}\delta_{ij}\,\,.\label{eq:D5}
\end{gathered}
```
To compute the viscoelastic tangent material matrix relating stress and strain, we need to compute the first term in {math:numref}`eq:58`.
From {math:numref}`eq:D4`, we have
```{math}
:label: eq:D12
\begin{gathered}
\frac{\partial\phantom{}^{t+\Delta t}S_{i}}{\partial\phantom{}^{t+\Delta t}e_{k}^{\prime}}=\frac{\delta_{ik}}{a_{E}+\frac{\alpha\Delta t}{2\eta}}\,\,.\label{eq:D12}
\end{gathered}
```
Using this, along with Equations {math:numref}`eq:58`, {math:numref}`eq:59`, and {math:numref}`eq:60`, the final material matrix relating stress and tensor strain is
```{math}
:label: eq:D13
\begin{gathered}
C_{ij}^{VE}=\frac{1}{3a_{m}}\left[\begin{array}{cccccc}
1 & 1 & 1 & 0 & 0 & 0\\
1 & 1 & 1 & 0 & 0 & 0\\
1 & 1 & 1 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0
\end{array}\right]+\frac{1}{3\left(a_{E}+\frac{\alpha\Delta t}{2\eta}\right)}\left[\begin{array}{cccccc}
2 & -1 & -1 & 0 & 0 & 0\\
-1 & 2 & -1 & 0 & 0 & 0\\
-1 & -1 & 2 & 0 & 0 & 0\\
0 & 0 & 0 & 3 & 0 & 0\\
0 & 0 & 0 & 0 & 3 & 0\\
0 & 0 & 0 & 0 & 0 & 3
\end{array}\right]\,.
\end{gathered}
```
Note that the coefficient of the second matrix approaches $E/3(1+\nu)=1/3a_{E}$ as $\eta$ goes to infinity.
To check the results we make sure that the regular elastic constitutive matrix is obtained for selected terms in the case where $\eta$ goes to infinity.
%
```{math}
:label: eq:D14
\begin{gathered}
C_{11}^{E}=\frac{E(1-\nu)}{(1+\nu)(1-2\nu)}\,\,\nonumber \\
C_{12}^{E}=\frac{E\nu}{(1+\nu)(1-2\nu)}\,.\\
C_{44}^{E}=\frac{E}{1+\nu}\,\,\nonumber
\end{gathered}
```
%
This is consistent with the regular elasticity matrix, and {math:numref}`eq:D13` should thus be used when forming the stiffness matrix.
We do not presently use this formulation, but it may be included in future versions.
