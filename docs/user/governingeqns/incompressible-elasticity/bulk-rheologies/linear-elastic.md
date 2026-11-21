% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# Linear Isotropic Incompressible Elastic Models

For isotropic, linear incompressible elasticity, the deviatoric elastic constants are
%
\begin{align}
  C^\mathit{dev}_{1111} &= C^\mathit{dev}_{2222} = C^\mathit{dev}_{3333} = +\frac{4}{3} \mu \\
  C^\mathit{dev}_{1122} &= C^\mathit{dev}_{1133} = C^\mathit{dev}_{2233} = -\frac{2}{3} \mu \\
  C^\mathit{dev}_{1212} &= C^\mathit{dev}_{1313} = C^\mathit{dev}_{2323} = \mu
\end{align}

% End of file
