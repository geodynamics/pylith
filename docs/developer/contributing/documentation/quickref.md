# MyST Quick reference

## Style guide

1. Use Markedly Structured Text (MyST), not reStructured Text (rST).
2. Place each sentence on its own single line.

## Headings

Use `# Heading 1` only at the top of each page.

(sec-quickref-heading2)=
## Heading 2

### Heading 3

#### Heading 4

Refer to {ref}`sec-quickref-heading2`.

## Admonitions

:::{admonition} General admonition as warning
:class: warning

Text goes here.
:::

:::{attention}
This is an `attention` admonition.
:::

:::{danger}
This is a `danger` admonition.
:::

:::{error}
This is an `error` admonition.
:::

:::{important}
This is an `important` admonition.
:::

:::{note}
This is a `note` admonition.
:::

:::{tip}
This is a `tip` admonition.
:::

:::{warning}
This is a `warning` admonition.
:::

:::{seealso}
This is a `seealso` admonition.
:::

:::{admonition} TODO
:class: error

This is a custom `TODO` admonition.
:::

## Lists

### Itemized lists

* Level 1a
  * Level 2a
    * Level 3a
    * Level 3b
  * Level 2b
* Level 1b
* Level 1c
  
### Definition lists

Term 1
: Definition of term 1

Term 2
: Definition of term 2

### Field lists

:field 1: Description of field 1
:field 2: Description of field 2

## Code blocks

```{code-block} c++
---
caption: C++ code block.
emphasize-lines: 3-4
---
int
main(int argc, char* argv[]) {
    // Emphasized lines corresponding to body of main().
    return 0;
}
```

```{code-block} python
---
caption: Python code block.
---
def square(x):
    return x**2
```

```{code-block} console
---
caption: Interactive shell.
---
$ ls
a b c
```

```{code-block} bash
---
caption: Bash code block.
---
# Comment
for i in "a b c"; do
  print $i
done
```

```{code-block} cfg
---
caption: Config code block.
---
# Comment
[pylithapp]
journal.info.problem = 1

[pylithapp.petsc]
ksp_rtol = 1.0e-3
```

## Tables

Please see {numref}`tab:quickref:example`.
Table labels cannot contain more than one dash, so we use colons to separate the words in the label.
This is likely a bug.

```{table} Table caption
:name: tab:quickref:example
|             Header 1 |   Header 2    | Header 3       |
|---------------------:|:-------------:|:---------------|
|        right aligned |   centered    | left aligned   |
| more data, more data | yet more data | even more data |
```

## Figures

Please see {numref}`fig:quickref:example`.
Figure labels cannot contain more than one dash, so we use colons instead.
This is likely a bug.

:::{figure-md} fig:quickref:example
<img src="../../../_static/images/cig_short_nolabel.*" alt="Screenshot"  width="100px"/>

This is the figure caption. Vector graphics should be provided in **both** PDF (latex output) and SVG (html output) formats.
Raster graphics should be provided in either PNG or JPG formats (whichever is more compact).
:::

## Math

For unlabled equations you can use AMSTex environments, such as `equation`, `gather` and `align`.
Unfortunately, Sphinx does not currently support labeled equations using AMSTex environments.
Labeled equations must use the Markdown syntax as in equation {math:numref}`eqn:quickref:example`.

:::{warning}
A group of equations only has a single label, so if you need to label individual equations, place each one in its own `math` block.
:::

Here is a single equation:
%
\begin{equation}
  F = m a
\end{equation}
%
Here is a group of equations:
\begin{gather}
  f_1(x,y) = a_1 x + b_1 y + c_1 \\
  f_2(x,y) = a_2 x + b_2 y + c_2
\end{gather}
%
Here is a group of aligned equations:
%
\begin{align}
  f_1(x,y) &= a_1 x + b_1 y + c_1 \\
  f_2(x,y) &= a_{xx} x^2 + a_{xy} x y + a_{yy} y^2 + a_x x + a_y y + a_0
\end{align}

Here is a labeled equation (equation {math:numref}`eqn:quickref:example`):
```{math}
:label: eqn:quickref:example
F = m a
```

## Citations

Traditional citation {cite}`Aagaard:etal:2007`.
Citation as noun {cite:t}`Bathe:1995`.

## Table of contents

```
# Table of contents tree with a maximum depth of 1
:::{toctree}
---
maxdepth: 1
---
one.md
two.md
:::
```

```
# Table of contents tree with no maximum depth specified.
:::{toctree}
one.md
two.md
:::
```
