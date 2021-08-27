---
orphan: true
---
# MyST Quick reference

(sec-quickref)=
# Heading 1
## Heading 2
### Heading 3
#### Heading 4

Refer to Section {ref}`sec-quickref`.

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

* Level 1
  * Level 2
    * Level 3
  
### Definition lists

Term 1
: Definition of term 1

Term 2
: Definition of term 2

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

Please see {numref}`tab-quickref`.

```{table} Table caption
:name: tab-quickref
|             Header 1 |   Header 2    | Header 3       |
| -------------------: | :-----------: | :------------- |
|        right aligned |   centered    | left aligned   |
| more data, more data | yet more data | even more data |
```

## Figures

Please see {numref}`fig-quickref`.


:::{figure-md} fig-quickref
<img src="_static/images/cig_short_nolabel.*" alt="Screenshot"  width="100px"/>

This is the figure caption.
:::

## Citations

Traditional citation {cite}`Aagaard:etal:2007`.
Citation as noun {cite:t}`Bathe:1995`.
