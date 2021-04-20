---
orphan: true
---
# MyST Quick reference

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
emphasize-lines: 3-3
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

```{code-block} bash
---
caption: Bash code block.
---
# Comment
for i in "a b c"; do
  print $i
done
```
