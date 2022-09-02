A tool for calculation of simplicial homotopy groups.

Input format: square-bracketed list of square-bracked lists
of points, like this:

```
[[1, 2, 3],
 [1, 4], [2, 4]]
```

Output:

```
H_0: Z^0
H_1: Z^1
```

This uses Rust nightly version. For help with setting up Rust,
see https://www.rust-lang.org/tools/install.