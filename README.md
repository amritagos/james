# James

<img src="https://github.com/amritagos/james/blob/main/resources/hbond_fig.png?raw=true" width="300" />

A small C++ library to create Network objects (with hydrogen bonds, distance-based cutoffs), and System (similar to the ASE Atoms) objects. Image rendered using [`solvis`](https://github.com/amritagos/solvis).

# Compilation and Running 

If you want to use [`micromamba`](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) as the package manager, create and activate the environment.

```bash
micromamba create -f environment.yml
micromamba activate seldonenv
```

We use `meson` to compile and build `bondfinder`. 

```bash
meson setup build
meson compile -C build
```

## Running Tests

To run the tests, run the following from the top-level directory

```bash
meson test -C build --verbose
```