# James

<img src="https://github.com/amritagos/james/blob/main/resources/hbond_fig.png?raw=true" width="300" />

A small C++ library to create Network objects (with hydrogen bonds, distance-based cutoffs), and System (similar to the ASE Atoms) objects. All images shown here have been rendered using [`solvis`](https://github.com/amritagos/solvis).

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

## Features

### Generating undirected networks 

Using the [`graph_lib`](https://github.com/seldon-code/graph_lib) library, undirected networks can be used
to represent bonds in an atomic system (represented by the `System` class).

### Finding ion pairs 

<p float="left">
    <img src="https://github.com/amritagos/james/blob/develop/resources/oct_with_hydrogens.png?raw=true" width="600" />
    <img src="https://github.com/amritagos/james/blob/develop/resources/oct_no_hydrogens.png?raw=true" width="600" />
</p>

Ion pairs can be found, provided that a network has been given in conjunction with a `System` object. Note that if you delete atoms in the `System` class or mess up the ordering, you will get errors unless you build the `UndirectedNetwork` object again. The ion pairs are found by searching for all shortest paths between a source and destination. Ion pairs should only contain water molcules in the middle and should have only ions at the end points. The images above show configurations where hydrogens are kept, and ignored, respectively. The procedure to check and find ion pairs would change slightly based on this. 