# Acoustic2d

Solver for a 2-dimensional acoustic system with reflective, MUR and PML absorbing boundary conditions.

## Getting started

### Prerequisites

```
CMake (ver. 2.8 or higher)
ParaView (ver. 5.6.0) (visualization)
```

### Building and running

First, create a local folder (e.g. `build`), where the results and techical data will be stored:

```
mkdir build
cd build
```

Then, build, using cmake and make, and run:

```
cmake ../ && make && ./solver
```

### Output and visualization
The output files are stored in the folder you created (e.g. `build`) in the subdirectory `data`.
To visualize, open `build/data/test.csv` in ParaView. Add FilterToPoints filter, select `x coord` as X Column, `y coord` as y Column, `z coord` as Z Column. Choose `scalar` Coloring. On top add Delaunay2D filter with `scalar` coloring.