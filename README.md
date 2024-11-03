Diffusion in brain
-------------------

We model tracer transport in the brain tissue over 72 hours. The total gadolinium concentration, $c(x, t)$ (amount per brain volume in mmol/l) at position $x$, and time $t \in [0, 259200]$ (in s), satisfies

```math
\begin{subequations}
\begin{align}
    \frac{\partial c}{\partial t} - \operatorname{div}\left( D^{\text{eff}}\nabla c \right)  & = - r \phi^{-1} c, \\
\intertext{subject to the boundary and initial conditions,}
    -D^{\text{eff}}\nabla c \cdot \boldsymbol{n} &= k (\phi^{-1}c - \hat c(x,t)),\\
    c(x, 0) &= 0,
\end{align}
\end{subequations}
```

where $D^\text{eff}$ is the effective diffusion tensor of gadolinium (in \si{\square\mm\per\s}), $r$ is a local clearance rate (in \si{\per\s}) due to tracer clearance to blood, $k$ is the brain surface conductivity (in \si{\mm\per\s}) and $\hat c (x,t)$ is the solute concentration in the cerebrospinal fluid just outside of the pial surface of the brain, and $\phi$ is the extra-cellular volume fraction of the brain tissue which is occupied by interstitial fluid. We assume that only the extra-cellular space is accessible to the tracer. Dividing the total concentration by $\phi$ computes the concentration per interstitial fluid volume.

We take the values for $D^\text{eff}$ and $\hat{c}(x,t)$ directly from the provided data set. Moreover we set $\phi = 0.2$, $k = \SI{1e-4}{\mm\per\s}$, $r = \SI{1e-4}{\per\s}$. A comparison of the simulation field $c(x,t)$ and the corresponding field provided in the data set at 5 time points in shown in the figure below.

Software requirements and installation
----------------------------------------

* cmake >= 3.18
* C++17 compliant compiler (e.g g++ >= 8 or clang >= 6)
* pkg-config
* MPI (e.g. openmpi) for parallel execution

Make a new folder (e.g. `dumux`) that will contain all modules.
Inside of this folder, clone this repo.
Folder structure would look like this

```
dumux
└───dumux-braindiffsion-miniapp
```

Then run from the top folder (`dumux`):

* `./dumux-braindiffsion-miniapp/setup.sh`

to download Dune/DuMux dependencies and configure and build the project.
The script takes care of this but if you are manually cloning the dependencies
make sure to use that branch.

After that folder structure will look like this:

```
dumux
├───dumux-braindiffsion-miniapp
│    ├───build-cmake/appl
│    ├───appl
│    ├───CMakeLists.txt
│    ....
├───dune-common
├───dune-geometry
├───dune-grid
├───dune-istl
├───dune-localfunctions
└───dumux
```

You can compile the application by running

* `cd dumux-braindiffsion-miniapp/build-cmake/app && make braindiffusion`

and run it with

* `mpirun -np 4 ./braindiffusion`

Runtime parameters can be configures via the parameter file `params.input`
that is located in the app folder but also linked to the build folder.
