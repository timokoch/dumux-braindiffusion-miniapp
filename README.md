![Simulation](./doc/img/mapped_mri_small.png)

Tracer diffusion in the brain
-------------------------------

We model tracer transport in the brain tissue over 72 hours. The total gadolinium concentration, $c(x, t)$ (amount per brain volume in mmol/l) at position $x$, and time $t \in [0, 259200]$ (in s), satisfies

```math
\begin{align}
    \frac{\partial c}{\partial t} - \mathrm{div}\left( D^{\mathrm{eff}}\nabla c \right) &= - r \phi^{-1} c
\end{align}
```
subject to the boundary and initial conditions

```math
\begin{align}
    -D^{\mathrm{eff}}\nabla c \cdot \mathbf{n} &= k (\phi^{-1}c - \hat{\phi}^{-1} \hat c(x,t)),\\
    c(x, 0) &= 0,
\end{align}
```

where $D^\mathrm{eff}$ is the effective diffusion tensor of gadolinium, $r$ is a local clearance rate due to tracer clearance to blood, $k$ is the brain surface conductivity and $\hat c (x,t)$ is the solute concentration in the cerebrospinal fluid just outside of the pial surface of the brain, and $\phi$, $\hat{\phi}$ are the extra-cellular volume fraction of the brain tissue and the SAS, respectively, occupied by interstitial fluid. We assume that only the extra-cellular space is accessible to the tracer. Dividing the total concentration by $\phi$ (resp. $\hat{\phi}$) computes the concentration per interstitial fluid volume.

We take the values for $D^\mathrm{eff}$ and $\hat{c}(x,t)$ directly from the provided data set. As $\hat{c}(x,t)$ is only provided at five discrete time points,
we provide two strategies to fill data inbetween: linear interpolation and a curve fit procedure where we fit an exponential of the form $f(t, (a, b)(x)) = a(x) t \mathrm{exp}(-b(x)t)$ to the data for every degree of freedom on the boundary. A comparison of the simulation field $c(x,t)$ and the corresponding field provided in the data set at 5 time points in shown in the figure above.

Software requirements and installation
----------------------------------------

* cmake >= 3.18
* C++17 compliant compiler (e.g g++ >= 8 or clang >= 6)
* pkg-config
* MPI (e.g. openmpi) for parallel execution
* Python 3.10 for postprocessing
* Suitesparse (at least CHOLMOD and UMFPack)

Make a new folder (e.g. `dumux`) that will contain all modules.
Inside of this folder, clone this repo

```
mkdir dumux
cd dumux
git clone https://github.com/timokoch/dumux-braindiffusion-miniapp.git
````

The folder structure will look like this

```
dumux
└───dumux-braindiffsion-miniapp
```

Then run from the top folder (`dumux`):

```
./dumux-braindiffsion-miniapp/setup.sh
```

to download Dune/DuMux dependencies and configure and build the project.
The script takes care of this but if you are manually cloning the dependencies
make sure to use that branch.

After that folder structure should look like this:

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

Usage
----------

You can compile the application by running

```
cd dumux-braindiffsion-miniapp/build-cmake/app
make braindiffusion
```

The first time this is executed, it also automatically downloads
the required dataset from Zenodo which can take a while.

You can run it the program in parallel with MPI

```
mpirun -np 4 ./braindiffusion
```

or in serial

```
./braindiffusion
```

Runtime parameters can be configures via the parameter file `params.input`
that is located in the app folder but also linked to the build folder.

For postprocessing, switch to the folder `dumux-braindiffsion-miniapp`.
Create a Python virtual environment (tested Python version 3.10) and install requirements:

```
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

Then switch to `dumux-braindiffsion-miniapp/post` and run the postprocessing
scripts to create result visualizations:

* `cd post`
* `python plot.py` (plots concentration in white and gray matter over time)
* `python map_to_mri.py` (maps the concentration field to the MRI reference image, may take a while)
* `python plot_mapped_mri.py` (plots the concentration field on the MRI reference image)

Troubleshooting
----------------

In case you had an error about a missing dependency, then after installing the dependency
you can reconfigure and rebuild the software by running from the topfolder (`dumux`):

```
./dumux-braindiffsion-miniapp/reconfigure.sh
```

After that continue building the `braindiffusion` app as described above.


Other info
----------

There is an app to test the curve fitting algorithm for the boudary data. Build

```
cd dumux-braindiffsion-miniapp/build-cmake/app
make curvefit
./curvefit
cd ../../post && python plot_curvefit.py
```

to inspect the results. This is thought to be used for debugging purposes only.

If you need an unpublished version of the dataset for testing during development,
you need to supply a Zenodo access token for `make`

```
ZENODO_ACCESS_TOKEN=<token> make braindiffusion
```

Removal
-------

If you don't need the app anymore, simply delete the created
folder `dumux` with all the content.


License
-------

The code is licensed under the GPLv3 (or any later version at your option) license.
See the LICENSES/ file for details.


Publication and reuse
-----------------------

This repository contains code in connection with the Gonzo dataset and the data publication
"MRI Data of CSF-Tracer Evolution Over 72h in Human Brain For Data-Integrated Simulations"
by Jørgen Riseth, Timo Koch, Sofie Lian, Tryggve Holck Storas, Ludmil Zikatanov,
Lars Magnus Valnes, Kaja Nordengen, and Kent-Andre Mardal.

Please refer to this publication for reuse of the data. When using the code in this repository
in a scientific context or otherwise, please cite the above publication and respect
the software license.
