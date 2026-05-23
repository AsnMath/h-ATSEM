# Prerequisites

All code has been tested on **Ubuntu 24.04** with **g++ 13.3.0**.

Other systems and compilers may also work, but **C++ 20** support is required.

The code is based on the Eigen (https://libeigen.gitlab.io/) and autodiff (https://autodiff.github.io/), which are included in ``3rdparty``.

# Build

Simply use ``make`` in the ``test/Poison`` and ``test/KohnSham``, then the executable file will be generated as ``./build/apps/main.exe``.

Otherwise, directly use ``g++`` command or other compiler. Since the library is header-only, the only source file for each test is ``main.cpp``, there is no extra operator needed.

# Run

Use ``./build/apps/main.exe`` to run the program. The test sample for Poisson equation is
$$
u = \exp(-10 \Vert x \Vert), x \in [-1, 1]^{3}
$$
and the test sample for Kohn-Sham equation is Helium atom.

The program will output informations such as the mesh size, error, number of iterations, where ``sample-output.txt`` gives the results by my computer (Intel(R) Core(TM) i7-8700 CPU with 32 GB memory).

