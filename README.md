# Running The Evaluations:

The code runs properly on the course VMs. Please use the following commit hashes to validate the runtimes listed and remember to `make clean && make` between runs.

Baseline: fc4f9a1d76fed9d5716d5bf6b7418ab829b00add

Optimization 1 - Compiler Flags: c150f8edfb894d81660bdb2e5a0ec4fe61d6324d

Optimization 2 - Localized Energy Calculations: 682c1713e8f35b1071f0aa015e30aafac21d63a2

Optimization 3 - Optimize RNG: cb2270b4a41c0ab855666669ff8aa1d6741ec969

Optimization 4 - Cache Exponentiation: e44c6cdf964b5856f4a11ff6640044444511d595


# 598APE-HW3

This repository contains code for homework 3 of 598APE.

This assignment is relatively simple in comparison to HW1 and HW2 to ensure you have enough time to work on the course project.

In particular, this repository implements a 2D Ising model Monte Carlo simulator (with Metropolis–Hastings algorithm) on an L×L lattice with periodic boundary conditions.

To compile the program run:
```bash
make -j
```

To clean existing build artifacts run:
```bash
make clean
```

This program assumes the following are installed on your machine:
* A working C compiler (gcc is assumed in the Makefile)
* make
* ImageMagick `convert` for PNG output

Usage (after compiling):
```bash
./main.exe <L> <temperature> <steps>
```

In particular, consider speeding up simple run like the following (which runs ~5 seconds on my local laptop under the default setup):
```bash
./main.exe 100 2.269 100000
```

Exact bitwise reproducibility is not required; sanity checks on energy/magnetization bounds must pass. In addition, at the critical temperature (T ≈ 2.269) the energy per spin should approach -1.414 in equilibrium.
