# Topology optimization of nonlinear structural dynamics with invariant manifold-based reduced order models

This code provides a MATLAB implementation of the density-based topology optimization algorithm for tailoring the hardening/softening dynamic response of nonlinear mechanical systems.

The coefficient $\gamma$ that controls this behavior is computed analytically using the third-order normal-form parametrization of the Lyapunov subcenter manifold.
The method leverages the adjoint method for efficiently computing sensitivities of the objective function and constraints, while the explicit formulation of nonlinear internal elastic forces through tensor notation simplifies these evaluations.

The optimization problem is solved using the Method of Moving Asymptotes (MMA, [Svanberg, 1987](https://doi.org/10.1002/nme.1620240207)).

## Installation

To run the code, clone the [YetAnotherFEcode](https://github.com/jain-shobhit/YetAnotherFEcode) repository and run the `startup.m` script to add the necessary paths.

## Examples

A detailed demo on the optimization of the $\gamma$ coefficient is provided both as MATLAB script and as live script. For convenience, the live script is also available as a PDF file.

Additionally, the following examples are included:

1. MBB beam: the optimization problem aims at maximizing the first natural frequency $\omega$ while imposing constraints on the total volume and on the gamma coefficient $\gamma$.
2. MEMS resonator: the optimization problem aims at maximizing/minimizing the coefficient $\gamma$ while imposing constraints on the total volume and on the first natural frequency $\omega$.

The folder `ssmtool_validation` contains the code used to validate the results of the optimization. To run the code, the package [SSMTool](https://github.com/jain-shobhit/SSMTool) is required.

## Citation
To cite this package, please cite the following:

M. Pozzi, J. Marconi, S. Jain, M. Li, and F. Braghin. "Topology optimization of nonlinear structural dynamics with invariant manifold-based reduced order models". *Structural and Multidisciplinary Optimization* (2025), DOI: ???.

## Contact
If you have any questions, create an issue or email <jacopo.marconi@polimi.it> or <matteo1.pozzi@polimi.it>.
