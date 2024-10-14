# Backbone curve tailoring via Lyapunov subcenter manifold optimization

This code provides a MATLAB implementation of the optimization algorithm for tailoring the backbone curve of nonlinear systems via Lyapunov subcenter manifold optimization. The backbone is obtained by targeting the desired response at a set of points in the amplitude-frequency domain. The optimization problem is solved by using the *interior-point* algorithm in MATLAB's `fmincon` function. This requires the [MATLAB Optimization Toolbox](https://it.mathworks.com/products/optimization.html).

The code follows the nomenclature and the methodology described in the paper:

Jain, S. & Haller, G. How to compute invariant manifolds and their reduced dynamics in high-dimensional finite element models. *Nonlinear Dyn* (2021). <https://doi.org/10.1007/s11071-021-06957-4>

A more efficient implementation can be pursued using, for instance, the [SSMTool](https://github.com/jain-shobhit/SSMTool/tree/master) package.

## Examples

The following examples are provided:

1. Duffing oscillator.
2. Uniform spring-mass chain.
3. Non-uniform spring-mass chain.
4. von Karman beam.

The last example requires the FE package [YetAnotherFEcode](https://github.com/jain-shobhit/YetAnotherFEcode/tree/master) to define the nonlinear model of the von Karman beam, while the other examples are self-contained.

## Citation
To cite this package, please cite the following:

Pozzi, M., Marconi, J., Jain, S. & Braghin, F. Backbone curve tailoring via Lyapunov subcenter manifold optimization. *Nonlinear Dyn* (2024). <https://doi.org/10.1007/s11071-024-09881-5>

## Contact
If you have any questions, create an issue or email <jacopo.marconi@polimi.it> or <matteo1.pozzi@polimi.it>.
