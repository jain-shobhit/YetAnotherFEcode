# Topology Optimization with YetAnotherFEcode
This folder contains basic examples to demonstrate the use of [YetAnotherFEcode](https://github.com/jain-shobhit/YetAnotherFEcode) for Topology Optimization problems. The source files are located in the `src/TopologyOptimization` folder.

The optimization problem is solved using the Method of Moving Asymptotes (MMA, [Svanberg, 1987](https://doi.org/10.1002/nme.1620240207)).

## Examples
To run the examples, clone the [YetAnotherFEcode](https://github.com/jain-shobhit/YetAnotherFEcode) repository and run the `startup.m` script to add the necessary paths.

The following examples are provided:

1. `ComplianceMinimization.m`: compliance minimization of the MBB beam structure subject to a volume constraint.

2. `CompliantMechanism.m`: compliant mechanism design subject to a volume constraint.

3. `ElectrostaticOptimization.m`: electrostatic compliance minimization (agave problem) subject to a volume constraint.

4. `FrequencyMaximization.m`: frequency maximization of a clamped-clamped resonator subject to a volume constraint.

## Citation
To cite this package, please cite the followings:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4011281.svg)](https://doi.org/10.5281/zenodo.4011281)

<!-- K. Svanberg. "The method of moving asymptotes - a new method for structural optimization". *Int J Numer Methods Engng* (1987). DOI: [10.1002/nme.1620240207](https://doi.org/10.1002/nme.1620240207). -->

## Contact
If you have any questions, create an issue or email <jacopo.marconi@polimi.it> or <matteo1.pozzi@polimi.it>.
