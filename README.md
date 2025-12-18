# QIC-IDADE
In this repository, the IDADE protocol presented in [Y. H. Teoh, M. Sajjan, Z. Sun, F. Rajabi, and R. Islam,
*Phys. Rev. A 104, 022420 (2021)*](https://doi.org/10.1103/PhysRevA.104.022420) is applied to a simple
1D ion-trap architecture.

# How to execute
- **Perform stability analysis on a conventional trap**:
Run stability_conv.mlx file.
- **Perform stability analysis on a trap + tweezers configuration**:
Run stability_conv.mlx first and then stability_t.mlx. Alternatively, run
stabilityN.mlx .
- **Execute IDADE protocol**:
Run idade.mlx. Note that it needs the following functions from the homonym 
.m files: AB_conv, full_IDA, DE.