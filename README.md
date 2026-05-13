# finite-element-method-2d
2D finite element solver for elliptic PDEs in C++

# Finite Element Method for 2D Elliptic PDEs

Scientific computing project developed during the M1 Applied Algebra program at Université Paris-Saclay.

## Overview

This project consists in the implementation of a 2D finite element solver in C++ for elliptic partial differential equations.

The implementation was developed from scratch without external FEM libraries.

## Main Features

- triangular mesh generation
- P1 finite elements
- barycentric basis functions
- numerical quadrature on reference triangles
- local diffusion and reaction matrices
- global matrix assembly
- matrix-vector products
- discrete variational formulation
- L2 and H1 norm computations

## Mathematical Background

The project is based on:
- variational formulations of elliptic PDEs
- finite element discretization
- numerical integration
- affine transformations on triangles
- sparse linear systems

## Project Structure

- `src/` : source code
- `include/` : headers
- `rapport_projet_fem.pdf` : full project report
- `main.cpp` : execution entry point

## Technologies

- C++
- Numerical analysis
- Finite Element Method (FEM)

## Author

Charbel Yannick ECLOU
