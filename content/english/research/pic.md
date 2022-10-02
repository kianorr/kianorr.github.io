---
title: "particle in cell from scratch"
math:
    enable: true
draft: false
---
## What is particle in cell?
Particle in cell (PIC) is a method to simulate the movement of particles under a force of some kind, which in my case is the electric force, essentially simulating a plasma. At it's core, a PIC code uses a mesh grid that has bins/cells. The edges of these bins are described by two $x$ coordinates, $x_j$ and $x_{j+1}$, and within that bin is a particle $i$ with location $r_i$. And from these particle locations, a denisty can be calculated at the edges as

$$
\rho_j = \sum_{i=0}^{N_p}\frac{x_{j+1} - r_i}{\Delta x}, \quad \rho_{j+1}=\sum_{i=0}^{N_p}\frac{r_i - x_{j}}{\Delta x},
$$

where $N_p$ is the number of particles within a bin and $\Delta x$ is the length of the cell. The particle in cell code that I made is the simplest case possible (I say that but it was very difficult and I needed a lot of guidance from my TA), as in it is 1D and electrostatic (time-independent) such that Maxwell's equations are

$$
\nabla \cdot E = \rho, \quad \nabla \cdot B = 0 
\\\\
\nabla \times E = 0, \quad \nabla \times B = 0.
$$

## Two stream instability