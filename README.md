# Deep Parallel Fluids
A repository to explore the benefits of parallelisable deep learning for high-fidelity fluid simulations along with error quantification.

Fluid simulation techniques implemented in this work are heavily based on the 2nd edition of the book [Fluid Simulation for Computer Graphics](https://www.cs.ubc.ca/~rbridson/fluidsimulation/) by Robert Bridson. Deep Learning approaches mostly derive from the works of the [Thuerey Group](https://ge.in.tum.de), notably their book ["Physics-based Deep Learning"](https://physicsbaseddeeplearning.org/intro.html). 

Our work touches and improves several areas:
1. HPC: highly scalable parallel simulations, models, and data
2. CFD: robust solver for various fluid phases and scales
3. DL: fully differentiable end-to-end solver for various deep learning tasks
4. UQ: error quantification on both simulation and the learning tasks

Furthermore, the techniques introduced here are generalisable, and have been applied to improve [this paper](https://www.sciencedirect.com/science/article/pii/S0010465522002466) with features 1 (HPC) and 4 (UQ).
