# NLA final project
Randomized SVD for Hilbert-Schmidt operators low-rank approximation

## Concept
____
The idea of the project is to implement the algorithm for Hilbert-Schmidt kernels' low-rank approximation found in [2022 ICLR article](https://arxiv.org/pdf/2105.13052.pdf) by N. Boulle and A. Townsend. Algorithm represents generalized Randomized SVD. Generalization is understood in the sense of application to infinite dimensional operators and use of Gaussian Process with non-standart covariance kernel for sampling. We have also considered application of Randomized SVD with nonstandart covariance kernel for sampling in finite case.
## Implementation
___
### Files and their meaning

*Implementation of RSVD for Hilbert-Schmidt kernels is inside `HS_RSVD_approx.m` function* 

*Experimnets with it is inside `main.m` or `main.mlx` or `main.pdf`*

*Implementation of RSVD with nonstandart kernel is inside `RSVD.py`*

### Algorithm describtion
**Input**: Kernel `G` for Hilbert-Schmidt operator, covariance kernel `K` for Gaussian Process, number `k` of functions to sample, left bound `a` and right bound `b` of functions' domain.

**Output**: Low-rank approximation `G_k` of integral kernel `G`.

**Steps**:
1. Sample the Gaussian Process `k` times, i.e. apply continous Cholesky factorization `L = chol(K)` and then get `U = LC`, where `C` is standart-normally distributed vectors.
2. Evaluate the integral operator on sampled functions : `Y = F * U`.
3. Orthonormalize functions from Y: `[Q,~] = qr(Y)`.
4. Compute approximation using adjoint operator: `G_k += q(i) * F_adj * q(i)`.

### Requirements and specifications
* Algorithm and experiments are implemented in **MATLAB**, so to reproduce the results one has to have MATLAB installed.
* To clone this repository use:

`git clone https://github.com/Evgurov/NLA_final_project.git`
*  **Chebfun** package can be installed via typing the following command inside MATLAB command window:

`unzip('https://github.com/chebfun/chebfun/archive/master.zip')
movefile('chebfun-master', 'chebfun'), addpath(fullfile(cd,'chebfun')), savepath`

* Files `main.m` and `main.mlx` are similar, but in different MATLAB formats: 'script' and 'livescript' respectively. If you want to see code from livescript without opening it in MATLAB you can check `main.pdf`.

## Contributors
___
* **Arkadiy Vladimirov** - idea of the project, article analysis, realization of the algorithm for infinite case  
* **Evgeny Gurov** - teamwork organization, article analysis, realization of the algorithm for infinite case, github handling
* **Emil Alkin** - perfoming experiments, presentation
* **Aleksandr Tolmachev** - RSVD with nonstandart kernel for finite case algorithm realization

## Resources
___
* [Original article: A generalization of the Randomized Singular Value Decomposition, Nicolas Boulle Alex Townsend.](https://arxiv.org/pdf/2105.13052.pdf)
* [MATLAB Chebfun mainpage](https://www.chebfun.org/)
* [MATLAB Chebfun repository](https://github.com/chebfun/chebfun)