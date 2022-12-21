# NLA final project
Randomized SVD for Hilbert-Schmidt kernels' low-rank approximation

## Concept
As a final project we decided to implement the algorithm for Hilbert-Schmidt kernels' low-rank approximation described [here](https://arxiv.org/pdf/2105.13052.pdf). Algorithm represents generalized Randomized SVD. Generalization is understood in the sense of application to 2-D functions and using Gaussian Process with non-standart covariance kernel for sampling.
## Implementation
### Algorithm describtion
*Implementation inside `HS_RSVD_approx` function* 

**Input**: Kernel `G` for Hilbert-Schmidt operator, covariance kernel `K` for Gaussian Process, number `k` of functions to sample, left bound `a` and right bound `b` of functions' domain.
**Output**: Low-rank approximation `G_k` of integral kernel `G`.
**Steps**:
1. Sample the Gaussian Process `k` times, i.e. apply continous Cholesky factorization $$ K = LL^* $$
2.  

## Contrubutors