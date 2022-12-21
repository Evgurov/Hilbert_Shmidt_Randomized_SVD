function [G_k, rel_errors] = HS_RSVD_approx(G, K, k, a, b)
% Function takes 5 parameters strictly:
% G : MATLAB function for Hibert Schmidt operator kernel
% K : MATLAB function for covariance kernel of Gaussian Process to sample
% from
% k : number of functions to be sampled
% a, b : left and right bounds for functions' domain
%   1.Generate chebfun operators for HS operator and its adjoint
%   using given operator kernels.
%   2.Sample k functions from GP with given kernel K using contionous 
%   Cholesky factorization.
%   3.Evaluate operator on sampled fucntions.
%   4.Make contionous QR factorization for obtained evaluation.
%   5.Construct QQ^*A representation

    K_ker = chebfun2(K);
    G_ker = chebfun2(G);

    F = chebop(a, b);
    F.op = @(x,f) fred(G_ker, f);
    
    G_adj = @(x,y) G_ker(y,x);
    G_adj_ker = chebfun2(G_adj);
    
    F_adg = chebop(a,b);
    F_adg.op = @(x,f) fred(G_adj_ker, f);
    
    L = chol(K_ker);
    n = length(L(:,1));
    
    C = randn(n,k);
    U = L' * C;
    Y = F * U;
    [Q,~] = qr(Y');

    rel_errors = ones(1, k);
    norm_G = norm(G_ker);

    G_k = chebfun2(@(x,y) 0 );
    for i = 1 : k
        q_i = Q(:, i);
        q_i_adj = (F_adg * q_i)'; 
        G_k = G_k + chebfun2(@(x,y) q_i(x) .* q_i_adj(y));
        rel_errors(i) = norm((G_ker - G_k)) / norm_G;
    end  
end

