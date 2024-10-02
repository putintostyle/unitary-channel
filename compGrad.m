function gradf = compGrad(U, Sigma, Rho)
global N M R


gradf = 2*Sigma*U*Rho;