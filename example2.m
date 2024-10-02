clear all
close all
global N M R rho sigma J
%% Initial  Setting 
%%% R: data pair \{\sigma_j, \rho_j\}_j=1^R
%%% N: size of U_k\in S_N*N
%%% M: number of unitary matrix
%%% min 1/2| U \rho_j U^*|_F^2

N = 10;
J = 1j;
itnumb = 10000;
rng(123)
%% Defingin \rho
rho = rand(N, N)+J*rand(N, N);
rho = (rho*rho');
rho = rho./trace(rho);

%% Defining U0 sigma
matrx = rand(N, N)+J*rand(N, N);
[U_true, S, V] = svd(matrx);
sigma = U_true*rho*U_true';
%% Initial Uk 
UHist = {};
Udiff = [];
Udiff_min = [];
ResAll_clean = [];
matrx = rand(N, N)+J*rand(N, N);
[Uinit, S, V] = svd(matrx);
U0 = Uinit;
TOL = 1e-30;
ResAll = [];
iteP = 0;
while iteP < itnumb
    hGrad = 2*sigma*U0*rho; % compute the gradient dh/dU_i
    [U_polor, P_polor] = poldec_new(hGrad); % do polar decomp
    
    
    U0 = U_polor;
    iteP = iteP+1;        
end
U_clean = U0;
%%
UHist = {};
Udiff = [];
Udiff_min = [];

for trial = 1:500
    
    noise = rand(N, N)+J*rand(N, N);
    noise = noise./norm(noise,'fro');
    noise_sigma = (1e-6)*noise+sigma;
    matrx = rand(N, N)+J*rand(N, N);
    U0 = Uinit;
    TOL = 1e-30;
    
    iteP = 0;
    while iteP < itnumb
        hGrad = 2*noise_sigma*U0*rho; % compute the gradient dh/dU_i
        [U_polor, P_polor] = poldec_new(hGrad); % do polar decomp
        
        
        U0 = U_polor;
    
        iteP = iteP+1;        
    end
    Udiff = [Udiff, 1/2*norm(U0-U_clean,'fro')^2];
    % Udiff = [Udiff, 1/2*norm(U_clean*rho*U_clean-U0*rho*U0','fro')^2];
    ResAll = [ResAll, 1/2*norm(noise_sigma-U0*rho*U0','fro')^2];
    ResAll_clean = [ResAll, 1/2*norm(sigma-U0*rho*U0','fro')^2];

end
%%
figure(1);
histogram(ResAll, 10)
xlabel('Optimal values')
ylabel('Counts')
title('Collection of Optimal Objective Values from Each Trial')
figure(2);
histogram(Udiff, 10)
xlabel('$\|\hat{U}-U\|_F^2$', 'FontSize',18, 'Interpreter','latex');
ylabel('Counts')
title('Collection of the Difference Between Noisy and Noise-Free Results:');
