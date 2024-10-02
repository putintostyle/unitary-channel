clear all
close all
global N M R rho sigma J
%% Initial  Setting 
%%% R: data pair \{\sigma_j, \rho_j\}_j=1^R
%%% N: size of U_k\in S_N*N
%%% M: number of unitary matrix
%%% min 1/2| U \rho_j U^*|_F^2

N = 4;
J = 1j;
itnumb = 1000;
rng(123)
%% Defingin \rho
rho = rand(N, N);%+J*rand(N, N);
rho = (rho*rho');
rho = rho./trace(rho);

%% Defining U0 sigma
matrx = rand(N, N)+J*rand(N, N);
U_true = 1/sqrt(2).*kron([1 1;1, -1], eye(2));
sigma = U_true*rho*U_true';
%% Initial Uk 
UHist = {};
Udiff = [];
Udiff_min = [];
matrx = rand(N, N);%J*rand(N, N);
[U0, S, V] = svd(matrx);
TOL = 1e-30;
ResAll = [];
iteP = 0;
while iteP < itnumb

    
    
    hGrad = 2*sigma*U0*rho; % compute the gradient dh/dU_i
    [U_polor, P_polor] = poldec_new(hGrad); % do polar decomp
    UHist{end+1} = U_polor;
    dU = U_polor-U0;
    Udiff = [Udiff, 1/2*norm(dU,'fro')^2];
    U0 = U_polor;

    ResAll = [ResAll, 1/2*norm(sigma-U0*rho*U0','fro')^2];
    % UItes{end+1} = UTrail;
    if ResAll(end)<TOL
        break
    end
    iteP = iteP+1;        
end
%%
figure(1);
loglog(Udiff);
title('Difference Between Iterations', 'FontSize',14);
xlabel('Iteration');
ylabel('$\frac{1}{2}\|U^{(s)}-U^{(s+1)}\|_F^2$', 'FontSize',18, 'Interpreter','latex');
% savefig('example1_diffU');
% saveas(gca, 'example1_diffU', 'eps')
%%
figure(2);
loglog(ResAll);
xlabel('Iteration');
ylabel('Objective function');
title('$\frac{1}{2}\|\sigma-U^{(s)}\rho (U^{(s)})^\dagger\|_F^2$', 'FontSize',18, 'Interpreter','latex');
% savefig('example1_obj');
% saveas(gca, 'example1_obj', 'eps')