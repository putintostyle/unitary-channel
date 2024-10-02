clear all
close all
global N J
%% Initial  Setting 
%%% R: data pair \{\sigma_j, \rho_j\}_j=1^R
%%% N: size of U_k\in S_N*N
%%% M: number of unitary matrix
%%% min 1/2| U \rho_j U^*|_F^2

N = 10;
J = 1j;
Nstate = 20;
itnumb = 1000;
rng(123)
%% Defingin \rho
rhoList = {};
for r_n = 1:Nstate
    rho = rand(N, N)+J*rand(N, N);
    rho = (rho*rho');
    rho = rho./trace(rho);
    rhoList{end+1} = rho;
end
%% Defining U0 sigma
matrx = rand(N, N)+J*rand(N, N);
[U_true, S, V] = svd(matrx);

sigmaList = {};
for r_n = 1:Nstate
    sigmaList{end+1} = U_true*rhoList{r_n}*U_true';
end
%% Initial Uk 
UHist = {};
Udiff = [];
Udiff_min = [];
matrx = rand(N, N)+J*rand(N, N);
[U0, S, V] = svd(matrx);
TOL = 1e-30;
ResAll = [];
iteP = 0;
while iteP < itnumb

    
    
    hGrad = 0; 
    for i = 1:Nstate
        hGrad = hGrad+2*sigmaList{i}*U0*rhoList{i};
    end
    
    [U_polor, P_polor] = poldec_new(hGrad); % do polar decomp
    UHist{end+1} = U_polor;
    dU = U_polor-U0;
    Udiff = [Udiff, 1/2*norm(dU,'fro')^2];
    U0 = U_polor;
    
    res = 0;
    for i = 1:Nstate
        res = res+1/2*norm(sigmaList{i}-U0*rhoList{i}*U0','fro')^2;
    end
    ResAll = [ResAll, res];
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
savefig('example1b_diffU');
saveas(gca, 'example1b_diffU', 'eps')
%%
figure(2);
loglog(ResAll);
xlabel('Iteration');
ylabel('Objective function');
title('$\frac{1}{2}\sum_{i=1}^{20}\|\sigma_i-U^{(s)}\rho_i (U^{(s)})^\dagger\|_F^2$', 'FontSize',18, 'Interpreter','latex');
savefig('example1b_obj');
saveas(gca, 'example1b_obj', 'eps')