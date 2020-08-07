function [ GetBreak, OptimalValue, OptimalValue_preStep, DownlinkRate_PerGroupPerUser, UplinkRate_PerGroupPerUser, RDown_next, RThDown_next, RUp_next, RTh_next, W_next, p_next, phi_next, time_next, w_tilde_next, p_bar_next, alpha_bar_next, alpha, beta_bar_next, beta ] = GetInitialization( Pbs, P, H, G_channel, G_SI, G_hat, W_current, p_current, phi_current, rho, Fixed_timegroup_assignment )
%GETINITIALIZATION Summary of this function goes here
%   Detailed explanation goes here

% global N_tx
% global G
% global K
% global L
N_tx  =4;
G = 3;
K = 4;
L = 4;

OptimalValue_preStep = 0;
OptimalValue = 0;

alpha_current = 1/2*ones(K, G);%1/2*ones(K, G);
beta_current = 1/2*ones(L, G);%1/2*ones(L, G);
if (Fixed_timegroup_assignment)
    time_current = 1/G*ones(1,G);
else
    time_current = rand(1,G);
    time_current = time_current/sum(time_current);
end

RDown_current = 0;
RThDown_current = zeros(K,G);
alpha_bar_current = zeros(K,G);
mu_current = zeros(K,G);

for g = 1:1:G
    for k = 1:1:K
        PerGroupPerUserRate = log2(1 + real(H(:,k)'*W_current((k-1)*N_tx+1:k*N_tx,g))^2/phi_current(k,g)^2);
        alpha_bar_current(k,g) = sqrt(alpha_current(k,g)*PerGroupPerUserRate);
%         mu_current(k,g) = norm(W_current((k-1)*N_tx+1:k*N_tx,g))^2/alpha_current(k,g);
        mu_current(k,g) = Pbs/K;
        RDown_current = RDown_current + time_current(g)*alpha_bar_current(k,g)^2;
%         RDown_current = RDown_current + time_current(g)*PerGroupPerUserRate;
        RThDown_current(k,g) = sqrt(time_current(g)*alpha_bar_current(k,g)^2);
    end
end

RUp_current = 0;
RTh_current = zeros(L,G);
beta_bar_current = zeros(L,G);

for g = 1:1:G
    for l = 1:1:L
        PerGroupPerUserRate = log2(1 + GetSINR(G_channel(:,l:L), G_SI, W_current(:,g), p_current(l:L,g), rho));
        beta_bar_current(l,g) = sqrt(beta_current(l,g)*PerGroupPerUserRate);
        RUp_current = RUp_current + time_current(g)*beta_bar_current(l,g)^2;
        RTh_current(l,g) = sqrt(time_current(g)*beta_bar_current(l,g)^2);
    end
end

w_tilde_current  = sum(W_current.^2);
p_bar_current = p_current.^2;%repmat(P, 1, G)/G*inv(diag(time_current));

n=0;
GetBreak = 0;

while (~(OptimalValue>0))
    
    if (n>0)
        OptimalValue_preStep = OptimalValue;
    end

    [ OptimalValue, DownlinkRate_PerGroupPerUser, UplinkRate_PerGroupPerUser, RDown_current, RThDown_current, RUp_current, RTh_current, W_current, p_current, phi_current, time_current, w_tilde_current, p_bar_current, alpha_bar_current, alpha, beta_bar_current, beta] = Get_optSolutionPerIteration4(Pbs, P, H, G_channel, G_SI, G_hat, W_current, p_current, phi_current, time_current, w_tilde_current, p_bar_current, RDown_current, RThDown_current, RUp_current, RTh_current, alpha_bar_current, beta_bar_current, rho, Fixed_timegroup_assignment);

    n = n+1
    
    if (n>50)
        GetBreak = 1;
        break;
    end
    
end

disp(['Iterations: ' num2str(n) ]);
OptimalValue_preStep
OptimalValue 
DownlinkRate_PerGroupPerUser
UplinkRate_PerGroupPerUser
RDown_next = RDown_current;
RThDown_next = RThDown_current;
RUp_next = RUp_current;
RTh_next = RTh_current;
W_next = W_current;
p_next = p_current;
phi_next = phi_current;
time_next = time_current;
w_tilde_next = w_tilde_current;
p_bar_next = p_bar_current;
alpha_bar_next = alpha_bar_current;
alpha
beta_bar_next = beta_bar_current;
beta

end

