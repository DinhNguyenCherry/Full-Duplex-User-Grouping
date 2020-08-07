function [vartheta, psi, lambda] = GetParameter_Uplink(G_channel, G_SI, W_current, p_current, W_next, p_next, rho, HD)

% global N_tx
N_tx=4;
if (nargin<8 || HD==0)
    HD = 0;
else
%     global N_rx
N_rx=4;
end
% global K
K=4;
% global rho
% global sigma
sigma = 0.01;

gamma = GetSINR(G_channel, G_SI, W_current, p_current, rho, HD);

vartheta = log(1 + gamma) - gamma;

psi = 2*gamma/p_current(1);



if (HD)
    W_current_mat = reshape(W_current, (N_tx+N_rx), K);
else
    W_current_mat = reshape(W_current, N_tx, K);
end

IN1 = GetNoiseplusInterf(G_channel(:,2:end), G_SI, W_current_mat, p_current(2:end), rho, HD);
IN2 = GetNoiseplusInterf(G_channel, G_SI, W_current_mat, p_current, rho, HD);

Theta = inv(IN1) - inv(IN2);

%lambda = sigma^2*trace(Theta);
lambda = [sigma*sqrt(trace(Theta))];

for j = 1:1:size(G_channel,2)
    %lambda = lambda + p_next(j)^2*G_channel(:,j)'*Theta*G_channel(:,j);
	lambda = [lambda p_next(j)*sqrt(G_channel(:,j)'*Theta*G_channel(:,j))];
end	
	
if (~HD)

    W_next_mat = reshape(W_next, N_tx, K);

    for k = 1:1:K

        %lambda = lambda + rho*W_next_mat(:,k)'*G_SI*Theta*G_SI'*W_next_mat(:,k);
        lambda = [lambda  sqrt(rho)*W_next_mat(:,k)'*G_SI*Theta^(0.5)];

    end
end

lambda = lambda';

end
