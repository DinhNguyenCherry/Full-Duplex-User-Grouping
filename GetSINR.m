function [SINR] = GetSINR( G_channel, G_SI, W_current, p_current, rho, HD)
%GETSINR Summary of this function goes here
%   Detailed explanation goes here

% This function is to get SINR for uplink
% global N_tx
% global K
N_tx = 4;
K = 4;
if (nargin<6)
    HD = 0;
end

if (HD)
%     global N_rx
    N_rx = 4;
    W_current_mat = reshape(W_current, N_tx+N_rx, K);
else
    W_current_mat = reshape(W_current, N_tx, K);
end

IN = GetNoiseplusInterf(G_channel(:,2:end), G_SI, W_current_mat, p_current(2:end), rho, HD);

SINR = p_current(1)^2*G_channel(:,1)'*inv(IN)*G_channel(:,1);

end