function [ Power ] = DownlinkPower_Analysis( W, time, HD )
%DOWNLINKPOWER_ANALYSIS Summary of this function goes here
%   Detailed explanation goes here


global K G N_tx
if (nargin<3 || HD==0)
    HD = 0;
else
    global N_rx
end

Power = [];

for g = 1:1:G
    
    if (HD)
        W_mat = reshape(W(:,g), N_tx+N_rx, K);
    else
        W_mat = reshape(W(:,g), N_tx, K);
    end
    
    Power = [Power real(time(g)*diag(W_mat'*W_mat))];
    
end

end

