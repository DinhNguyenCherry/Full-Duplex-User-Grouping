% function [W, p, phi, z, z_bar, upsilon, upsilon_bar, time, x, y, OptimalValue] = Get_optSolutionPerIteration3( Pbs, P, H, G_channel, G_SI, G_hat, W_current, p_current, phi_current, z_current, z_bar_current, upsilon_current, upsilon_bar_current, alpha_current, beta_current, time_current, x_current, y_current, mu_current, theta_current)
function [OptimalValue, DownlinkRate_PerUser, UplinkRate_PerUser, RDown, RThDown, RUp, RTh, W, p, phi, time, w_tilde, p_bar, alpha_bar, alpha, beta_bar, beta] = Get_optSolutionPerIteration4( Pbs, P, H, G_channel, G_SI, G_hat, W_current, p_current, phi_current, time_current, w_tilde_current, p_bar_current, RDown_current, RThDown_current, RUp_current, RTh_current, alpha_bar_current, beta_bar_current, rho, Fixed_timegroup_assignment )
%GET_OPTSOLUTIONPERITERATION3 Summary of this function goes here
%   Detailed explanation goes here

% This function is to get the optimal solution per iteration

% global N_tx
% global G
% global K
% global L
% global sigma_K
% global Rate_Threshold
N_tx = 4;
G = 3;
K  = 4;
L = 4;
sigma_K = 0.01*ones(1, K);
Rate_Threshold = 1;

BarrierValue = 100;

if (Fixed_timegroup_assignment)
    
    alpha_next = ones(K,G); %
    beta_next = ones(L,G); %
    time_next = 1/G*ones(1,G);
    
else
    
    alpha_next = sdpvar(K,G,'full','real'); 
    beta_next = sdpvar(L,G,'full','real'); 
    time_next = sdpvar(1,G,'full','real');

end
% tau_next = ones(1,G);%sdpvar(1,G,'full','real');
% mu_next = sdpvar(1,G,'full','real');
% theta_next = sdpvar(L,G,'full','real');
    

W_next = sdpvar(N_tx*K,G,'full','complex');
p_next = sdpvar(L, G,'full','real');
% p_next = sdpvar(L, G);
phi_next = sdpvar(K, G,'full','real');
% eta_next = sdpvar(1, G);
% z_next = sdpvar(K, G,'full','real');
% kapa_next = sdpvar(1, G);
% upsilon_next = sdpvar(L, G,'full','real');


obj = 0;
cons = [];

% add objective function and constraints for downlink rate in Eq. (33)

% z_bar_next = sdpvar(K, G);  % 1/z_next <= z_bar_next
% W_bar_next = sdpvar(K,G); % sqrt( real(h'*W_next) ) <= W_bar_next
z_tilde_next = sdpvar(K, G); % (phi_next^2 + real(h'*W_next)^2)/z_next <= z_tilde_next

% cons = [cons, z_next>=0];
% cons = [cons, z_bar_next>=0];

GammaDown_next = sdpvar(1,G); % = sum{1->K} of log2(1+gamma(k,g)(w,p))
RDown_next = sdpvar(1,G); % = time(g) * GammaDown(g)

RThDown_next = sdpvar(K,G);
RThDown_bar_next = sdpvar(K,G);

alpha_bar_next = sdpvar(K,G); % alpha_next*log2()>= alpha_bar_next^2
alpha_tilde_next = sdpvar(K,G); % alpha_bar_next^2 >= alpha_tilde_next
PerGroupDownlinkRate = sdpvar(K,G);

for g = 1:1:G
    
%     PerGroupDownlinkRate = 0;
    
    for k = 1:1:K

        
		
		% get parameters in Eq. (16)
		[varphi, chi, varpi] = GetParameter_Downlink(H(:,k), W_current(((k-1)*N_tx+1):(k*N_tx),g), phi_current(k,g));
		
		% add first term in Eq. (33)
		
% 		obj = obj + real(varphi);
%         PerGroupDownlinkRate = PerGroupDownlinkRate + time_next(g)*real(varphi);
        PerGroupDownlinkRate(k,g) = real(varphi);
%         Downlink_Rate = Downlink_Rate + real(varphi);
% 		
% 		cons = [cons, cone([1, 0.5*(z_next(k,g)-z_bar_next(k,g))], 0.5*(z_next(k,g)+z_bar_next(k,g))) ];
%         ratio_zz = sqrt(z_current(k,g)/z_bar_current(k,g));
%         cons = [cons, cone([sqrt(1/2)*z_next(k,g)/ratio_zz, sqrt(1/2)*z_bar_next(k,g)*ratio_zz], 1)];
		
		
		% add second term in Eq. (33)
		
% 		obj = obj + real(chi)*real(H(:,k)'*W_next(((k-1)*N_tx+1):(k*N_tx),g));
%         PerGroupDownlinkRate = PerGroupDownlinkRate + time_next(g)*real(chi)*real(H(:,k)'*W_next(((k-1)*N_tx+1):(k*N_tx),g));
        PerGroupDownlinkRate(k,g) = PerGroupDownlinkRate(k,g) + real(chi)*real(H(:,k)'*W_next(((k-1)*N_tx+1):(k*N_tx),g));
%         Downlink_Rate = Downlink_Rate + real(chi)*real(H(:,k)'*W_next(((k-1)*N_tx+1):(k*N_tx),g));
		
% 		cons = [cons, real(H(:,k)'*W_next(((k-1)*N_tx+1):(k*N_tx),g)) <= (W_bar_next(k,g))^2];
% 		cons = [cons, cone(sqrt(real(H(:,k)'*W_next(((k-1)*N_tx+1):(k*N_tx),g))), W_bar_next(k,g) )];
        
%         cons = [cons, real(H(:,k)'*W_next(((k-1)*N_tx+1):(k*N_tx),g))>=real(W_bar_next(k,g)^2)];

		% add third term in Eq. (33)
		
% 		obj = obj - real(varpi)*z_tilde_next(k,g);
%         PerGroupDownlinkRate = PerGroupDownlinkRate - time_next(g)*real(varpi)*z_tilde_next(k,g);
        PerGroupDownlinkRate(k,g) = PerGroupDownlinkRate(k,g) - real(varpi)*z_tilde_next(k,g);
%         Downlink_Rate = Downlink_Rate - real(varpi)*z_tilde_next(k,g);

%         cons = [cons, z_tilde_next(k,g)^2<=z_bar_next(k,g)];

        cons = [cons, cone([phi_next(k,g), real(H(:,k)'*W_next(((k-1)*N_tx+1):(k*N_tx),g)), 0.5*(z_tilde_next(k,g)-1)], 0.5*(z_tilde_next(k,g)+1) ) ];

        cons = [cons, alpha_next(k,g)>=0];
        
        cons = [cons, alpha_next(k,g)<=1];
        
        cons = [cons, alpha_next(k,g)<=BarrierValue*PerGroupDownlinkRate(k,g)];
        
        cons = [cons, PerGroupDownlinkRate(k,g)>=0];
        
%         cons = [cons, alpha_next(k,g)<=real(H(:,k)'*W_next(((k-1)*N_tx+1):(k*N_tx),g))];

        cons = [cons, cone([alpha_bar_next(k,g), 0.5*(alpha_next(k,g)-PerGroupDownlinkRate(k,g))], 0.5*(alpha_next(k,g)+PerGroupDownlinkRate(k,g))) ];
        
        cons = [cons, alpha_bar_current(k,g)^2+2*alpha_bar_current(k,g)*(alpha_bar_next(k,g)-alpha_bar_current(k,g))>=alpha_tilde_next(k,g)];
        
%         cons = [cons, cone([RThDown_next(k,g) 0.5*(time_next(g)-PerGroupDownlinkRate(k,g))],0.5*(time_next(g)+PerGroupDownlinkRate(k,g)))];
        cons = [cons, cone([RThDown_next(k,g) 0.5*(time_next(g)-alpha_tilde_next(k,g))],0.5*(time_next(g)+alpha_tilde_next(k,g)))];
        
        cons = [cons, RThDown_current(k,g)^2+2*RThDown_current(k,g)*(RThDown_next(k,g)-RThDown_current(k,g))>= RThDown_bar_next(k,g)];
        
        cons = [cons, RThDown_bar_next(k,g)>=0];
		        
    end

%     cons = [cons, real(PerGroupDownlinkRate) >= GammaDown_next(g)];
%     cons = [cons, real(sum(alpha_bar_next(:,g))) >= GammaDown_next(g)];
    cons = [cons, real(sum(alpha_tilde_next(:,g))) >= GammaDown_next(g)];
    
    cons = [cons, cone([RDown_next(g), 0.5*(time_next(g)-GammaDown_next(g))],0.5*(time_next(g)+GammaDown_next(g)))];
    
%     cons = [cons, RDown_next(g)>=0];
    
end

cons = [cons, sum(RThDown_bar_next,2)>=log(2)*Rate_Threshold];

cons = [cons, sum(alpha_next,2)>=1];

% obj = obj + Downlink_Rate;
% obj = obj + sum(time_next.*GammaDown_next);

RDown_bar_next = sdpvar(1,G);

cons = [cons, RDown_current.^2+2*RDown_current.*(RDown_next - RDown_current)>=RDown_bar_next];

% X = sdpvar(1,1);

% cons = [cons, X <= sum(RDown_bar_next)];

obj = obj + sum(RDown_bar_next);
% obj = obj + X;



% add objective function and constraints for UPLINK rate in Eq. (36)

% upsilon_bar_next = sdpvar(L, G); % 1/upsilon_next <= upsilon_bar_next
p_bar_next = sdpvar(L, G); % sqrt(p_next) <= p_bar_next
s_next = sdpvar(L, G); % lambda_vec/upsilon_next <= s_next


GammaUp_next = sdpvar(1,G); % = sum{1->L} of log2(1+gamma(l,g)(w,p))
RUp_next = sdpvar(1,G); % = time(g) * GammaUp(g)
RTh_next = sdpvar(L,G); % threshold for uplink rate
RTh_bar_next = sdpvar(L,G);

beta_bar_next = sdpvar(L,G); % alpha_next*log2()>= alpha_bar_next^2
beta_tilde_next = sdpvar(L,G); % alpha_bar_next^2 >= alpha_tilde_next


Uplink_Rate = 0;
UplinkRate_PerGroupPerUser = sdpvar(L,G);
PerGroupUplinkRate = sdpvar(L,G);

for g = 1:1:G

%     PerGroupUplinkRate = 0;
    
    for l = 1:1:L
			
		% get parameters in Eq. (21)
		[vartheta, psi, lambda_vec] = GetParameter_Uplink(G_channel(:,l:L), G_SI, W_current(:,g), p_current(l:L,g), W_next(:,g), p_next(l:L,g), rho);
        
		% add first term in Eq. (36)
		
% 		PerGroupUplinkRate = PerGroupUplinkRate + time_next(g)*real(vartheta);
        PerGroupUplinkRate(l,g) = real(vartheta);
%         PerGroupUplinkRate = PerGroupUplinkRate + real(vartheta);
%         Uplink_Rate = Uplink_Rate + real(vartheta);
% 		
% 		cons = [cons, cone([1, 0.5*(upsilon_next(l,g)-upsilon_bar_next(l,g))], 0.5*(upsilon_next(l,g)+upsilon_bar_next(l,g))) ];
        
%         ratio_uu = sqrt(upsilon_current(l,g)/upsilon_bar_current(l,g));
%         cons = [cons, cone([sqrt(1/2)*upsilon_next(l,g)/ratio_uu, sqrt(1/2)*upsilon_bar_next(l,g)*ratio_uu], 1)];
		
		
		% add second term in Eq. (36)
		
% 		PerGroupUplinkRate = PerGroupUplinkRate + time_next(g)*real(psi)*p_next(l,g);
        PerGroupUplinkRate(l,g) = PerGroupUplinkRate(l,g) + real(psi)*p_next(l,g);
%         PerGroupUplinkRate = PerGroupUplinkRate + real(psi)*p_next(l,g);
%         Uplink_Rate = Uplink_Rate + real(psi)*p_next(l,g);
		
% 		cons = [cons, sqrt(real(p_next(l,g))) <= p_bar_next(l,g)];
%         cons = [cons, cone(sqrt(real(p_next(l,g))),p_bar_next(l,g))];
		
        
%         cons = [cons, real(p_next(l,g))>=real(p_bar_next(l,g)^2)];

		% add third term in Eq. (36)
		
% 		PerGroupUplinkRate = PerGroupUplinkRate - time_next(g)*s_next(l,g);
        PerGroupUplinkRate(l,g) = PerGroupUplinkRate(l,g) - s_next(l,g);
%         PerGroupUplinkRate = PerGroupUplinkRate - s_next(l,g);
%         Uplink_Rate = Uplink_Rate - s_next(l,g);

        UplinkRate_PerGroupPerUser(l,g) = PerGroupUplinkRate(l,g);
		
		cons = [cons, cone([lambda_vec; 0.5*(s_next(l,g)-1)], 0.5*(s_next(l,g)+1) ) ];
        
        cons = [cons, beta_next(l,g)>=0];
        
        cons = [cons, beta_next(l,g)<=1];
        
        cons = [cons, beta_next(l,g)<=BarrierValue*PerGroupUplinkRate(l,g)];
        
        cons = [cons, PerGroupUplinkRate(l,g)>=0];

        cons = [cons, cone([beta_bar_next(l,g), 0.5*(beta_next(l,g)-PerGroupUplinkRate(l,g))], 0.5*(beta_next(l,g)+PerGroupUplinkRate(l,g))) ];
        
        cons = [cons, beta_bar_current(l,g)^2+2*beta_bar_current(l,g)*(beta_bar_next(l,g)-beta_bar_current(l,g))>=beta_tilde_next(l,g)];
        
%         cons = [cons, cone([RTh_next(l,g) 0.5*(time_next(g)-PerGroupUplinkRate(l,g))],0.5*(time_next(g)+PerGroupUplinkRate(l,g)))];

        cons = [cons, cone([RTh_next(l,g) 0.5*(time_next(g)-beta_tilde_next(l,g))],0.5*(time_next(g)+beta_tilde_next(l,g)))];
        
        cons = [cons, RTh_current(l,g)^2+2*RTh_current(l,g)*(RTh_next(l,g)-RTh_current(l,g))>= RTh_bar_next(l,g)];
        
        cons = [cons, RTh_bar_next(l,g)>=0];
        
    end
    
%     cons = [cons, real(PerGroupUplinkRate) >= GammaUp_next(g)];
    cons = [cons, real(sum(beta_tilde_next(:,g))) >= GammaUp_next(g)];
    
    cons = [cons, cone([RUp_next(g), 0.5*(time_next(g)-GammaUp_next(g))],0.5*(time_next(g)+GammaUp_next(g)))];
    
%     cons = [cons, RUp_next(g)>=0];
	
end

cons = [cons, sum(RTh_bar_next,2)>=log(2)*Rate_Threshold];

% obj = obj + Uplink_Rate;
% obj = obj + sum(time_next.*GammaUp_next);

RUp_bar_next = sdpvar(1,G);

cons = [cons, RUp_current.^2+2*RUp_current.*(RUp_next - RUp_current)>=RUp_bar_next];

% Y = sdpvar(1,1);

% cons = [cons, Y <= sum(RUp_bar_next)];

obj = obj + sum(RUp_bar_next);
% obj = obj + Y;

% cons = [cons, Uplink_Rate>=1];


%% Add power constraints

cons = [cons, vec(p_next)>=0]; %8d 10d

% cons = [cons, sum(p_next.^2*diag(time_next),2)<=P]; %8d 10c

p_bar_next = sdpvar(L,G);

pl_perUser_next = sdpvar(L,G);


cons = [cons, p_next.^2 <= p_bar_next];


for l = 1:1:L
    
    for g = 1:1:G
                 
        ratio_time_pbar = sqrt(time_current(g)/p_bar_current(l,g));
        
        cons = [cons, cone([sqrt(1/2)*time_next(g)/ratio_time_pbar, sqrt(1/2)*p_bar_next(l,g)*ratio_time_pbar, 0.5*(pl_perUser_next(l,g)-1)],0.5*(pl_perUser_next(l,g)+1))];
        
    end
    
end

cons = [cons, sum(pl_perUser_next,2) <= P];


% cons = [cons, vec(p_next)<=P.^(0.5)]; %8d 10c

cons = [cons, time_next > 0 ]; %8g

cons = [cons, sum(time_next)<=1]; %8g

% cons = [cons, tau_next > 0];


% 12a 12b
Beam_intime = [];

w_bar_next = sdpvar(1,G);

% mu_next = sdpvar(K,G);
% w_bar_next = sdpvar(K,G);

% mu_bar_next = sdpvar(K,G);
w_tilde_next = sdpvar(1,G);
pd_onegroup_next = sdpvar(1,G);

for g = 1:1:G
    
    
    W_next_mat_all = reshape(W_next(:,g), N_tx, K);
    
    cons = [cons, cone(W_next(:,g), w_bar_next(g))];
%     
%     cons = [cons, real(w_bar_next(g)^2) <= w_tilde_next(g)];
    cons = [cons, cone([w_bar_next(g), 0.5*(w_tilde_next(g)-1)], 0.5*(w_tilde_next(g)+1)) ];
    
    

    for k = 1:1:K
		
		W_next_mat = W_next_mat_all(:,[1:(k-1) (k+1):K]);
		
% 		SumCons = sum((H(:,k)'*W_next_mat).^2) + sum((p_next(:,g).^2).*(G_hat(:,k).^2)) + sigma_K(k)^2;
        Sum_vec = [(H(:,k)'*W_next_mat) (p_next(:,g).*G_hat(:,k))' sigma_K(k)];
		
% 		cons = [cons, sqrt(real(SumCons))<=phi_next(k,g)];
        cons = [cons, cone(Sum_vec, phi_next(k,g))]; % 12a
        
        cons = [cons, real(H(:,k)'*W_next(((k-1)*N_tx+1):k*N_tx,g))>=0]; % 12b
        
%         cons = [cons, mu_next(k,g)>=0];
        
%         cons = [cons, cone(W_next_mat_all(:,k), w_bar_next(k,g))];
        
%         cons = [cons, cone([w_bar_next(k,g) 0.5*(alpha_next(k,g)-mu_next(k,g))],0.5*(alpha_next(k,g)+mu_next(k,g)))];
        
%         ratio_alpha_mu = sqrt(alpha_current(k,g)/mu_current(k,g));
        
%         cons = [cons, cone([sqrt(1/2)*alpha_next(k,g)/ratio_alpha_mu, sqrt(1/2)*mu_next(k,g)*ratio_alpha_mu], mu_bar_next(k,g) )];
		
    end
    
%     cons = [cons, cone([mu_bar_next(:,g); 0.5*(w_tilde_next(g)-1)], 0.5*(w_tilde_next(g)+1))];
    
    ratio_time_wtilde = sqrt(time_current(g)/w_tilde_current(g));
    
    cons = [cons, cone([sqrt(1/2)*time_next(g)/ratio_time_wtilde, sqrt(1/2)*w_tilde_next(g)*ratio_time_wtilde, 0.5*(pd_onegroup_next(g)-1)], 0.5*(pd_onegroup_next(g)+1))];
    
% 	cons = [cons, cone([sqrt(time_next(g))*W_next(:,g)],sqrt(Pbs))]; % 10b
%     Beam_intime = [Beam_intime; sqrt(time_next(g))*W_next(:,g)];
%     Beam_intime = [Beam_intime; sqrt(time_next(g))*w_tilde_next(g)];
%     cons = [cons, cone([W_next(:,g)],sqrt(Pbs))]; % 10b
end

% cons = [cons, cone(Beam_intime,sqrt(Pbs))]; % 10b
% cons = [cons, sum(time_next.*w_tilde_next) <=Pbs]; % 10b
cons = [cons, sum(pd_onegroup_next) <=Pbs]; % 10b


% 27b 27c
% 
% cons = [cons, alpha_next>=0];
% cons = [cons, alpha_next<=1];

% for k = 1:1:K
%     for g = 1:1:G
%         cons = [cons, -alpha_next(k,g)<=0];
%         cons = [cons, alpha_next(k,g)<=1];
%     end
% end

% cons = [cons, beta_next>=0];
% cons = [cons, beta_next<=1];

% for l = 1:1:K
%     for g = 1:1:G
%         cons = [cons, -beta_next(l,g)<=0];
%         cons = [cons, beta_next(l,g)<=1];
%     end
% end

% x_next = sdpvar(K, G, 'full', 'real');
% x_bar_next = sdpvar(K, G, 'full', 'real');
% y_next = sdpvar(L, G, 'full', 'real');
% y_bar_next = sdpvar(L, G, 'full', 'real');

% for k = 1:1:K
%     for g = 1:1:G
% %         cons = [cons, cone([x_next(k,g), 0.5*(alpha_next(k,g)-time_next(g))], 0.5*(alpha_next(k,g)+time_next(g)) ) ];
% %         cons = [cons, x_current(k,g)*(2*x_next(k,g) - x_current(k,g))>=x_bar_next(k,g)];
% %         cons = [cons, cone([1, 0.5*(z_next(k,g) - x_bar_next(k,g))],0.5*(z_next(k,g) + x_bar_next(k,g)))];
%         cons = [cons, cone([1, 0.5*(z_next(k,g) - time_next(g))],0.5*(z_next(k,g) + time_next(g)))];
%     end
% end



% for l = 1:1:L
%     for g = 1:1:G
% %         cons = [cons, cone([y_next(l,g), 0.5*(beta_next(l,g)-time_next(g))], 0.5*(beta_next(l,g)+time_next(g)) ) ];
% %         cons = [cons, y_current(l,g)*(2*y_next(l,g) - y_current(l,g))>=y_bar_next(l,g)];
% %         cons = [cons, cone([1, 0.5*(upsilon_next(l,g) - y_bar_next(l,g))],0.5*(upsilon_next(l,g) + y_bar_next(l,g)))];
%         cons = [cons, cone([1, 0.5*(upsilon_next(l,g) - time_next(g))],0.5*(upsilon_next(l,g) + time_next(g)))];
%     end
% end



% 28
% tau_tilde_next = sdpvar(1,G,'full','real');
% for g = 1:1:G
% % 	cons = [cons, cone([1, 0.5*(tau_next(g)-time_next(g))], 0.5*(tau_next(g)+time_next(g)) )];
% %     cons = [cons, cone([tau_next(g) time_next(g)], sqrt(2)) ];
% 
% %     a_bar = sqrt(time_current(g)/tau_current(g));
% %     cons = [cons, cone([sqrt(1/2)*time_next(g)/a_bar, sqrt(1/2)*tau_next(g)*a_bar], 1 ) ];
%     cons = [cons, cone([1, 0.5*(tau_next(g)-tau_tilde_next(g))], 0.5*(tau_next(g)+tau_tilde_next(g)) ) ];
% end
% 
% cons = [cons, sum(tau_tilde_next)<=1];

% 30 31

    
% for g = 1:1:G
% 
%     for k = 1:1:K
%         
%         ratio_alpha_z = sqrt(alpha_current(k,g)/z_current(k,g));
%         
%         cons = [cons, cone([sqrt(1/2)*alpha_next(k,g)/ratio_alpha_z, sqrt(1/2)*z_next(g)*ratio_alpha_z, 0.5*(tau_next(g)-1)], 0.5*(tau_next(g)+1) ) ];
%         
%     end
%     
%     for l = 1:1:L
%         
%         ratio_beta_upsilon = sqrt(beta_current(l,g)/upsilon_current(l,g));
%         
%         cons = [cons, cone([sqrt(1/2)*beta_next(l,g)/ratio_beta_upsilon, sqrt(1/2)*upsilon_next(g)*ratio_beta_upsilon, 0.5*(tau_next(g)-1)], 0.5*(tau_next(g)+1) ) ];
%         
%     end
%     
% end

% 39a 39b

% cons = [cons, sum(mu_next)<=Pbs];

% cons = [cons, mu_next > 0];

% for g = 1:1:G
% % 	cons = [cons, cone([W_next(:,g); 0.5*(mu_next(g)-tau_next(g))], 0.5*(mu_next(g)+tau_next(g)) )];
%     cons = [cons, cone([W_next(:,g); 0.5*(mu_next(g)-1)], 0.5*(mu_next(g)+1) )];
% end

% 40a 40b

% cons = [cons, vec(theta_next) > 0];

% for l = 1:1:L
% 
% 	cons = [cons, sum(theta_next(l,:)) <= P(l)]; %40a
% 
% 	for g = 1:1:G
% 	
% % 		cons = [cons, cone([p_next(l,g), 0.5*(theta_next(l,g)-tau_next(g))],0.5*(theta_next(l,g)+tau_next(g)) )]; % 40b
%         cons = [cons, cone([p_next(l,g), 0.5*(theta_next(l,g)-1)],0.5*(theta_next(l,g)+1) )]; % 40b
% 		
% 	end
% 	
% end

%41a

% for k = 1:1:K
% 	
% 	for g = 1:1:G
% 	
% 		cons = [cons, cone([eta_next(g), 0.5*(alpha_next(k,g)-z_next(k,g))], 0.5*(alpha_next(k,g)-z_next(k,g)))];
% 	
% 	end
% 	
% end


%43

% cons = [cons, tau_next <= eta_current.^2 + 2*eta_current.*(eta_next-eta_current)];

%44a 44b

% for l = 1:1:L
% 
% 	for g = 1:1:G
% 	
% 		cons = [cons, cone([kapa_next(g), 0.5*(beta_next(l,g)-upsilon_next(l,g))], 0.5*(beta_next(l,g)+upsilon_next(l,g)))];
% 		
% 	end
% 	
% end
% 
% 
% cons = [cons, tau_next <= kapa_current.^2 + 2*kapa_current.*(kapa_next-kapa_current)];


myops = sdpsettings('solver','sdpt3','verbose',0);

% obj_t = sdpvar(1,1);
% 
% cons = [cons, obj_t<=obj];

diagnotics = solvesdp(cons, -obj, myops);
	
OptimalValue = double(obj);

RDown = real(double(RDown_next));

RThDown = real(double(RThDown_next));

DownlinkRate_PerUser = real(double(PerGroupDownlinkRate));

alpha = real(double(alpha_next));

alpha_bar = real(double(alpha_bar_next));

% mu = real(double(mu_next));

RUp = real(double(RUp_next));

RTh = real(double(RTh_next));

UplinkRate_PerUser = real(double(PerGroupUplinkRate));

beta = real(double(beta_next));

beta_bar = real(double(beta_bar_next));

% UplinkRate = double(Uplink_Rate)

W = double(W_next);
p = double(p_next);
phi = double(phi_next);
% eta = double(eta_next);
% z = double(z_next);
% z_bar = double(z_bar_next);
% kapa = double(kapa_next);
% upsilon = double(upsilon_next);
% upsilon_bar = double(upsilon_bar_next);

% alpha = double(alpha_next)
% beta = double(beta_next)
time = real(double(time_next));
w_tilde = real(double(w_tilde_next));
p_bar = real(double(p_bar_next));

% x = double(x_next);
% y = double(y_next);

% tau = double(tau_next);
% mu = double(mu_next);
% theta = double(theta_next);
% z_tilde = double(z_tilde_next);

end

