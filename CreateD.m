function [ D_H, D_G_channel, D_G_hat ] = CreateD( K, L, Parameters, AllCells, Order )
%CREATED Summary of this function goes here
%   Detailed explanation goes here

RadiusOfCell = Parameters(1);
RadiusOfNearestUser = Parameters(2);
StandardDeviation = 10^(Parameters(3)/10);
ploss = Parameters(4);


%% large-fading for K downlink users    
    Zvector = StandardDeviation*randn(1,K);
    rvector = RadiusOfNearestUser*ones(1,K) + (RadiusOfCell-RadiusOfNearestUser)*rand(1,K);
    anglevector = 2*pi*rand(1,K);
    positionDownlinkUsers = [(rvector.*cos(anglevector))' (rvector.*sin(anglevector))'];
%     positionBSj = AllCells(1,:);
        
%     positionBSl = positionBSj;
        
%     vector_BSl_To_UsersInj = repmat(positionBSj,K,1) + positionUsersInCellj - repmat(positionBSl,K,1);
        
    distanceBS_To_DownlinkUsers = sqrt(sum(positionDownlinkUsers'.^2));
        
    betavector_downlink = (10.^(Zvector/10))./((distanceBS_To_DownlinkUsers/RadiusOfNearestUser).^ploss);

        
    D_H = (diag(betavector_downlink)).^(0.5)
    
    
%% large-fading for L uplink users  

    Zvector = StandardDeviation*randn(1,L);
    rvector = RadiusOfNearestUser*ones(1,L) + (RadiusOfCell-RadiusOfNearestUser)*rand(1,L);
    anglevector = 2*pi*rand(1,L);
    positionUplinkUsers = [(rvector.*cos(anglevector))' (rvector.*sin(anglevector))'];
%     positionBSj = AllCells(1,:);
        
%     positionBSl = positionBSj;
        
%     vector_BSl_To_UsersInj = repmat(positionBSj,K,1) + positionUsersInCellj - repmat(positionBSl,K,1);
        
    distanceBS_To_UplinkUsers = sqrt(sum(positionUplinkUsers'.^2));
        
    betavector_uplink = (10.^(Zvector/10))./((distanceBS_To_UplinkUsers/RadiusOfNearestUser).^ploss);

        
    D_G_channel = (diag(betavector_uplink)).^(0.5)
    
    
%% large-fading between uplink and downlink users

    Zvector = StandardDeviation*randn(L*K,1);
    positionDownUsers_Ltimes = kron(positionDownlinkUsers,ones(L,1));
    positionUpUsers_Ktimes = repmat(positionUplinkUsers, K, 1);
    
    distance_UpDownUsers = sqrt(sum((positionDownUsers_Ltimes-positionUpUsers_Ktimes).^2 ,2));
    
    betavector_updownusers = (10.^(Zvector/10))./((distance_UpDownUsers/RadiusOfNearestUser).^ploss);
    
    D_G_hat = (reshape(betavector_updownusers, L, K)).^(0.5)
        
        
%         if ((j==1))
%             j
%             l
%             D(:,:,j,l)
%         end
%         


%         posl = AllCells(Order==l,:);
%         distancejl = norm(positionBSj - posl);
%         
%         if (distancejl<1000)
%             D(:,:,j,l) = eye(K);
%         elseif (distancejl<2500)
%             D(:,:,j,l) = 0.8*eye(K);
%         else
%             D(:,:,j,l) = 0.2*eye(K);
%         end

end

