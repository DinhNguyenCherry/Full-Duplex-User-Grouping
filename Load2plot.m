
load('Rate_all_4.mat', 'rho_dB')



Rate_AllChannel1 = []; % 2 Groups Alg1
Rate_AllChannel2 = []; % 2 Groups Alg2
Rate_AllChannel3 = []; % 1 Group

load('Rate_all_4.mat', 'Alg1_Rate_all')
load('Rate_all_4.mat', 'Alg2_Rate_all')

Rate_AllChannel1 = [Rate_AllChannel1; Alg1_Rate_all];
Rate_AllChannel2 = [Rate_AllChannel2; Alg2_Rate_all];

load('Rate_all_5.mat', 'Alg1_Rate_all')
load('Rate_all_5.mat', 'Alg2_Rate_all')

SelectedChannel = [1 3:5];

Rate_AllChannel1 = [Rate_AllChannel1; Alg1_Rate_all(SelectedChannel,:)];
Rate_AllChannel2 = [Rate_AllChannel2; Alg2_Rate_all(SelectedChannel,:)];

load('Rate_all_6.mat', 'Alg1_Rate_all')
load('Rate_all_6.mat', 'Alg2_Rate_all')

SelectedChannel = [1:3 5];

Rate_AllChannel1 = [Rate_AllChannel1; Alg1_Rate_all(SelectedChannel,:)];
Rate_AllChannel2 = [Rate_AllChannel2; Alg2_Rate_all(SelectedChannel,:)];

load('Rate_all_7.mat', 'Alg2_Rate_all')
SelectedChannel = [2 4:7 9:10];
Rate_AllChannel2 = [Rate_AllChannel2; Alg2_Rate_all(SelectedChannel,:)];

load('Rate_all_8.mat', 'Alg2_Rate_all')
SelectedChannel = [1 2 3 4 7 9:11 13 14 16 18:22 24 26 27 29 30];
Rate_AllChannel2 = [Rate_AllChannel2; Alg2_Rate_all(SelectedChannel,:)];

Rate_AllChannel2post = [];

load('Rate_all_9.mat', 'Alg2_Rate_all')
SelectedChannel = [1:12 14 15];
Rate_AllChannel2post = [Rate_AllChannel2post; Alg2_Rate_all(SelectedChannel,:)];

load('Rate_all_10.mat', 'Alg2_Rate_all')
SelectedChannel = [1 2 4 6:15 17 19 20 22 25];
Rate_AllChannel2post = [Rate_AllChannel2post; Alg2_Rate_all(SelectedChannel,:)];

load('Rate_all_11.mat', 'Alg2_Rate_all')
SelectedChannel = [1:4 7:10];
Rate_AllChannel2post = [Rate_AllChannel2post; Alg2_Rate_all(SelectedChannel,:)];

load('Rate_all_12.mat', 'Alg2_Rate_all')
SelectedChannel = [2 3 6 7 10];
Rate_AllChannel2post = [Rate_AllChannel2post; Alg2_Rate_all(SelectedChannel,:)];

load('Rate_all_13.mat', 'Alg2_Rate_all')
SelectedChannel = [1 2 6 7 9:27 29 30];
Rate_AllChannel2post = [Rate_AllChannel2post; Alg2_Rate_all(SelectedChannel,:)];

load('Rate_all_14.mat', 'Alg2_Rate_all')
SelectedChannel = [1:3 5 7:9];
Rate_AllChannel2post = [Rate_AllChannel2post; Alg2_Rate_all(SelectedChannel,:)];

load('Rate_all_15_1Group.mat', 'Alg1_Rate')
SelectedChannel = [1:30];
size(Rate_AllChannel3)
% Rate_AllChannel3 = [Rate_AllChannel3; Alg1_Rate_all(SelectedChannel,:)]
Rate_AllChannel3 = Alg1_Rate

Rate_AllChannel1
size(Rate_AllChannel2)
Rate_AllChannel2
size(Rate_AllChannel2post)
Rate_AllChannel2post

Rate_AllChannel2post = mean(Rate_AllChannel2post);


Alg1_Rate = mean(Rate_AllChannel1)
Alg2_Rate = mean(Rate_AllChannel2)
Alg2_Rate = mean([Alg2_Rate; Alg2_Rate(1:5) Rate_AllChannel2post])
Alg3_Rate = Rate_AllChannel3

hold on
plot(rho_dB, Alg1_Rate, 'rs--', 'linewidth', 2, 'markersize',9);
plot(rho_dB, Alg2_Rate, 'r+-', 'linewidth', 2, 'markersize',9);
plot(rho_dB, Alg3_Rate, 'bo-', 'linewidth', 2, 'markersize',9);