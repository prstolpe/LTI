%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     author: Mario Senden (mario.senden@maastrichtuniversity.nl)     %%%

% This is an implementation of the simulation experiments described in
% Lange, G., Senden, M., Radermacher, A., & De Weerd, P. (2020). 
% Interfering with a memory without erasing its trace. 
% Neural Networks, 121, 339–355. 
% https://doi.org/10.1016/j.neunet.2019.09.027

clear all;close all;clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             settings                                %%%

OD_0     =   4.5;           % initial orientation difference
Sessions =   8;             % number of sessions
Reps     =  2;             % number of times each experiment is repeated
Trials   = 50;             % number of trials per session


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             parameters                              %%%

N        = 512;             % number of neurons
alpha    =  10;             % gain of spike encoder
sigma_ff =  45;             % width of feedforward bias
J_ff     =   0.5;           % forward connection strength
J_rec    =   1;             % recurrent connection strength
a_e      =   2.2;           % exponent exc. connections
a_i      =   1.4;           % exponent inh. connections
c_e      =   1.2025e-3;     % normalization exc. connection
c_i      =   1.6875e-3;     % normalization inh. connection
k        =   1.05;          % scaling of variance
C        =   0.53;          % decision criterion
eta      =   1.4e-11;       % learning rate
mu       =   0;             % exponent of power law weight dependence
t_sim    =   0.5;           % simulation time (seconds)
tau      =   1.5e-2;        % membrane time constant (seconds)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                 setup                               %%%

Q = cell(3,1);
Exp = cell(3,1);    Int = cell(3, 1);
qp = cell(3,1);     qr = cell(3,1);     rp = cell(3,1);     rr = cell(3,1);
for i=1:3
    Q{i} = RM(...   % create a model for each of three quadrants
        N,...               % (i.e. experiments)
        alpha,...
        sigma_ff,...
        J_ff,...
        J_rec,...
        a_e,...
        a_i,...
        c_e,...
        c_i,...
        k,...
        C,...
        eta,...
        mu,...
        t_sim,...
        tau,...
        Trials,...
        OD_0);
    Exp{i}.Ab = zeros(Reps,Sessions);
    Exp{i}.At = zeros(Reps,Sessions);
    Int{i}.left = zeros(Reps,Sessions);
    Int{i}.right = zeros(Reps,Sessions);
    qp{i}.Ab        = NaN(Reps, Sessions, Trials);
    qp{i}.At        = NaN(Reps, Sessions, Trials);
    qp{i}.left      = NaN(Reps, Sessions, Trials);
    qp{i}.right     = NaN(Reps, Sessions, Trials);
    
    qr{i}.Ab = NaN(Reps, Sessions, Trials);
    qr{i}.At = NaN(Reps, Sessions, Trials);
    qr{i}.left = NaN(Reps, Sessions, Trials);
    qr{i}.right = NaN(Reps, Sessions, Trials);
    
    rp{i}.Ab = NaN(Reps, Sessions, Trials, N);
    rp{i}.At = NaN(Reps, Sessions, Trials, N);
    rp{i}.left = NaN(Reps, Sessions, Trials, N);
    rp{i}.right = NaN(Reps, Sessions, Trials, N);
    
    rr{i}.Ab = NaN(Reps, Sessions, Trials, N);
    rr{i}.At = NaN(Reps, Sessions, Trials, N);
    rr{i}.left = NaN(Reps, Sessions, Trials, N);
    rr{i}.right = NaN(Reps, Sessions, Trials, N);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             experiments                             %%%
                                    % 45 = left
% Exp1  (green: 135°    ->          45°          -> 135°)
% Exp2  (blue:  135°    ->      105° & 165°      -> 135°)
% Exp3  (red:    //     ->      105° & 165°      -> 135°)
% W       = NaN(Sessions, Trials, N, N);
% W_exc   = NaN(Sessions, Trials, N, N);
% W_inh   = NaN(Sessions, Trials, N, N);
% V_ff    = NaN(Sessions, Trials, N);
% v       = NaN(Sessions, Trials, 100, N);
% correct = NaN(Sessions, Trials);
% Phi_probe = NaN(Sessions, Trials);
% r_ref   = NaN(Sessions, Trials, 100, N);
% r_pro   = NaN(Sessions, Trials, 100, N);
% % dW_inh  = NaN(Sessions, Trials, N);
% m_JND   = NaN(Sessions, Trials);
% JND     = NaN(Sessions, Trials);


for r=1:Reps
    fprintf('\n - participant %.2d',r)
    % part 1 (135° - baseline)
%     Q{1}.set_PHI(135);
    Q{2}.set_PHI(135);
    for s=1:Sessions
%         [W(s, :, :, :), V_ff(s, :, :), v(s, :, :, :), correct(s, :), W_exc(s, :, :, :), W_inh(s, :, :, :), Phi_probe(s, :), r_ref(s, :, :, :), r_pro(s, :, :, :), m_JND(s, :), JND(s, :)] = Q{1}.session();
        [qp{2}.Ab(r, s, :), qr{2}.Ab(r, s, :), rp{2}.Ab(r, s, :, :), rr{2}.Ab(r, s, :, :)] = Q{2}.session();
%         Exp{1}.Ab(r,s) = Q{1}.get_JND * 2;
        Exp{2}.Ab(r,s) = Q{2}.get_JND * 2;
    end
    
    % part 2a (105° & 45° - interference)
%     Q{1}.set_PHI(45);
    Q{2}.set_PHI(105);
%     Q{3}.set_PHI(105);
%     Q{1}.set_OD();
    Q{2}.set_OD();
%     Q{3}.set_OD();
    for s=1:Sessions
%         Q{1}.session();
        [qp{2}.left(r, s, :), qr{2}.left(r, s, :), rp{2}.left(r, s, :, :), rr{2}.left(r, s, :, :)] = Q{2}.session();
%         Q{3}.session();
%         Int{1}.left(r,s) = Q{1}.get_JND * 2;
        Int{2}.left(r,s) = Q{2}.get_JND * 2;
%         Int{3}.left(r,s) = Q{3}.get_JND * 2;
    end
    
%     % part 2b (165° - interference)
    Q{2}.set_PHI(165);
%     Q{3}.set_PHI(165);
    Q{2}.set_OD();
%     Q{3}.set_OD();
    for s=1:Sessions
        [qp{2}.right(r, s, :), qr{2}.right(r, s, :), rp{2}.right(r, s, :, :), rr{2}.right(r, s, :, :)] = Q{2}.session();
%         Q{3}.session();
        Int{2}.right(r,s) = Q{2}.get_JND * 2;
%         Int{3}.right(r,s) = Q{3}.get_JND * 2;
%         
    end
    
    % part 3 (135° - test)
%     Q{1}.set_PHI(135);
    Q{2}.set_PHI(135);
%     Q{3}.set_PHI(135);
%     Q{1}.set_OD(Exp{1}.Ab(r,end));
    Q{2}.set_OD(Exp{2}.Ab(r,end));
%     Q{3}.set_OD();
    for s=1:Sessions
%         Q{1}.session();
        [qp{2}.At(r, s, :), qr{2}.At(r, s, :), rp{2}.At(r, s, :, :), rr{2}.At(r, s, :, :)] = Q{2}.session();
%         Q{3}.session();
%         Exp{1}.At(r,s) = Q{1}.get_JND * 2;
        Exp{2}.At(r,s) = Q{2}.get_JND * 2;
%         Exp{3}.At(r,s) = Q{3}.get_JND * 2;
    end
    
%     Q{1}.reset();
    Q{2}.reset();
%     Q{3}.reset();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             plotting                                %%%

% Pos = [200 200  950 350];
% figure('Color','w','Position' ,Pos)

% experiment 1
% subplot(1,3,1,'Fontsize',9)
% hold all
% figure(1)
% hold all
% plot(mean(Exp{1}.Ab),'color',[0 .75 0],'linestyle','--','linewidth',2.5)
% plot(mean(Exp{1}.At),'color',[0 .75 0],'linewidth',2.5)
% set(gca, 'XTick', 1:8)
% set(gca, 'YScale', 'log')
% xlim([0.5 8.5])
% ylim([1.5 8.5])
% xlabel('session')
% title('Experiment 1 (ACA)')
% legend('A_B','A_T')
% legend('boxoff')

% experiment 2
figure(2)
subplot(121)
hold all, grid on
plot(mean(Exp{2}.Ab),'color',[0 0 .75],'linestyle','--','linewidth',2.5)
plot(mean(Exp{2}.At),'color',[0 0 .75],'linewidth',2.5)
set(gca, 'XTick', 1:8)
set(gca, 'YScale', 'log')
xlim([0.5 8.5])
ylim([1.5 8.5])
xlabel('session')
ylabel('JND [deg]')
title('Experiment 2 (ABA)')
legend('A_B','A_T')
legend('boxoff')

subplot(122)
hold all, grid on 
plot(mean(Int{2}.left),'r-.', 'linewidth',2.5)
plot(mean(Int{2}.right),'g--', 'linewidth',2.5)
set(gca, 'XTick', 1:8)
set(gca, 'YScale', 'log')
xlim([0.5 8.5])
ylim([1.5 8.5])
xlabel('session')
ylabel('JND [deg]')
title('Experiment 2 (ABA)')
legend('$105^\circ$','$165^\circ$', 'interpreter', 'latex')
legend('boxoff')

figure(3)
subplot(121)
plot(mean(squeeze(qp{2}.Ab(1, :, :)), 2), 'color',[0 0 .75],'linestyle','--','linewidth',2.5), hold on
plot(mean(squeeze(qp{2}.At(1, :, :)), 2), 'color',[0 0 .75],'linewidth',2.5), grid on
xlabel('session')
xlim([0.5 8.5])

subplot(122)
plot(mean(squeeze(qp{2}.left(1, :, :)), 2), 'r-.', 'linewidth',2.5), grid on, hold on
plot(mean(squeeze(qp{2}.right(1, :, :)), 2),'g--', 'linewidth',2.5), hold off
xlabel('session')
xlim([0.5 8.5])
legend('$105^\circ$','$165^\circ$', 'interpreter', 'latex')

figure(4)
subplot(121)
plot(mean(squeeze(qr{2}.Ab(1, :, :)), 2), 'color',[0 0 .75],'linestyle','--','linewidth',2.5), hold on 
plot(mean(squeeze(qr{2}.At(1, :, :)), 2), 'color',[0 0 .75],'linewidth',2.5), grid on
xlabel('session')
xlim([0.5 8.5])

subplot(122)
plot(mean(squeeze(qr{2}.left(1, :, :)), 2), 'r-.', 'linewidth',2.5), grid on, hold on
plot(mean(squeeze(qr{2}.right(1, :, :)), 2),'g--', 'linewidth',2.5), hold off
xlabel('session')
xlim([0.5 8.5])
legend('$105^\circ$','$165^\circ$', 'interpreter', 'latex')

figure(5)
[coeff,score,latent] = pca(squeeze(rp{2}.Ab(1, 2, :, :)));
% 
% % experiment 3
% subplot(1,3,3,'Fontsize',9)
% hold all
% plot(mean(Exp{2}.Ab),'color',[0 0 .75],'linestyle','--','linewidth',2.5)
% plot(mean(Exp{3}.At),'color',[.75 0 0],'linewidth',2.5)
% set(gca, 'XTick', 1:8)
% set(gca, 'YScale', 'log')
% xlim([0.5 8.5])
% ylim([1.5 8.5])
% xlabel('session')
% title('Experiment 3 (BA)')
% legend('A_B','A_T')
% legend('boxoff')