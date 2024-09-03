clc
clear all
close all

addpath("functions/");
para = para_init();

theta = 45*pi/180; % user direction
r = 10; % user distance

para.N_T = 16; % number of TTDs
para.M = 256; % number of subcarriers

figure; 

%% Bandwidth B = 10 GHz
B = 1e10; % bandwidth
m = 1:para.M;
para.fm_all =  para.fc + B*(2*m-1-para.M) / (2*para.M); % subcarrier frequencies

% calculate array gain achieved by different methods
[P_prop, P_prop_robust, P_conv_CF, P_con_MCCM, P_conv_MCM] = beampattern(para, theta, r);
subplot(3,1,1); hold on; box on;
plot(para.fm_all/1e9, 10*log10(P_prop/para.N), '-', 'LineWidth', 1.5);
plot(para.fm_all/1e9, 10*log10(P_prop_robust/para.N), '-.', 'LineWidth', 1.5);
plot(para.fm_all/1e9, 10*log10(P_con_MCCM/para.N), '--', 'LineWidth', 1.5);
plot(para.fm_all/1e9, 10*log10(P_conv_MCM/para.N), '--', 'LineWidth', 1.5);
plot(para.fm_all/1e9, 10*log10(P_conv_CF/para.N), ':', 'LineWidth', 1.5);
set(groot,'defaultAxesTickLabelInterpreter','latex');
legend("TTD-BF, Proposed method","TTD-BF, Proposed robust method",...
    "Conventional BF, MCCM", "Conventional BF, MCM", "Conventional BF, CF", 'Interpreter', 'Latex');
xlabel('Frequency (GHz)', 'Interpreter', 'Latex');
title('$B = 10$ GHz', 'Interpreter', 'Latex');

%% Bandwidth B = 20 GHz
B = 2e10; % bandwidth
m = 1:para.M;
para.fm_all =  para.fc + B*(2*m-1-para.M) / (2*para.M); % subcarrier frequencies

% calculate array gain achieved by different methods
[P_prop, P_prop_robust, P_conv_CF, P_con_MCCM, P_conv_MCM] = beampattern(para, theta, r);
subplot(3,1,2); hold on; box on;
plot(para.fm_all/1e9, 10*log10(P_prop/para.N), '-', 'LineWidth', 1.5);
plot(para.fm_all/1e9, 10*log10(P_prop_robust/para.N), '-.', 'LineWidth', 1.5);
plot(para.fm_all/1e9, 10*log10(P_con_MCCM/para.N), '--', 'LineWidth', 1.5);
plot(para.fm_all/1e9, 10*log10(P_conv_MCM/para.N), '--', 'LineWidth', 1.5);
plot(para.fm_all/1e9, 10*log10(P_conv_CF/para.N), ':', 'LineWidth', 1.5);
xlabel('Frequency (GHz)', 'Interpreter', 'Latex');
title('$B = 20$ GHz', 'Interpreter', 'Latex');

%% Bandwidth B = 30 GHz
B = 3e10; % bandwidth
m = 1:para.M;
para.fm_all =  para.fc + B*(2*m-1-para.M) / (2*para.M); % subcarrier frequencies

% calculate array gain achieved by different methods
[P_prop, P_prop_robust, P_conv_CF, P_con_MCCM, P_conv_MCM] = beampattern(para, theta, r);
subplot(3,1,3); hold on; box on;
plot(para.fm_all/1e9, 10*log10(P_prop/para.N), '-', 'LineWidth', 1.5);
plot(para.fm_all/1e9, 10*log10(P_prop_robust/para.N), '-.', 'LineWidth', 1.5);
plot(para.fm_all/1e9, 10*log10(P_con_MCCM/para.N), '--', 'LineWidth', 1.5);
plot(para.fm_all/1e9, 10*log10(P_conv_MCM/para.N), '--', 'LineWidth', 1.5);
plot(para.fm_all/1e9, 10*log10(P_conv_CF/para.N), ':', 'LineWidth', 1.5);
xlabel('Frequency (GHz)', 'Interpreter', 'Latex'); 
title('$B = 30$ GHz', 'Interpreter', 'Latex');

