clc
clear all
close all

addpath("functions/");
set(groot,'defaultAxesTickLabelInterpreter','latex');

%% parameters
para = para_init();
user_r = rand(para.K, 1) * 10 + 5; % user distances 5 ~ 15 m
user_theta = sort(rand(para.K, 1) * pi); % user directions 0 ~ 180 degree

%% generate channel matrix
[H] = generate_channel(para, user_r, user_theta);

%% Fully-digital approximation (FDA) approach

[R_convergence, penalty_convergence, A, D] = algorithm_FDA_penalty(para, H, user_r, user_theta);

% plot convergence results
figure; 
subplot(1,2,1); 
plot(R_convergence, '-b', 'LineWidth', 1.5);
xlabel('Number of outer-loop iterations', 'Interpreter', 'Latex');
ylabel('Spectral efficiency (bit/s/Hz)', 'Interpreter', 'Latex');
box on; grid on;

subplot(1,2,2);
semilogy(penalty_convergence, '-b', 'LineWidth', 1.5);
xlabel('Number of outer-loop iterations', 'Interpreter', 'Latex');
ylabel('Penalty value', 'Interpreter', 'Latex');
box on; grid on;

R_FDA = R_convergence(end);

%% Heuristic two-stage (HTS) approach
[R_HTS_PNF] = algorithm_HTS_PNF(para, H, user_r, user_theta);
[R_HTS_robust] = algorithm_HTS_robust(para, H, user_r, user_theta);
