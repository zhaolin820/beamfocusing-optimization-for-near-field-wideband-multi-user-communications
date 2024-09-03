function [para] = para_init()
%Construct a struct of the initial values for all the parameters 
%  [para] = para_init()
%Inputs:
%   None
%Outputs:
%   values: a struct
%Date: 22/07/2024
%Author: Zhaolin Wang

para.fc = 1e11; % carrier frequency (Hz)
para.B = 1e10; % system bandwidth (Hz)
c = 3e8; % speed of light (m/s)
para.lambda = c/para.fc; % carrier wavelength (m)
para.d = para.lambda/2; % antenna spacing (m)

para.Pt = 10^(20/10); % overall transmit power (dBm)
para.K = 4; % user number

para.N = 512; % number of antennas
para.N_RF = para.K; % number of RF chains
para.N_sub = para.N/para.N_RF; % the number of antennas connected to each RF chain
para.N_T = 16; % number of TTDs connected to each RF chain
para.t_max = para.N/(2*para.fc); % maximum delay of TTDs (s)

para.D = (para.N-1)*para.d; % aperture of antenna array (m)
para.Rayl = 2*para.D^2/para.lambda; % Rayleigh distance (m)

para.M = 10; % number of subcarriers
para.Lcp = 4; % length of cyclic prefix

m = 1:para.M;
para.fm_all =  para.fc + para.B*(2*m-1-para.M) / (2*para.M); % subcarrier frequency (Hz)

para.noise_dB = -174; % noise power density in dBm/Hz
para.noise_dB = para.noise_dB + 10*log10(para.B/10); % noise power in dBm

para.Gt = 15; % dBi, transmit antenna gain
para.Gr = 5; % dBi, receive antenna gain


end

