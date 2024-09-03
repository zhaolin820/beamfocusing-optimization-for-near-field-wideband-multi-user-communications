function [H] = generate_channel(para, user_r, user_theta)
%Generate the BS-user, RIS-user and BS-RIS channels 
%  [h, Hc, Hr, dc, dr, LOS, NLOS] = generate_channel(para, angle, path_loss, ar, ac)
%Inputs:
%   para: structure of the initial parameters
%   angle: struture of the angles
%   path_loss: structure of the path loss
%   ar: steering vector of the radar antennas
%   ac: steering vector of the communication antennas
%Outputs:
%   h: RIS-user channels
%   H: BS-RIS channel
%Date: 01/04/2022
%Author: Zhaolin Wang

HITRANparams = importdata('data_freq_abscoe.txt');

L = 4;
r_NLoS = rand(para.K, L) * 10 + 5; % 5 ~ 15 m
theta_NLoS = rand(para.K, L) * pi; % 0 ~ 180 degree


H = zeros(para.N, para.K, para.M);
for k = 1:para.K

    reflection_factor = sqrt(10^(-15/10))*(randn(L,1) + 1i*randn(L,1));
    for m = 1:para.M
        fm = para.fm_all(m);
        path_loss = getSpreadLoss(fm, user_r(k)) + getAbsLoss(fm, user_r(k), HITRANparams );
        path_loss = 10.^((-path_loss - para.noise_dB + para.Gt + para.Gr)/10);
        % LoS channel
        H(:,k,m) = sqrt(path_loss) * array_response_vector(user_r(k), user_theta(k), para.N, para.d, fm);
        
        % NLoS channel
        for l = 1:L
            H(:,k,m) = H(:,k,m) + sqrt(path_loss) * reflection_factor(l) * array_response_vector(r_NLoS(k,l), theta_NLoS(k,l), para.N, para.d, fm);
        end
    end
end
H = conj(H);

end