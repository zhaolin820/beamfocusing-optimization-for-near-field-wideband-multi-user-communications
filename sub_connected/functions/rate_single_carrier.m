function [R_sum, R] = rate_single_carrier(para, W, H)
%Calculate sum rate of a single subcarrier
%  [R_sum, R] = rate_fully_digital(para, P, h)
%Inputs:
%   para: structure of the initial parameters
%   W: beamformer
%   h: channel for all users
%Outputs:
%   R_sum: spectral efficiency
%   R: rate of each user
%Date: 04/04/2024
%Author: Zhaolin Wang

R = zeros(para.K, 1);
for k = 1:para.K
    hk = H(:,k);
    wk = W(:,k); 
    P_I = W; P_I(:,k) = [];
    Ik = norm(hk'*P_I)^2 + norm(W, 'fro')^2/para.Pt; 
    R(k) = log2( 1 + abs(hk'*wk)^2/Ik );
end
R_sum = sum(R);



end

