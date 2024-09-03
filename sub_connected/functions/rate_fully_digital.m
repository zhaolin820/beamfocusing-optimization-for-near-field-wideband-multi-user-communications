function [R_sum, R] = rate_fully_digital(para, W, h)
%Calculate sum rate
%  [R_sum, R] = rate_fully_digital(para, P, h)
%Inputs:
%   para: structure of the initial parameters
%   W: beamformer
%   h: channel for all users
%Outputs:
%   R_sum: spectral efficiency
%   R: rate of each user at each subcarrier
%Date: 04/04/2024
%Author: Zhaolin Wang


R = zeros(para.K, para.M);
for m = 1:para.M
    Wm = W(:,:,m);
    for k = 1:para.K
        hmk = h(:,k,m);
        wmk = Wm(:,k); 
        P_I = Wm; P_I(:,k) = [];
        Imk = norm(hmk'*P_I)^2 + norm(Wm, 'fro')^2/para.Pt; 
        R(k,m) = log2( 1 + abs(hmk'*wmk)^2/Imk );
    end
end
R_sum = sum(sum(R));


end

