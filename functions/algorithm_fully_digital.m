function [R, P] = algorithm_fully_digital(para, h, P_initial)
%WMMSE algorithm for fully digital systems
%  [R_convergence] = algorithm_TTD_parallel(para, h, P, user_r, user_theta)
%Inputs:
%   para: structure of the initial parameters
%   h: channel for all users
%   P_initial: initial beamforming vectors
%Outputs:
%   R: achievable rates
%   P: optimal fully digital beamformers
%Date: 04/04/2024
%Author: Zhaolin Wang

P = zeros(para.N, para.K, para.M);
for m = 1:para.M
    [Rm, Pm] = RWMMSE(para, h(:,:,m), P_initial);
    % disp(['Fully digital, subcarrier - ' num2str(m) ', rate - ' num2str(Rm)]);
    P(:,:,m) = Pm;
end
[R] = rate_fully_digital(para, P, h);
R = R/(para.M+para.Lcp);

end

function [R, P] = RWMMSE(para, h, P)
R_pre = 0;
for i = 1:20
    [P] = update_P(para, h, P);

    % check convergence
    [R] = rate_single(para, P, h);
    if abs(R - R_pre)/R <= 1e-4
        break;
    end
    R_pre = R;
end

end

%% closed-form WMMSE updates
function [P] = update_P(para, h, P)

E = eye(para.K);
Phi = 0; Upsilon = 0;
for k = 1:para.K
    hk = h(:,k);
    pk = P(:,k); 
    I = norm(hk'*P)^2 + norm(P, 'fro')^2/para.Pt; 
    w_k = 1 + abs(hk'*pk)^2 / (I - abs(hk'*pk)^2);
    v_k = hk'*pk / I;

    Phi = Phi + w_k*abs(v_k)^2 * ( hk*hk' + eye(para.N)/para.Pt );
    Upsilon = Upsilon + w_k*conj(v_k)*E(:,k)*hk';

end

P = Phi'\Upsilon';

end

%% achievable rate
function [R_sum, R] = rate_single(para, P, h)

R = zeros(para.K, 1);
for k = 1:para.K
    hk = h(:,k);
    pk = P(:,k); 
    P_I = P; P_I(:,k) = [];
    Ik = norm(hk'*P_I)^2 + norm(P, 'fro')^2/para.Pt; 
    R(k) = log2( 1 + abs(hk'*pk)^2/Ik );
end
R_sum = sum(R);



end


