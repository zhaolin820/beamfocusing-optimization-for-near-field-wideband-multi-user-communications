function [R, A, D, t] = algorithm_HTS_robust(para, H, user_r, user_theta, W_initial)
%The robust heuristic two-stage approach
%  [R, A, D, t] = algorithm_HTS_robust(para, H, user_r, user_theta, W_initial)
%Inputs:
%   para: structure of the initial parameters
%   H: channel for all users
%   user_r: distance for all users
%   user_theta: angle for all users
%   W_initial: initialized digital beamformers (for FDA approach)
%Outputs:
%   R: optimized spectral efficiency
%   A: optimized analog beamforming matrix
%   D: optimized digital beamforming matrix
%   t: optimized time delays of TTDs
%Date: 22/07/2024
%Author: Zhaolin Wang

%% Initialization
switch nargin
    case 4
    W_initial = randn(para.N, para.K) + 1i * randn(para.N, para.K);
    W_initial = W_initial / norm(W_initial, 'fro') * sqrt(para.Pt);
end

c = 3e8; % speed of light
t_search = 0:para.t_max/1e3:para.t_max; % search space of TTDs' time delay

N_sub = para.N/para.N_T; % number of antennas connected to each TTD
e = ones(N_sub, 1);

%% Analog beamformer design
t = zeros(para.N_T, para.N_RF); % time delay of TTDs
A_PS = zeros(para.N, para.N_RF); % PS based analog beamformer

% design the analog beamformer for each RF chain
% the n-th beamformer is designed for the n-th user
for n = 1:para.N_RF
    theta = user_theta(n); r = user_r(n);
    
    % calculate the exact array response at all subcarriers
    array_response = zeros(para.N, para.M);
    for m = 1:para.M
        fm = para.fm_all(m);
        bm = array_response_vector(r, theta, para.N, para.d, fm);
        array_response(:,m) = conj(bm);
    end
    
    % initialization using the piecewise-near-field approximation
    r_n = zeros(para.N_T, 1);
    t_n = zeros(para.N_T, 1);
    a_n = zeros(para.N, 1);
    for l = 1:para.N_T
        xi_l = (l-1-(para.N_T-1)/2)*N_sub;
        r_l = sqrt(r^2 + xi_l^2*para.d^2 - 2*r*xi_l*para.d*cos(theta)); % Equation (50)
        theta_l = acos( (r*cos(theta) - xi_l*para.d)/r_l ); % Equation (49)
        r_n(l) = r_l;
        t_n(l) = - (r_l - r)/c; % delay difference between the center of the subarray and the center of the entire array

        % PS coefficients
        q = (0:(N_sub-1))';
        delta_q = (q-(N_sub-1)/2) * para.d;
        a_n((l-1)*N_sub+1 : l*N_sub) = exp( 1i * 2 * pi * para.fc/c...
            * (sqrt(r_l^2 + delta_q.^2 - 2*r_l*delta_q*cos(theta_l)) - r_l) ); % Equation (51)
    end 
    t_n = t_n - min(t_n);
    t_n(t_n>para.t_max) = para.t_max;
          
    % iterative optimization using the proposed robust design
    obj_value_max_pre = 0;
    iter_max = 40; % maximum iteration number
    for i = 1:iter_max
        % update PS coefficients
        a_n_pre = a_n; % the PS coefficients obtained in the previous iteration
        q_n = 0;
        for m = 1:para.M
            fm = para.fm_all(m); % frequency of the m-th subcarrier
            bm = array_response(:,m); % array response vector of the m-th subcarrier
            eta_m = bm .* exp(1i*2*pi*fm*kron(t_n, e)); % Equation (66)
            q_n = q_n + eta_m*eta_m'*a_n_pre/abs(eta_m'*a_n_pre); % Equation (70)
        end
        a_n = q_n./abs(q_n); % PS-based analog beamformer, Equation (69)
    
    
        % update TTD coefficients
        gamma = zeros(para.N_T, para.M);
        for l = 1:para.N_T
            phi_n_l = a_n(((l-1)*N_sub+1):l*N_sub); % PS-based analog beamformer for the l-th sub-array
            for m = 1:para.M
                gamma(l,m) = array_response( ((l-1)*N_sub+1):l*N_sub ,m)' * phi_n_l;
            end
        end
     

        for l = 1:para.N_T
            t_null = t_n; t_null(l) = []; % remove the l-th entry
            
            obj_value = 0;
            for m = 1:para.M
                gamma_m = gamma(:,m);
                gamma_m_null = gamma_m; gamma_m_null(l) = [];
                fm = para.fm_all(m);
                fixed_term = sum(gamma_m_null.*exp(-1i*2*pi*fm*t_null));
                search_term = gamma_m(l)*exp(-1i*2*pi*fm*t_search );
                obj_value = obj_value + abs(fixed_term + search_term); % Objective value in Equation (63)
            end
            [~,I] = max(obj_value); % one-dimensional search
            t_n(l) = t_search(I);
        end
        
        % check convergence of the PS and TTD optimization
        obj_value_max = obj_value(I);
        if abs((obj_value_max-obj_value_max_pre)/obj_value_max) < 1e-4
            break;
        end
        obj_value_max_pre = obj_value_max;
    end

    t(:,n) = t_n; % Time delay for TTDs connected to the n-th RF chain
    A_PS(:,n) = a_n; % PS-based analog beamformer connected to the n-th RF chain
end

% calculate the overall analog beamformer and the equivalent channel
A = zeros(para.N, para.N_RF, para.M);
H_equal = zeros(para.N_RF, para.K, para.M);
for m = 1:para.M
    A(:,:,m) = analog_bamformer(para, A_PS, t, para.fm_all(m)); % overall analog beamformer
    H_equal(:,:,m) = A(:,:,m)'*H(:,:,m); % equivalent channel
end

%% Digital beamformer design
D = zeros(para.N_RF, para.K, para.M);
for m = 1:para.M
    Dm = pinv(A(:,:,m))*W_initial;
    [~, Dm] = RWMMSE(para, H(:,:,m), H_equal(:,:,m), Dm, A(:,:,m));
    D(:,:,m) = Dm;
end

%% Calculate the spectral efficiency
W = zeros(para.N, para.K, para.M);
for m = 1:para.M
    W(:,:,m) = A(:,:,m)*D(:,:,m);
end

[R] = rate_fully_digital(para, W, H);
R = R/(para.M+para.Lcp);
end


%% Calculate the overall analog beamformer at frequency f
function [A] = analog_bamformer(para, A_PS, t, f)
    e = ones(para.N/para.N_T,1);
    T = exp(-1i*2*pi*f*t);
    A = A_PS .* kron(T, e);
end


%% RWMMSE method for optimizing the digital beamformer
function [R, D] = RWMMSE(para, H, H_equal, D, A)
R_pre = 0;
for i = 1:20
    E = eye(para.K);
    Phi = 0; Upsilon = 0;
    for k = 1:para.K
        hk = H_equal(:,k);
        dk = D(:,k); 
        I = norm(hk'*D)^2 + norm(A*D, 'fro')^2/para.Pt; 
        w_k = 1 + abs(hk'*dk)^2 / (I - abs(hk'*dk)^2);
        v_k = hk'*dk / I;
    
        Phi = Phi + w_k*abs(v_k)^2 * ( hk*hk' + eye(para.N_RF)/para.Pt );
        Upsilon = Upsilon + w_k*conj(v_k)*E(:,k)*hk';
    
    end
    
    D = Phi\Upsilon';

    % check convergence
    [R] = rate_single_carrier(para, A*D, H);
    if abs(R - R_pre)/R <= 1e-4
        break;
    end
    R_pre = R;
end

end

