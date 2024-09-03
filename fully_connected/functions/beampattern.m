function [P_PNF,P_robust,P_conv_CF, P_con_MCCM, P_conv_MCM] = beampattern(para, theta, r)
%Calculate the beam pattern achieved by different analog beamforming method
%  [P_PNF,P_robust,P_con, P_con_mean_cov, P_con_mean_vec] = beampattern(para, theta, r)
%Inputs:
%   para: structure of the initial parameters
%   r: distance
%   theta: angle
%Outputs:
%   P_PNF: beampattern achieved by the PNF method
%   P_robust: beampattern achieved by the robust method
%   P_conv_CF: beampattern achieved by the conventional CF method
%   P_con_MCCM: beampattern achieved by the MCCM method
%   P_conv_MCM: beampattern achieved by the MCM method
%Date: 22/07/2024

c = 3e8; % speed of light
N_sub = para.N/para.N_T; % the number of antennas connected to each TTD
t_search = 0:para.t_max/1e3:para.t_max; % search space of TTDs' time delay
e = ones(N_sub, 1);

%% Proposed PNF method
r_n = zeros(para.N_T, 1);
t_PNF = zeros(para.N_T, 1);
a_PNF = zeros(para.N, para.N_T);
for l = 1:para.N_T
    xi_l = (l-1-(para.N_T-1)/2)*N_sub;
    r_l = sqrt(r^2 + xi_l^2*para.d^2 - 2*r*xi_l*para.d*cos(theta)); % Equation (50)
    theta_l = acos( (r*cos(theta) - xi_l*para.d)/r_l ); % Equation (49)
    r_n(l) = r_l;
    t_PNF(l) = - (r_l - r)/c; % delay difference between the center of the subarray and the center of the entire array

    % Calculate PS coefficients
    q = (0:(N_sub-1))';
    delta_q = (q-(N_sub-1)/2) * para.d;
    a_PNF((l-1)*N_sub+1 : l*N_sub, l) = exp( 1i * 2 * pi * para.fc/c...
        * (sqrt(r_l^2 + delta_q.^2 - 2*r_l*delta_q*cos(theta_l)) - r_l) ); % Equation (51)
end

t_PNF = t_PNF - min(t_PNF);
t_PNF(t_PNF>para.t_max) = para.t_max;

% Optimize TTD coefficients
obj_value_max_pre = 0;
for step = 1:40
    for l = 1:para.N_T
        t_n_null = t_PNF; t_n_null(l) = []; % remove the l-th entry
        r_n_null = r_n; r_n_null(l) = []; % remove the l-th entry
        obj_value = 0;
        for m = 1:para.M
            fm = para.fm_all(m);
            fixed_term = sum(exp(-1i*2*pi*fm*( (r_n_null - r)/c +  t_n_null  )));
            search_term = exp(-1i*2*pi*fm* ( (r_n(l) - r )/c + t_search) );
            obj_value = obj_value + abs(fixed_term + search_term);
        end
        [~,I] = max(obj_value); % one-dimensional search
        t_PNF(l) = t_search(I);
    end
    obj_value_max = obj_value(I);
    if abs((obj_value_max-obj_value_max_pre)/obj_value_max) < 1e-4
        break;
    end
    obj_value_max_pre = obj_value_max;
end

% Calculate beam pattern
P_PNF = zeros(para.M, 1);
for m = 1:para.M
    fm = para.fm_all(m);
    a_PNF_m = a_PNF*exp(-1i*2*pi*fm*t_PNF);

    bm = array_response_vector(r, theta, para.N, para.d, fm);
    P_PNF(m) = abs(bm.' * a_PNF_m);
end

%% Proposed robust method
array_response = zeros(para.N, para.M);
for m = 1:para.M
    fm = para.fm_all(m);
    bm = array_response_vector(r, theta, para.N, para.d, fm);
    array_response(:,m) = conj(bm);
end

% initialization using the PNF approximation
r_robust = zeros(para.N_T, 1);
t_robust = zeros(para.N_T, 1);
a_robust = zeros(para.N, 1);
for l = 1:para.N_T
    xi_l = (l-1-(para.N_T-1)/2)*N_sub;
    r_l = sqrt(r^2 + xi_l^2*para.d^2 - 2*r*xi_l*para.d*cos(theta)); % Equation (50)
    theta_l = acos( (r*cos(theta) - xi_l*para.d)/r_l ); % Equation (49)
    r_robust(l) = r_l;
    t_robust(l) = - (r_l - r)/c; % delay difference between the center of the subarray and the center of the entire array

    % PS coefficients
    q = (0:(N_sub-1))';
    delta_q = (q-(N_sub-1)/2) * para.d;
    a_robust((l-1)*N_sub+1 : l*N_sub) = exp( 1i * 2 * pi * para.fc/c...
        * (sqrt(r_l^2 + delta_q.^2 - 2*r_l*delta_q*cos(theta_l)) - r_l) ); % Equation (51)
end
t_robust = t_robust - min(t_robust);
t_robust(t_robust>para.t_max) = para.t_max;

% iterative optimization using the proposed robust design
obj_value_max_pre = 0;
for i = 1:40
    % update PS coefficients
    a_robust_pre = a_robust; % the PS coefficients obtained in the previous iteration
    q = 0;
    for m = 1:para.M
        fm = para.fm_all(m); % frequency of the m-th subcarrier
        bm = array_response(:,m); % array response vector of the m-th subcarrier
        eta_m = bm .* exp(1i*2*pi*fm*kron(t_robust, e)); % Equation (66)
        q = q + eta_m*eta_m'*a_robust_pre/abs(eta_m'*a_robust_pre); % Equation (70)
    end
    a_robust = q./abs(q);


    % update TTD coefficients
    gamma = zeros(para.N_T, para.M);
    for l = 1:para.N_T
        phi_l = a_robust(((l-1)*N_sub+1):l*N_sub); % PS-based analog beamformer for the l-th sub-array
        for m = 1:para.M
            gamma(l,m) = array_response( ((l-1)*N_sub+1):l*N_sub ,m)' * phi_l;
        end
    end
 
    for l = 1:para.N_T
        t_null = t_robust; t_null(l) = []; % remove the l-th entry
        obj_value = 0;
        for m = 1:para.M
            gamma_m = gamma(:,m);
            gamma_m_null = gamma_m; gamma_m_null(l) = [];
            fm = para.fm_all(m);
            fixed_term = sum(gamma_m_null.*exp(-1i*2*pi*fm*t_null));
            search_term = gamma_m(l)*exp(-1i*2*pi*fm*t_search );
            obj_value = obj_value + abs(fixed_term + search_term);
        end
        [~,I] = max(obj_value);
        t_robust(l) = t_search(I);
    end
    
    % check convergence of the PS and TTD optimization
    obj_value_max = obj_value(I);
    if abs((obj_value_max-obj_value_max_pre)/obj_value_max) < 1e-4
        break;
    end
    obj_value_max_pre = obj_value_max;
end

% Calculate beam pattern
P_robust = zeros(para.M, 1);
for m = 1:para.M
    fm = para.fm_all(m);
    bm = array_response_vector(r, theta, para.N, para.d, fm);
    am = a_robust.* exp(-1i*2*pi*fm*kron(t_robust, e));
    P_robust(m) = abs(bm.'*am);
end

%% Conventional method - Mean channel covariance matrix

% Calculate the mean of channel covariance matrix
cov_mean = zeros(para.N, para.N, para.M);
for m = 1:para.M
    fm = para.fm_all(m);
    bm = array_response_vector(r, theta, para.N, para.d, fm);
    cov_mean(:,:,m) = (conj(bm)*bm.');
end
cov_mean = sum(cov_mean,3)/para.M;
[V,D] = eig(cov_mean); 
D = diag(D);

% Calculate analog beamformer
a = V(:,1)*sqrt(D(1));
a = a ./ abs(a);

% Calculate beam pattern
P_con_MCCM = zeros(para.M, 1);
for m = 1:para.M
    fm = para.fm_all(m);
    bm = array_response_vector(r, theta, para.N, para.d, fm);
    P_con_MCCM(m) = abs(bm.' * a);
end


%% Conventional method - Mean channel matrix

% Calculate vector mean
vec_mean = 0;
for m = 1:para.M
    fm = para.fm_all(m);
    bm = array_response_vector(r, theta, para.N, para.d, fm);
    vec_mean = vec_mean + conj(bm);
end
vec_mean = vec_mean ./ abs(vec_mean);

% Calculate analog beamformer
a = vec_mean;

% Calculate beam pattern
P_conv_MCM = zeros(para.M, 1);
for m = 1:para.M
    fm = para.fm_all(m);
    bm = array_response_vector(r, theta, para.N, para.d, fm);
    P_conv_MCM(m) = abs(bm.' * a);
end


%% Conventional method - Central frequency
a = conj(array_response_vector(r, theta, para.N, para.d, para.fc));

P_conv_CF = zeros(para.M, 1);
for m = 1:para.M
    fm = para.fm_all(m);
    bm = array_response_vector(r, theta, para.N, para.d, fm);
    P_conv_CF(m) = abs(bm.' * a);
end


end

