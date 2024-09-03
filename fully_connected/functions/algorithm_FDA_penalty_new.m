function [R_convergence, penalty_convergence, A, D] = algorithm_FDA_penalty_new(para, H, user_r, user_theta)
%The new version of the panalty-based fully-digital approximation (FDA) method based on the method in the following paper
%Compared to the method in this paper, the new version is more stable and efficient.
% [1] Z. Wang, et al., "TTD Configurations for Near-Field Beamforming: Parallel, Serial, or Hybrid?," 
%     IEEE Transactions on Communications, vol. 72, no. 6, pp. 3783-3799, June 2024, doi: 10.1109/TCOMM.2024.3370832.
%
%   [R_convergence, penalty_convergence, A, D] = algorithm_FDA_penalty(para, H, W, user_r, user_theta, P_initial)
%Inputs:
%   para: structure of the initial parameters
%   H: channel for all users
%   user_r: distance for all users (for initialization)
%   user_theta: angle for all users (for initialization)
%Outputs:
%   R_convergence: achievable rates at each iteration
%   penalty_convergence: penalty value at each iteration
%   A: optimized analog beamforming matrix
%   D: optimized digital beamforming matrix
%Date: 22/07/2024

R_convergence = [];
penalty_convergence = [];

% Initialization using the fully-digital solution and the PNF-based HTS approach
W_initial = randn(para.N, para.K) + 1i * randn(para.N, para.K);
W_initial = W_initial / norm(W_initial, 'fro') * sqrt(para.Pt);
[~, W] = algorithm_fully_digital(para, H, W_initial);
[~, A, D, t] = algorithm_HTS_PNF(para, H, user_r, user_theta, W_initial);

% Optimization
penalty_factor = 1e5; % initial value of the penalty factor
iter_max = 40; % maximum iteration number
for outer_step = 1:iter_max
    obj_pre = 0;
    for inner_step = 1:iter_max
        [W] = update_fully_digital(para, H, W, A, D, penalty_factor);
        [A, t] = update_analog_beamformer(para, W, D, t);
        [D] = update_digital_beamformer(para, W, A, D);
        
        % calculate objective value
        [R_sum_FD] = rate_fully_digital(para, W, H);
        [penalty_value, ~] = penalty_value_calculator(para, W, A, D);
        obj = R_sum_FD - 1/penalty_factor*penalty_value;

        % calculate the rate achieved by the real hybrid beamformers
        W_hybrid = zeros(para.N, para.K, para.M);
        for m = 1:para.M
            W_hybrid(:,:,m) = A(:,:,m)*D(:,:,m);
        end
        [R_sum] = rate_fully_digital(para, W_hybrid, H);
        
        % display the output of the inner loop
        disp(['Inner loop - ' num2str(inner_step) ', obj - ' num2str(obj) ...
            ', rate_FD - ' num2str(R_sum_FD/(para.M+para.Lcp)) ', rate_hybrid - ' num2str(R_sum/(para.M+para.Lcp))...
            ', penalty_value - ' num2str(penalty_value)]);
        % check the convergence of the inner loop
        if abs((obj-obj_pre)/obj) < 1e-4
            break;
        end
        obj_pre = obj;
    end
    
    % update the penalty factor
    penalty_factor = 0.5*penalty_factor;

    % display the output of the outer loop
    [~, penalty_value_max] = penalty_value_calculator(para, W, A, D);    
    disp(['Outer loop - ' num2str(outer_step)...
        ', penalty_value_max - ' num2str(penalty_value_max)...
        ', penalty_factor - ' num2str(penalty_factor)]);
    disp('-----------------------------------------------------------------------------------------------------');
    
    % check the convergence of the outer loop
    if penalty_value_max < 1e-4
        break;
    end

    R_convergence = [R_convergence, R_sum/(para.M+para.Lcp)];
    penalty_convergence = [penalty_convergence, penalty_value];
end

end

%% Update auxiliary fully-digital beamformer
function [W] = update_fully_digital(para, H, W, A, D, penalty_factor)
    E = eye(para.K);
    for m = 1:para.M
        hm = H(:,:,m); Dm = D(:,:,m); Am = A(:,:,m); Wm = W(:,:,m);

        Phi = 0; Upsilon = 0;
        for k = 1:para.K
            hk = hm(:,k);
            wk = Wm(:,k); 
            I = norm(hk'*Wm)^2 + norm(Wm, 'fro')^2/para.Pt; 
            mu_k = abs(hk'*wk)^2 / (I - abs(hk'*wk)^2); % Equation (35)
            lambda_k = sqrt(1+mu_k)*hk'*wk / I; % Equation (36)

            Phi = Phi + abs(lambda_k)^2 * ( hk*hk' + eye(para.N)/para.Pt );
            Upsilon = Upsilon + sqrt(1+mu_k)*conj(lambda_k)*E(:,k)*hk';
        end
        Upsilon = Upsilon + 1/penalty_factor * pinv(Dm)*Am';
        Psi = 1/penalty_factor * pinv(Dm)*pinv(Dm)';
        Wm = sylvester(Psi, Phi, Upsilon); % solve the Sylvester equation, Equation (69) in reference [1]
        W(:,:,m) = Wm';
    end
end

%% Update digital beamformer
function [D] = update_digital_beamformer(para, W, A, D)
    for m = 1:para.M
        D(:,:,m) = pinv(A(:,:,m))*W(:,:,m); % Equation (32)
    end
end

%% Update analog beamformer
function [A, t] = update_analog_beamformer(para, P, D, t)

    t_search = 0:para.t_max/1e3:para.t_max; % search space of TTDs' time delay

    PD = zeros(para.N, para.N_RF, para.M);
    for m = 1:para.M
        PD(:,:,m) = P(:,:,m)*pinv(D(:,:,m));
    end

    % update PS coefficients
    N_sub = para.N/para.N_T;
    A_PS = zeros(para.N, para.N_RF);
    for n = 1:para.N_RF
        for l = 1:para.N_T
            p_l = PD((l-1)*N_sub+1:l*N_sub, n, :);  
            p_l = squeeze(p_l);
            t_l = t(l,n);
            a_l = sum(p_l*diag(exp(1i*2*pi*para.fm_all*t_l)), 2);
            A_PS((l-1)*N_sub+1:l*N_sub, n) = a_l./abs(a_l);
        end
    end

    % update TTD coefficients
    t = zeros(para.N_T, para.N_RF);
    for n = 1:para.N_RF
        for l = 1:para.N_T
            p_l = PD((l-1)*N_sub+1:l*N_sub, n, :);  
            p_l = squeeze(p_l);
            a_l = A_PS((l-1)*N_sub+1:l*N_sub, n);
            psi_nq = p_l'*a_l;
            obj_value = real(psi_nq.'*exp(-1i*2*pi*para.fm_all'*t_search));
            [~,I] = max(obj_value); % one-dimensional search
            t(l,n) = t_search(I);
        end
    end

    % calculate overall analog beamformer
    A = zeros(para.N, para.N_RF, para.M);
    for m = 1:para.M
        A(:,:,m) = analog_bamformer(para, A_PS, t, para.fm_all(m));
    end
end



%% Calculate the overall TTD-based analog beamformer
function [Am] = analog_bamformer(para, A, t, fm)
    e = ones(para.N/para.N_T,1);
    Tm = exp(-1i*2*pi*fm*t);
    Am = A .* kron(Tm, e);
end


%% Calculate the penalty value
function [penalty_value, penalty_value_max] = penalty_value_calculator(para, W, A, D)
    penalty_value = 0; % Overall penalty value
    penalty_value_max = zeros(para.M,1); % Maximum entry of the penalty matrix
    for m = 1:para.M
        Wm = W(:,:,m); Dm = D(:,:,m); Am = A(:,:,m);
        penalty_value = penalty_value + norm(Wm*pinv(Dm) - Am, 'fro')^2;
        penalty_value_max(m) = norm(Wm - Am*Dm, 'inf');
    end

    penalty_value_max = max(penalty_value_max);
end





