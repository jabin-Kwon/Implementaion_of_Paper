function [eta, c_user] = solve_feasible_load_PRB(g_ik, P_RB_W, b_i, p)
%SOLVE_FEASIBLE_LOAD_PRB  Fixed-point solver for feasible load (paper Eq.(6)).
%
% Inputs:
%   g_ik   : (Nuser x M) linear link gains
%   P_RB_W : (M x 1) per-RB transmit power (Pk)
%   b_i    : (Nuser x 1) serving BS index per user
%   p      : parameters
%
% Outputs:
%   eta    : (M x 1) cell load in [0,1]
%   c_user : (Nuser x 1) user rate c_i = N_RB * r_i  [bps]

    [Nuser, M] = size(g_ik);

    % init
    eta = ones(M,1) * p.eta_init;

    % demand per user in bps (Ω MB/hour -> bps)
    Omega_bits = p.Omega_MB_per_hour * 8e6;
    demand_bps = Omega_bits / 3600;

    % noise power per RB: kTB * NF
    sigma2 = thermal_noise_W(p.W_RB) * db2lin(p.NF_dB);

    c_user = zeros(Nuser,1);

    for it = 1:p.eta_iter_max
        eta_prev = eta;

        % Use eta as "interference probability" in [0,1]
        eta_eff = min(max(eta,0), 1);

        % total received load-scaled (signal+interference) term
        I_all = sum(g_ik .* (P_RB_W'.* eta_eff'), 2);  % (Nuser x 1)

        idx  = (1:Nuser)';
        serv = b_i;

        % desired signal (paper Eq.(1) numerator has no eta)
        S = g_ik(sub2ind([Nuser,M], idx, serv)) .* P_RB_W(serv);

        % remove serving BS from interference sum (k != serving)
        I_serv = g_ik(sub2ind([Nuser,M], idx, serv)) .* P_RB_W(serv) .* eta_eff(serv);

        denom = (I_all - I_serv) + sigma2;
        denom = max(denom, 1e-30);

        gamma = S ./ denom;

        % rate per RB (Eq.(2)) and total rate c_i
        se = min(log2(1 + gamma), p.nu_max);
        r_RB = p.W_RB * se;          % bps per RB
        c_user = p.N_RB * r_RB;      % bps total

        % load update (Eq.(6) style)
        eta_new = zeros(M,1);
        for k = 1:M
            uidx = (b_i == k);
            if any(uidx)
                eta_new(k) = sum(demand_bps ./ max(c_user(uidx), 1e-30));
            end
        end

        % enforce definition eta in [0,1]
        eta = min(max(eta_new, 0), 1);

        if max(abs(eta - eta_prev)) < p.eta_eps
            break;
        end
    end
end
