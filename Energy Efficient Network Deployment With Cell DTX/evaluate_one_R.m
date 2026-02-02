function out = evaluate_one_R(R_m, delta, p)
%EVALUATE_ONE_R  Evaluate objective and QoS statistics for a given cell range R.
%
% Output:
%   out.daily_avg_Parea_W_per_km2 : daily avg area power [W/km^2] using alpha_24 (operation)
%   out.busy_cchi_bps             : busy-hour c_{chi%} [bps] under WORST-CASE alpha=1 (QoS sizing)
%   out.busy_eta_mean             : busy-hour mean load under OPERATION alpha=alpha_24(busy_t)

    % 19-site layout
    bs_pos = generate_19site_hex_positions(R_m);

    % per-RB TX power from coverage
    P_RB_W = compute_PRB_from_Pmin(R_m, p);

    % network area (19 cells)
    A_km2 = area_km2(R_m, p.M);

    % ---------- (1) Daily avg energy under operation alpha(t) ----------
    Parea_t = zeros(24,1);

    % store operation busy-hour load (for Fig.2)
    eta_bh_drop_op = nan(p.nDrops,1);

    for t = 1:24
        Parea_drop = zeros(p.nDrops,1);

        for dd = 1:p.nDrops
            seed = p.rng_seed + 100000*round(R_m) + 100*t + dd;
            rng(seed);

            % OPERATION: active users follow alpha_24(t)
            N_act = floor(p.rho_users_km2 * p.alpha_24(t) * A_km2);
            N_act = min(N_act, p.users_cap);

            if N_act <= 0
                eta = zeros(p.M,1);
            else
                users = sample_users_uniform_union_hex(N_act, R_m, bs_pos, p);
                g_ik  = link_gain_linear(bs_pos, users, p);

                rx = g_ik .* (P_RB_W');
                [~, b_i] = max(rx, [], 2);

                [eta, c_user] = solve_feasible_load_PRB(g_ik, P_RB_W, b_i, p);
            end

            % area power for this hour (operation)
            Parea_drop(dd) = compute_Parea_t(eta, P_RB_W, delta, R_m, p);

            % operation busy-hour load statistic (Fig.2)
            if t == p.busy_t
                eta_bh_drop_op(dd) = mean(eta);  % avg over 19 cells
            end
        end

        Parea_t(t) = mean(Parea_drop);
    end

    out.daily_avg_Parea_W_per_km2 = mean(Parea_t);
    out.busy_eta_mean = mean(eta_bh_drop_op(~isnan(eta_bh_drop_op)));

    % ---------- (2) Busy-hour QoS under worst-case alpha = 1 ----------
    % Only used for QoS constraint c_{chi%} >= rmin (Rhat, Algorithm 1)
    c_bh_all = [];

    t = p.busy_t;
    for dd = 1:p.nDrops
        seed = p.rng_seed + 100000*round(R_m) + 100*t + dd + 9999; % offset to decouple from op drop
        rng(seed);

        % WORST-CASE: assume alpha = 1 at busy hour for QoS sizing
        N_act = floor(p.rho_users_km2 * 1.0 * A_km2);
        N_act = min(N_act, p.users_cap);

        if N_act <= 0
            c_user = [];
        else
            users = sample_users_uniform_union_hex(N_act, R_m, bs_pos, p);
            g_ik  = link_gain_linear(bs_pos, users, p);

            rx = g_ik .* (P_RB_W');
            [~, b_i] = max(rx, [], 2);

            [~, c_user] = solve_feasible_load_PRB(g_ik, P_RB_W, b_i, p);
        end

        c_bh_all = [c_bh_all; c_user]; %#ok<AGROW>
    end

    if isempty(c_bh_all)
        out.busy_cchi_bps = NaN;
    else
        out.busy_cchi_bps = prctile(c_bh_all, p.chi_percent);
    end
end
