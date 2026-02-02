function Parea_t = evaluate_hourly_profile(R_m, delta_op, p)
%EVALUATE_HOURLY_PROFILE  Hourly area power profile for Fig.3.

    bs_pos = generate_19site_hex_positions(R_m);

    P_RB_W = compute_PRB_from_Pmin(R_m, p);
    A_km2  = area_km2(R_m, p.M);

    Parea_t = zeros(24,1);

    for t = 1:24
        Parea_drop = zeros(p.nDrops,1);

        for dd = 1:p.nDrops
            N_act = floor(p.rho_users_km2 * p.alpha_24(t) * A_km2);
            N_act = min(N_act, p.users_cap);

            if N_act <= 0
                eta = zeros(p.M,1);
            else
                users = sample_users_uniform_union_hex(N_act, R_m, bs_pos, p);
                g_ik  = link_gain_linear(bs_pos, users, p);

                rx = g_ik .* (P_RB_W');
                [~, b_i] = max(rx, [], 2);

                [eta, ~] = solve_feasible_load_PRB(g_ik, P_RB_W, b_i, p);
            end

            Parea_drop(dd) = compute_Parea_t(eta, P_RB_W, delta_op, R_m, p);
        end

        Parea_t(t) = mean(Parea_drop);
    end
end
