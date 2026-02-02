function P_RB = compute_PRB_from_Pmin(R_m, p)
%COMPUTE_PRB_FROM_PMIN  Per-RB TX power from coverage constraint.
%
% Paper: choose minimum Pk to satisfy g_edge * Pk = Pmin.
% Here Pk is per-RB power and Pmin is interpreted per RB (consistent).
%
% Output:
%   P_RB : (M x 1) per-RB TX power [W]

    % received threshold per reference BW
    Pmin_ref_W = dbm2w(p.Pmin_dBm);

    % convert to per-RB received threshold
    Pmin_perRB_W = Pmin_ref_W * (p.W_RB / p.Pmin_ref_BW_hz);

    % edge distance in km
    d_km = max(R_m/1000, 1e-6);

    % path loss at edge
    PL = hata_urban_PL(p.fc_mhz, d_km, p.hb_m, p.hm_m, p.large_city);

    % antenna gain
    Gtot_dB = p.Gtx_dBi + 0;
    g_edge = 10.^(-(PL - Gtot_dB)/10);

    % invert
    P_RB_scalar = Pmin_perRB_W / max(g_edge, 1e-30);
    P_RB = ones(p.M,1) * P_RB_scalar;
end
