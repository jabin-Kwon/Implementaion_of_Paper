function g_ik = link_gain_linear(bs_pos, users, p)
%LINK_GAIN_LINEAR  Compute linear link gains g_ik from BSs to users.
%
%   g_ik = link_gain_linear(bs_pos, users, p)
%
% Inputs:
%   bs_pos : (M x 2) BS coordinates [m]
%   users  : (N x 2) UE coordinates [m]
%   p      : parameter struct
%
% Output:
%   g_ik   : (N x M) linear link gains (including antenna gain)
%
% Steps:
%   1) compute distances UE-BS
%   2) compute path loss using Okumura-Hata urban model
%   3) apply antenna gain and convert PL(dB) to linear gain

    Nuser = size(users,1);
    M = size(bs_pos,1);

    d_m = zeros(Nuser, M); % 행: user, 열: 기지국
    for k = 1:M
        dx = users(:,1) - bs_pos(k,1);
        dy = users(:,2) - bs_pos(k,2);
        d_m(:,k) = sqrt(dx.^2 + dy.^2);
    end

    d_km = max(d_m/1000, 1e-6);

    PL = hata_urban_PL(p.fc_mhz, d_km, p.hb_m, p.hm_m, p.large_city);
    Gtot_dB = p.Gtx_dBi + 0; % UE gain assumed 0 dBi
    g_ik = 10.^(-(PL - Gtot_dB)/10);
end