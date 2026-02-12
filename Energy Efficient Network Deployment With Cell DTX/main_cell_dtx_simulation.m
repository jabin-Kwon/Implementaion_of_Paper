%% ============================================================
%  Energy Efficient Network Deployment With Cell DTX
%  MATLAB implementation aligned with paper Eq.(1)~(13) + Algorithm 1
%
%  Outputs:
%   - Fig.1: Daily avg area power vs R (Case1 δ=0.1 vs Case2 δ=1) + R* + Rhat
%   - Fig.2: R*(δ) and E[η^bh(R*)] vs δ
%   - Fig.3: Hourly Parea(t) 3-case bar chart (as described in paper)
%
%  Key modeling points (paper-consistent):
%   - Pk: transmit power per RB (power spectral density per RB)
%   - SINR (Eq.1): load-dependent interference (η-scales interference)
%   - Rate (Eq.2): r_i = W_RB * min(log2(1+γ), νmax)  [bit/s per RB]
%   - User total rate: c_i = N_RB * r_i
%   - Load (Eq.6): η_k = sum_{i∈βk} demand / c_i
%   - Area power: [W/km^2]
%
%  IMPORTANT: Several parameters are not explicitly specified in the paper.
%  Those are marked as [ASSUMPTION].
%
%  This version includes robustness/consistency fixes:
%   (1) Pmin reference bandwidth handling (Pmin can be defined on BW_hz or W_RB)
%   (2) Geometry "shaking" removal: robust convex-polygon point-in-hex test + sane tolerance
%   (3) Clarified load solver: function name / comments match actual inputs (per-RB power)
%% ============================================================

clear; clc; close all;

%% =========================
% 1) Parameters (paper + assumptions)
% =========================
p = struct();

% ---- Paper simulation settings ----
p.M = 19;                  % 19 sites (2-ring hex layout)
p.fc_mhz = 2000;           % 2 GHz
p.BW_hz  = 10e6;           % 10 MHz
p.Gtx_dBi = 15;            % BS antenna gain 15 dBi
p.NF_dB = 8;               % UE noise figure 8 dB
p.rho_users_km2 = 800;    % user density (users/km^2)
p.Omega_MB_per_hour = 18;  % Ω = 18 MB in 1 hour (traffic demand per user)
p.nu_max = 6;              % νmax = 6 bps/Hz
p.zeta = 4.7;              % power model coefficient
p.P0_W = 130;              % fixed power offset per BS (W)

p.rmin_bps = 8e6;          % rmin = 8 Mbps

% ---- Cell range sweep ----
p.R_min = 100;
p.R_max = 1000;
p.R_step = 10;             % [ASSUMPTION] not specified in the paper

% ---- QoS percentile χ ----
p.chi_percent = 20;         % [ASSUMPTION] paper does not explicitly state χ

% ---- RB settings (paper defines RB concept but not numeric mapping) ----
p.N_RB = 50;        % RB 개수는 그대로
p.W_RB = 180e3;     % RB bandwidth를 180 kHz로 고정
p.Pmin_dBm = -70;          % paper: Pmin = -70 dBm
p.Pmin_ref_BW_hz = p.W_RB;   % -70 dBm per-RB (RB 기준)

% ---- Okumura-Hata parameters ----
p.hb_m = 20;               % [ASSUMPTION] BS height (m)
p.hm_m = 1.5;              % [ASSUMPTION] UE height (m)
p.large_city = false;       % [ASSUMPTION] large city urban model

% ---- User placement / Monte-Carlo ----
p.rng_seed = 7;
p.nDrops = 1;              % [ASSUMPTION] paper does not specify MC repetitions
p.users_cap = 2000000;       % [ASSUMPTION] memory safety cap

% ---- Load-coupling fixed-point solver settings ----
p.eta_init = 0.2;          % [ASSUMPTION]
p.eta_eps = 1e-6;          % [ASSUMPTION]
p.eta_iter_max = 200;      % [ASSUMPTION]

% Reference active user percentage based on downlink traffic analysis
pct = [13 10 7 4 3 2 2 3 4 6 8 9 10 10 11 12 12 12 13 14 15 16 16 15]; % 0~23h
p.alpha_24 = pct/100;      % active fraction (peak=0.16)
[~, p.busy_t] = max(p.alpha_24);

[~, p.busy_t] = max(p.alpha_24);  % busy hour index

% ---- IMPORTANT FIX #1: Pmin reference bandwidth ----
% The paper often lists only "Pmin = -70 dBm" without explicit bandwidth.
% Two common interpretations:
%  (A) Pmin defined over the *system bandwidth* (10 MHz)
%  (B) Pmin defined over the *RB bandwidth* (W_RB)
%
% Our model requires per-RB Pk. We convert as:
%   Pmin_perRB = Pmin_ref * (W_RB / Pmin_ref_BW)
%
% Choose ONE:
% p.Pmin_ref_BW_hz = p.BW_hz;  % RECOMMENDED DEFAULT: assume -70 dBm over 10 MHz

% ---- Geometry robustness tolerance ----
p.hex_in_tol = 1e-9;       % tolerance scale for point-in-hex test (relative to R)

% ---- Cases ----
delta_case1 = 0.1;
delta_case2 = 1.0;

%% ===== Rhat only (QoS max feasible R) =====
R_list_tmp = p.R_min:p.R_step:p.R_max;
cchi_tmp = nan(size(R_list_tmp));

for i = 1:numel(R_list_tmp)
    R = R_list_tmp(i);
    out = evaluate_one_R(R, delta_case1, p); % QoS는 내부에서 alpha=1로 계산됨
    cchi_tmp(i) = out.busy_cchi_bps;
end

feasible_mask = (cchi_tmp >= p.rmin_bps);
if any(feasible_mask)
    Rhat_only = max(R_list_tmp(feasible_mask));
else
    Rhat_only = NaN;
end

fprintf('Rhat (QoS-only, alpha=1 in sizing) = %.1f m\n', Rhat_only);
%% ============================================================
% 2) Fig.1: sweep R → curves + R* + Rhat
%% ============================================================
rng(p.rng_seed);

R_list = p.R_min:p.R_step:p.R_max;

Parea_case1 = nan(size(R_list));
Parea_case2 = nan(size(R_list));
cchi_case1  = nan(size(R_list));
cchi_case2  = nan(size(R_list));

for i = 1:numel(R_list)
    R = R_list(i);

    out1 = evaluate_one_R(R, delta_case1, p);
    out2 = evaluate_one_R(R, delta_case2, p);

    Parea_case1(i) = out1.daily_avg_Parea_W_per_km2;
    Parea_case2(i) = out2.daily_avg_Parea_W_per_km2;
    cchi_case1(i)  = out1.busy_cchi_bps;
    cchi_case2(i)  = out2.busy_cchi_bps;
end

% Rhat: maximum R satisfying QoS constraint (busy-hour cχ >= rmin)
feasible_mask = (cchi_case1 >= p.rmin_bps); % using δ=0.1 feasibility (paper does similar)
if any(feasible_mask)
    Rhat = max(R_list(feasible_mask));
else
    Rhat = NaN;
end

% Algorithm 1 to find R*
Rstar_case1 = algorithm1_find_Rstar(delta_case1, p);
Rstar_case2 = algorithm1_find_Rstar(delta_case2, p);

figure; grid on; hold on;
set(gca,'YScale','log');
plot(R_list, Parea_case1, 'b*');
plot(R_list, Parea_case2, 'r*');

if ~isnan(Rhat), xline(Rhat, 'g', 'LineWidth', 2); end
xline(Rstar_case1, 'b--', 'LineWidth', 2);
xline(Rstar_case2, 'r--', 'LineWidth', 2);

ylim([1e1 1e4]); % [ASSUMPTION] approximate range seen in paper figures
xlabel('Cell Range [meters]');
ylabel('Daily Average Area Power Consumption [W/km^2]');
legend('Case 1 (\delta=0.1)','Case 2 (\delta=1)', 'Max R satisfying QoS', 'Location','best');
title('Daily average area power consumption vs cell range (r_{min}=8 Mbps)');

%% ============================================================
% 3) Fig.2: sweep δ → R*(δ) and E[η^{bh}(R*)]
%% ============================================================
delta_list = 0:0.04:1.0;
Rstar_list = nan(size(delta_list));
eta_bh_list = nan(size(delta_list));
p2 = p;
p2.nDrops = 1;  

for i = 1:numel(delta_list)
    d = delta_list(i);

    Rstar = algorithm1_find_Rstar(d, p2);
    Rstar_list(i) = Rstar;

    out = evaluate_one_R(Rstar, d, p2);
    eta_bh_list(i) = out.busy_eta_mean;
end

figure; grid on; hold on;
plot(delta_list, Rstar_list, 'b*-');
xlabel('\delta'); ylabel('R^* [meters]');
title('Impact of cell DTX performance on optimum cell range and busy-hour load');

% Paper-style: show E[η^{bh}] as text at selected δ points
markD = [0.1 0.5 1.0];
for md = markD
    [~, ii] = min(abs(delta_list - md)); % 매트랩 관용 표현
    if ~isnan(Rstar_list(ii))
        text(delta_list(ii), Rstar_list(ii), sprintf('  E[\\eta^{bh}]=%.0f%%', 100*eta_bh_list(ii)), ...
            'VerticalAlignment','bottom');
    end
end
yline(Rstar_list(end), 'r--'); % δ=1 baseline (visual cue)

%% ============================================================
% 4) Fig.3: hourly Parea(t) for δ=0.1 (3 cases)
%% ============================================================
Rplan_case2 = Rstar_case2; % planned with δ=1
Rplan_case1 = Rstar_case1; % planned with δ=0.1

Parea_blue  = evaluate_hourly_profile(Rplan_case2, 1.0, p);  % operation δ=1
Parea_green = evaluate_hourly_profile(Rplan_case2, 0.1, p);  % operation δ=0.1 (same plan)
Parea_red   = evaluate_hourly_profile(Rplan_case1, 0.1, p);  % plan+op δ=0.1

figure; grid on;
bar(1:24, [Parea_blue(:), Parea_green(:), Parea_red(:)]);
xlabel('T (Hours)');
ylabel('Area Power Consumption (W/km^2)');
legend('R* for Case 2, without Cell DTX', ...
       'R* for Case 2, with Cell DTX', ...
       'R* for Case 1, with Cell DTX', ...
       'Location','northwest');

title('Daily area power consumption variation for \delta=0.1');
