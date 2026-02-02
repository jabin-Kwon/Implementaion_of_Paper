function Rstar = algorithm1_find_Rstar(delta, p)
%ALGORITHM1_FIND_RSTAR  Paper Algorithm 1 (early-stopping search for R*).
%
%   Rstar = algorithm1_find_Rstar(delta, p)
%
% Inputs:
%   delta : cell DTX performance parameter δ
%   p     : parameter struct
%
% Output:
%   Rstar : optimum cell range R* found by early stopping rule
%
% Logic (paper-like):
%   Sweep R from R_min to R_max:
%     - compute objective E_t[Parea(R)] (daily average area power)
%     - compute QoS metric at busy hour: c_{χ%}(R)
%   Continue while:
%     - objective strictly decreases AND QoS constraint satisfied
%   Stop and return last passing R.

    rng(p.rng_seed);

    R_list = p.R_min:p.R_step:p.R_max;
    last_pass_R = R_list(1);
    last_pass_P = inf;

    for i = 1:numel(R_list)
        R = R_list(i);
        out = evaluate_one_R(R, delta, p);

        P = out.daily_avg_Parea_W_per_km2;
        cchi = out.busy_cchi_bps;

        if i == 1
            if cchi >= p.rmin_bps
                last_pass_R = R;
                last_pass_P = P;
            else
                Rstar = R; % if QoS fails already at R_min
                return;
            end
        else
            tolP = 1e-6 * max(1, last_pass_P);
            if (P < last_pass_P - tolP) && (cchi >= p.rmin_bps)
                last_pass_R = R;
                last_pass_P = P;
            else
                Rstar = last_pass_R;
                return;
            end
        end
    end

    Rstar = last_pass_R;
end