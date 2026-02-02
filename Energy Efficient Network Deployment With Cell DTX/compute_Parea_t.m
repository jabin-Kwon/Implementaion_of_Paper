function Parea = compute_Parea_t(eta, P_RB_W, delta, R_m, p)
%COMPUTE_PAREA_T  Area power (paper Eq.(7),(11)).
%
% Ek = ζ Pk ηk + (1-δ) P0 ηk + δ P0
% Pk here is per-RB power (as defined in paper text).
% Parea = sum_k Ek / A(R),  A(R)=M*(3*sqrt(3)/2)*R^2

    Ek = p.zeta .* P_RB_W .* eta ...
       + (1 - delta) * p.P0_W .* eta ...
       + delta * p.P0_W;

    A_km2 = area_km2(R_m, p.M);
    Parea = sum(Ek) / A_km2;
end
