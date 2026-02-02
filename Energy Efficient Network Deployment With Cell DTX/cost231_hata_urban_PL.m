function PL = cost231_hata_urban_PL(fc_mhz, d_km, hb_m, hm_m, metropolitan)
%COST231_HATA_URBAN_PL  COST-231 Hata model for urban (1500~2000 MHz).
%
% Inputs:
%   metropolitan: true면 C=3 dB(대도시), 아니면 C=0 dB(중소도시)
%
% PL[dB] = 46.3 + 33.9 log10(fc) - 13.82 log10(hb) - a(hm)
%          + (44.9 - 6.55 log10(hb)) log10(d) + C

    fc = fc_mhz;
    d  = max(d_km, 1e-6);

    % a(hm): urban correction (COST231 commonly uses this form)
    a_hm = (1.1*log10(fc) - 0.7)*hm_m - (1.56*log10(fc) - 0.8);

    C = 0;
    if metropolitan
        C = 3;
    end

    PL = 46.3 + 33.9*log10(fc) - 13.82*log10(hb_m) ...
       - a_hm + (44.9 - 6.55*log10(hb_m)).*log10(d) + C;
end
