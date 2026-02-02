function PL = hata_urban_PL(fc_mhz, d_km, hb_m, hm_m, large_city)
%HATA_URBAN_PL  Okumura-Hata path loss model (urban).
%
%   PL = hata_urban_PL(fc_mhz, d_km, hb_m, hm_m, large_city)
%
% Inputs:
%   fc_mhz     : carrier frequency [MHz]
%   d_km       : distance [km] (can be vector/matrix)d
%   hb_m       : BS height [m]
%   hm_m       : UE height [m]
%   large_city : boolean (large city correction)
%
% Output:
%   PL : path loss [dB]

    if large_city
        if fc_mhz <= 300
            a_hm = 8.29*(log10(1.54*hm_m))^2 - 1.1;
        else
            a_hm = 3.2*(log10(11.75*hm_m))^2 - 4.97;
        end
    else
        a_hm = (1.1*log10(fc_mhz) - 0.7)*hm_m - (1.56*log10(fc_mhz) - 0.8);
    end

    PL = 69.55 + 26.16*log10(fc_mhz) - 13.82*log10(hb_m) ...
       - a_hm + (44.9 - 6.55*log10(hb_m)).*log10(d_km);
end