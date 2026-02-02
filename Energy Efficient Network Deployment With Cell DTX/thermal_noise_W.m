function Nw = thermal_noise_W(B_hz)
%THERMAL_NOISE_W  Thermal noise power kTB (W).
%
%   Nw = thermal_noise_W(B_hz)
%
% Input:
%   B_hz : bandwidth [Hz]
%
% Output:
%   Nw   : noise power [W]

    k = 1.38064852e-23;
    T = 290;
    Nw = k*T*B_hz;
end