function A_km2 = area_km2(R_m, M)
%AREA_KM2  Network area using the paper's hex-cell approximation.
%
%   A_km2 = area_km2(R_m, M)
%
% Inputs:
%   R_m : cell range (m), interpreted as center-to-vertex distance
%   M   : number of cells (sites)
%
% Output:
%   A_km2 : total area [km^2]
%
% Paper area formula:
%   A = M * (3*sqrt(3)/2) * R^2

    A_m2 = M * (3*sqrt(3)/2) * R_m^2;
    A_km2 = A_m2 / 1e6;
end