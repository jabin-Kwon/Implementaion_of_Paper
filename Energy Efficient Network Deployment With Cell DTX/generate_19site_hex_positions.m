function pos = generate_19site_hex_positions(R_m)
%GENERATE_19SITE_HEX_POSITIONS  Generate 19-site (2-ring) hexagonal layout.
%
%   pos = generate_19site_hex_positions(R_m)
%
% Input:
%   R_m : cell range (m), center-to-vertex distance
%
% Output:
%   pos : (19 x 2) BS coordinates [m]
%
% Geometry:
%   Use axial coordinates for a pointy-top hex grid.
%   Neighbor center distance (ISD) becomes sqrt(3)*R_m for pointy-top geometry.

    radius = 2; % 2-ring → 1 + 6 + 12 = 19 cells
    axial = [];

    for q = -radius:radius
        r1 = max(-radius, -q - radius);
        r2 = min(radius, -q + radius);
        for r = r1:r2
            axial = [axial; q r]; %#ok<AGROW>
        end
    end

    q = axial(:,1);
    r = axial(:,2);

    % axial -> xy (pointy-top)
    x = R_m * sqrt(3) * (q + r/2);
    y = R_m * (3/2) * r;

    pos = [x y];
end