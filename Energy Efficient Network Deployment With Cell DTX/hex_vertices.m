function V = hex_vertices(R)
%HEX_VERTICES  Regular hexagon vertices (pointy-top) in CCW order.
%
%   V = hex_vertices(R)
%
% Input:
%   R : center-to-vertex distance
%
% Output:
%   V : (6 x 2) vertices in CCW order

    ang = (0:5)' * (pi/3) + pi/6;  % +30° => pointy-top
    V = [R*cos(ang), R*sin(ang)];
end