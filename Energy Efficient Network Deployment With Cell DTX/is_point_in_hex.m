function inside = is_point_in_hex(P, R, p)
%IS_POINT_IN_HEX  Robust point-in-convex-polygon test for a regular hexagon.
%
%   inside = is_point_in_hex(P, R, p)
%
% Inputs:
%   P : (N x 2) points in coordinates where hex center is at origin
%   R : hex radius (center-to-vertex)
%   p : parameter struct (uses p.hex_in_tol)
%
% Output:
%   inside : (N x 1) boolean (inside or on boundary)
%
% Method:
%   - Use CCW vertices from hex_vertices()
%   - A point is inside a convex CCW polygon if it is always on the left side
%     of each directed edge (cross product >= 0 with tolerance)
%
% IMPORTANT FIX #2:
%   Use a *reasonable* tolerance. Too-small tolerance can cause boundary jitter.

    V = hex_vertices(R);     % (6 x 2) CCW vertices
    Vn = V([2:6 1],:);
    E  = Vn - V;

    eps0 = p.hex_in_tol * max(R,1);  % e.g., 1e-9 * R

    inside = true(size(P,1),1);
    for i = 1:6
        c = E(i,1).*(P(:,2)-V(i,2)) - E(i,2).*(P(:,1)-V(i,1));
        inside = inside & (c >= -eps0);
        if ~any(inside), break; end
    end
end