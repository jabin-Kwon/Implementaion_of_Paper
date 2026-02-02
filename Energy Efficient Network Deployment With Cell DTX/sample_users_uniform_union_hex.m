function users = sample_users_uniform_union_hex(N, R_m, bs_pos, p)
%SAMPLE_USERS_UNIFORM_UNION_HEX  Uniform sampling over union of hex-cell polygons.
%
%   users = sample_users_uniform_union_hex(N, R_m, bs_pos, p)
%
% Inputs:
%   N      : number of users
%   R_m    : hex radius (center-to-vertex) in meters
%   bs_pos : (M x 2) cell centers
%   p      : parameter struct (uses p.hex_in_tol)
%
% Output:
%   users  : (N x 2) user coordinates [m]
%
% Method:
%   - Build a bounding box that covers all hex vertices in the 19-cell layout
%   - Draw candidate points uniformly in the bounding box
%   - Accept those that fall inside ANY of the hexagons
%   - Repeat until N points are accepted

    V0 = hex_vertices(R_m); % (6 x 2) base hex at origin

    allV = [];
    for k = 1:size(bs_pos,1)
        allV = [allV; V0 + bs_pos(k,:)]; %#ok<AGROW>
    end

    xmin = min(allV(:,1)); xmax = max(allV(:,1));
    ymin = min(allV(:,2)); ymax = max(allV(:,2));

    users = zeros(N,2);
    cnt = 0;

    batch = max(5000, 5*N);   % [ASSUMPTION]
    safety_iter = 4000;       % [ASSUMPTION]
    iter = 0;

    while cnt < N
        iter = iter + 1;
        if iter > safety_iter
            error('Rejection sampling too slow. Check geometry or increase batch.');
        end

        x = xmin + (xmax - xmin) * rand(batch,1);
        y = ymin + (ymax - ymin) * rand(batch,1);
        cand = [x y];

        in = points_in_union_hex(cand, bs_pos, R_m, p);

        acc = cand(in,:);
        nacc = size(acc,1);

        if nacc > 0
            take = min(nacc, N - cnt);
            users(cnt+1:cnt+take,:) = acc(1:take,:);
            cnt = cnt + take;
        end
    end
end