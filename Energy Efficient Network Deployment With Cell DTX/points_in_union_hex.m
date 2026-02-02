function in = points_in_union_hex(P, centers, R, p)
% 코드가 잘 이해가 안됨
%POINTS_IN_UNION_HEX  Check if points belong to union o
% f regular hexagons.
%
%   in = points_in_union_hex(P, centers, R, p)
%
% Inputs:
%   P       : (N x 2) points
%   centers : (M x 2) hex centers
%   R       : hex radius
%   p       : parameter struct (uses p.hex_in_tol)
%
% Output:
%   in      : (N x 1) boolean (true if inside any hex)

    N = size(P,1); % 행렬의 행 개수
    in = false(N,1); % 행 개수 만큼 벡터 생성

    for k = 1:size(centers,1) %center의 행 개수만큼
        idx = ~in; %새로운 리스트
        if ~any(idx), break; end % 남아 있는 점이 없으면 break

        Pin = is_point_in_hex(P(idx,:) - centers(k,:), R, p);
        tmp = in;
        tmp(idx) = Pin;
        in = tmp;
    end
end