function P = compute_camera_from_two_views(m_1, m_2)
%COMPUTE_CAMERA_FROM_TWO_VIEWS Summary of this function goes here
%   Detailed explanation goes here

% Custom function
ssm = @(x) [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];

% Compute the Fundamental matrix from two correspondences.
common_points = ~isnan(m_1(1, :) .* m_2(1, :));

x_1 = [m_1(:, common_points); ones(1, size(m_1(1, common_points), 2))];
x_2 = [m_2(:, common_points); ones(1, size(m_2(1, common_points), 2))];

[F, ~, e_2] = fundmatrix(x_1, x_2);
F = F / norm(F, 'fro');
e_2 = e_2 / norm(e_2);

P = [ssm(e_2) * F, e_2];
P = P / norm(P, 'fro');

end

