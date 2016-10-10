function C = compute_camera_centres(P)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

f = size(P, 1) / 3;
C = nan(4, f);
for i = 1 : f
  C(:, i) = null(P([3 * i - 2, 3 * i - 1, 3 * i], :));
  C(:, i) = C(:, i) / C(4, i);
end

end

