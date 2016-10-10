function JX = compute_jacobian_over_points(P, X, opts, store, point_num)
%COMPUTE_JACOBIAN_OVER_POINTS Summary of this function goes here
%   Detailed explanation goes here

if nargin < 5,
  % Check whether we want the Jacobian over the entire points or only one.
  rng_j = 1 : store.dim.num_points;
  JX = cell(store.dim.num_points, 1);
else
  rng_j = point_num(1);
end

for j = rng_j
  % s = p13i * x_j
  if strcmpi(opts.camera_model, 'affine')
    % If affine, compute the corresponding derivatives.
    s = P(store.index.visible_camera_r3{j}, 4) * X(4, j);
    p3 = [sparse(store.dim.num_visible_frames_per_point(j), 3), ...
      P(store.index.visible_camera_r3{j}, 4)];
  elseif opts.mu == 1
    % Otherwise, assume the projective model and compute the derivatives.
    p3 = P(store.index.visible_camera_r3{j}, :);
    s = p3 * X(:, j);
  else
    % Otherwise, assume the projective model and compute the derivatives.
    p3 = (repmat([opts.mu, opts.mu, opts.mu, 1], ...
      store.dim.num_visible_frames_per_point(j), 1) .* ...
      P(store.index.visible_camera_r3{j}, :));
    s = p3 * X(:, j);
  end
  
  % P12i * x_j / s_ij
  P12xj = P(store.index.visible_camera_r12{j}, :) * X(:, j) ./ repelem(s, 2, 1);
  
  if nargin < 5,
    % If the point number is not given, compute the entire Jacobian.
    JX{j} = (P(store.index.visible_camera_r12{j}, :) ...
      - repmat(P12xj, 1, store.dim.num_point_params) ...
      .* repelem(p3, 2, 1)) ...
      ./ repelem(s, 2, store.dim.num_point_params);
  else
    % Otherwise, compute only the Jacobian with respect to point j.
    JX = (P(store.index.visible_camera_r12{j}, :) ...
      - repmat(P12xj, 1, store.dim.num_point_params) ...
      .* repelem(p3, 2, 1)) ...
      ./ repelem(s, 2, store.dim.num_point_params);
  end
end

end
