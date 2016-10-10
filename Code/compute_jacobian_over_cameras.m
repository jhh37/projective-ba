function JP = compute_jacobian_over_cameras(P, X, opts, store, camera_num)
%COMPUTE_JACOBIAN_OVER_CAMERAS Summary of this function goes here
%   Detailed explanation goes here

if nargin < 5,
  % Check whether we want the Jacobian over the entire points or only one.
  rng_i = 1 : store.dim.num_frames;
  JP = cell(store.dim.num_frames, 1);
else
  rng_i = camera_num(1);
end

for i = rng_i
  % s = p13i * x_j
  if strcmpi(opts.camera_model, 'affine')
    % If affine, compute the corresponding derivatives.
    s = (P(3 * i, 4) * X(4, store.index.visible_points{i}))';
    x3 = [sparse(store.dim.num_visible_points_per_frame(i), 3), ...
      X(4, store.index.visible_points{i})'];
  elseif opts.mu == 1
    % Otherwise, assume the projective model and compute the derivatives.
    x3 = X(:, store.index.visible_points{i})';
    s = x3 * P(3 * i, :)';
  else
    x3 = repmat([opts.mu, opts.mu, opts.mu, 1], ...
      store.dim.num_visible_points_per_frame(i), 1) .* ...
      X(:, store.index.visible_points{i})';
    s = x3 * P(3 * i, :)';  
  end
  
  % P12i * x_j / s_ij
  P12xj = reshape(P([3 * i - 2, 3 * i - 1], :) * ...
    X(:, store.index.visible_points{i}), ...
    2 * store.dim.num_visible_points_per_frame(i), 1) ./ repelem(s, 2, 1);
  
  if nargin < 5,
    % If the camera number is not given, compute the entire Jacobian.
    JP{i} = [repmat(repelem(speye(2), 1, 4), ...
      store.dim.num_visible_points_per_frame(i), 1) ...
      .* repmat(repelem(X(:,store.index.visible_points{i})', 2, 1), 1, 2), ...
      - repmat(P12xj, 1, 4) .* repelem(x3, 2, 1)] ...
      ./ repelem(s, 2, store.dim.num_camera_params);
  else
    % Otherwise, compute only the Jacobian with respect to camera i.
    JP = [repmat(repelem(speye(2), 1, 4), ...
      store.dim.num_visible_points_per_frame(i), 1) ...
      .* repmat(repelem(X(:,store.index.visible_points{i})', 2, 1), 1, 2), ...
      - repmat(P12xj, 1, 4) .* repelem(x3, 2, 1)] ...
      ./ repelem(s, 2, store.dim.num_camera_params);
  end
end

end
