function err = compute_residual(P, X, opts, store, point_num)
%COMPUTE_RESIDUAL Summary of this function goes here
%   Detailed explanation goes here

if nargin < 5,
  % Check whether we want the Jacobian over the entire points or only one.
  rng_j = 1 : store.dim.num_points;
  err = nan(store.dim.nnz_frames * 2, 1);
else
  rng_j = point_num(1);
end

for j = rng_j
  % s = p13i * x_j
  if strcmpi(opts.camera_model, 'affine')
    s = P(store.index.visible_camera_r3{j}, 4) * X(4, j);
  elseif opts.mu == 1
    s = P(store.index.visible_camera_r3{j}, :) * X(:, j);
  else
    s = (repmat([opts.mu, opts.mu, opts.mu, 1], ...
      store.dim.num_visible_frames_per_point(j), 1) .* ...
      P(store.index.visible_camera_r3{j}, :)) * X(:, j);
  end
  
  % P12i * x_j
  if nargin < 5
    % If the point number is not given, compute the entire error vector.
    err(store.rng.proj_per_point{j}) = P(store.index.visible_camera_r12{j}, :) * X(:, j) ./ repelem(s, 2, 1);
  else
    % Otherwise, compute the error vector for point j.
    err = P(store.index.visible_camera_r12{j}, :) * X(:, j) ./ repelem(s, 2, 1);
  end
end

if nargin < 5
  err = err - store.vec.mt;
else
  err = err - store.vec.mt(store.rng.proj_per_point{j});
end

end
