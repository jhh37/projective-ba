function err = compute_residual_framewise(P, X, opts, store, camera_num)
%COMPUTE_RESIDUAL_FRAMEWISE Summary of this function goes here
%   Detailed explanation goes here

if nargin < 5,
  % Check whether we want the Jacobian over the entire points or only one.
  rng_i = 1 : store.dim.num_frames;
  err = nan(store.dim.nnz_frames * 2, 1);
else
  rng_i = camera_num(1);
end

for i = rng_i
  % s = p13i * x_j
  if strcmpi(opts.camera_model, 'affine')
    s = (P(3 * i, 4) * X(4, store.index.visible_points{i}))';
  elseif opts.mu == 1
    s = (P(3 * i, :) * X(:, store.index.visible_points{i}))';
  else
    s = ([opts.mu, opts.mu, opts.mu, 1] .* P(3 * i, :) * ...
      X(:, store.index.visible_points{i}))';
  end
  
  % P12i * x_j
  if nargin < 5
    % If the point number is not given, compute the entire error vector.
    err(store.rng.proj_per_frame{i}) = ...
      reshape(P([3 * i - 2, 3 * i - 1], :) ...
      * X(:, store.index.visible_points{i}), ...
      2 * store.dim.num_visible_points_per_frame(i), 1) ./ repelem(s, 2, 1);
  else
    % Otherwise, compute the error vector for point j.
    err = ...
      reshape(P([3 * i - 2, 3 * i - 1], :) ...
      * X(:, store.index.visible_points{i}), ...
      2 * store.dim.num_visible_points_per_frame(i), 1) ./ repelem(s, 2, 1);
  end
end

if nargin < 5
  err = err - store.vec.mt_framewise;
else
  err = err - store.vec.mt(store.rng.proj_per_frame{i});
end

end
