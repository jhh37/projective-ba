function JP = compute_inhom_jacobian_over_cameras(X, opts, store, camera_num)
%COMPUTE_INHOM_JACOBIAN_OVER_CAMERAS Summary of this function goes here
%   Detailed explanation goes here

if nargin < 4,
  % Check whether we want the Jacobian over the entire points or only one.
  rng_i = 1 : store.dim.num_frames;
  JP = cell(store.dim.num_frames, 1);
else
  rng_i = camera_num(1);
end

for i = rng_i
  if strcmpi(opts.camera_model, 'affine')
    % Compute the Jacobian for inhomogeneous affine camera.
    if nargin < 4,
      JP{i} = repmat(repelem(speye(2), 1, 4), ...
      store.dim.num_visible_points_per_frame(i), 1) ...
      .* repmat(repelem(X(:,store.index.visible_points{i})', 2, 1), 1, 2);
    else
      JP = repmat(repelem(speye(2), 1, 4), ...
      store.dim.num_visible_points_per_frame(i), 1) ...
      .* repmat(repelem(X(:,store.index.visible_points{i})', 2, 1), 1, 2);
    end
  end

end
