function JX = compute_inhom_jacobian_over_points(P, opts, store, point_num)
%COMPUTE_INHOM_JACOBIAN_OVER_POINTS Summary of this function goes here
%   Detailed explanation goes here

if nargin < 4,
  % Check whether we want the Jacobian over the entire points or only one.
  rng_j = 1 : store.dim.num_points;
  JX = cell(store.dim.num_points, 1);
else
  rng_j = point_num(1);
end

for j = rng_j
  if strcmpi(opts.camera_model, 'affine')
    % Compute the Jacobian for inhomogeneous affine camera.
    if nargin < 4,
      JX{j} = P(store.index.visible_camera_r12{j}, 1:3);
    else
      JX = P(store.index.visible_camera_r12{j}, 1:3);
    end
  end

end
