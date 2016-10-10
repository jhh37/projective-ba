function P = compute_algebraic_cameras(X, opts, store)
%COMPUTE_ALGEBRAIC_POINTS Summary of this function goes here
%   Detailed explanation goes here
% for j = 1 : store.dim.num_points
%   A = store.blocks.mcp{j} * P(store.index.visible_camera_rows{j}, 1 : 3);
%   b = store.blocks.mcp{j} * X(4, j) * P(store.index.visible_camera_rows{j}, 4);
%   X(1 : 3, j) = - A \ b;
% end

P = nan(3 * store.dim.num_frames, 4);
for i = 1 : store.dim.num_frames
  
  % Compute sparse matrix values.
  vals = reshape(store.blocks.measurements_framewise{i}, 2, store.dim.num_visible_points_per_frame(i));
%   vals = [ repmat([1; -1], 1, store.dim.num_visible_points_per_frame(i)); ...
%     vals(2, :); - vals(1, :) ];
  vals_left = [X(:, store.index.visible_points{i})', ...
    - X(:, store.index.visible_points{i})'];
  vals_right = reshape([vals(2, :); -vals(1, :)], ...
    2 * store.dim.num_visible_points_per_frame(i), 1);
  vals_right = repelem(vals_right, 1, 4) .* ...
  repelem(X(:, store.index.visible_points{i})', 2, 1);
  % Compute row and column indexes.
  row_idx_left = bsxfun(@plus, ...
    [2 2 2 2 1 1 1 1], 2 * (0 : store.dim.num_visible_points_per_frame(i) - 1)');
  row_idx_right = repmat( ...
    (1 : 2 * store.dim.num_visible_points_per_frame(i))', 1, 4);
  col_idx_left = repmat([1 2 3 4 5 6 7 8], ...
    store.dim.num_visible_points_per_frame(i), 1);
  col_idx_right = repmat(9 : 12, ...
    2 * store.dim.num_visible_points_per_frame(i), 1);
  A = au_sparse(...
    int32([row_idx_left(:); row_idx_right(:)]), ...
    int32([col_idx_left(:); col_idx_right(:)]), ...
    [vals_left(:); vals_right(:)]); 
  
  [~, ~, V] = svd(full(A), 0);
  if opts.homogeneous,
    p = V(:, end) / norm(V(:, end));
  else
    p = V(:, end) / V(end);
  end
  P([3 * i - 2, 3 * i - 1, 3 * i], :) = ...
    reshape(p, 4, 3)';
end

end

