function X = compute_algebraic_points(P, opts, store)
%COMPUTE_ALGEBRAIC_POINTS Summary of this function goes here
%   Detailed explanation goes here
% for j = 1 : store.dim.num_points
%   A = store.blocks.mcp{j} * P(store.index.visible_camera_rows{j}, 1 : 3);
%   b = store.blocks.mcp{j} * X(4, j) * P(store.index.visible_camera_rows{j}, 4);
%   X(1 : 3, j) = - A \ b;
% end

X = nan(4, store.dim.num_points);
for j = 1 : store.dim.num_points
  A = store.blocks.mcp{j} * P(store.index.visible_camera_rows{j}, :);
  [~, ~, V] = svd(A, 0);
  if opts.homogeneous,
    X(:, j) = V(:, end) / norm(V(:, end));
  else
    X(:, j) = V(:, end) / V(end);
  end
end

end

