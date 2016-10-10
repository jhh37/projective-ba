function store = preprocess_dataset(M, W)
%PREPROCESS_DATASET Summary of this function goes here
%   Detailed explanation goes here

% Custom function
vec = @(X) X(:);

% Create a reduced frame matrix from the weight matrix.
W_frame = floor(W(1:2:end, :));

%% DIMENSIONS
% Get number of points and frames.
[store.dim.num_frames, store.dim.num_points] = size(W_frame);

store.dim.num_camera_params = 12;
store.dim.num_point_params = 4;
store.dim.p = store.dim.num_camera_params * store.dim.num_frames;
store.dim.x = store.dim.num_point_params * store.dim.num_points;

store.dim.hess = store.dim.x + store.dim.p;
store.index.starting_dx = store.dim.p + 1;

%% INDEXES
store.index.affine_camera = vec(bsxfun(@plus, 3 * (0 : 1 : (store.dim.num_frames - 1)), [1 2]'));

% Build various indexes for visible frames.
store.index.visible_frames = cell(store.dim.num_points, 1);
store.index.visible_proj_rows = cell(store.dim.num_points, 1);
store.index.visible_camera_rows = cell(store.dim.num_points, 1);
store.index.visible_camera_r3 = cell(store.dim.num_points, 1);
store.index.visible_camera_r12 = cell(store.dim.num_points, 1);
for j = 1 : store.dim.num_points,
  store.index.visible_frames{j} = find(W_frame(:, j));
  store.index.visible_proj_rows{j} = ...
    vec(bsxfun(@plus, 2 * store.index.visible_frames{j}, [-1 0])');
  store.index.visible_camera_rows{j} = ...
    vec(bsxfun(@plus, 3 * store.index.visible_frames{j}, [-2 -1 0])');
  store.index.visible_camera_r3{j} = ...
    store.index.visible_camera_rows{j}(3 : 3 : end);
  store.index.visible_camera_r12{j} = ...
    vec(bsxfun(@plus, 3 * store.index.visible_frames{j}, [-2 -1 ])');
end
store.dim.num_visible_frames_per_point = sum(W_frame)';
store.index.starting_frame = cumsum([1; store.dim.num_visible_frames_per_point]);
store.index.starting_proj = 2 * store.index.starting_frame - 1;
store.index.starting_elem = 3 * store.index.starting_frame - 2;
store.dim.nnz_frames = sum(store.dim.num_visible_frames_per_point);
% [store.index.ii, store.index.jj] = find(W_frame);

% Build various indexes for visible points.
store.index.visible_points = cell(store.dim.num_frames, 1);
for i = 1 : store.dim.num_frames
  store.index.visible_points{i} = find(W_frame(i, :)');
end
store.dim.num_visible_points_per_frame = sum(W_frame, 2);
store.index.starting_point = cumsum([1; store.dim.num_visible_points_per_frame]);
store.index.starting_proj_framewise = 2 * store.index.starting_point - 1;
store.index.starting_elem_framewise = 3 * store.index.starting_point - 2;

%% MEASUREMENT BLOCKS
store.blocks.measurements = cell(store.dim.num_points, 1);
for j = 1 : store.dim.num_points
  store.blocks.measurements{j} = M(store.index.visible_proj_rows{j}, j);
end

store.blocks.measurements_framewise = cell(store.dim.num_frames, 1);
for i = 1 : store.dim.num_frames
  store.blocks.measurements_framewise{i} = ...
    reshape(M([2 * i - 1, 2 * i], store.index.visible_points{i}), ...
    2 * store.dim.num_visible_points_per_frame(i), 1);
end

% Projective form
store.vec.mt = vec(reshape(M(W~=0), 2, store.dim.nnz_frames));
store.vec.mt_framewise = nan(size(store.vec.mt));

% Trilinear form
store.vec.mt2 = vec([reshape(M(W~=0), 2, store.dim.nnz_frames); ...
  ones(1, store.dim.nnz_frames)]);

% tic
% store.blocks.measurements = cell(store.dim.nnz_frames, 1);
% for kk = 1 : store.dim.nnz_frames
%     store.blocks.measurements{kk} = M(bsxfun(@plus, ...
%         2 * repmat(store.index.ii(kk), 1, 2), [-1 0]), ...
%         store.index.jj(kk));
% end
% toc

%% RANGE PATTERNS
% Projections and elements
store.rng.frames_per_point = cell(store.dim.num_points, 1);
store.rng.proj_per_point = cell(store.dim.num_points, 1);
store.rng.elem_per_point = cell(store.dim.num_points, 1);
for j = 1 : store.dim.num_points
  store.rng.frames_per_point{j} = ...
    store.index.starting_frame(j) : store.index.starting_frame(j + 1) - 1;
  store.rng.proj_per_point{j} = ...
    store.index.starting_proj(j) : store.index.starting_proj(j + 1) - 1;
  store.rng.elem_per_point{j} = ...
    store.index.starting_elem(j) : store.index.starting_elem(j + 1) - 1;
end

store.rng.points_per_frame = cell(store.dim.num_frames, 1);
store.rng.proj_per_frame = cell(store.dim.num_frames, 1);
store.rng.elem_per_frame = cell(store.dim.num_frames, 1);
for i = 1 : store.dim.num_frames
  store.rng.points_per_frame{i} = ...
    store.index.starting_point(i) : store.index.starting_point(i + 1) - 1;
  store.rng.proj_per_frame{i} = ...
    store.index.starting_proj_framewise(i) : store.index.starting_proj_framewise(i + 1) - 1;
  store.rng.elem_per_frame{i} = ...
    store.index.starting_elem_framewise(i) : store.index.starting_elem_framewise(i + 1) - 1;
end

% RW1 patterns
store.rng.RW1.ii = cell(store.dim.num_frames, 1);
store.rng.RW1.jj = cell(store.dim.num_frames, 1);
for i = 1 : store.dim.num_frames
  store.rng.RW1.ii{i} = ...
    vec(bsxfun(@plus, 4 * store.index.visible_points{i}, [-3 -2 -1 0])');
  store.rng.RW1.ii{i} = ...
    repmat(store.rng.RW1.ii{i}, store.dim.num_camera_params, 1);
  store.rng.RW1.jj{i} = ...
    vec(repelem(12 * i + (-11 : 0), 4 * store.dim.num_visible_points_per_frame(i), 1));
end
  store.rng.RW1.iii = int32(vertcat(store.rng.RW1.ii{:}));
  store.rng.RW1.jjj = int32(vertcat(store.rng.RW1.jj{:}));

% Add regularizer to projective
store.rng.aff_reg = int32(vec( ...
  bsxfun(@plus, [9 10 11]', 0 : 12 : 12 * (store.dim.num_frames - 1)) ...
  ));
store.rng.aff_reg = [ store.rng.aff_reg; ...
  store.dim.num_camera_params * store.dim.num_frames ];
  
% Jp, Jq (dimension-reduced)
store.rng.Jp = int32(vec(bsxfun(@plus, ...
  [1 2 3 4 5 6 7 11]', ...
  0 : 11 : (store.dim.num_frames - 1) * 11)));

store.rng.Jq = int32(vec(bsxfun(@plus, ...
  [8 9 10]', ...
  0 : 11 : (store.dim.num_frames - 1) * 11)));

%% ALGEBRAIC COMPUTATION
% [M]x * P
store.blocks.mcp = cell(store.dim.num_points, 1);
for j = 1 : store.dim.num_points
  % Compute sparse matrix values.
  vals = reshape(store.blocks.measurements{j}, 2, store.dim.num_visible_frames_per_point(j));
  vals = [ repmat([1; -1], 1, store.dim.num_visible_frames_per_point(j)); ...
    vals(2, :); - vals(1, :) ];
  
  % Compute row and column indexes.
  row_idx = bsxfun(@plus, [2 1 1 2]', 2 * (0 : store.dim.num_visible_frames_per_point(j) - 1));
  col_idx = bsxfun(@plus, [1 2 3 3]', 3 * (0 : store.dim.num_visible_frames_per_point(j) - 1));
  store.blocks.mcp{j} = au_sparse(int32(row_idx(:)), int32(col_idx(:)), vals(:));
end
store.mat.Mxt = blkdiag(store.blocks.mcp{:});

%% Wt matrix
col_idx = cell(store.dim.num_points, 1);
for j = 1 : store.dim.num_points,
  col_idx{j} = store.index.visible_camera_rows{j} + 3 * store.dim.num_frames * (j - 1);
end
col_idx = vertcat(col_idx{:});
row_idx = (1 : 3 * store.dim.nnz_frames)';
[col_idx, sort_idx] = sort(col_idx);
row_idx = row_idx(sort_idx);
store.mat.Wt = au_sparse(int32(row_idx), int32(col_idx), ones(3 * store.dim.nnz_frames, 1));

%% PERMUTATIONS
store.mat.K_P = sparse(1:12 * store.dim.num_frames, vec(reshape(1: 12 * store.dim.num_frames, 3 * store.dim.num_frames, 4)'), 1);

%% mt framewise
for i = 1 : store.dim.num_frames
  tmp_mt = M([2 * i - 1, 2 * i], :);
  store.vec.mt_framewise(store.rng.proj_per_frame{i}) = tmp_mt(~isnan(tmp_mt));
end

W2 = W(1:2:2*store.dim.num_frames, :);
kk = find(W2);
[~, jj] = sort(mod(kk - 1, store.dim.num_frames));
store.index.pointwise_to_framewise = vec(bsxfun(@plus, 2 * (jj - 1), [1 2])');
store.mat.K_p2f = au_sparse(int32(store.index.pointwise_to_framewise), ...
  int32(1 : 2 * store.dim.nnz_frames), ...
  ones(2 * store.dim.nnz_frames, 1))';

end

