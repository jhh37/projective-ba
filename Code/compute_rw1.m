function RW1 = compute_rw1(P, X, err, opts, store, camera_num)
%COMPUTE_RW1 Summary of this function goes here
%   Detailed explanation goes here

% dJxTe_dp_proj_ij = @(p, xt_j, e_ij) ...
%   [ kron(e_ij', eye(4)), zeros(4,4) ] / (p(9:12)' * xt_j) ...
%   - p(9:12) * [ e_ij(1) * xt_j', e_ij(2) * xt_j', zeros(1,4)] / (p(9:12)' * xt_j) ^ 2 ... 
%   - [ zeros(4,8), [p(1:4), p(5:8)] * e_ij * xt_j'] / (p(9:12)' * xt_j) ^ 2 ...
%   + 2 * [ zeros(4,8), p(9:12) * [ p(1:4)' * xt_j, p(5:8)' * xt_j ] * e_ij * xt_j'] / (p(9:12)' * xt_j) ^ 3 ...
%   - [ zeros(4,8), ([p(1:4)' * xt_j, p(5:8)' * xt_j] * e_ij) * eye(4)] / (p(9:12)' * xt_j) ^ 2;
% 
% dJxTe_dp_aff_ij = @(p, xt_j, e_ij) ...
%   [ kron(e_ij', eye(4)), zeros(4,4) ] / (p(12) * xt_j(4)) ...
%   - [zeros(3,1); p(12)] * [ e_ij(1) * xt_j', e_ij(2) * xt_j', zeros(1,4)] / (p(12) * xt_j(4)) ^ 2 ... 
%   - [ zeros(4,8), zeros(4,3), [p(1:4), p(5:8)] * e_ij * xt_j(4)] / (p(12) * xt_j(4)) ^ 2 ...
%   + 2 * [ zeros(4,8), zeros(4,3), [zeros(3,1); p(12)] * [ p(1:4)' * xt_j, p(5:8)' * xt_j ] * e_ij * xt_j(4)] / (p(12)' * xt_j(4)) ^ 3 ...
%   - [ zeros(4,8), zeros(4,3), [zeros(3, 1); [p(1:4)' * xt_j, p(5:8)' * xt_j] * e_ij] ] / (p(12) * xt_j(4)) ^ 2;


if nargin < 6,
  % Check whether we want the Jacobian over the entire points or only one.
  rng_i = 1 : store.dim.num_frames;
  RW1 = cell(store.dim.num_frames, 1);
else
  rng_i = camera_num(1);
end

for i = rng_i
  % s = p13i * x_j
  if strcmpi(opts.camera_model, 'affine')
    % If affine, compute the corresponding derivatives.
    e_i = err(store.rng.proj_per_frame{i});
    s = (P(3 * i, 4) * X(4, store.index.visible_points{i}))';
    p3 = [sparse(1, 3), P(3 * i, 4)]';
    x3 = [sparse(store.dim.num_visible_points_per_frame(i), 3), ...
      X(4, store.index.visible_points{i})'];
    e3 = [sparse(3, 1); 1];
  elseif opts.mu == 1
    % Otherwise, assume the projective model and compute the derivatives.
    e_i = err(store.rng.proj_per_frame{i});   
    p3 = P(3 * i, :)';
    x3 = X(:, store.index.visible_points{i})';
    s = x3 * p3;
    e3 = ones(4, 1);
  else
    % Otherwise, assume the projective model and compute the derivatives.
    e_i = err(store.rng.proj_per_frame{i});   
    p3 = P(3 * i, :)';
    x3 = repmat([opts.mu, opts.mu, opts.mu, 1], ...
      store.dim.num_visible_points_per_frame(i), 1) .* ...
      X(:, store.index.visible_points{i})';
    s = x3 * p3;
    p3 = [opts.mu, opts.mu, opts.mu, 1]' .* p3;
    e3 = [opts.mu, opts.mu, opts.mu, 1];
  end
  
  % P12i * x_j / s_ij
  P12xj = (P([3 * i - 2, 3 * i - 1], :) * ...
    X(:, store.index.visible_points{i}))' ./ ...
    repmat(s, 1, 2);
  
  % Create e_ij / s_ij.
  T0 = reshape(e_i, 2, store.dim.num_visible_points_per_frame(i))' ./ ...
    repmat(s, 1, 2);
  
  % Create (P12i * x_j / s_ij) * (e_ij / s_iij).
  T1 = sum(P12xj .* T0, 2);
  
  % Create term 2 and 3 (T0 = reshaped e_i for camera i).
  % T1 = kron(T0, speye(4));
  T2 = repelem(T0 ./ repmat(s, 1, 2), 1, 4) .* ...
    repmat(X(:, store.index.visible_points{i})', 1, 2);
  T2 = kron(T0, speye(4)) - repmat(p3, store.dim.num_visible_points_per_frame(i), 8) .* ...
    repelem(T2, 4, 1);
  
  T3 = - reshape(P([3*i-2, 3*i-1] ,:)' * (T0 ./ repmat(s, 1, 2))', ...
    4 * store.dim.num_visible_points_per_frame(i), 1);
  T3 = repmat(T3, 1, 4) .* repelem(x3, 4, 1);
  
  T4 = 2 * repmat(p3, store.dim.num_visible_points_per_frame(i), 4) .* ...
    repelem(repmat(T1 ./ s, 1, 4) .* x3, 4, 1) ...
    - kron(T1, diag(e3));
  
  RW1_unsorted = [T2, T3 + T4];

  if nargin < 6,
    % If the camera number is not given, compute the entire Jacobian.
    RW1{i} = full(RW1_unsorted(:));
  else
    % Otherwise, compute only the Jacobian with respect to camera i.
    RW1 = ...
      au_sparse(int32(store.rng.RW1.ii{i}), ...
      int32(store.rng.RW1.jj{i}), RW1_unsorted(:));
    return
  end
end

RW1 = au_sparse(store.rng.RW1.iii, store.rng.RW1.jjj, vertcat(RW1{:}));

end
