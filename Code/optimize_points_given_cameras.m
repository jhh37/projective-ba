function [P, X, cost, iter] = optimize_points_given_cameras(M, W, P, X, store, varargin)
%OPTIMIZE_POINTS_GIVEN_CAMERAS Summary of this function goes here
%   Detailed explanation goes here

%% Settings
opts = au_opts( ...
  'camera_model=projective', ...
  'check_derivatives=0', ...
  'max_iter=300', ...
  'max_trial=50', ...
  'func_tol=1e-9', ...
  'max_lambda=1e+14', ...
  'min_lambda=1e-14', ...
  'lambda_inc_factor=10', ...
  'lambda=1e-1', ...
  'homogeneous=1', ...
  'display=1', ...
  'alpha=0', ...
  'init_from_alg_points=1', ...
  'mu=1', ...
  varargin{:});

% Preprocess dataset if not already done.
if nargin < 5 || ~isstruct(store)
  store = preprocess_dataset(M, W);
end

%% Derivative check
if opts.check_derivatives,
  % Check the Jacobian with respect to points.
  JX_blocks = compute_jacobian_over_points(P, X, opts, store);
  JX = blkdiag(JX_blocks{:});
  au_check_derivatives(@(x) ...
    compute_residual(P, ...
    reshape(x, store.dim.num_point_params, store.dim.num_points), ...
    opts, store), X(:), JX, 'tol=1e-4');
  return
end

%% Initialization
% Normalize each point for numerical stability.
if opts.init_from_alg_points,
  X = compute_algebraic_points(P, opts, store);
else
  for j = 1 : store.dim.num_points,
    X(:, j) = X(:, j) / norm(X(:, j));
  end
end

%% Second-order solver
num_max_iters = 0;
for j = 1 : store.dim.num_points
  % Solver parameters
  eval = 0;
  lambda = opts.lambda;
  err = compute_residual(P, X, opts, store, j);
  cost = norm(err) / sqrt(2 * store.dim.num_visible_frames_per_point(j));
  
  % Solve iteratively for each point.
  for iter = 1 : opts.max_iter
    %     if ~opts.homogeneous,
    %       J = compute_inhom_jacobian_over_points(P, opts, store, j);
    %       JTJ = J' * J;
    %       g = J' * err;
    %     else
    % Compute the Jacobian with respect to the points.
    JX = compute_jacobian_over_points(P, X, opts, store, j);
    
    % Project the Jacobian and the gradient to the tangent space of the
    % points.
    Pr = null(X(:, j)');
    J = JX * Pr;
    JTJ = J' * J;
    g = J' * err;
    %     end
    
    for trial = 1 : opts.max_trial
      eval = eval + 1;
      % Compute the update
      dw = - (JTJ + lambda * speye(3)) \ g;
      X_eval = X;
      %       if ~opts.homogeneous,
      %         X_eval(1:3, j) = X(1:3, j) + dw;
      %       else
      X_eval(:, j) = X(:, j) + Pr * dw;
      X_eval(:, j) = X_eval(:, j) / norm(X_eval(:, j));
      %       end
      err_eval = compute_residual(P, X_eval, opts, store, j);
      
      if norm(err_eval) - norm(err) < 0 ...
          || abs(norm(err_eval) - norm(err)) / ...
          sqrt(2 * store.dim.num_visible_frames_per_point(j)) ...
          < opts.func_tol,
        X = X_eval;
        err = err_eval;
        lambda = max(lambda / opts.lambda_inc_factor, opts.min_lambda);
        break
      elseif trial == opts.max_trial
        break
      else
        lambda = min(lambda * opts.lambda_inc_factor, opts.max_lambda);
      end
    end
    
    % Print the current status.
    cost_iter = norm(err) / sqrt(2 * store.dim.num_visible_frames_per_point(j));
    if opts.display
      fprintf('[Point %04d][%04d][%04d] %.6e\n', j, iter, eval, cost_iter);
    end
    
    if abs(cost - cost_iter) < opts.func_tol,
      % If the function tolerance is reached, quit.
      break
    else
      % Otherwise, update the best cost value.
      cost = cost_iter;
    end
  end
  
  if iter == opts.max_iter
    num_max_iters = num_max_iters + 1;
  end
end

% Set the scales to 1.
if ~opts.homogeneous,
  X = X ./ repmat(X(4, :), 4, 1);
end

% Output the number of incomplete convergence.
if opts.display,
  fprintf('Final error: %.6e\n', ...
    norm(compute_residual(P, X, opts, store)) / sqrt(store.dim.nnz_frames * 2));
  fprintf('Total number of points with incomplete number of convergence: %d\n', num_max_iters);
end

end
