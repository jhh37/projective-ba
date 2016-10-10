function [P, X, cost, iter] = optimize_cameras_given_points(M, W, P, X, store, varargin)
%OPTIMIZE_CAMERAS_GIVEN_POINTS Summary of this function goes here
%   Detailed explanation goes here

%% Settings
% Custom functions
vec = @(X) X(:);

% Options
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
  varargin{:});

% Preprocess dataset if not already done.
if nargin < 5 || ~isstruct(store)
  store = preprocess_dataset(M, W);
end

%% Derivative check
if opts.check_derivatives,
  % Check the Jacobian with respect to cameras.
  JP_blocks = compute_jacobian_over_cameras(P, X, opts, store);
  JP = blkdiag(JP_blocks{:});
  au_check_derivatives(@(p) ...
    compute_residual_framewise(...
    reshape(p, 4, 3 * store.dim.num_frames)', ...
    X, opts, store), vec(P'), JP, 'tol=1e-4');
  return
end

%% Initialization
% Normalize each camera for numerical stability (retraction).
if opts.homogeneous,
  for i = 1 : store.dim.num_frames,
    P([3 * i - 2, 3 * i - 1, 3 * i], :) = ...
      P([3 * i - 2, 3 * i - 1, 3 * i], :) ...
      / norm(vec(P([3 * i - 2, 3 * i - 1, 3 * i], :)));
  end
end

%% Second-order solver
num_max_iters = 0;
for i = 1 : store.dim.num_frames
  % Solver parameters
  eval = 0;
  lambda = opts.lambda;
  cost = inf;
  err = compute_residual_framewise(P, X, opts, store, i);
  
  % Solve iteratively for each point.
  for iter = 1 : opts.max_iter
    if ~opts.homogeneous,
      J = compute_inhom_jacobian_over_cameras(X, opts, store, i);
      JTJ = J' * J;
      g = J' * err;
    else
      % Compute the Jacobian with respect to the points.
      JP = compute_jacobian_over_cameras(P, X, opts, store, i);
      
      % Project the Jacobian and the gradient to the tangent space of the
      % points.
      Pr = null(vec(P([3 * i - 2, 3 * i - 1, 3 * i], :)')');
      J = JP * Pr;
      JTJ = J' * J;
      g = J' * err;
    end
    
    for trial = 1 : opts.max_trial
      eval = eval + 1;
      % Compute the update
      dw = - (JTJ + lambda * speye(size(JTJ, 1))) \ g;
      P_eval = P;
      if ~opts.homogeneous,
        P_eval([3 * i - 2, 3 * i - 1], :) ...
          = P_eval([3 * i - 2, 3 * i - 1], :) + reshape(dw, 4, 2)';
      else
        p_eval = vec(P_eval([3 * i - 2, 3 * i - 1, 3 * i], :)') + Pr * dw;
        P_eval([3 * i - 2, 3 * i - 1, 3 * i], :) ...
          = reshape(p_eval / norm(p_eval), 4, 3)';
      end
      err_eval = compute_residual_framewise(P_eval, X, opts, store, i);
      
      if norm(err_eval) - norm(err) < 0 ...
          || abs(norm(err_eval) - norm(err)) / ...
          sqrt(2 * store.dim.num_visible_points_per_frame(i)) ...
          < opts.func_tol,
        P = P_eval;
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
    cost_iter = norm(err) / sqrt(2 * store.dim.num_visible_points_per_frame(i));
    if opts.display
      fprintf('[Camera %04d][%03d][%03d] %.6e\n', i, iter, eval, cost_iter);
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

% Output the number of incomplete convergence.
if opts.display,
  fprintf('Final error: %.6e\n', ...
    norm(compute_residual(P, X, opts, store)) / sqrt(store.dim.nnz_frames * 2));
  fprintf('Total number of cameras with incomplete number of convergence: %d\n', num_max_iters);
end

end
