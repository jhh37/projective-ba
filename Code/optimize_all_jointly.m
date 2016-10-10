function [P, X, cost, iter, msg, store] = optimize_all_jointly(M, W, P, X, store, varargin)
%OPTIMIZE_ALL_JOINTLY Summary of this function goes here
%   Only 3 dimensions used.

%% PRE-PROCESSING
% Fetch all the options.
opts = au_opts( ...
  'init_algebraic_points=0', ...
  'epsilon=1', ...
  'camera_model=projective', ...
  'start_from_optimal_x=1', ...
  'homogeneous=1', ...
  'lambda=1e-1', ... % damping parameter
  'lambda_inc_factor=10', ... % damping increase rate
  'lambda_dec_factor=10', ... % damping decrease rate
  'max_lambda=1e+14', ... % maximum value of lambda
  'min_lambda=1e-14', ... % minimum value of lambda
  'max_iter=300', ... % maximum number of iterations
  'max_trial=50', ... % maximum number of trials for each iteration
  'func_tol=1e-9', ... % function tolerance
  'display=1', ... % whether to display progress
  'mu=1', ...
  varargin{:} );

% Preprocess dataset if not already done.
if nargin < 5 || ~isstruct(store)
  store = preprocess_dataset(M, W);
end

%% INITIALIZATION
% Constrain the norm of p.
if opts.homogeneous,
  for i = 1 : store.dim.num_frames
    p = reshape(P([3 * i - 2, 3 * i - 1, 3 * i], :)', store.dim.num_camera_params, 1);
    if strcmpi(opts.camera_model, 'affine'),
      a = p([9 10 11]);
      p([9 10 11]) = 0;
      if opts.homogeneous,
        p = p / norm(p);
      else
        p = p / p(end);
      end
      p([9 10 11]) = a;
    else
      if opts.homogeneous,
        p = p / norm(p);
      else
        p = p / p(end);
      end
    end
    P([3 * i - 2, 3 * i - 1, 3 * i], :) ...
      = reshape(p, 4, 3)';
  end
end

% Compute necessary components from the initial set of P & X.
err = compute_residual_framewise(P, X, opts, store);
cost = norm(err) / sqrt(2 * store.dim.nnz_frames);
% [p0x0.err_vec, p0x0.cost, p0x0.pi, p0x0.z, p0x0.z0] = compute_residual(P, X, store, opts);
% If display is on, show the initial cost.
if opts.display,
  fprintf('[%04d][%04d] %.2e %.6e\n', 0, 0, 0, cost);
end

if opts.init_algebraic_points,
  % Initialize points with algebraic error if necessary.
  if opts.display,
    fprintf('Computing algebraic solution for points...\n');
  end
  X = compute_algebraic_points(P, opts, store);
  err = compute_residual_framewise(P, X, opts, store);
  cost = norm(err) / sqrt(2 * store.dim.nnz_frames);
  % If display is on, show the initial cost.
  if opts.display,
    fprintf('[%04d][%04d] %.2e %.6e\n', 0, 0, 0, cost);
  end
end

if opts.start_from_optimal_x
  % Compute optimal initial x*
  if opts.display,
    fprintf('Computing optimal points for initial cameras ...\n');
  end
  
  [~, X] = optimize_points_given_cameras(M, W, P, X, nan, varargin{:}, ...
    'display=0', ['homogeneous=', opts.homogeneous]);
  err = compute_residual_framewise(P, X, opts, store);
  cost = norm(err) / sqrt(2 * store.dim.nnz_frames);
  % If display is on, show the initial cost.
  if opts.display,
    fprintf('[%04d][%04d] %.2e %.6e\n', 0, 0, 0, cost);
  end
end

%% SECOND-ORDER SUBPROBLEM SOLVER
% Compute Joint from p0, x0
if opts.display,
  fprintf('Running joint optimization...\n');
end
eval = 0;
lambda = opts.lambda;
Xp = cell(store.dim.num_points, 1);
Pr_blocks = cell(store.dim.num_frames, 1);
for iter = 1 : opts.max_iter
  
  if ~opts.homogeneous,
    JX = compute_inhom_jacobian_over_points(P, opts, store);
    % If not only optmizing over point, compute Jacobian over cameras.
    JP = compute_inhom_jacobian_over_cameras(X, opts, store);
  else
    JX = compute_jacobian_over_points(P, X, opts, store);
    % If not only optmizing over point, compute Jacobian over cameras.
    JP = compute_jacobian_over_cameras(P, X, opts, store);
    
    % Project cameras to the space orthogonal to scaling.
    for i = 1 : store.dim.num_frames
      p = reshape(P([3 * i - 2, 3 * i - 1, 3 * i], :)', 12, 1);
      Pr_blocks{i} = null(p');
      JP{i} = sparse(JP{i} * Pr_blocks{i});
    end
    
    % Project points to the space orthogonal to scaling.
    for j = 1 : store.dim.num_points
      Xp{j} = sparse(null(X(:, j)'));
      JX{j} = sparse(JX{j} * Xp{j});
    end
  end
  
  J = [sparse(blkdiag(JP{:})), store.mat.K_p2f * blkdiag(JX{:})];
  
  g = J' * err;
  JTJ = J' * J;
  
  % Evaluate until the cost decreases or reaches the max number of fail.
  for trial = 1 : opts.max_trial
    eval = eval + 1;
    
    % Compute new P and X evaluations.
    dw = - (JTJ + lambda * speye(size(JTJ, 1))) \ g;
    
    % Check if the condition number is too high.
    if ~isempty(strfind(lastwarn, 'NaN'))
      lastwarn('');
      msg = 3;
      return
    end

    if strcmpi(opts.camera_model, 'affine') && ~opts.homogeneous
        P_eval = P;
        dp = dw(1 : 8 * store.dim.num_frames);
        P_eval(store.index.affine_camera, :) = ...
        P_eval(store.index.affine_camera, :) + ...
        reshape(dp, 4, 2 * store.dim.num_frames)';
        dx = dw(8 * store.dim.num_frames + 1 : end);
        X_eval = X;
        X_eval(1 : 3, :) = reshape(dx, 3, store.dim.num_points);
    else
      dp = blkdiag(Pr_blocks{:}) * dw(1 : store.dim.num_frames * 11);
      dx = blkdiag(Xp{:}) * dw(store.dim.num_frames * 11 + 1 : end);
    end
    
    if opts.homogeneous,
      P_eval = P + reshape(dp, 4, 3 * store.dim.num_frames)';
      for i = 1 : store.dim.num_frames
        P_eval([3 * i - 2, 3 * i - 1, 3 * i], :) = ...
          P_eval([3 * i - 2, 3 * i - 1, 3 * i], :) / ...
          norm(P_eval([3 * i - 2, 3 * i - 1, 3 * i], :), 'fro');
      end

      X_eval = X + reshape(dx, 4, store.dim.num_points);
      for j = 1 : store.dim.num_points
        X_eval(:, j) = X_eval(:, j) / norm(X_eval(:, j));
      end
    end
    
    % Compute the cost: this will be replaced with a more efficient code
    % later on.
    err_eval = compute_residual_framewise(P_eval, X_eval, opts, store);
    
    if norm(err_eval) - norm(err) < 0 ...
        || abs(norm(err_eval) - norm(err)) / ...
        sqrt(2 * store.dim.nnz_frames) ...
        < opts.func_tol,
      X = X_eval;
      %       X = X ./ repmat(X(4, :), 4, 1); % normalize
      P = P_eval;
      %       for i = 1 : store.dim.num_frames
      %        P([3 * i - 2, 3 * i - 1, 3 * i], :) = ...
      %          P([3 * i - 2, 3 * i - 1, 3 * i], :) / P(3 * i, 4);
      %       end
      err = err_eval;
      lambda = max(lambda / opts.lambda_inc_factor, opts.min_lambda);
      break
    elseif trial == opts.max_trial
      break
    else
      lambda = min(lambda * opts.lambda_inc_factor, opts.max_lambda);
    end
  end
  
  if trial == opts.max_trial
    msg = 1; % reached max. no. of inner trials
    return
  end
  
  % Print the current status.
  cost_iter = norm(err) / sqrt(2 * store.dim.nnz_frames);
  if opts.display
    fprintf('[%03d][%03d] %.6e %.1e\n', iter, eval, cost_iter, lambda);
  end
  
  if abs(cost - cost_iter) < opts.func_tol,
    % If the function tolerance is reached, quit.
    msg = 0; % reached func. tol.
    return
  else
    % Otherwise, update the best cost value.
    cost = cost_iter;
  end
end

if ~exist('msg', 'var'),
  msg = 2; % reached max no. of iterations.
end

end
