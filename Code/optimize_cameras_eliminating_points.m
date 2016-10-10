function [P, X, cost, iter, msg, store] = optimize_cameras_eliminating_points(M, W, P, X, store, varargin)
%OPTIMIZE_CAMERAS_ELIMINATING_POINTS Summary of this function goes here
%   Detailed explanation goes here

%% Settings
% Custom functions
vec = @(X) X(:);

% Options
opts = au_opts( ...
  'camera_model=projective', ...
  'rw=2', ...
  'nu=0', ...
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
  'penalize_gauge_freedom=1', ...
  'mu=1', ...
  'mu_=0', ...
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


%% Second-order solver

% Solver parameters
eval = 0;
lambda = opts.lambda;
[~, X] = optimize_points_given_cameras(M, W, P, X, store, ...
  ['camera_model=', opts.camera_model], ...
  ['homogeneous=', opts.homogeneous], ...
  ['mu=', opts.mu], ...
  'display=0');
err = compute_residual_framewise(P, X, opts, store);
cost = norm(err) / sqrt(2 * store.dim.nnz_frames);

% Solve iteratively for each point.
JXR_blocks = cell(store.dim.num_points, 1);
JXQ_blocks = cell(store.dim.num_points, 1);
Pr_blocks = cell(store.dim.num_frames, 1);
K = sparse(1 : store.dim.num_camera_params * ...
  store.dim.num_frames, ...
  vec(reshape(1:store.dim.num_camera_params * ...
  store.dim.num_frames, 3 * store.dim.num_frames, 4)'), 1);
Xp = cell(store.dim.num_points, 1);
dw = nan(11 * store.dim.num_frames, 1);
if ~opts.homogeneous && strcmpi(opts.camera_model, 'affine')
  % Indices for filling the inhomogeneous affine rank
  a = bsxfun(@plus, [-3 -2 -1]', 4 : 4 : 8 * store.dim.num_frames);
  b = 4 : 4 : 8 * store.dim.num_frames;
end
for iter = 1 : opts.max_iter
  if ~opts.homogeneous,
    % Only RW2 implemented for inhomogeneous affine camera.
    JP = compute_inhom_jacobian_over_cameras(X, opts, store);
    JX = compute_inhom_jacobian_over_points(P, opts, store);
    
    for j = 1 : store.dim.num_points
      [JXQ_blocks{j}, ~] = qr(blkdiag(JX{j}), 0);
      JXQ_blocks{j} = sparse(JXQ_blocks{j});
    end
    JXQ2 = store.mat.K_p2f * blkdiag(JXQ_blocks{:});
    JP2 = blkdiag(JP{:});
    JTJ = JP2' * JP2 - (JXQ2' * JP2)' * (JXQ2' * JP2);
    
    if opts.penalize_gauge_freedom,
      P13 = P(store.index.affine_camera, 1:3);
      p4 = P(store.index.affine_camera, 4);
      JTJ(a, a) = JTJ(a, a) + kron(P13 * P13', speye(3));
      JTJ(b, b) = JTJ(b, b) + p4 * p4';
    end
    g = JP2' * err;
  else
    % Compute the Jacobian with respect to the points.
    JP = compute_jacobian_over_cameras(P, X, opts, store);
    JX = compute_jacobian_over_points(P, X, opts, store);
    for j = 1 : store.dim.num_points
      Xp{j} = sparse(null(X(:, j)'));
      JX{j} = sparse(JX{j} * Xp{j});
    end
    if opts.rw == 1,
      for j = 1 : store.dim.num_points
        [JXQ_blocks{j}, JXR_blocks{j}] = qr(JX{j}, 0);
      end
    else
      for j = 1 : store.dim.num_points
        [JXQ_blocks{j}, ~] = qr(JX{j}, 0);
      end
    end
    JXQ2 = store.mat.K_p2f * blkdiag(JXQ_blocks{:});
    JP2 = blkdiag(JP{:});
    
    % Project the Jacobian and the gradient to the tangent space of the
    % points.
    for i = 1 : store.dim.num_frames
      p = vec(P([3 * i - 2, 3 * i - 1, 3 * i], :)');
      if strcmpi(opts.camera_model, 'affine'),
        % If using an affine model, set p31, p32 and p33 to 0 and then
        % truncate the Jacobian to 9 parameters per camera (from 12
        % parameters).
        p([9 10 11]) = 0;
        Pr_blocks{i} = null(p');
        Pr_blocks{i} = Pr_blocks{i}(:, [1 : 7, 11]);
      else
        Pr_blocks{i} = null(p');
      end
    end
    Pr = blkdiag(Pr_blocks{:});
    JP2 = JP2 * Pr;
    
    % Create a RW2 J.
    % J = (speye(2 * store.dim.nnz_frames) - JXQ2 * JXQ2') * JP2;
    if opts.rw == 1,
      % If using RW1, add the corresponding term.
      % J = J - JXQ2 * (JXR' \ ...
      %   (blkdiag(Xp{:})' * compute_rw1(P, X, err, opts, store) * Pr) ...
      %  );
      RW1 = JXQ2 * ( ...
        blkdiag(JXR_blocks{:})' \ (blkdiag(Xp{:})' * compute_rw1(P, X, err, opts, store) * Pr));
    end
    
    J = (speye(2 * store.dim.nnz_frames) - JXQ2 * JXQ2') * JP2;
    if opts.rw == 1,
      J = J - RW1;
    end
    
    if opts.penalize_gauge_freedom,
      PA = P;
      if strcmpi(opts.camera_model, 'affine')
        % If affine, zero out p31, p32 and p33 elements.
        PA(3 : 3 : 3 * store.dim.num_frames, 1:3) = ...
          zeros(store.dim.num_frames, 3);
      end
      PPT = kron(speye(4), PA') * K' * Pr;
      JTJ = J' * J + PPT' * PPT;
    else
      JTJ = J' * J;
    end

    g = J' * err;
    %     JTJ = JP2' * JP2 - JP2' * (JXQ2 * JXQ2') * JP2 + PPT' * PPT;
    %     g = JP2' * err;
    
    if opts.mu_ > 0 && ~strcmpi(opts.camera_model, 'affine'),
      p3_13 = vec(P(3:3:3*store.dim.num_frames, 1:3)');
      JTJ3 = au_sparse(...
        store.rng.aff_reg, ...
        store.rng.aff_reg, ...
        [opts.mu_ * ones(3 * store.dim.num_frames, 1); 0]);
      JTJ = JTJ + Pr' * JTJ3 * Pr;
      g3 = au_sparse(...
        store.rng.aff_reg, ...
        int32(ones(3 * store.dim.num_frames + 1, 1)), ...
        [opts.mu_ * p3_13; 0] ...
        );
      g = g + Pr' * g3;
      %       JTJ2 = au_sparse(...
      %         int32(12 : 12 : 12 * store.dim.num_frames), ...
      %         int32(12 : 12 : 12 * store.dim.num_frames), ...
      %         opts.mu_ * ones(store.dim.num_frames, 1));
      %       JTJ = JTJ + Pr' * JTJ2 * Pr;
      %       g2 = kron( ...
      %         opts.mu_ * (P(3 : 3 : 3 * store.dim.num_frames, 4) - 1), ...
      %         [sparse(11, 1); 1]);
      %       g = g + Pr' * g2;
    end
    
  end
  
  for trial = 1 : opts.max_trial
    eval = eval + 1;
    % Compute the update
    if strcmpi(opts.camera_model, 'projective') && opts.nu > 0,
      lambda_mat = lambda * ones(store.dim.num_camera_params, 1);
      lambda_mat([9 10 11]) = opts.nu;
      %       lambda_mat = lambda * ones(store.dim.num_camera_params - 1, 1);
      %       lambda_mat([8 9 10]) = opts.nu;
      lambda_mat = ...
        Pr' * diag(sparse(repmat(lambda_mat, store.dim.num_frames, 1))) * Pr;
      % diag(sparse(repmat(lambda_mat, store.dim.num_frames, 1)));
      H = JTJ + lambda_mat;
      % dw = - (JTJ + lambda_mat) \ g;
      A = H(store.rng.Jq, store.rng.Jq);
      B = H(store.rng.Jp, store.rng.Jq);
      C = H(store.rng.Jp, store.rng.Jp);
      RC = chol(C);
      %       D = inv(C);
      dq = (A - B' * (RC \ (RC' \ B))) \ ...
        (B' * (RC \ (RC' \ g(store.rng.Jp))) - g(store.rng.Jq));
      % dq = (abs(dq) > 1e-6) .* dq;
      dp = - RC \ (RC' \ (g(store.rng.Jp) + B * dq));
      dw(store.rng.Jp) = dp;
      dw(store.rng.Jq) = dq;
      %       dw = - (JTJ + lambda_mat) \ g;
    else
      dw = - (JTJ + lambda * speye(size(JTJ, 1))) \ g;
    end
    
    if ~isempty(strfind(lastwarn, 'NaN'))
      lastwarn('');
      msg = 3;
      return
    end
    
    % disp(rank(full(JTJ)));
    P_eval = P;
    if ~opts.homogeneous,
      P_eval(store.index.affine_camera, :) = ...
        P_eval(store.index.affine_camera, :) + ...
        reshape(dw, 4, 2 * store.dim.num_frames)';
    else
      P_eval = P_eval + reshape(Pr * dw, 4, 3 * store.dim.num_frames)';
      for i = 1 : store.dim.num_frames
        p_eval = reshape(P_eval([3 * i - 2, 3 * i - 1, 3 * i], :)', store.dim.num_camera_params, 1);
        if strcmpi(opts.camera_model, 'affine'),
          a = p_eval([9 10 11]);
          p_eval([9 10 11]) = 0;
          %           p_eval = p_eval / p_eval(12);
          p_eval = p_eval / norm(p_eval);
          p_eval([9 10 11]) = a;
        else
          p_eval = p_eval / norm(p_eval);
          %           p_eval = p_eval / p_eval(12);
        end
        P_eval([3 * i - 2, 3 * i - 1, 3 * i], :) ...
          = reshape(p_eval, 4, 3)';
      end
    end
    [~, X_eval] = optimize_points_given_cameras(...
      M, W, P_eval, X, store, ...
      varargin, ...
      ['homogeneous=', opts.homogeneous], ...
      ['camera_model=', opts.camera_model], ...
      ['mu=', opts.mu], ...
      'display=0');

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
