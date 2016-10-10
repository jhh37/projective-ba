function M_est = je_extrapolate_points(P, X, varargin)
%EXTRAPOLATE_POINTS Summary of this function goes here
%   Detailed explanation goes here

% Custom function.
vec = @(X) X(:);

opts = au_opts( ...
  'epsilon=1.0', ...
  'projection_function=pinhole', ...
  varargin{:});

id_xy = vec(bsxfun(@plus, [0; 1], 1 : 3 : max(size(P))));
PX = P * X;
Z = PX(3 : 3 : size(PX, 1), :);

if strcmpi(opts.projection_function, 'affine')
  Z = ones(size(Z));
elseif strcmpi(opts.projection_function, 'soft')
  Z = opts.epsilon * log(exp(1 / opts.epsilon) + exp(Z / opts.epsilon));
elseif strcmpi(opts.projection_function, 'relu')
  Z = max(Z, ones(size(Z)));
end

M_est = PX(id_xy, :) ./ kron(Z, ones(2, 1));

end

