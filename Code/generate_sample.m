function [P, X] = generate_sample(sample_id, f, n, varargin)
%GENERATE_SAMPLE Generate a set of samples with input dimension.
%

%% Default options
opts = au_opts(...
  'base=1348035934', ...
  'sample_width=10', ...
  'affine=0', ...
  'homogeneous=1', ...
  varargin{:});

%% Generate a random seed
rng(opts.base + sample_id * opts.sample_width, 'twister');

%% Generate the cameras
P = randn(3 * f, 4);
if opts.affine,
  % If call for affine, remove the 3 entries.
  P(3 : 3 : 3*f, 1 : 3) = 0;
end
if opts.homogeneous
  % If not homogeneous, divide by the last element.
  for i = 1 : f
    P([3 * i - 2, 3 * i - 1, 3 * i], :) = ...
      P([3 * i - 2, 3 * i - 1, 3 * i], :) / ...
      norm(P([3 * i - 2, 3 * i - 1, 3 * i], :), 'fro');
  end
else
  for i = 1 : f
    P([3 * i - 2, 3 * i - 1, 3 * i], :) = ...
      P([3 * i - 2, 3 * i - 1, 3 * i], :) / ...
      P(3 * i, 4);
  end
end

%% Generate points
X = randn(4, n);
if ~opts.homogeneous,
  X = X ./ repmat(X(4, :), 4, 1);
end

end
