function run_experiment(dataset, varargin)
% RUN_EXPERIMENT

%% OPTIONS
% Extend varargin to include tolerance, maximum number of iterations and
% display.
warning('off', 'MATLAB:nearlySingularMatrix'); % turn off warning for better visual
varargin = [ ...
  { 'tol=1e-9', ...
  'max_iter=500', ...
  'display=0', ...
  'output_filename=0', ...
  }, ...
  varargin, ...
  ];

% Convert varargin to opts.
opts = au_opts(...
  'sample_start=1;sample_end=0;num_samples=100', ...
  'overwrite=0', ...
  'overwrite_first_stage=1', ...
  varargin{:});

% Make the first letter of the dataset name uppercase.
dataset = regexprep(lower(dataset),'(\<[a-z])','${upper($1)}');
if ~opts.output_filename,
  file = ['../Results/', dataset, '/', lower(dataset)];
else
  file = ['../Results/', dataset, '/', lower(opts.output_filename)];
end

% If a previous version of the file exists, load it.
clearvars -global pba_results
if exist([file, '.mat'], 'file') == 2, load(file);
end

% Set pba_results as a global variable.
global pba_results

%% Dataset
if ~opts.sample_end
  opts.sample_end = opts.sample_start + opts.num_samples - 1;
end

data = load(['../Datasets/' dataset '/' lower(dataset)]); % Load the measurement matrix.

% Load dataset info.
store = preprocess_dataset(data.M, data.W);

%% List of algorithms
algorithms = {
  % AHRW1P: Affine homogeneous RW1 with manifold Projection
  'AHRW1P', @(P, X) optimize_cameras_eliminating_points(...
  data.M, data.W, P, X, store, varargin{:}, ...
  'camera_model=affine', 'homogeneous=1', 'rw=1'), '';
  % AHRW2P: Affine homogeneous RW2 with manifold Projection
  'AHRW2P', @(P, X) optimize_cameras_eliminating_points(...
  data.M, data.W, P, X, store, varargin{:}, ...
  'camera_model=affine', 'homogeneous=1', 'rw=2'), '';
  % AIRW2P: Affine Inhomogeneous RW2 with manifold Projection
  'AIRW2P', @(P, X) optimize_cameras_eliminating_points(...
  data.M, data.W, P, X, store, varargin{:}, ...
  'camera_model=affine', 'homogeneous=0', 'rw=2'), '';
  
  % AHRW1: Affine homogeneous RW1 with manifold Projection
  'AHRW1', @(P, X) optimize_cameras_eliminating_points(...
  data.M, data.W, P, X, store, varargin{:}, ...
  'camera_model=affine', 'homogeneous=1', 'rw=1', ...
  'penalize_gauge_freedom=0'), '';
  % AHRW2: Affine homogeneous RW2 with manifold Projection
  'AHRW2', @(P, X) optimize_cameras_eliminating_points(...
  data.M, data.W, P, X, store, varargin{:}, ...
  'camera_model=affine', 'homogeneous=1', 'rw=2', ...
  'penalize_gauge_freedom=0'), '';
  % AIRW2: Affine Inhomogeneous RW2 with manifold Projection
  'AIRW2', @(P, X) optimize_cameras_eliminating_points(...
  data.M, data.W, P, X, store, varargin{:}, ...
  'camera_model=affine', 'homogeneous=0', 'rw=2', ...
  'penalize_gauge_freedom=0'), '';
  
  % PHRW1P: Projective homogeneous RW1 with manifold Projection
  'PHRW1P', @(P, X) optimize_cameras_eliminating_points(...
  data.M, data.W, P, X, store, varargin{:}, ...
  'camera_model=projective', 'homogeneous=1', 'rw=1'), '';
  % PHRW2P: Projective homogeneous RW2 with manifold Projection
  'PHRW2P', @(P, X) optimize_cameras_eliminating_points(...
  data.M, data.W, P, X, store, varargin{:}, ...
  'camera_model=projective', 'homogeneous=1', 'rw=2'), '';
  % PHJP: Projective homogeneous Joint with manifold Projection
  'PHJP', @(P, X) optimize_all_jointly(...
  data.M, data.W, P, X, store, varargin{:}, ...
  'camera_model=projective', 'homogeneous=1', 'rw=2'), '';
  
  % AHPHRW1P: Affine(homogeneous RW2)-then-Projective homogeneous RW1 with
  % manifold Projection
  'AHPHRW1P', @(P, X) optimize_cameras_eliminating_points(...
  data.M, data.W, P, X, store, varargin{:}, ...
  'camera_model=projective', 'homogeneous=1', 'rw=1'), 'AHRW2P'
  % AHPHRW2P: Affine(homogeneous RW2)-then-Projective homogeneous RW2 with
  % manifold Projection
  'AHPHRW2P', @(P, X) optimize_cameras_eliminating_points(...
  data.M, data.W, P, X, store, varargin{:}, ...
  'camera_model=projective', 'homogeneous=1', 'rw=2'), 'AHRW2P'
  % AIPHRW1P: Affine(inhomogeneous RW2)-then-Projective homogeneous RW1 with
  % manifold Projection
  'AIPHRW1P', @(P, X) optimize_cameras_eliminating_points(...
  data.M, data.W, P, X, store, varargin{:}, ...
  'camera_model=projective', 'homogeneous=1', 'rw=1'), 'AIRW2P'
  % AIPHRW2P: Affine(inhomogeneous RW2)-then-Projective homogeneous RW2 with
  % manifold Projection
  'AIPHRW2P', @(P, X) optimize_cameras_eliminating_points(...
  data.M, data.W, P, X, store, varargin{:}, ...
  'camera_model=projective', 'homogeneous=1', 'rw=2'), 'AIRW2P'
  % AHPHJP: Affine(homogeneous RW2)-then-Projective joint with
  % manifold Projection
  'AHPHJP', @(P, X) optimize_all_jointly(...
  data.M, data.W, P, X, store, varargin{:}, ...
  'camera_model=projective', 'homogeneous=1'), 'AHRW2P'
  % AIPHJP: Affine(inhomogeneous RW2)-then-Projective joint with
  % manifold Projection
  'AIPHJP', @(P, X) optimize_all_jointly(...
  data.M, data.W, P, X, store, varargin{:}, ...
  'camera_model=projective', 'homogeneous=1'), 'AIRW2P'
  
  % AHPHRW1P: Affine(homogeneous RW2)-then-Projective homogeneous RW1 with
  % manifold Projection
  'AHPHRW1P_', @(P, X) optimize_cameras_eliminating_points(...
  data.M, data.W, P, X, store, varargin{:}, ...
  'camera_model=projective', 'homogeneous=1', 'rw=1'), 'AHRW2'
  % AHPHRW2P: Affine(homogeneous RW2)-then-Projective homogeneous RW2 with
  % manifold Projection
  'AHPHRW2P_', @(P, X) optimize_cameras_eliminating_points(...
  data.M, data.W, P, X, store, varargin{:}, ...
  'camera_model=projective', 'homogeneous=1', 'rw=2'), 'AHRW2'
  % AIPHRW1P: Affine(inhomogeneous RW2)-then-Projective homogeneous RW1 with
  % manifold Projection
  'AIPHRW1P_', @(P, X) optimize_cameras_eliminating_points(...
  data.M, data.W, P, X, store, varargin{:}, ...
  'camera_model=projective', 'homogeneous=1', 'rw=1'), 'AIRW2'
  % AIPHRW2P: Affine(inhomogeneous RW2)-then-Projective homogeneous RW2 with
  % manifold Projection
  'AIPHRW2P_', @(P, X) optimize_cameras_eliminating_points(...
  data.M, data.W, P, X, store, varargin{:}, ...
  'camera_model=projective', 'homogeneous=1', 'rw=2'), 'AIRW2'
  % AHPHJP: Affine(homogeneous RW2)-then-Projective joint with
  % manifold Projection
  'AHPHJP_', @(P, X) optimize_all_jointly(...
  data.M, data.W, P, X, store, varargin{:}, ...
  'camera_model=projective', 'homogeneous=1'), 'AHRW2'
  % AIPHJP: Affine(inhomogeneous RW2)-then-Projective joint with
  % manifold Projection
  'AIPHJP_', @(P, X) optimize_all_jointly(...
  data.M, data.W, P, X, store, varargin{:}, ...
  'camera_model=projective', 'homogeneous=1'), 'AIRW2'
  
  };

%% EXPERIMENT
for sample_num = opts.sample_start : opts.sample_end
  
  fprintf('Experiments running on [%s] sample [%d] \n', dataset, sample_num);
  pba_results.(dataset)(sample_num).sample_num = sample_num;
  
  % Generate sample (affine inhomogeneous by default)
  [P0, X0] = generate_sample(sample_num, ...
    store.dim.num_frames, ...
    store.dim.num_points, ...
    'affine=1', ...
    'homogeneous=0');
  
  for k = 1 : length(algorithms)
    
    alg = algorithms{k,1};
    alg_caller = algorithms{k,2};
    alg_prev = algorithms{k, 3};
    if isfield(opts, (alg))
      if opts.(alg)
        if ~opts.overwrite && ...
            isfield(pba_results.(dataset)(sample_num), alg) && ...
            isfield(pba_results.(dataset)(sample_num).(alg), 'cost')
          fprintf('Algorithm [%s] already run on [%s] sample [%d]: %.6e\n', ...
            alg, dataset, sample_num, pba_results.(dataset)(sample_num).(alg).cost);
        else
          if ~isempty(alg_prev)
            % Warm start
            % either check for the previous result or run the previous
            % algorithm.
            run_first_stage_algorithm = 1;
            if isfield(pba_results.(dataset)(sample_num), alg_prev)
              if ~isempty(pba_results.(dataset)(sample_num).(alg_prev))
                run_first_stage_algorithm = 0 | opts.overwrite_first_stage;
              end
            end
            
            if run_first_stage_algorithm
              % If the first stage has not been run, run it first.
              alg_prev_idx = strcmpi({algorithms{:, 1}}, alg_prev);
              alg_prev_caller = algorithms{alg_prev_idx, 2};
              
              fprintf('Algorithm [%s] (the 1st stage of [%s]) running on [%s] sample [%d]', alg_prev, alg, dataset, sample_num);
              
              time_start = tic;
              [P, X, cost, iter, msg] = alg_prev_caller(P0, X0);
              runtime = toc(time_start);
              
              pba_results.(dataset)(sample_num).(alg_prev).iters = iter;
              pba_results.(dataset)(sample_num).(alg_prev).cost = cost;
              pba_results.(dataset)(sample_num).(alg_prev).runtime = runtime;
              pba_results.(dataset)(sample_num).(alg_prev).P = P;
              pba_results.(dataset)(sample_num).(alg_prev).X = X;
              pba_results.(dataset)(sample_num).(alg_prev).msg = msg;
              
              
              fprintf(': %.6e, %03d iters, %.3f iters per sec\n', ...
                cost, iter, iter / runtime);
              
              % Check if the folder and/or the file exists.
              if ~exist('../Results', 'dir'),
                mkdir('../Results');
              end
              if ~exist(['../Results/', dataset], 'dir'),
                mkdir(['../Results/', dataset]);
              end
              
              save(file, 'pba_results');
              
            end
            
            % Now, initialize cameras and points for  warm start.
            P_init = pba_results.(dataset)(sample_num).(alg_prev).P;
            X_init = pba_results.(dataset)(sample_num).(alg_prev).X;
          else
            P_init = P0;
            X_init = X0;
          end
          
          fprintf('Algorithm [%s] running on [%s] sample [%d]', alg, dataset, sample_num);
          time_start = tic;
          [P, X, cost, iter, msg] = alg_caller(P_init, X_init);
          runtime = toc(time_start);
          
          if ~isempty(alg_prev)
            % Adjust the number of iterations and runtime for warm-starts.
            iter = iter + ...
              pba_results.(dataset)(sample_num).(alg_prev).iters;
            runtime = runtime + ...
              pba_results.(dataset)(sample_num).(alg_prev).runtime;
          end
          
          pba_results.(dataset)(sample_num).(alg).iters = iter;
          pba_results.(dataset)(sample_num).(alg).cost = cost;
          pba_results.(dataset)(sample_num).(alg).runtime = runtime;
          pba_results.(dataset)(sample_num).(alg).P = P;
          pba_results.(dataset)(sample_num).(alg).X = X;
          pba_results.(dataset)(sample_num).(alg).msg = msg;
          
          fprintf(': %.6e, %03d iters, %.3f iters per sec\n', ...
            cost, iter, iter / runtime);
          
          % Check if the folder and/or the file exists.
          if ~exist('../Results', 'dir'),
            mkdir('../Results');
          end
          if ~exist(['../Results/', dataset], 'dir'),
            mkdir(['../Results/', dataset]);
          end
          
          save(file, 'pba_results');
        end
      end
    end
  end
end

clearvars -global pba_results

end
