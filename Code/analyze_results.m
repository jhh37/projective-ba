function analyze_results(dataset, varargin)
%ANALYZE_RESULTS Analyzes useful statistics from the results.

% These flags say which methods to run
%% OPTIONS
opts = au_opts(...
  'sample_start=1;sample_end=0;num_samples=100', ...
  'tol=1e-6', ...
  'csv=0', ...
  'display=1', ...
  'output_filename=0', ...
  'folder_suffix=0', ...
  'display_reconstructions=0', ...
  'save_reconstructions=0', ...
  'AIPHRW1P=0', ...
  'AHRW1P=0', ...
  'AHRW2P=0', ...
  'AIRW2P=0', ...
  'AHRW1=0', ...
  'AHRW2=0', ...
  'AIRW2=0', ...
  'PHRW2P=0', ...
  'PHRW1P=0', ...
  'PHJP=0', ...
  'AHPHRW1P=1', ...
  'AIPHRW1P=1', ...
  'AHPHRW2P=0', ...
  'AIPHRW2P=0', ...
  'AHPHJP=1', ...
  'AIPHJP=1', ...
  'AHPHRW1P_=0', ...
  'AIPHRW1P_=0', ...
  'AHPHRW2P_=0', ...
  'AIPHRW2P_=0', ...
  'AHPHJP_=0', ...
  'AIPHJP_=0', ...
  varargin{:});

alg_names = {
  'AHRW1P' ...
  'AHRW2P' ...
  'AIRW2P' ...
  'PHRW2P' ...
  'PHRW1P' ...
  'PHJP' ...
  'AHPHRW1P' ...
  'AIPHRW1P' ...
  'AHPHRW2P' ...
  'AIPHRW2P' ...
  'AHPHJP' ...
  'AIPHJP' ...
  }';

% Make the first letter of the dataset name uppercase.
dataset = regexprep(lower(dataset),'(\<[a-z])','${upper($1)}');

% Add folder suffix
file = '../Results';
if opts.folder_suffix ~= 0, file = [file, '_', opts.folder_suffix];
end

if ~opts.output_filename,
  file = [file, '/', dataset, '/', lower(dataset)];
else
  file = [file, '/', dataset, '/', lower(opts.output_filename)];
end

% If a previous version of the file exists, load it.
if exist([file, '.mat'], 'file') == 2, load(file);
else
  error('The specified results file does not exist.');
end

best_min = nan;
N = size(alg_names, 1);

% Create a list of used algorithms and search for the best optimum.
j = 1;
alg_list.idx = zeros(N, 1);
for i=1:N
  alg_name = alg_names{i,1};
  if isfield(opts, (alg_name))
    if opts.(alg_name)
      
      % If the specified algorithm does not exist, continue with the rest
      % of the for loop.
      if ~isfield([pba_results.(dataset)], (alg_name)), continue
      end
      
      if opts.sample_end == 0,
        alg = [pba_results.(dataset)(opts.sample_start:end).(alg_name)];
      else
        alg = [pba_results.(dataset)(opts.sample_start:opts.sample_end).(alg_name)];
      end
      if isempty(alg)
        % If no information is available, skip to the next algorithm.
        continue
      end
      
      alg_list.idx(j) = i;
      j = j + 1;
      
      if (isnan(best_min)) || (best_min > min([alg.cost]))
        best_min = min([alg.cost]);
      end
    end
  end
end
alg_list.idx = alg_list.idx(alg_list.idx ~= 0);
N = length(alg_list.idx);

% Initialize arrays list.algorithms = algorithms{alg_idx,1};

alg_list.runs.all = nan(N,1);
alg_list.runs.conv = nan(N,1);
alg_list.iters.all.mean = nan(N,1);
alg_list.iters.all.median = nan(N,1);
alg_list.iters.success.mean = nan(N,1);
alg_list.iters.success.median = nan(N,1);
alg_list.iters.fail.mean = nan(N,1);
alg_list.iters.fail.median = nan(N,1);
alg_list.runtime.all.mean = nan(N,1);
alg_list.runtime.all.median = nan(N,1);
alg_list.runtime.success.mean = nan(N,1);
alg_list.runtime.success.median = nan(N,1);
alg_list.runtime.fail.mean = nan(N,1);
alg_list.runtime.fail.median = nan(N,1);
alg_list.ips.all.mean = nan(N,1);
alg_list.ips.all.median = nan(N,1);

for i=1:N
  alg_name = alg_names{alg_list.idx(i),1};
  
  % If the specified algorithm exists, do the calculations. Otherwise,
  % continue with the rest of the for loop.
  if opts.(alg_name) && isfield([pba_results.(dataset)], (alg_name))
    if opts.sample_end == 0
      alg = [pba_results.(dataset)(opts.sample_start:end).(alg_name)];
    else
      alg = [pba_results.(dataset)(opts.sample_start:opts.sample_end).(alg_name)];
    end
    alg_conv = [alg.cost] - best_min < opts.tol  *best_min;
    alg_all_iters = [alg.iters];
    alg_success_iters = alg_all_iters(alg_conv);    alg_fail_iters = alg_all_iters(~alg_conv);
    alg_all_runtime = [alg.runtime];
    alg_success_runtime = alg_all_runtime(alg_conv);
    alg_fail_runtime = alg_all_runtime(~alg_conv);
    
    alg_list.runs.all(i) = sum([alg.cost]>0);
    alg_list.runs.conv(i) = sum(alg_conv);
    alg_list.iters.all.mean(i) = mean(alg_all_iters);
    alg_list.iters.all.median(i) = median(alg_all_iters);
    alg_list.iters.success.mean(i) = mean(alg_success_iters);
    alg_list.iters.success.median(i) = median(alg_success_iters);
    alg_list.iters.fail.mean(i) = mean(alg_fail_iters);
    alg_list.iters.fail.median(i) = median(alg_fail_iters);
    alg_list.runtime.all.mean(i) = mean(alg_all_runtime);
    alg_list.runtime.all.median(i) = median(alg_all_runtime);
    alg_list.runtime.success.mean(i) = mean(alg_success_runtime);
    alg_list.runtime.success.median(i) = median(alg_success_runtime);
    alg_list.runtime.fail.mean(i) = mean(alg_fail_runtime);
    alg_list.runtime.fail.median(i) = median(alg_fail_runtime);
    alg_list.ips.all.mean(i) = mean(alg_all_iters ./ alg_all_runtime);
    alg_list.ips.all.median(i) = median(alg_all_iters ./ alg_all_runtime);
    
    p = alg_list.runs.conv(i) / alg_list.runs.all(i);
    
    % If no success run is found, estimate the lower bound of time to 2
    % successes.
    if p == 0,
      p = 1 / (alg_list.runs.all(i) + 1);
      alg_list.runtime.success.mean(i) = 0;
    end
  end
end

% a = zeros(2, 1);
% a(1) = alg_list.runtime.t2s.mean(i);
% a(2) = alg_list.runtime.t2s.std(i);

if opts.display
  fprintf(['\n--- [ ', dataset,' ] --- \n']);
end

T = table(...
  alg_names(alg_list.idx, 1), ...
  [ alg_list.runs.conv alg_list.runs.all], ...
  alg_list.runs.conv ./ alg_list.runs.all * 100, ...
  alg_list.runtime.all.median, ...
  alg_list.iters.success.median, ...
  alg_list.iters.fail.median, ...
  alg_list.iters.all.median, ...
  alg_list.runtime.success.median, ...
  alg_list.runtime.fail.median, ...
  alg_list.ips.all.median, ...
  alg_list.iters.success.mean, ...
  alg_list.iters.fail.mean, ...
  alg_list.iters.all.mean, ...
  alg_list.runtime.success.mean, ...
  alg_list.runtime.fail.mean, ...
  alg_list.runtime.all.mean, ...
  alg_list.ips.all.mean, ...
  'VariableNames',{'ALGORITHM', 'CONV_RUNS', 'CONV_RATE', 'TIME_MED', 'SI_MED', 'FI_MED', 'ITERS_MED', 'ST_MED', 'FT_MED', 'IPS_MED', 'SI_AVR', 'FI_AVR', 'ITERS_AVR', 'ST_AVR', 'FT_AVR', 'TIME_AVR', 'IPS_AVR'});

if opts.display, disp(T(:,[1:4 7]));
end

if opts.csv
  writetable(T, [file, '.csv']);
end

if opts.display
  fprintf('Best minimum: \t\t%.7e\n', best_min);
  fprintf('Point of convergence: \t%.6e\n\n', best_min + opts.tol * best_min);
end

end
