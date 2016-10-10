% Run the setup file prior to performing any experiment.
run setup

% Perform experiment on one algorithm with the following settings.
% [ max_iter = 1000, func_tol = 1e-9, display = 0, sample_end = 1 ]
% We will limit sample_end to 1.
run_experiment('dinosaur_trimmed', 'AIPHRW1P=1', ...
  'max_iter=1000', 'display=1', 'sample_end=1');

% Run multiple algorithms.
% NB: Display is OFF by default.
% - Setting overwrite to 1 makes the script to re-run the experiment for
% that particular seed.
% - NB: overwrite_first_stage is ON by default.
run_experiment('dinosaur_trimmed', 'AIPHRW1P=1', 'AHPHRW1P=1', ...
  'max_iter=1000', 'sample_end=1', 'overwrite=1', 'overwrite_first_stage=1');

% Analyze algorithm performance
% [ csv = 1 will write the csv file to
% Results/<Dataset_name>/<Dataset_name_along with rank>.csv. ]
je_analyze_results('dinosaur_trimmed', 'AIPHRW1P=1', 'AHPHRW1P=1', 'csv=1');


%% LIST OF ALGORITHMS
% Please refer to the paper
% AIRW1P: Affine Inhomogeneous RW1 with manifold Projection
% AIRW2P: Affine Inhomogeneous RW2 with manifold Projection
% AHRW1P: Affine Homogeneous RW1 with manifold Projection
% AHRW2P: Affine Homogeneous RW2 with manifold Projection
% AHJP: Affine Homogeneous Joint with manifold Projection
% PHRW1P: Projective Homogeneous RW1 with manifold Projection
% PHRW2P: Projective Homogeneous RW2 with manifold Projection
% PHJP: Projective Homogeneous Joint with manifold Projection
% AIPHRW1P: AIRW2P -> PHRW1P
% AIPHRW2P: AIRW2P -> PHRW2P
% AHPHRW1P: AHRW2P -> PHRW1P
% AHPHRW2P: AHRW2P -> PHRW2P
