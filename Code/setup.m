function setup()

%% SETUP EXTERNAL LIBRARIES
% AWFUL
try
  au_sparse(int32([1 2 3 4]), int32([1 3 5 7]), randn(4,1));
catch
  disp('Setting up au_sparse...');
  cd('External/awful/matlab/');
  addpath(genpath(pwd));
  mex('-largeArrayDims', 'au_sparse.cxx');
  cd('../../..');
end

disp('Setup completed.');

end
