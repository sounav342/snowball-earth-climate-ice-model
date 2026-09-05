function root = setup_snowball()
%SETUP_SNOWBALL Add the publication release folders to the MATLAB path.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'src'));
addpath(fullfile(root,'src','EBMInput'));
addpath(fullfile(root,'examples'));
addpath(fullfile(root,'tests'));
end
