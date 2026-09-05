function run_single_experiment(userConfig)
%RUN_SINGLE_EXPERIMENT Run one asynchronously coupled EBM/FIS experiment.
%
% Example:
%   run_single_experiment(struct('FCO2',10,'solarConstant',1292));

if nargin < 1
    userConfig = struct();
end

setup_snowball();
cfg = struct( ...
    'atmosphericEmissivity',0.70, ...
    'flowMode','flow', ...
    'startCycle',1, ...
    'solarConstant',1285, ...
    'experimentNumber',1, ...
    'FCO2',35);

if ~isscalar(userConfig) || ~isstruct(userConfig)
    error('run_single_experiment:invalidConfig', ...
        'Configuration must be a scalar structure.');
end

names = fieldnames(userConfig);
for k = 1:numel(names)
    name = names{k};
    if ~isfield(cfg,name)
        error('run_single_experiment:unknownOption', ...
            'Unknown configuration option: %s',name);
    end
    cfg.(name) = userConfig.(name);
end

FISEBM(cfg.atmosphericEmissivity,cfg.flowMode,cfg.startCycle, ...
    cfg.solarConstant,cfg.experimentNumber,cfg.FCO2);
end
