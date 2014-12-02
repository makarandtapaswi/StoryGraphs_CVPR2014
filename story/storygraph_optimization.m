function [SG, final_fval] = storygraph_optimization(SG, TolFun)
%STORYGRAPH_OPTIMIZATION Actually call fmincon and run the optimization
%
% Author: Makarand Tapaswi
% Created: 23-10-2013

% IMPORTANT: TolFun is being set a little high to prevent excessive iterations.
% Put this to 1e-6 for the final runs. 
if ~exist('TolFun', 'var'), TolFun = 1e-6; end

lin_cooc = SG.lin_cooc;
all_presence = SG.presence;
init_indices = SG.init_indices;
objparams = SG.objparams;
castlist = SG.castlist;
block_times = SG.block_times;

%%% Setup lower/upper bounds
% 1 <= x_{i,t} <= num_cast
num_cast = length(castlist);
lower_bound = ones(size(init_indices));
upper_bound = num_cast * ones(size(init_indices));

%% Create optimization options and call fmincon
if isfield(SG, 'optimopts')
    options = SG.optimopts;
else

% * * * Final Optimization options calling * * * %
% options = optimoptions(@fmincon, 'Algorithm', 'interior-point', 'Display', 'iter', 'GradObj', 'on', 'TolFun', TolFun);

% * * * Options with Gradient checking * * * %
options = optimoptions(@fmincon, 'Algorithm', 'interior-point', 'Display', 'iter', 'TolFun', TolFun, 'MaxFunEvals', 3000, ...
                   'GradObj', 'on', 'GradConstr', 'off', 'DerivativeCheck', 'off', 'FinDiffType', 'central');

% * * * Options using Gradient and drawing StoryGraph at every 5th iteration * * * %
if SG.inputopts.debug
    options = optimoptions(@fmincon, 'Algorithm', 'interior-point', 'Display', 'iter', 'GradObj', 'on', 'TolFun', TolFun, 'MaxFunEvals', 3000, ...
              'GradObj', 'on', 'GradConstr', 'off', 'DerivativeCheck', 'off', 'FinDiffType', 'central', ...
              'OutputFcn', @(x, optimValues, state) outfun(x, optimValues, state, lin_cooc, block_times, objparams, all_presence, castlist));
end

end

% * * * Single run of FMINCON * * * %
% [indices, final_fval] = fmincon(@(indices) storygraph_objective(indices, lin_cooc, all_presence, objparams), init_indices, ...
%                               [], [], [], [], lower_bound, upper_bound, [], options);

% * * * FMINCON with manual perturbation by limiting MaxIter * * * %
for k = 1:SG.inputopts.iter.outer
    if k == SG.inputopts.iter.outer % if in final iteration
        options.MaxIter = 1000;  % reset max iterations to default of 1000
        if strcmp(SG.VS.series, 'game_of_thrones')
            options.MaxIter = 300; % stop early for game of thrones preventing massive collapse
        end
    else
        options.MaxIter = SG.inputopts.iter.inner;  % (optionally) reduce this as outer iter's increase
    end
    %%%% Perform FMINCON
    indices = fmincon(@(indices) storygraph_objective(indices, lin_cooc, all_presence, objparams), init_indices, ...
                      [], [], [], [], lower_bound, upper_bound, [], options);
    %%%% Perform Perturbation
    SG.indices = indices;
    pSG = perturb_sg_indices(SG);
    % prepare initial indices for next round
    init_indices = pSG.indices;
end
[indices, final_fval] = fmincon(@(indices) storygraph_objective(indices, lin_cooc, all_presence, objparams), init_indices, ...
                          [], [], [], [], lower_bound, upper_bound, [], options);

%% Return
SG.indices = indices;
SG.optimopts = options;

end


% Optimization OutputFcn
function stop = outfun(x, optimValues, state, lin_cooc, block_times, objparams, all_presence, castlist)
stop = false;
% draw the things every 5th iteration
if ~mod(optimValues.iteration, 5)
    figure(1); clf;
    draw_storygraph(x, block_times, 'names', castlist, 'presence', all_presence.orig, ...
                    'gradients', reshape(optimValues.gradient, size(x)));
    % figure(2), imagesc(x); colorbar;
    storygraph_objective(x, lin_cooc, all_presence, objparams, 1);
    drawnow;
    % keyboard;
    % pause
end
end
