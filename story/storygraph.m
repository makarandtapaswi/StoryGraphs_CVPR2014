function SG = storygraph(VideoStruct, varargin)
%STORYGRAPH Global optimization on ordering of characters over time
%
% Author: Makarand Tapaswi
% Created: 14-09-2013

if ~exist('VideoStruct', 'var'), VideoStruct = BBT(1,1); end
% if ~exist('VideoStruct', 'var'), VideoStruct = BUFFY(5,1); end
% if ~exist('VideoStruct', 'var'), VideoStruct = GOT(1,1); varargin = [varargin, 'facetracks_type', 'pf2']; end
close all;

%% Prepare and load all data
dopts.autopersonid = 0;         % use automatic personid results?
dopts.inloop = false;           % will save final image to cache for comparison
dopts.ignore_tracks = {'false_positive', 'trackswitch', 'unknown'}; % ignore tracks with groundTruthIdentity among these
dopts.facetracks_type = 'pf';   % face tracks type
dopts.xcoord_method = 'scenes'; % xcoordinates are obtained by 'scenes' | 'blocky' time intervals
dopts.blocky.block_time = 60;   % default discrete block time
dopts.blocky.window_size = 1;   % default block-window size for co-occurrence
dopts.scenes.min_scene_dur = 5; % minimum scene duration below which it will be collapsed is 5s
dopts.scenes.window_size = 1;   % default block-window size for co-occurrence
dopts.minsep_val = 0.09;         % minimum separation between lines
dopts.init_method = 'olo_special'; % valid options: 'alphabetical', 'olo', 'olo_special', 'olo_random', 'all_random'
dopts.norm_cooc = 'hysteresis'; % valid options: 'all', 'each', 'hysteresis', 'binarize'
% objective function weights
dopts.weights.straight = 1;     % likes to keep lines as straight (same y-coord) as possible
dopts.weights.cross = 1;        % line crossing penalty
dopts.weights.sep = 1;          % ensures separation between lines
dopts.weights.prox = 1;         % pulls together lines when cooc ~= 0
dopts.weights.proxpush = 0.1;   % pushes lines apart when cooc == 0
dopts.debug = 0;                % debug mode
% optimization iterations
dopts.iter.outer = 5;           % main iterations of [fmincon & perturbation] xN
dopts.iter.inner = 50;          % iterations for fmincon (default 1k)
opts = cvhci_process_options(varargin, dopts);


%%% Get co-occurence matrix
if length(VideoStruct) == 1 % do the simple co-occurrence function calling thing
    switch opts.xcoord_method
        case 'blocky'
            [cooc, presence, castlist] = co_occurrence_blocky(VideoStruct, opts.facetracks_type, false, opts.blocky.block_time, opts.blocky.window_size, opts.ignore_tracks);
            bt = opts.blocky.block_time * ones(size(presence, 2), 1);
            block_times = [cumsum(bt)-opts.blocky.block_time, cumsum(bt)];
        case 'scenes'
                [cooc, presence, castlist, block_times] = co_occurrence_scenes(VideoStruct, ...
                    opts.facetracks_type, false, opts.blocky.window_size, opts.scenes.min_scene_dur, opts.autopersonid, opts.ignore_tracks);
        otherwise
            error('Unknown xcoordinate choosing method!');
    end
    big_block_times = [];
else % pass it to the big guns!
    [cooc, presence, castlist, block_times, big_block_times] = mega_storygraph_data(VideoStruct, opts.facetracks_type, ...
                            false, opts.scenes.window_size, opts.scenes.min_scene_dur, opts.autopersonid, opts.ignore_tracks);
end

%%% Filter character list based on whether the person is never present
remove_characters = find(~any(presence, 2));
if ~isempty(remove_characters)
    castlist(remove_characters) = [];
    presence(remove_characters, :) = [];
    cooc(remove_characters, :, :) = [];
    cooc(:, remove_characters, :) = [];
end


num_cast = length(castlist);
num_blocks = size(cooc, 3);
assert(size(cooc, 1) == num_cast);

%%% Index the upper triangular part of the cooc matrices
% generate linear index
linind = zeros(num_cast);
for k = 1:num_cast, linind(k, 1:k-1) = 1; end
linind = find(linind);
% pick that part from cooc
lin_cooc = zeros(length(linind), num_blocks);
for k = 1:num_blocks
    c = cooc(:, :, k)';
    lin_cooc(:, k) = c(linind);
end

%%% Chart title
if length(VideoStruct) == 1
    chart_title = strrep(upper(VideoStruct.name), '_', '-');
else
    chart_title = sprintf('%s Episodes: %02d - %02d', strrep(upper(VideoStruct(1).series), '_', ' '), ...
                            VideoStruct(1).episode, VideoStruct(end).episode);
end

%% Prepare stuff for optimization
%%% initialize the indices
switch opts.init_method
    case 'alphabetical' % based on simple A-Z cast-list names
        init_indices = (1:num_cast)' * ones(1, num_blocks);
    case 'olo' % based on global co-occurrence
        D = 1./squareform(sum(lin_cooc, 2));
        D(isinf(D)) = 0;
        olo = optimalleaforder(linkage(D), D)';
        [~, order] = sort(olo);
        init_indices = order * ones(1, num_blocks);
        init_indices(init_indices == 1) = 1 + 0.001;
        init_indices(init_indices == num_cast) = num_cast - 0.001;
    case 'olo_special' % special co-occurrence XOR based function
        % Filtering to clean-up intermediate holes
        se = strel('line', round(num_blocks/7), 0);
%         olopresence = imclose(presence, se);
        olopresence = presence;
        dist_mat = zeros(num_cast);
        for ii = 1:num_cast
            for jj = ii+1:length(castlist)
                sx1 = sum(olopresence(ii, :));
                sx2 = sum(olopresence(jj, :));
                x1_and_x2 = olopresence(ii, :) .* olopresence(jj, :);
                x1_xor_x2 = xor(olopresence(ii, :), olopresence(jj, :));
                dist_mat(ii, jj) = (1 - ( sum(x1_and_x2)/sqrt(sx1*sx2) )) + ...
                                   (2*sum(x1_xor_x2)/(sx1+sx2));
            end
        end
        Y = squareform(dist_mat');
        links = linkage(Y, 'single');
%         [~, ~, olo] = dendrogram(links);
        olo = optimalleaforder(links, Y);
        [~, order] = sort(olo);
        init_indices = order' * ones(1, num_blocks);
        init_indices(init_indices == 1) = 1 + 0.001;
        init_indices(init_indices == num_cast) = num_cast - 0.001;
    case 'olo_random'
        init_indices = randperm(num_cast)' * ones(1, num_blocks);
    case 'all_random'
        init_indices = rand(size(presence)) * (num_cast - 1) + 1;
    otherwise
        error('Invalid indices initialization method');
end

if ~opts.inloop
    figure(3);
    draw_storygraph(init_indices, block_times, 'names', castlist, 'presence', presence, ...
                    'fig.title', ['Narrative Chart: ' chart_title]);
end

%%% Normalize the co-occurrence values
switch opts.norm_cooc
    case 'each' % normalize each column to 1
        lin_cooc = bsxfun(@rdivide, lin_cooc, max(1, max(lin_cooc, [], 1)));
    case 'all' % normalize max value to 1
        lin_cooc = lin_cooc / max(max(lin_cooc));
    case 'hysteresis'
        c = numel(lin_cooc) * lin_cooc / sum(sum(lin_cooc));
        c2 = c; c2(c2 > 1) = 1;
        c3 = numel(c2) * c2 / sum(sum(c2));
        lin_cooc = c3;
    case 'binarize'
        lin_cooc = double(logical(lin_cooc));
    otherwise
        warning('Co-occurrence will NOT be normalized.');
end
lin_cooc = lin_cooc / max(lin_cooc(:));

%%% Prepare objective function parameters
objparams.wstraight = opts.weights.straight;
objparams.wcross = opts.weights.cross;
objparams.wsep = opts.weights.sep;
objparams.wprox = opts.weights.prox;
if isstr(opts.weights.proxpush), objparams.wppush = eval(opts.weights.proxpush);
else                             objparams.wppush = opts.weights.proxpush;
end

objparams.lossfun.minsep = @zeroed_inverse_loss;
objparams.lossfun.crossing = @softhinge_loss;

%%% Prepare matrix which will compute column element differences
% NOTE: This matrix will be used for lots of stuff! :)
col_diff_mat = [];
for k = 1:num_cast
    addrows = num_cast - k;
    col_diff_mat = [col_diff_mat; [zeros(addrows, k-1), ones(addrows, 1), diag(-1*ones(addrows, 1))]];
end
objparams.col_diff_mat = col_diff_mat;
objparams.minsep_val = opts.minsep_val;

%%% Create cast pairs list that is used to in the variables with "diff" (diff_indices, diff_presence, etc.)
for k = 1:size(col_diff_mat, 1)
    castpairs(k, :) = castlist(logical(col_diff_mat(k, :)));
end

%%% Create the appropriate presence matrix
se_presence = zeros(size(presence));
for k = 1:num_cast
    se_presence(k, find(presence(k, :), 1, 'first') : find(presence(k, :), 1, 'last')) = 1;
end
diff_se_presence = (logical(col_diff_mat) * se_presence) == 2;
diff_presence = (logical(col_diff_mat) * double(presence)) == 2;
all_presence.diff_se = diff_se_presence;
all_presence.diff = diff_presence;
all_presence.se = se_presence;
all_presence.orig = presence;

%%% Call the function minimization
% Test run with initial objective
init_fval = storygraph_objective(init_indices, lin_cooc, all_presence, objparams, 1);

%%% Collect stuff to return and store
SG.VS = VideoStruct;
SG.castlist = castlist;
SG.castpairs = castpairs;
SG.cooc = cooc;
SG.presence = all_presence;
SG.lin_cooc = lin_cooc;
SG.init_indices = init_indices;
SG.block_times = block_times;
SG.big_block_times = big_block_times;
SG.inputopts = opts;
SG.objparams = objparams;

%%% Run "fmincon"
[SG, final_fval] = storygraph_optimization(SG);
fprintf('Optimization finished. Objective %.2f --> %.2f\n', init_fval, final_fval);
SG.init_fval = init_fval;
SG.final_fval = final_fval;

%% Preview results
if ~opts.inloop
    finalfig = figure(100);
    draw_storygraph(SG.indices, block_times, 'names', castlist, 'presence', presence, ...
                    'fig.title', ['Narrative Chart: ' chart_title]);
end

end

