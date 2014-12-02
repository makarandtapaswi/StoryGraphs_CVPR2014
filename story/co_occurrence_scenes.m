function [cooc, presence, castlist, block_times, frame_counts] = co_occurrence_scenes(VideoStruct, ...
     facetracks_type, visualize, window_size, merge_scenes_below_x_sec, autopersonid, ignore_tracks)
%CO_OCCURRENCE_SCENES
% Loads the scenes (saved in from VideoEvents using write_out_scenes_for_storygraph.m)
% and does the same thing as co_occurrence_blocky.
%
% This is called by StoryGraph, and shouldn't need to be called separately.
% Anyway, we still provide all necessary default arguments.
%
% When personid_method is not empty, it will try to load the PersonID results
% from cache and use them. Else it will just use ground-truth results.
%       Default: groundtruth.
%
% Author: Makarand Tapaswi
% Created: 16-10-2013

if ~exist('VideoStruct', 'var'), VideoStruct = BUFFY(5,1); end
if ~exist('facetracks_type', 'var'), facetracks_type = 'pf'; end
if ~exist('autopersonid', 'var'), autopersonid = false; end
if ~exist('visualize', 'var'), visualize = false; end
if ~exist('window_size', 'var'), window_size = 1; end % co-occurrence windowed grouping
if ~exist('merge_scenes_below_x_sec', 'var'), merge_scenes_below_x_sec = 5; end
if ~exist('ignore_tracks', 'var'), ignore_tracks = {'false_positive', 'trackswitch', 'unknown'}; end

if autopersonid
    personid_method = 'fmlr-frac20';
else
    personid_method = '';
end

if nargout == 0, visualize = true; end

%% Initial Setup
% Load the Scenes!
cache_fname = VideoStruct.data.sg_scenes;
try
    fprintf('Loading scenes from cache... ');
    load(cache_fname);
    Scenes_all = Scenes;
    clear Scenes;
    fprintf('Success\n');
catch
    fprintf('Failed\n');
    error('Scene mat file does not exist. First run VideoEvents/scenes/write_out_scenes_for_storygraph.m first! Or contact Makarand ;)');
end

% Merge very short scenes into the previous scene to simplify drawing ;)
scene_durations = arrayfun(@(x) x.time.finish - x.time.start, Scenes_all);
small_scenes = scene_durations < merge_scenes_below_x_sec;
Scenes(1) = Scenes_all(1);
for k = 2:length(Scenes_all)
    if small_scenes(k) || (k == 2 && small_scenes(1)) % merge k'th scene into previous
        Scenes(end).shots.finish = Scenes_all(k).shots.finish;
        Scenes(end).frames.finish = Scenes_all(k).frames.finish;
        Scenes(end).time.finish = Scenes_all(k).time.finish;
        Scenes(end).shots.duration = Scenes(end).shots.finish - Scenes(end).shots.start + 1;
        Scenes(end).frames.duration = Scenes(end).frames.finish - Scenes(end).frames.start + 1;
        Scenes(end).time.duration = Scenes(end).time.finish - Scenes(end).time.start;
    else
        Scenes(end+1) = Scenes_all(k);
    end
end
assert(abs(length(Scenes) + sum(small_scenes) - length(Scenes_all)) <= 1);

% Normalize character names (if capital and space separated like in GOT)
castlist = cellfun(@(x) lower(strrep(x, ' ', '_')), VideoStruct.characters, 'UniformOutput', false);
castlist = setdiff(castlist, 'unknown');

% Get tracks and results
if ~isempty(personid_method)
    % Check if they can be loaded directly!
    cache_fname = VideoStruct.data.personid;
    try
        fprintf('Loading PersonID results from cache... ');
        pid = load(cache_fname);
        ft = pid.FaceTracks;
        [ft.assigned] = deal(pid.FaceResults.assigned{:});
        if any(strcmp(unique({ft.assigned}), 'gabelhauser'))
            all_names = {ft.assigned};
            all_names(strcmp(all_names, 'gabelhauser')) = {'gablehauser'};
            [ft.assigned] = deal(all_names{:});
        end
        ft = remove_ignore_tracks(ft, ignore_tracks);
        unkidx = strcmp({ft.assigned}, 'unknown');
        if any(unkidx) % This should not occur with frac20 method now. Unknowns were not trained.
            % WARNING: hack here! all main characters assigned to unknown are fixed
            % manually by putting them back to ground truth. too complicated to do
            % anything else.
            keyboard;
            [ft(unkidx).assigned] = deal(ft(unkidx).groundTruthIdentity);
        end
        % fix castlist issues
        castlist = intersect(castlist, union(unique({ft.groundTruthIdentity}), unique({ft.assigned})));
        fprintf('Success\n');
    catch
        fprintf('Failed\n');
        error('Some error occured here OR PersonID mat file does not exist');
    end
    
else
    load(VideoStruct.data.facetracks);
    % remove tracks whose name is unlisted in the castlist
    remove_names = setdiff(unique({ft.groundTruthIdentity}), castlist);
    cleanup = [];
    for f = 1:length(remove_names)
        cleanup = [cleanup, find(strcmp({ft.groundTruthIdentity}, remove_names{f}))];
    end
    ft(cleanup) = [];
    % remove names which don't have any tracks
    castlist = intersect(castlist, unique({ft.groundTruthIdentity}));
end
num_cast = length(castlist);

% Extract important info from tracks
start_end_times = cell2mat(arrayfun(@(x) [x.timestamps(1), x.timestamps(end)], ft, 'UniformOutput', false)');
track_lengths = arrayfun(@(x) length(x.frames), ft);

if isempty(personid_method) % using ground truth ids
    names = {ft.groundTruthIdentity};
else
    names = {ft.assigned};
end

extra_characters = setdiff(unique(names), castlist);
if ~isempty(extra_characters), keyboard; end

%% Go through first time-blocks
fprintf('Binning facetracks into %d scene blocks\n', length(Scenes));
numblocks = length(Scenes);
block_times = zeros(numblocks, 2);
track_counts = zeros(length(unique(names)), numblocks);
frame_counts = zeros(length(unique(names)), numblocks);
% go through all scenes
for k = 1:numblocks
    block_time_start = Scenes(k).time.start;
    block_time_end = Scenes(k).time.finish;
    block_times(k, :) = [block_time_start, block_time_end];
    tracks_in_interval = start_end_times(:, 1) >= block_time_start & start_end_times(:, 2) <= block_time_end;
    names_in_interval = names(tracks_in_interval);
    tlengths_in_interval = track_lengths(tracks_in_interval);
    n_inblock = unique(names_in_interval);
    for n = 1:length(n_inblock)
        this_name_idx = strcmp(n_inblock{n}, names_in_interval);
        row_idx = find(strcmp(n_inblock{n}, castlist));
        frame_counts(row_idx, k) = sum(tlengths_in_interval(this_name_idx));
        track_counts(row_idx, k) = sum(this_name_idx);
    end
end

if isempty(personid_method) % using ground truth ids
    presence = logical(frame_counts);
else % use the result inferred from cooc_differences_gt_auto_personid.m: gives the 1.4 magic number
%     presence = frame_counts > 100; % stupid way to do it. at least should depend on scene size
    block_time_durations = diff(block_times, 1, 2);
    presence = bsxfun(@gt, frame_counts, 1.4*block_time_durations');
    frame_counts(~presence) = 0;
    track_counts(~presence) = 0;
end

if visualize
%     figure(1);
    subplot(211); imagesc(frame_counts); set(gca, 'YTick', 1:num_cast, 'YTickLabel', castlist); colorbar;
    subplot(212); imagesc(track_counts); set(gca, 'YTick', 1:num_cast, 'YTickLabel', castlist); colorbar;
end
% keyboard;
%% Window'ed co-occurrence
fprintf('Computing co-occurrence on window over blocks of size %d\n', window_size);
cooc_idx = nchoosek(1:num_cast, 2);
cooc_count = zeros(num_cast, numblocks);
cooc = zeros(num_cast, num_cast, numblocks);
winadd = floor(window_size/2);
for k = 1:numblocks
    % get window legnth
    idx = max(1, k-winadd):min(numblocks, k+winadd);
    % compute co-occurrence scores based on frame counts
    cooc_count(:, k) = sum(frame_counts(:, idx), 2);
    cooc_scores = geomean(nchoosek(cooc_count(:, k), 2), 2);
    % insert them into a matrix
    a = zeros(num_cast);
    a(sub2ind(size(a), cooc_idx(:, 1), cooc_idx(:, 2))) = cooc_scores;
    cooc(:, :, k) = a;
end

if visualize
    figure(2);
    for k = 1:size(cooc, 3)
        imagesc(cooc(:, :, k)); set(gca, 'YTick', 1:num_cast, 'YTickLabel', castlist);
        colorbar; title(num2str(k)); pause(0.5);
    end
end

end
