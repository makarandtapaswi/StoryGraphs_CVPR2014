% This script essentially generates all the results finally used in the
% paper and supplementary material

startup;

%% Character colors:
%%% BBT
bbt_colors = {...
'leonard',      [ 24    88   242] / 255; ...
'sheldon',      [151   242    24] / 255; ...
'penny',        [242    24   215] / 255; ...
'howard',       [242    24    88] / 255; ...
'raj',          [242   142    24] / 255; ...
'kurt',         [121   219   242] / 255; ...
'doug',         [121   219   242] / 255; ...
'leslie',       [121   219   142] / 255; ...
'summer',       [121   219   242] / 255; ...
'mary',         [121   178   242] / 255; ...
'gablehauser',  [159   121   242] / 255};

figure(1); set(gcf, 'Name', 'BBT Character colors');
col = cat(1, bbt_colors{:, 2});
col = shiftdim(col, -1);
imshow(uint8(255*col)); set(gcf, 'Position', [100 700 1600 300], 'Color', 'w');
set(gca, 'Position', [0, 0, 1, 1.5]);
xticklabel_rotate(1:length(col), 90, cellfun(@(x) strrep(x, '_', ' '), bbt_colors(:, 1), 'UniformOutput', false));

%%% BUFFY
buffy_colors = {...
'buffy',        [242     0   238] / 255; ...
'willow',       [  0   242   176] / 255; ...
'xander',       [168   242     0] / 255; ...
'riley',        [242   105     0] / 255; ...
'giles',        [ 34     0   242] / 255; ...
'spike',        [242     0   107] / 255; ...
'anya',         [ 37   242     0] / 255; ...
'tara',         [242    73   148] / 255; ...
'dawn',         [ 73    73    73] / 255; ...
'joyce',        [ 73   242   196] / 255; ...
'dracula',      [ 97    73   242] / 255; ...
'harmony',      [242   121   240] / 255; ...
'mort',         [242   174   121] / 255; ...
'manager',      [242   121   240] / 255; ...
'toth',         [138   121   242] / 255; ...
'xander2',      [168   242     0] / 255; ...
'ben',          [138   121   242] / 255; ...
'graham',       [205   242   121] / 255; ...
'overheiser',   [242   174   121] / 255; ...
'glory',        [130   100   100] / 255; ...
'monk',         [242   201   170] / 255; ...
'watchman',     [242   174   121] / 255; ...
'beth',         [139   242   121] / 255; ...
'maclay',       [138   121   242] / 255; ...
'leiach',       [242   201   170] / 255; ...
'sandy',        [139   242   121] / 255; ...
'donny',        [242   174   121] / 255;};

figure(2); set(gcf, 'Name', 'BUFFY Character colors');
col = cat(1, buffy_colors{:, 2});
col = shiftdim(col, -1);
imshow(uint8(255*col)); set(gcf, 'Position', [100 400 1600 300], 'Color', 'w');
set(gca, 'Position', [0, 0, 1, 1.2]);
xticklabel_rotate(1:length(col), 90, cellfun(@(x) strrep(x, '_', ' '), buffy_colors(:, 1), 'UniformOutput', false));

%%% SAVE FOLDER
save_folder = 'sg/';

%%
close all; drawnow;
runexp = 2;

% expids:
%     1. BBT GT PersonID
%     2. BUFFY GT PersonID

%% 1. BBT -- Annotated PersonID
if any(runexp == 1)
% check for cache
wcross = 1; wprox = 1; wsep = 100; wstraight = 1; wppush = '1/(length(castlist).^2)';
cache_fname = sprintf('%smat/SG_BBT_ANNOT_wcross-%d.wprox-%d.wsep-%d.wstraight-%d.mat', save_folder, wcross, wprox, wsep, wstraight);
if exist(cache_fname, 'file'), load(cache_fname); end
% go through all episodes
for k = 1:6
    % run through each episode 3 times to get different starting positions
    if ~exist('SG_BBT', 'var')  || length(SG_BBT) < k
        SG_BBT(k) = storygraph(BBT(1, k), 'init_method', 'olo_special', 'weights.straight', wstraight, 'weights.sep', wsep, ...
                           'weights.cross', wcross, 'weights.prox', wprox, 'weights.proxpush', wppush, 'inloop', 1);
    end
    % prepare to draw the storygraph and save to file
    figure(100 + k);
    SG = SG_BBT(k);
    colors = [];
    for c = 1:length(SG.castlist)
        colors = [colors; bbt_colors{strcmp(bbt_colors(:, 1), SG.castlist{c}), 2}];
    end
    draw_storygraph(SG.indices, SG.block_times, 'names', SG.castlist, 'presence', SG.presence.orig, ...
                'fig.title', ['Narrative Chart: ' strrep(upper(SG.VS.name), '_', ' ')], 'colors', colors);

    % export figure if toolbox is available
    try
        export_fig(sprintf('%sfig/%s_personid-gt.pdf', save_folder, SG.VS.name));
    catch
        fprintf('Not saving SG as figure. Please install export-fig toolbox\n');
    end
end
% save to mat file if it doesn't exist
if ~exist(cache_fname, 'file'), save(cache_fname, 'SG_BBT', 'bbt_colors'); end
end

%% 2. BUFFY -- Annotated PersonID
if any(runexp == 2)
% check for cache
wcross = 0.1; wprox = 5; wsep = 100; wstraight = 1; wppush = 0.005;
cache_fname = sprintf('%smat/SG_BUFFY_ANNOT_wcross-%d.wprox-%d.wsep-%d.wstraight-%d.mat', save_folder, wcross, wprox, wsep, wstraight);
if exist(cache_fname, 'file'), load(cache_fname); end
% go through all episodes
for k = 1:6
    % if SG and perturbed SG is loaded, then skip. else compute them.
    if ~exist('SG_BUFFY', 'var')  || length(SG_BUFFY) < k
        SG_BUFFY(k) = storygraph(BUFFY(5, k), 'init_method', 'olo_special', 'weights.straight', wstraight, 'weights.sep', wsep, ...
                           'weights.cross', wcross, 'weights.prox', wprox, 'weights.proxpush', wppush, 'inloop', 1);
    end
    % prepare to draw the storygraph and save to file
    figure(200 + k);
    SG = SG_BUFFY(k);
    colors = [];
    for c = 1:length(SG.castlist)
        colors = [colors; buffy_colors{strcmp(buffy_colors(:, 1), SG.castlist{c}), 2}];
    end
    draw_storygraph(SG.indices, SG.block_times, 'names', SG.castlist, 'presence', SG.presence.orig, ...
                'fig.title', ['Narrative Chart: ' strrep(upper(SG.VS.name), '_', ' ')], 'colors', colors);

    % export figure if toolbox is available
    try
        export_fig(sprintf('%sfig/%s_personid-gt.pdf', save_folder, SG.VS.name));
    catch
        fprintf('Not saving SG as figure. Please install export-fig toolbox\n');
    end
end
% save to mat file if it doesn't exist
if ~exist(cache_fname, 'file'), save(cache_fname, 'SG_BUFFY', 'buffy_colors'); end
end
