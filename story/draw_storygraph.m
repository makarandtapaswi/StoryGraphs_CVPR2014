function draw_storygraph(indices, block_times, varargin)
%DRAW_STORYGRAPH Creates the story-graph given which line goes where (indices)
%
%   indices: One character per row, value indicating position in y-axis
%            An index value of 0 means the person is not visible at that time.
%   block_times: [start, end]
%
%   varargin:
%         names:    specify name (string) for each line/character
%         colors:   specify color (r,g,b) for each line/character
%         presence: bool matrix, size(indices), indicates character visible or not
%         offset:   offset, when character appears in block, the line is drawn
%                   from x - offset : x + offset (default, 0.2)
%
%         fig.title:  VideoStruct.name for example (string)
%         fig.xlabel: XLabel (minutes, 30seconds, etc.)
%
%         appear.first:  first appearance (a row of numbers with column index)
%         appear.last:   last appearance (a row of numbers with column index)
%
%         graphics.backlines: (true|false) should i draw gray background lines?
%         graphics.xkcdify:   (true|false) whether to output in XKCD mode
%         graphics.linewidth: (integer) default 3
%         graphics.fontsize:  (integer) default 16
%
%
% Author: Makarand Tapaswi
% Created: 13-09-2013

% close all;

%% Generate defaults, parse arguments
if ~exist('indices', 'var')
    % sample data to test the drawing skills
    indices = [1 1 1 1 1 1 1 2 2 2 2 1 1 1 1 1
               2 2 2 2 2 2 2 1 1 1 1 2 2 2 2 2
               3 3 3 3 3 3 3 3 3 3 3 1 1 1 1 1
               5 5 5 5 5 5 5 4 4 4 4 3 3 3 3 3
               6 6 6 6 6 6 6 5 5 5 5 5 6 6 6 6
               4 4 4 4 4 4 4 6 6 6 6 6 5 5 5 5];
   presence = [0 0 1 1 1 1 1 1 1 1 1 0 0 0 0 0
               1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
               1 1 1 1 1 1 0 0 0 1 1 1 1 1 1 1
               1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1
               0 0 0 0 0 1 1 1 1 0 0 1 1 1 1 1
               0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1];
%     block_times = [000 030; 030 060; 060 090; 090 120; 120 150; 150 180; 180 210; 210 240;
%                    240 270; 270 300; 300 330; 330 360; 360 390; 390 420; 420 450; 450 480];
    % create fake block_times to test drawing
    times = 60*(rand(1, size(indices, 2))+0.1);
    block_times = [cumsum([0, times(1:end-1)])', cumsum(times)'];
    varargin = [varargin, 'presence', presence];
    % a random example
%     indices = rand(10, 25);
end
numcharacter = size(indices, 1);
numcols = size(indices, 2);

%%% Generate colors (uses the golden number color generation scheme for now)
% colors will repeat for the same saturation value, so confine it somehow later
colors = colormap(jet(numcharacter));

%%% Generate names in case they don't exist
names = cell(numcharacter, 1);
for k = 1:numcharacter
    names{k} = sprintf('name%02d', k);
end

%%% Varargin
dopts.colors = colors;
dopts.names = names;
dopts.presence = ones(size(indices));
dopts.offset = 0.25;
% Default figure labels
dopts.fig.title = 'Narrative Chart';
dopts.fig.xlabel = 'Time in Minutes';
dopts.fig.xtickspace = 5;
dopts.fig.gui = false; % is the drawing part of a gui?
% Default graphics
dopts.graphics.backlines = true;
dopts.graphics.xkcdify = false;
dopts.graphics.linewidth = 3;
dopts.graphics.fontsize = 16;
% Default 0 gradients
dopts.gradients = zeros(size(indices));
opts = cvhci_process_options(varargin, dopts);
colors = opts.colors;
% Default first/last appearance
opts.se_presence = zeros(size(opts.presence));
for k = 1:numcharacter
    if all(~opts.presence(k, :))
        opts.appear.first(k) = -1;
        opts.appear.last(k) = -1;
        opts.se_presence(k, :) = 0;
    else
        opts.appear.first(k) = find(opts.presence(k, :) == 1, 1, 'first');
        opts.appear.last(k) = find(opts.presence(k, :) == 1, 1, 'last');
        opts.se_presence(k, opts.appear.first(k):opts.appear.last(k)) = 1;
    end
end

%% Pre-drawing Graphics Options
hold off; cla; hold on;
ylim_whitespace = (max(indices(:)) - min(indices(:))) * 0.05;
minylim = min(min(nonzeros(opts.se_presence .* indices)))-ylim_whitespace;
maxylim = max(max(nonzeros(opts.se_presence .* indices)))+ylim_whitespace;

%%% Vertical lines option
x_coordinates = [block_times(:, 1)', block_times(end, end)];
if opts.graphics.backlines
    for jj = x_coordinates
        line([jj, jj], [minylim, maxylim], 'Color', [0.8, 0.8, 0.8]);
    end
end

%% Prepare line segments (1 for each character, draw later!)
% y-coordinate is indicated by the value in the indices column
% x-coordinate column spacing is now related to timing (in seconds)
off = opts.offset;

%%% Convert block_times to column indices
index_duration = block_times(:, 2) - block_times(:, 1);
xcoord_startend = bsxfun(@plus, [off * index_duration, (1 - off) * index_duration], block_times(:, 1));
xcoord_forward_fillers = [(1-off) * index_duration(1:end-1) + block_times(1:end-1, 1), ...
                              off * index_duration(2:end)   + block_times(1:end-1, 2), ...
                          block_times(1:end-1, 2)];
xcoord_forward_fillers(end+1, :) = [xcoord_startend(end, end), block_times(end, end), NaN];

%%% Create line structures with x, y points mapped for each identity
NewLine = struct('x', [], 'y', [], 'id', [], 'type', []);
for ii = 1:numcharacter         % ii is character color code
    % Lines for each character extend only from the first appearance to last.
    % Everything else is empty and will not be plotted anyways.
    NewLine(end+1).id = ii;
    % skip if the character never appears
    if all(~opts.presence(ii, :)), continue; end
    NewLine(end).x = block_times(opts.appear.first(ii), 1);
    NewLine(end).y = indices(ii, opts.appear.first(ii));
    for jj = opts.appear.first(ii):opts.appear.last(ii) % jj is x-axis location
        % indices(ii, jj) is y-axis position
        NewLine(end).x = [NewLine(end).x, xcoord_startend(jj, :)];
        NewLine(end).y = [NewLine(end).y, indices(ii, jj), indices(ii, jj)];
        NewLine = add_transition_segment(NewLine, indices, ii, jj, xcoord_forward_fillers, opts);
    end
    NewLine(end).x = [NewLine(end).x, block_times(opts.appear.last(ii), 2)];
    NewLine(end).y = [NewLine(end).y, NewLine(end).y(end)];
end
NewLine = NewLine(2:end);

%% Draw lines, first-last appearance, names, etc.
%%% Actually draw the lines
for ii = 1:numcharacter
    % skip if the character never appears
    if all(~opts.presence(ii, :)), continue; end
    for jj = 1:numcols
        if ~opts.se_presence(ii, jj)
            continue; % just skip over. Out of first/last appearance bounds
        end
        [~, start_draw_idx] = min(abs(NewLine(ii).x - block_times(jj, 1)));
        [~, end_draw_idx] = min(abs(NewLine(ii).x - block_times(jj, 2)));
        draw_idx = start_draw_idx:end_draw_idx;
        % choose marker type based on whether person is visible
        switch opts.presence(ii, jj)
            case 0, marker = ':'; linewidth = opts.graphics.linewidth - 2;
            case 1, marker = '-'; linewidth = opts.graphics.linewidth;
        end
        % draw the line
        line(NewLine(ii).x(draw_idx), NewLine(ii).y(draw_idx), 'Color', colors(NewLine(ii).id, :), ...
             'LineStyle', marker, 'LineWidth', linewidth);
    end
end

%%% First / Last Appearance drawings
for k = 1:numcharacter
    % skip if the character never appears
    if all(~opts.presence(k, :)), continue; end
    % first appearance
    plot(block_times(opts.appear.first(k), 1), indices(k, opts.appear.first(k)), 'Marker', '>', ...
         'MarkerSize', 10, 'MarkerFaceColor', colors(k, :), 'MarkerEdgeColor', colors(k, :));
    % last appearance
    plot(block_times(opts.appear.last(k), 2), indices(k, opts.appear.last(k)), 'Marker', 'o', ...
         'MarkerSize', 10, 'MarkerFaceColor', colors(k, :), 'MarkerEdgeColor', colors(k, :));
end

%%% Names
for ii = 1:numcharacter
    % skip if the character never appears
    if all(~opts.presence(ii, :)), continue; end
    format_name = strrep(opts.names{ii}, '_', ' ');
    idx = [1, strfind(format_name, ' ')+1];
    format_name(idx) = upper(format_name(idx));
    % Text at the first appearance
    text(block_times(opts.appear.first(ii), 1)-30, indices(ii, opts.appear.first(ii)), format_name, ...
         'Color', colors(ii, :), 'FontSize', opts.graphics.fontsize, ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    % Text at the last appearance
    text(block_times(opts.appear.last(ii), 2)+30, indices(ii, opts.appear.last(ii)), format_name, ...
         'Color', colors(ii, :), 'FontSize', opts.graphics.fontsize, ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end

%% Draw the Gradients as tiny arrows on the lines
if any(any(opts.gradients))
% if any gradient is nonzero, basically given some gradient draw
for ii = 1:numcharacter
    for jj = 1:numcols
        if opts.se_presence(ii, jj) % if line is present, draw. else ignore
            block_center = mean(block_times(jj, :));
%             arrow([block_center, indices(ii, jj), 10], [block_center, indices(ii, jj) + opts.gradients(ii, jj), 10], ...
%                   'EdgeColor', [0.8, 0.8, 0.8], 'FaceColor', [0.8, 0.8, 0.8]);
            line([block_center, block_center], [indices(ii, jj), indices(ii, jj) + opts.gradients(ii, jj)], ...
                 'Color', colors(ii, :), 'LineWidth', 2);
        end
    end
end
end

%% Decorate
%%% Axes setup
set(gca, 'YLim', [minylim, maxylim]);
set(gca, 'XLim', [x_coordinates(1)-5, x_coordinates(end)+5]); % 5 seconds +/- x-coordinates
set(gca, 'YDir', 'reverse');
set(gca, 'YTick', [], 'YColor', [1 1 1]);
set(gca, 'FontSize', opts.graphics.fontsize-2);
% title(opts.fig.title, 'FontSize', opts.graphics.fontsize);
maxtime = block_times(end,end);
set(gca, 'XTick', 0:opts.fig.xtickspace*60:maxtime, 'XTickLabel', (0:opts.fig.xtickspace*60:maxtime)/60);
xlabel(opts.fig.xlabel, 'FontSize', opts.graphics.fontsize);
if ~opts.fig.gui
    set(gcf, 'Color', 'w');
    set(gcf, 'MenuBar', 'none');
    set(gcf, 'Position', [100, 100, maxtime, 200+(maxylim-minylim)*50]);
end

%%% XKCD effect
if opts.graphics.xkcdify
    warning off;
    xkcdify(gca);
    warning on;
end

end


function [NewLine, draw_balance] = add_transition_segment(NewLine, indices, ii, jj, xcoord_forward_fillers, opts)
% Creates the intermediate transition region between the lines
% If the transition is not flat, inserts a sigmoidal transition
sigmoid_npoints = 20;
if jj < size(indices, 2) && jj < opts.appear.last(ii) % add transition segment
    xc = xcoord_forward_fillers(jj, :);
    if opts.presence(ii, jj) == opts.presence(ii, jj+1)
        %%% character presence is same in current and next scene
        % spread sigmoid transition evenly
        draw_balance = 0;
        xcoords = linspace(xc(1), xc(2), sigmoid_npoints);
        ycoords = mySigmoid(xcoords, [], [], xcoords(1), xcoords(end));
        ycoords = (indices(ii, jj+1) - indices(ii, jj)) * ycoords;
    elseif opts.presence(ii, jj) && ~opts.presence(ii, jj+1)
        %%% character is present in current scene, but not in next scene
        % spread sigmoid transition in next scene
        draw_balance = -1;
        xcoords = [xc(1), linspace(xc(3), xc(2), sigmoid_npoints)];
        ycoords = [0,     mySigmoid(xcoords(2:end), [], [], xcoords(2), xcoords(end))];
        ycoords = (indices(ii, jj+1) - indices(ii, jj)) * ycoords;
    elseif opts.presence(ii, jj+1)
        %%% character is present in next scene, but not in current scene
        % spread sigmoid transition in current scene
        draw_balance = 1;
        xcoords = [linspace(xc(1), xc(3), sigmoid_npoints), xc(2)];
        ycoords = [mySigmoid(xcoords(1:end-1), [], [], xcoords(1), xcoords(end-1)), 1];
        ycoords = (indices(ii, jj+1) - indices(ii, jj)) * ycoords;
    end 
    NewLine(end).x = [NewLine(end).x, xcoords];
    NewLine(end).y = [NewLine(end).y, ycoords+indices(ii, jj)];
end

end

