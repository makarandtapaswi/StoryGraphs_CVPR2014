function SG = perturb_sg_indices(SG)
%PERTURB_SG_INDICES Perturb indices to resolve unnecessary crossings
%
% The optimization is highly nonlinear. This function attempts to solve a part
% of that problem by switching potential crossings if energy is reduced.
%
% Author: Makarand Tapaswi
% Created: 23-10-2013

origSG = SG;
init_obj = storygraph_objective(SG.indices, SG.lin_cooc, SG.presence, SG.objparams, 1);

if SG.inputopts.debug
    figure(11); set(gcf, 'Name', 'Original StoryGraph');
    draw_storygraph(SG.indices, SG.block_times, 'names', SG.castlist, 'presence', SG.presence.orig, ...
         'fig.title', ['Narrative Chart: ' strrep(upper(SG.VS.name), '_', ' ')]);
    posn = get(gcf, 'Position');
    set(gcf, 'Position', [100, 500, posn(3), posn(4)]);
end

%% Iterate (probably)...
bestobj = init_obj;
while true
    % Search for crossing pairs
    candidates = search_crossing_candidates(SG);
    % Run a swapping pass
    [SG, newobj] = run_swapping_pass(SG, candidates, bestobj);
    % Check breaking condition
    if abs(bestobj - newobj) < 1e-10
        break;
    else
        bestobj = newobj;
    end
end

fprintf('Overall improvement in objective: %e\n', init_obj - bestobj);

end


function candidates = search_crossing_candidates(SG)
% Finds the list of crossing positions (points which need to be flipped)

diff_indices = SG.objparams.col_diff_mat * SG.indices;
% crossings_at = diff_indices(:, 1:end-1) .* diff_indices(:, 2:end) ...
%                 .* SG.presence.diff_se(:, 1:end-1) .* SG.presence.diff_se(:, 2:end) < 0;
crossings_at = diff_indices(:, 1:end-1) .* diff_indices(:, 2:end) < 0;

% check whether flipping first column (present characters only) helps
[i, j] = find(SG.presence.se(:, 1) * SG.presence.se(:, 1)');
candidates = [i, j, ones(length(i), 1)];

% add other (t=2 and onwards) crossing locations
for t = 1:size(crossings_at, 2)
    [i, j] = find(triu(squareform(crossings_at(:, t))));
    candidates = [candidates; [i, j, (t+1)*ones(length(i), 1)]];
end

end


function [SG, bestobj] = run_swapping_pass(SG, candidates, bestobj)
% example of swapping using deal -- [b(1), b(2)] = deal(a(2), a(1))
for k = 1:size(candidates, 1)
    SGc = SG;
    [i, j, t] = deal(candidates(k, 1), candidates(k, 2), candidates(k, 3));
    [SG.indices(j, t), SG.indices(i, t)] = deal(SG.indices(i, t), SG.indices(j, t));
    curr_obj = storygraph_objective(SG.indices, SG.lin_cooc, SG.presence, SG.objparams, 0);
    
    if curr_obj >= bestobj % ignore this perturbation
        SG = SGc;
    else % keep this perturbation, draw it and show off :)
        if SG.inputopts.debug
            fprintf('Improved final objective by: %e\n', bestobj - curr_obj);
            figure(12); set(gcf, 'Name', 'Perturbed StoryGraph');
            draw_storygraph(SG.indices, SG.block_times, 'names', SG.castlist, 'presence', SG.presence.orig, ...
                 'fig.title', ['Narrative Chart: ' strrep(upper(SG.VS.name), '_', ' ')]);
            drawnow;
        end
        bestobj = curr_obj;
    end
    
end

end
