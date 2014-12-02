function [obj, grad] = storygraph_objective(indices, lin_cooc, presence, objparams, printscores)
%STORYGRAPH_OBJECTIVE Objective function to be minimized to find lines order
% Became too big to keep in the main storygraph file :)
%
% Author: Makarand Tapaswi
% Created: 02-10-2013

% Computes the objective of the StoryGraph optimization
% NOTE: Objective needs to be minimized (f"min"con)

if ~exist('printscores', 'var'), printscores = false; end
[nr, nc] = size(indices);

%% Basic preparations
calculate_gradient = false;
if nargout > 1, calculate_gradient = true; end

diff_se_presence = presence.diff_se;
se_presence = presence.se;
diff_presence = presence.diff;

%%% diff_indices = AX
diff_indices = objparams.col_diff_mat * indices;
nd = size(diff_indices, 1); % nd = nr * (nr-1)/2 % used for normalization

%%% gradient computation pre-prep
% get easy index into character pairs
diff_indexing = nchoosek(1:nr, 2);
char_indexing = zeros(nr, nr);
for k = 1:nr
    ind = find(any(diff_indexing == k, 2));
    char_indexing(:, k) = [ind(1:k-1); size(lin_cooc, 1); ind(k:end)];
end
% append a set of zeros to simplify indexing :)
g_lin_cooc = [lin_cooc; zeros(1, size(lin_cooc, 2))];
g_diff_presence = [diff_presence; zeros(1, size(diff_presence, 2))];
g_diff_se_presence = [diff_se_presence; zeros(1, size(diff_se_presence, 2))];


%% Local Straight Line loss %%
% % Prefer to keep lines straight
% % se presence is multiplied with next time slot for the diff to work correctly in one shot.
% straight_line_loss = sum(sum(se_presence(:, 1:end-1) .* se_presence(:, 2:end) .* diff(indices, 1, 2).^2));
% 
% %%%% NO PRESENCE LOSS
% % straight_line_loss = sum(sum(diff(indices, 1, 2).^2));
% 
% %%%% GRADIENT
% if calculate_gradient
%     % d L /dx_{it} = - 2 (x_{i,t-1} - x_{it}) + 2 (x_{it} - x_{i,t+1})
%     % % % Loop
%     g_straight_line_loop = zeros(size(indices));
%     for i = 1:nr
%         g_straight_line_loop(i, 1) = 2 * (se_presence(i, 1) * se_presence(i, 2)) * (indices(i, 1) - indices(i, 2));
%         for t = 2:(nc-1)
%             g_straight_line_loop(i, t) = -2 * (se_presence(i, t-1) * se_presence(i, t)) * (indices(i, t-1) - indices(i, t)) ...
%                                          +2 * (se_presence(i, t) * se_presence(i, t+1)) * (indices(i, t) - indices(i, t+1));
%         end
%         g_straight_line_loop(i, nc) = -2 * (se_presence(i, nc-1) * se_presence(i, nc)) * (indices(i, nc-1) - indices(i, nc));
%     end
% 
%     % % % Matrix manipulation
%     % g_straight_line_mat = imfilter(indices, [-2, 4, -2]);
%     % g_straight_line_mat(:, 1) = 2 * (indices(:, 1) - indices(:, 2));
%     % g_straight_line_mat(:, end) = 2 * (indices(:, end) - indices(:, end-1));
% end

%% Global Straight Line loss %%
% Prefer to keep lines straight as a whole line considered at once
% basically checks "variance" during the se_presence parts
% to simplify the gradient, we keep out the current value in the mean :)
straight_line_loss = 0;
for i = 1:nr
    this_line{i} = nonzeros(se_presence(i, :) .* indices(i, :))';  % transpose to make it a row-vector
    if length(this_line{i}) == 1, continue; end
    mean_line{i} = remove_one_element_mean(this_line{i});             % compute mean leaving that one element out
    straight_line_loss = straight_line_loss + sum((this_line{i} - mean_line{i}).^2) / length(this_line{i});
end

%%%% GRADIENT
if calculate_gradient
    % d L /dx_{it} = 2 / T (x_{i,t} - mean(x_{i,t'}))
    % % % Loop
    g_straight_line_loop = zeros(size(indices));
    for i = 1:nr
        if length(this_line{i}) == 1, continue; end
        g_straight_line_loop(i, find(se_presence(i, :) == 1, 1, 'first'):find(se_presence(i, :) == 1, 1, 'last')) = 2 * (this_line{i} - mean_line{i});
    end

    % % % Matrix manipulation
    % g_straight_line_mat = imfilter(indices, [-2, 4, -2]);
    % g_straight_line_mat(:, 1) = 2 * (indices(:, 1) - indices(:, 2));
    % g_straight_line_mat(:, end) = 2 * (indices(:, end) - indices(:, end-1));
end


%% Proximity Pull loss %%
% Co-occurrence proximity (Proximity pull)
% brings together lines that are present and co-occur
prox = sum(sum(diff_presence .* (diff_indices .^ 2) .* lin_cooc)); % minimize!

%%%% NO PRESENCE LOSS
% prox = sum(sum((diff_indices .^ 2) .* lin_cooc)); % minimize!

%%%% GRADIENT
if calculate_gradient
    % d L /dx_{it} = 2 (x_{it} - x_{jt}) * c_{ijt}; j = [1 ... n];
    % % % Loop
    g_prox_loop = zeros(size(indices));
    for i = 1:nr
        for t = 1:nc
            this_cooc = g_lin_cooc(char_indexing(:, i), t);
            g_prox_loop(i, t) = sum(2 * g_diff_presence(char_indexing(:, i), t) .* ...
                                   (indices(i, t) - indices(:, t)) .* (this_cooc));
        end
    end

    % % % Matrix manipulation
    % g_prox_mat = zeros(size(indices));
    % for t = 1:nc
    %     indices_alldiff = indices(:, t) * ones(1, nr) - (indices(:, t) * ones(1, nr))';
    %     this_time_cooc = reshape(lin_cooc(char_indexing(:), t), [nr, nr]);
    %     g_prox_mat(:, t) = diag(2 * indices_alldiff * this_time_cooc);
    % end
end


%% Proximity Push loss %%
% we push away lines which have no co-occ
% The loss is a no presence loss. Adding presence to this doesn't make
% sense, since lin_cooc and diffpresence are opposites
proxpush = sum(sum((diff_indices .^ 2) .* ~lin_cooc)); % maximize!

%%%% NO PRESENCE LOSS
% proxpush = sum(sum((diff_indices .^ 2) .* ~lin_cooc)); % maximize!

%%%% GRADIENT
if calculate_gradient
    % d L /dx_{it} = 2 (x_{it} - x_{jt}) * c_{ijt}; \forall c_{ijt} = 0 (else it does not exist)
    % % % Loop
    g_proxpush_loop = zeros(size(indices));
    for i = 1:nr
        for t = 1:nc
            this_cooc = g_lin_cooc(char_indexing(:, i), t);
            g_proxpush_loop(i, t) = sum(2 * (indices(i, t) - indices(:, t)) .* ~(this_cooc));
        end
    end

    % % % Matrix manipulation
    % g_proxpush_mat = zeros(size(indices));
end


%% Minimum Separation loss %%
% Satisfy minimum separation between lines
% we apply this loss only on those lines which are both present at same time
% \forall i, j, i \neq j, x_{i,t} - x_{j,t} >= minsep_val!
minsep_loss = sum(sum(diff_se_presence .* objparams.lossfun.minsep((diff_indices.^2+1e-6), objparams.minsep_val)));

%%%% NO PRESENCE LOSS
% minsep_loss = sum(sum(objparams.lossfun.minsep(diff_indices.^2, objparams.minsep_val)));

%%%% GRADIENT
if calculate_gradient
    % d L /dx_{it} = 2 (x_{it} - x_{jt}) * dh_{ijt}; j = [1 ... n];
    % dh_{ijt} = hubergrad((x_{it} - x_{jt})^2, mu)
    [~, zeroinvgrad] = zeroed_inverse_loss((diff_indices.^2+1e-6), objparams.minsep_val);
    zeroinvgrad = [zeroinvgrad; zeros(1, size(zeroinvgrad, 2))];
    % % % Loop
    g_minsep_zeroinv_loop = zeros(size(indices));
    for i = 1:nr
        for t = 1:nc
            this_zeroinvgrad = zeroinvgrad(char_indexing(:, i), t);
            g_minsep_zeroinv_loop(i, t) = sum(2 * g_diff_se_presence(char_indexing(:, i), t) .* ...
                                           (indices(i, t) - indices(:, t)) .* (this_zeroinvgrad));
            % g_minsep_zeroinv_loop(i, t) = sum(2 * (indices(i, t) - indices(:, t)) .* (this_zeroinvgrad)); % NO PRESENCE
        end
    end

    % % % Matrix manipulations
    % g_minsep_zeroinv_mat = 0;
end


%% Line Crossing loss %%
% Pick the diff_indices at current and next time stamp and multiply them
% Simultaneously, make sure both lines are present at both times to penalize the cross
valid_crossing_scores =   softhinge_loss( diff_indices(:, 1:end-1) .* diff_indices(:, 2:end), 0 ) ...
                                   .* diff_se_presence(:, 1:end-1) .* diff_se_presence(:, 2:end);
crossing_loss = sum(sum(valid_crossing_scores));

%%%% NO PRESENCE LOSS
% valid_crossing_scores = softhinge_loss(diff_indices(:, 1:end-1) .* diff_indices(:, 2:end), 0);
% crossing_loss = sum(sum(valid_crossing_scores));

%%%% GRADIENT
if calculate_gradient
    % d L /dx_{it} = (x_{it-1} - x_{1t-1}) dh_{(1it-1)(1it)} + ... + (x_{it-1} - x_{nt-1}) dh_{(int-1)(int)}
    %               +(x_{it+1} - x_{1t+1}) dh_{(1it+1)(1it)} + ... + (x_{it+1} - x_{nt+1}) dh_{(int+1)(int)}
    % % % Loop
    [~, crossgrad_timediff] = softhinge_loss(diff_indices(:, 1:end-1) .* diff_indices(:, 2:end), 0);
    crossgrad_timediff = crossgrad_timediff .* g_diff_se_presence(1:end-1, 1:end-1) .* g_diff_se_presence(1:end-1, 2:end);
    crossgrad_timediff = [crossgrad_timediff; zeros(1, size(crossgrad_timediff, 2))];
    g_linecross_loop = zeros(size(indices));
    for i = 1:nr
        % t == 1
        hg_tnext = crossgrad_timediff(char_indexing(:, i), 1);
        g_linecross_loop(i, 1) = sum((indices(i, 2) - indices(:, 2)) .* hg_tnext);
        % t = [2, nc-1]
        for t = 2:nc-1
            % prev = t-1, next = t due to the way diff was computed!
            hg_tprev = crossgrad_timediff(char_indexing(:, i), t-1);
            hg_tnext = crossgrad_timediff(char_indexing(:, i), t);
            g_linecross_loop(i, t) = sum((indices(i, t-1) - indices(:, t-1)) .* hg_tprev) + ...
                                     sum((indices(i, t+1) - indices(:, t+1)) .* hg_tnext);
        end
        % t == nc, use col == nc-1
        hg_tprev = crossgrad_timediff(char_indexing(:, i), nc-1);
        g_linecross_loop(i, nc) = sum((indices(i, nc-1) - indices(:, nc-1)) .* hg_tprev);
    end

    % % % Matrix manipulations
    % g_linecross_mat = 0;    
end


%% Combine %%

%%% Objective Function
norm_crossing_loss  = crossing_loss / (nr * (nc - 1));
norm_minsep_loss    = minsep_loss / (nd * nc);
norm_prox           = prox / (nd * nc);
norm_proxpush       = proxpush / (nd * nc);
norm_straight_line  = straight_line_loss / (nd * (nc - 1));
obj = objparams.wcross      * norm_crossing_loss ...
      + objparams.wsep      * norm_minsep_loss ...
      + objparams.wprox     * norm_prox ...
      - objparams.wppush    * norm_proxpush  ... % note the subtraction! :)
      + objparams.wstraight * norm_straight_line;

%%% Gradients
if calculate_gradient
    use_glinecross      = objparams.wcross    * g_linecross_loop / (nr * (nc - 1));
    use_gminsep_zeroinv = objparams.wsep      * g_minsep_zeroinv_loop / (nd * nc);
    use_gprox           = objparams.wprox     * g_prox_loop / (nd * nc);
    use_gproxpush       = objparams.wppush    * g_proxpush_loop / (nd * nc);
    use_gstraight       = objparams.wstraight * g_straight_line_loop / (nd * (nc - 1));
    grad = use_glinecross + use_gminsep_zeroinv + use_gprox - use_gproxpush + use_gstraight;
else
    grad = zeros(size(indices));
end

if isnan(obj), keyboard; end

%% Print %%
if printscores
    fprintf('\tCROSS: %8g | PROX: %8.2g | PPUSH: %8.2g | SEP: %8g | STR8: %8g | TOT: %8.4g\n', ...
             norm_crossing_loss, norm_prox, -norm_proxpush, norm_minsep_loss, norm_straight_line, obj);
end


end


function meanvec = remove_one_element_mean(invec)
% Returns the mean of the rest of the vector for each position
meanvec = zeros(1, length(invec));
for k = 1:length(invec)
    meanvec(k) = mean([invec(1:(k-1)) invec((k+1):end)]);
end
end

