function [f, g] = zeroed_inverse_loss(in, offset, delta)
%ZEROED_INVERSE_LOSS f(x=0) = Inf, and f(x=offset) = 0
%   Yet, differentiable!
%
% Loss function made up by Martin and Makarand :)
% 
% Input:
%     in: data
%     offset: x-coordinate to start rising, default = 0
%     delta: not used, just keep so that it looks similar to the softhinge_loss
%
% Author: Makarand & Martin
% Created: 09-10-2013

%% function in math form
% f(x) = 0, x > offset
%      = 1/x * (sqrt(1 + (x-offset)^2) - 1)
%      = (g(x) - 1) / x
% where, g(x) = sqrt(1 + (x-offset)^2)

% f'(x) = (x-offset)/(x*g(x)) - (g(x)-1)/x^2

if ~exist('offset', 'var'), offset = 1; end

%% objective
common = sqrt(1 + (in - offset).^2); % g(x)
f = zeros(size(in));
% 0 part
zeroidx = in > offset;
% other part
f(~zeroidx) = (common(~zeroidx) - 1) ./ in(~zeroidx);
% f = (common - 1) ./ in;

%% gradient
g = zeros(size(in));
% 0 part
zeroidx = in >= offset;
% other part
% g = (in - offset)./(in .* common) - (common - 1) ./ (in .^2);
% g(~zeroidx) = ((in(~zeroidx) .* (in(~zeroidx) - offset) ./ common(~zeroidx)) - common(~zeroidx) + 1) ./ (in(~zeroidx).^2);
g(~zeroidx) = (in(~zeroidx) - offset) ./ (in(~zeroidx) .* common(~zeroidx)) - (common(~zeroidx) - 1) ./ (in(~zeroidx).^2);

end
