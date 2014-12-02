function [fun, grad] = softhinge_loss(in, offset, delta)
%SOFTHINGE_LOSS SoftHinge loss breaking at "offset"
%   Differentiable!
%
% Formally known as the -- Pseudo Huber loss function (one-sided)
% http://en.wikipedia.org/wiki/Huber_loss_function
% 
% Input:
%     in: data
%     offset: x-coordinate to start rising, default = 0
%     delta: slope, default = 1
%
% Author: Makarand Tapaswi
% Created: 05-10-2013

%% function in math form
% f(x) = 0, x > offset
%      = d^2 * (sqrt(1 + ((x-offset) / d)^2) - 1)

% g(x) = f'(x)
%      = x ./ sqrt(1 + ((x-offset) / d)^2)

if ~exist('offset', 'var'), offset = 0; end
if ~exist('delta', 'var'), delta = 1; end

%% objective
fun = zeros(size(in));
% 0 part
zeroidx = in > offset;
% other part
fun(~zeroidx) = delta^2 * (sqrt(1 + ((in(~zeroidx)-offset)/delta).^2) - 1);

%% gradient
grad = zeros(size(in));
% 0 part
zeroidx = in > offset;
% other part
grad(~zeroidx) = (in(~zeroidx)-offset) ./ sqrt(1 + ((in(~zeroidx)-offset)/delta).^2);

end
