%% Initialize the StoryGraphs project
% This file will be called automatically on starting Matlab in this directory.

clear all;
clc;
close all;

% Make some directories
if ~exist('tmp/', 'dir'),           mkdir('tmp/'); end
if ~exist('sg/fig', 'dir'),         mkdir('sg/fig'); end
if ~exist('sg/mat', 'dir'),         mkdir('sg/mat'); end

global SGparams

%% Working directory
SGparams.base_dir = [fileparts(mfilename('fullpath')) '/'];
if ispc
    SGparams.base_dir = strrep(base_dir, '\', '/');
end

% Check first initialization
first_init;

%% Repository folders
addpath(genpath('utilities/'));
addpath('losses/');
addpath('story/');
addpath('initializers/');


%% Go go go :)
fprintf('=======================================================\n');
fprintf(['Initialized StoryGraphs repository. Supported videos:', ...
         '\n\t%20s : The Big Bang Theory', ...
         '\n\t%20s : Buffy the Vampire Slayer', ...
         '\n\t%20s : Game of Thrones', ...
     '\n'], ...
     'BBT(se, ep)', 'BUFFY(se, ep)', 'GOT(se, ep)');

