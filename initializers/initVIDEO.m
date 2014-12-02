function VideoStruct = initVIDEO(video_name, varargin)
%INITVIDEO - Returns VideoStructure array for list of input video names
% Usually called by other stuff, although it can also be called directly.
% Examples:
%   video_name = 'bbt_s01e01', 'buffy_s05e01', 'got_s01e01'
%
% length(varargin) == 2, {season, episode}
%
% This code adds or calls the adders of stuff to the video structure
% ADD: video_info, labels, cache

%% Process arguments
default_args.series = '';
default_args.season = [];
default_args.episode = [];
default_args.movie = [];
VideoStruct = cvhci_process_options(varargin, default_args);

%% Add all the information
VideoStruct.name = video_name;

%%% data file name templates
VideoStruct.data.castlist =    ['data/castlist/', video_name, '.cast'];
VideoStruct.data.sg_scenes =   ['data/scene_boundaries/', video_name, '.method-dp.scenes.mat'];
VideoStruct.data.facetracks =  ['data/facetracks/', video_name, '.facetracks.mat'];
VideoStruct.data.personid =    ['data/personid/', video_name, '.personid.mat'];


%% Add list of characters
try
    [VideoStruct.characters, VideoStruct.CharacterStruct] = get_character_names(VideoStruct);
catch
    VideoStruct.characters = {''};
    VideoStruct.CharacterStruct = empty_struct();
    warning('Failed to load list of characters for %s. Setting to empty.', VideoStruct.name);
end

end
