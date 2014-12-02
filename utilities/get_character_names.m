function [names, CharStruct] = get_character_names(VideoStruct, normalize_names)
%GET_CHARACTER_NAMES Provides a list of character names appearing in the series
% Reads from a JSON character file
% OUTPUTS:
%       names: a cell array of names
%       CharStruct: contains any detailed information for character
%           possible options: age, gender, attributes, nicknames, firstname-lastname, etc.
%
% Author: Makarand Tapaswi
% Last modified: 01-02-2013

if ~exist('normalize_names', 'var'), normalize_names = false; end

CharacterStruct = loadjson(VideoStruct.data.castlist);
CharStruct = [CharacterStruct.Characters{:}];
names = {CharStruct(:).name};

%% Convert all "aliases" field to cell arrays even if empty
if isfield(CharStruct, 'aliases')
    for k = 1:length(CharStruct)
        if ~iscell(CharStruct(k).aliases)
            CharStruct(k).aliases = {CharStruct(k).aliases};
        end
    end
end

if normalize_names
    names = cellfun(@(x) strrep(lower(x), ' ', '_'), names, 'UniformOutput', false);
end

end
