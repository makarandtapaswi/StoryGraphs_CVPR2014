function VideoStruct = GOT(season, episode)
% Initialize STORYGRAPH project for Game of Thrones
% Given season and episode
%

k = 1;
for ep = episode
    video_name = sprintf('got_s%02de%02d', season, ep);
    VideoStruct(k) = initVIDEO(video_name, 'series', 'game_of_thrones', 'season', season, 'episode', ep);
    VideoStruct(k).characters = cellfun(@(x) lower(strrep(x, ' ', '_')), VideoStruct(k).characters, 'UniformOutput', false);
    k = k + 1;
end

end
