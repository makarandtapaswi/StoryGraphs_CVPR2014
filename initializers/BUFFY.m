function VideoStruct = BUFFY(season, episode)
% Initialize STORYGRAPH project for Buffy the Vampire Slayer
% Given season and episode
%

k = 1;
for ep = episode
    video_name = sprintf('buffy_s%02de%02d', season, ep);
    VideoStruct(k) = initVIDEO(video_name, 'series', 'buffy', 'season', season, 'episode', ep);
    k = k + 1;
end

end
