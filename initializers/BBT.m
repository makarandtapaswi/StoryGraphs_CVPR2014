function VideoStruct = BBT(season, episode)
% Initialize STORYGRAPH project for The Big Bang Theory
% Given season and episode
%

k = 1;
for ep = episode
    video_name = sprintf('bbt_s%02de%02d', season, ep);
    VideoStruct(k) = initVIDEO(video_name, 'series', 'bbt', 'season', season, 'episode', ep);
    k = k + 1;
end

end
