% Performs some compilation / messaging operations on first initialization of the repository
global SGparams;

% Create a file in PROJECTROOT/tmp/ to know whether it has been accessed before
tmp_fname = 'tmp/first_init';
if exist(tmp_fname, 'file')
    return;
else
    if ~isdir('tmp'), mkdir('tmp'); end
    fid = fopen(tmp_fname, 'w');
    fprintf(fid, 'Created temporary file on %s\n', date);
    fclose(fid);
end

%% jsonlab
fprintf(2, 'Matlab - JSON interface.\n');
fprintf('http://mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files-in-matlab-octave\n');
fprintf('Please download and unpack the folder to "%s/ext/jsonlab" and make sure the file loadjson.m exists in the folder.\n', SGparams.base_dir);
json_fname = 'ext/jsonlab/loadjson.m';
while ~exist(json_fname, 'file')
    fprintf('Press any key to continue...\n');
    pause;
end
% Check that the loading works
addpath(genpath('ext/jsonlab/'));
try
    loadjson('data/castlist/bbt_s01e01.cast');
    fprintf('Success!\n\n');
catch
    fprintf('Error in reading JSON file. The jsonlab toolbox is not correctly installed.\n');
    delete(tmp_fname);
end


%% export_fig
fprintf(2, 'Export Figures\n');
fprintf('http://mathworks.com/matlabcentral/fileexchange/23629-export-fig\n');
fprintf('Please download and unpack the folder to "%s/ext/export_fig" and make sure the file export_fig.m exists in the folder.\n', SGparams.base_dir);
efig_fname = 'ext/export_fig/export_fig.m';
while ~exist(efig_fname, 'file')
    fprintf('Press any key to continue...\n');
    pause;
end
% Check that a figure can be saved properly
addpath(genpath('ext/export_fig/'));
try
    draw_storygraph;
    export_fig('tmp/check_drawing_and_saving.pdf');
    close all;
    fprintf('Success!\n\n');
catch
    fprintf('Error in saving Figure as pdf file. The export_fig toolbox is not correctly installed.\n');
    delete(tmp_fname);
end


%% xticklabel_rotate
fprintf(2, 'XTick Label Rotate\n');
fprintf('http://mathworks.com/matlabcentral/fileexchange/3486-xticklabel-rotate\n');
fprintf('Please download and unpack the folder to "%s/ext/xticklabel_rotate" and make sure the file xticklabel_rotate.m exists in the folder.\n', SGparams.base_dir);
efig_fname = 'ext/xticklabel_rotate/xticklabel_rotate.m';
while ~exist(efig_fname, 'file')
    fprintf('Press any key to continue...\n');
    pause;
end
% Check that a figure can be saved properly
addpath(genpath('ext/xticklabel_rotate/'));
try
    plot(1:10, 1:10);
    xticklabel_rotate(45);
    close all;
    fprintf('Success!\n\n');
catch
    fprintf('Error in saving XTick label rotate as pdf file. The xticklabel_rotate toolbox is not correctly installed.\n');
    delete(tmp_fname);
end





