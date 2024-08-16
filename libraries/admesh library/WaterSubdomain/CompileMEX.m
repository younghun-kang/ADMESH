function CompileMEX

% Find all .c files in admesh library
Subfolders = {
    'admesh library',...
    'InsidePoly',...
    ['topotoolbox-master',filesep,'@FLOWobj',filesep,'private']
    };

c_files = [];
for i = 1 : length(Subfolders)
    c_files = [c_files; dir(['libraries/',Subfolders{i},'/**/*.c'])];
    c_files = c_files(:);
end

% Identify extension of mex file
extension = ['.',mexext];

% Compile mex files over loop
PWD = pwd; % Keep current folder
for i = 1 : length(c_files)

    folder = c_files(i).folder;
    name = c_files(i).name;

    foldername = [folder,'/',name];
    foldername_mex = strrep(foldername,'.c',extension);
    
    if ~exist(foldername_mex,'file')
        cd(folder);
        mex(foldername);
    end
end

% Return to original folder
cd(PWD);