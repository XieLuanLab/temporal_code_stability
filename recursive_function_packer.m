function result = recursive_function_packer(scriptFile, searchDir, copyToDir, doCopy, skipIfExists)
% RECURSIVE_FUNCTION_PACKER
% Tracks custom .m file dependencies recursively, optionally copies them to a package folder.
%
% Inputs:
%   scriptFile     - path to the main .m script to analyze
%   searchDir      - root directory where all your custom functions (.m files) reside
%   copyToDir      - folder where dependencies will be copied
%   doCopy         - true/false: whether to copy files
%   skipIfExists   - true/false: whether to skip files if already exist in copyToDir or subfolders
%
% === CONFIGURATION ===
% scriptFile = 'main_tracking.m';                  % Script to analyze
% searchDir = fullfile('/home/hanlin/Downloads/codeFolderICML_mnist/codeFolderICML');          % Directory to search custom functions
% copyToDir = fullfile(pwd,'tracking'); % Optional: where to copy files
% doCopy = true;                                 % Set to false to skip copying
% skipExist=true;
% result = recursive_function_packer(scriptFile, searchDir, copyToDir, doCopy,skipExist)
% Output:
%   result         - Nx2 cell array: {function name, full path to .m file}

    if nargin < 5
        skipIfExists = false;
    end
    if nargin < 4
        doCopy = true;
    end

    % === Step 1: Map all .m files in the project ===
    mFiles = dir(fullfile(searchDir, '**', '*.m'));
    allFuncMap = containers.Map();
    for i = 1:length(mFiles)
        [~, fname] = fileparts(mFiles(i).name);
        allFuncMap(fname) = fullfile(mFiles(i).folder, mFiles(i).name);
    end

    % === Step 2: Initialize visited map and queue ===
    visited = containers.Map();  % funcName -> fullPath
    queue = {};

    % Add the initial script file
    [~, initName] = fileparts(scriptFile);
    visited(initName) = scriptFile;
    queue{end+1} = scriptFile;

    % === Step 3: Recursively process dependencies ===
    idx = 1;
    pattern = '(?<![\.\w])([a-zA-Z_]\w*)\s*\(';

    while idx <= length(queue)
        currentFile = queue{idx};
        [~, currentName] = fileparts(currentFile);

        if ~isfile(currentFile)
            warning('File not found: %s (skipped)', currentFile);
            idx = idx + 1;
            continue;
        end

        lines = readlines(currentFile);

        for i = 1:length(lines)
            line = strtrim(lines(i));
            if isempty(line) || startsWith(line, '%')
                continue;
            end
            matches = regexp(line, pattern, 'tokens');
            for j = 1:length(matches)
                funcName = matches{j}{1};
                % Skip if already visited or not found
                if isKey(visited, funcName)
                    continue;
                end
                if isKey(allFuncMap, funcName)
                    fullPath = allFuncMap(funcName);
                    visited(funcName) = fullPath;
                    queue{end+1} = fullPath;
                end
            end
        end

        idx = idx + 1;
    end

    % === Step 4: Filter out existing files in copyToDir if requested ===
    if skipIfExists && exist(copyToDir, 'dir')
        existingFiles = dir(fullfile(copyToDir, '**', '*.m'));
        existingNames = erase({existingFiles.name}, '.m');
    else
        existingNames = {};
    end

    filteredNames = {};
    filteredPaths = {};
    funcNames = keys(visited);
    funcPaths = values(visited);

    for i = 1:numel(funcNames)
        fname = funcNames{i};
        fpath = funcPaths{i};
        if skipIfExists && ismember(fname, existingNames)
            continue;
        end
        filteredNames{end+1} = fname;
        filteredPaths{end+1} = fpath;
    end

    % === Step 5: Copy files if enabled ===
    result = [filteredNames(:), filteredPaths(:)];

    if doCopy
        if ~exist(copyToDir, 'dir')
            mkdir(copyToDir);
        end
        for i = 1:size(result, 1)
            destFile = fullfile(copyToDir, [result{i, 1}, '.m']);
            copyfile(result{i, 2}, destFile);
        end
        fprintf('Copied %d new function(s) to %s\n', size(result, 1), copyToDir);
    end

    % === Display Summary ===
    disp('Final Custom Function Dependencies:');
    disp(cell2table(result, 'VariableNames', {'FunctionName', 'FullPath'}));
end
