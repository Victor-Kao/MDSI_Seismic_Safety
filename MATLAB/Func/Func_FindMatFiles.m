function mat_files = Func_FindMatFiles(folder_path)

    addpath("D:\MDSI_project\MATLAB\Func")

    % Get a list of all files and folders in the directory
    files_and_folders = dir(folder_path);

    % Initialize a cell array to store the names of CSV files
    mat_files = {};

    % Iterate through each file and folder
    for i = 1:numel(files_and_folders)
        % Skip '.' and '..'
        if strcmp(files_and_folders(i).name, '.') || strcmp(files_and_folders(i).name, '..')
            continue;
        end
        % Full path to the file or folder
        full_path = fullfile(folder_path, files_and_folders(i).name);

        % Check if it's a directory
        if files_and_folders(i).isdir
            % Recursively search for CSV files in the subfolder
            csv_files_in_subfolder = Func_FindCsvFiles(full_path);

            % Append CSV files found in the subfolder to the main list
            mat_files = [mat_files, csv_files_in_subfolder];
        elseif endsWith(files_and_folders(i).name, '.mat', 'IgnoreCase', true)
            % Add the CSV file to the list
            mat_files{end+1} = full_path;
        end
    end
end