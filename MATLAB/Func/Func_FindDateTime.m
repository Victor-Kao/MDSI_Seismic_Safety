function DateTimeMat = Func_FindDateTime(csv_files)
    % Initialize a cell array to store datetime information
    DateTimeMat = cell(numel(csv_files), 1);

    for i = 1:length(csv_files)
        % Extract the filename and remove the directory part
        [~, filename, ~] = fileparts(csv_files{i});

        % Split the string at the space character
        splitStr = strsplit(filename, ' ');
        
        % Concatenate the first two parts to get the desired substring
        extractedStr = [splitStr{1}, ' ', splitStr{2}];

        DateTimeMat{i} = extractedStr;   
    end
end