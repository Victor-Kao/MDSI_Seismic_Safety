clear;
clc;
close all;

addpath("D:\MDSI_project\MATLAB\Func");

All_path = Func_FindCsvFiles("D:\MDSI_project\DATA_GM_RawData\DATA_MENHIR_Alldata\Mp1-Bodenplatte\events\2023");

All_path_char = {};
for i = 1:length(All_path)
    All_path_char{i,1} = char(All_path{1,i});
end
save('AllEvent_path_EG.mat', 'All_path_char');



% Initialized 
All_event = cell(1, length(All_path));
%All_event = cell(1, 10);
for i = 1:length(All_path)
    if rem(i,100) ==0 
        disp(i)
    end
    All_event{i} = Func_ImportMenhirData2Tab(All_path{i});
end
save('CellArrayEvent_EG.mat', 'All_event');