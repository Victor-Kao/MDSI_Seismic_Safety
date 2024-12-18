clear
clc

addpath("D:\MDSI_project\MATLAB\Func\")


%% Set the number of samples

% PCE testing 1
%type_list = ["Uniform","Uniform","Uniform","Uniform","Uniform",...
%             "Uniform","Uniform","Uniform","Uniform","Uniform",...
%             "Uniform","Uniform","Uniform","Uniform","Uniform",...
%             "Uniform","Uniform","Uniform","Uniform","Uniform",...
%             "Uniform","Uniform","Uniform","Uniform","Uniform",...
%             "Uniform","Uniform","Uniform","Uniform"];

% PCE testing 2
type_list = ["Uniform","Uniform","Uniform","Uniform","Uniform",...
             "Uniform","Uniform","Uniform","Uniform"];


numSamples = 10000;
numVars = length(type_list);


%% DOE
seed = 2;
lhsSamples = lhsdesign(numSamples, numVars);
X = zeros(numSamples, numVars);
rng(seed)
for i = 1:numVars
    X_i = Func_InvTras_generator(lhsSamples(:, i),type_list(i));
    X(:,i) = X_i;
end

output_file_name = sprintf('RNDNUM_building_all_uni_VAR_%d_DOE_%d_SEED_%d.mat',numVars,numSamples,seed);
%output_file_name = sprintf('TEST%d_X_SBAGM_V%d_%s_DOE_%d_DIR_%s.mat',TEST_CASE,numVars,Learning_type,numSamples,dir);
save(output_file_name,'X','-mat');