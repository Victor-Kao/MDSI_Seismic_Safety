clear;
clc;
close all;

load("Acc_1G_Z_SSICOV_damp.mat");
load("Acc_1G_Z_SSICOV_freq.mat");

x_axis = [];
y_axis = [];
z_axis = [];

list_a = [10,11];
k = 1;
for i_file = 1:15
    if ismember(i_file, list_a)
        continue
    end
    x_axis = [x_axis , All_F_res{i_file}];
    y_axis = [y_axis , All_F_damp{i_file}];
    z_axis = [z_axis , k*ones(1,length(All_F_res{i_file}))];
    k = k+1;
end

scatter3(x_axis,y_axis,z_axis)