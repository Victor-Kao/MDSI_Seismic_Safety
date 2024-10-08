%Project Statistics
%Robert Gumuchian
%Stand:November 2023

clear
clc
close all

%% Data Import
load Data_Processing\Event_Statistics\Velocity\29112023\13_00\BP_291123_1301.mat
load Data_Processing\Event_Statistics\Velocity\29112023\13_00\OG1_291123_1301.mat
load Data_Processing\Event_Statistics\Velocity\29112023\13_00\OG2_291123_1301.mat



t1=28.51;          %Measurement start
t2=43.671;          %Measurement end
T1=32;              %Sampling start
T2=35;               %Sampling end
delta_t=0.001;      %Timestep of sensor
t=t1:0.001:t2;
T_samp=T1:0.001:T2;
index1=3491;
index2=6491;

%% INPUT FOR DENOISING
denoise_off=0; %0 denoise on, 1 denoise off
flag = 1; % 1 for working on velocity, 2 for working on accelerations          
n1=100; % Decomposition level (for acceleration), e.g.: 100            
n2=120; % Decomposition level (for velocity),e.g.: 120

V_X_BP=table2array(BP_291123_1301(:,2));
V_X_BP_samp(:,1)=T_samp;
V_X_BP_samp(:,2)=V_X_BP(index1:index2);
%[DSW, V_X_BP_samp(:,2), DSWB] = DeNoising(flag,V_X_BP_samp(:,2),1,'sym8',n1,'sym8',n2,delta_t);

V_Y_BP=table2array(BP_291123_1301(:,3));
V_Y_BP_samp(:,1)=T_samp;
V_Y_BP_samp(:,2)=V_Y_BP(index1:index2);
%[DSW, V_Y_BP_samp(:,2), DSWB] = DeNoising(flag,V_Y_BP_samp(:,2),1,'sym8',n1,'sym8',n2,delta_t);

V_Z_BP=table2array(BP_291123_1301(:,4));
V_Z_BP_samp(:,1)=T_samp;
V_Z_BP_samp(:,2)=V_Z_BP(index1:index2);
%[DSW, V_Z_BP_samp(:,2), DSWB] = DeNoising(flag,V_Z_BP_samp(:,2),1,'sym8',n1,'sym8',n2,delta_t);

V_X_OG1=table2array(OG1_291123_1301(:,2));
V_X_OG1_samp(:,1)=T_samp;
V_X_OG1_samp(:,2)=V_X_OG1(index1:index2);
%[DSW, V_X_OG1_samp(:,2), DSWB] = DeNoising(flag,V_X_OG1_samp(:,2),1,'sym8',n1,'sym8',n2,delta_t);

V_Y_OG1=table2array(OG1_291123_1301(:,3));
V_Y_OG1_samp(:,1)=T_samp;
V_Y_OG1_samp(:,2)=V_Y_OG1(index1:index2);
%[DSW, V_Y_OG1_samp(:,2), DSWB] = DeNoising(flag,V_Y_OG1_samp(:,2),1,'sym8',n1,'sym8',n2,delta_t);

V_Z_OG1=table2array(OG1_291123_1301(:,4));
V_Z_OG1_samp(:,1)=T_samp;
V_Z_OG1_samp(:,2)=V_Z_OG1(index1:index2);
%[DSW, V_Z_OG1_samp(:,2), DSWB] = DeNoising(flag,V_Z_OG1_samp(:,2),1,'sym8',n1,'sym8',n2,delta_t);

V_X_OG2=table2array(OG2_291123_1301(:,2));
V_X_OG2_samp(:,1)=T_samp;
V_X_OG2_samp(:,2)=V_X_OG2(index1:index2);
%[DSW, V_X_OG2_samp(:,2), DSWB] = DeNoising(flag,V_X_OG2_samp(:,2),1,'sym8',n1,'sym8',n2,delta_t);

V_Y_OG2=table2array(OG2_291123_1301(:,3));
V_Y_OG2_samp(:,1)=T_samp;
V_Y_OG2_samp(:,2)=V_Y_OG2(index1:index2);
%[DSW, V_Y_OG2_samp(:,2), DSWB] = DeNoising(flag,V_Y_OG2_samp(:,2),1,'sym8',n1,'sym8',n2,delta_t);

V_Z_OG2=table2array(OG2_291123_1301(:,4));
V_Z_OG2_samp(:,1)=T_samp;
V_Z_OG2_samp(:,2)=V_Z_OG2(index1:index2);
%[DSW, V_Z_OG2_samp(:,2), DSWB] = DeNoising(flag,V_Z_OG2_samp(:,2),1,'sym8',n1,'sym8',n2,delta_t);

%% BP measurments seem much smaller so they will be plotted in all directions
% figure
% plot(t,V_X_BP);
% xlabel('Time [s]')
% ylabel('Velocity [mm/s]')
% title('BP only')
% hold on
% plot(t,V_Y_BP);
% plot(t,V_Z_OG2);
% legend('X','Y','Z')
% 
figure
plot(T_samp,V_X_BP_samp(:,2));
xlabel('Time [s]')
ylabel('Velocity [mm/s]')
ylim([-0.5, 0.4])
title('X-Direction Slab Velocity for Event on November 29, 2023 at 13:01')
hold on
plot(T_samp,V_X_OG1_samp(:,2));
plot(T_samp,V_X_OG2_samp(:,2));
legend('BP','1OG Floor', '2OG Floor')


figure
plot(T_samp,V_Y_BP_samp(:,2));
xlabel('Time [s]')
ylabel('Velocity [mm/s]')
ylim([-0.5, 0.4])
title('Y-Direction Slab Velocity for Event on November 29, 2023 at 13:01')
hold on
plot(T_samp,V_Y_OG1_samp(:,2));
plot(T_samp,V_Y_OG2_samp(:,2));
legend('BP','1OG Floor', '2OG Floor')


figure
plot(T_samp,V_Z_BP_samp(:,2));
xlabel('Time [s]')
ylabel('Velocity [mm/s]')
ylim([-0.5, 0.4])
title('Z-Direction Slab Velocity for Event on November 29, 2023 at 13:01')
hold on
plot(T_samp,V_Z_OG1_samp(:,2));
plot(T_samp,V_Z_OG2_samp(:,2));
legend('BP','1OG Floor', '2OG Floor')

%Obtain the displacement values through numerical integration
D_X_BP=cumtrapz(t,V_X_BP);
D_X_BP_samp(:,1)=T_samp;
D_X_BP_samp(:,2)=D_X_BP(index1:index2);

D_Y_BP=cumtrapz(t,V_Y_BP);
D_Y_BP_samp(:,1)=T_samp;
D_Y_BP_samp(:,2)=D_Y_BP(index1:index2);

D_Z_BP=cumtrapz(t,V_Z_BP);
D_Z_BP_samp(:,1)=T_samp;
D_Z_BP_samp(:,2)=D_Z_BP(index1:index2);

D_X_OG1=cumtrapz(t,V_X_OG1);
D_X_OG1_samp(:,1)=T_samp;
D_X_OG1_samp(:,2)=D_X_OG1(index1:index2);

D_Y_OG1=cumtrapz(t,V_Y_OG1);
D_Y_OG1_samp(:,1)=T_samp;
D_Y_OG1_samp(:,2)=D_Y_OG1(index1:index2);

D_Z_OG1=cumtrapz(t,V_Z_OG1);
D_Z_OG1_samp(:,1)=T_samp;
D_Z_OG1_samp(:,2)=D_Z_OG1(index1:index2);

D_X_OG2=cumtrapz(t,V_X_OG2);
D_X_OG2_samp(:,1)=T_samp;
D_X_OG2_samp(:,2)=D_X_OG2(index1:index2);

D_Y_OG2=cumtrapz(t,V_Y_OG2);
D_Y_OG2_samp(:,1)=T_samp;
D_Y_OG2_samp(:,2)=D_X_OG2(index1:index2);

D_Z_OG2=cumtrapz(t,V_Z_OG2);
D_Z_OG2_samp(:,1)=T_samp;
D_Z_OG2_samp(:,2)=D_Z_OG2(index1:index2);


figure
plot(t,D_X_BP);
xlabel('Time [s]')
ylabel('Displacement [mm]')
ylim([-0.015, 0.015])
title('X-Direction Slab Displacement for Event on November 29, 2023 at 13:01')
hold on
plot(t,D_X_OG1);
plot(t,D_X_OG2);
legend('BP','1OG Floor', '2OG Floor')

figure
plot(t,D_Y_BP);
xlabel('Time [s]')
ylabel('Displacement [mm]')
ylim([-0.015, 0.015])
title('Y-Direction Slab Displacement for Event on November 29, 2023 at 13:01')
hold on
plot(t,D_Y_OG1);
plot(t,D_Y_OG2);
legend('BP','1OG Floor', '2OG Floor')

figure
plot(t,D_Z_BP);
xlabel('Time [s]')
ylabel('Displacement [mm]')
ylim([-0.015, 0.015])
title('Z-Direction Slab Displacement for Event on November 29, 2023 at 13:01')
hold on
plot(t,D_Z_OG1);
plot(t,D_Z_OG2);
legend('BP','1OG Floor', '2OG Floor')

% %X-Direction Fourier Analysis ALL
% Fs = 1/delta_t;                                                 % Sampling frequency       
% L_all = length(t);                                     % Length of signal
% t_vec = (0:L_all-1)*delta_t;                                        % Time vector
% X_BP_FFT_all = fft(V_X_BP);                  % fft
% X_OG1_FFT_all = fft(V_X_OG1);
% X_OG2_FFT_all = fft(V_X_OG2);
% 
% X_BP_scale_all = X_BP_FFT_all/Fs;                                 % scaling according to Parsival theorem
% X_OG1_scale_all = X_OG1_FFT_all/Fs;
% X_OG2_scale_all = X_OG2_FFT_all/Fs;
% 
% X_BP_FFT_out_all(:,2) = X_BP_scale_all(1:round(L_all/2+1));                  % single sided
% X_BP_FFT_out_all(2:end-1,2) = 2*X_BP_FFT_out_all(2:end-1,2);      % doubled
% X_BP_FFT_out_all(:,1)=Fs*(0:(L_all/2))/L_all;                         % frequency
% 
% X_OG1_FFT_out_all(:,2) = X_OG1_scale_all(1:round(L_all/2+1));                  % single sided
% X_OG1_FFT_out_all(2:end-1,2) = 2*X_OG1_FFT_out_all(2:end-1,2);      % doubled
% X_OG1_FFT_out_all(:,1)=Fs*(0:(L_all/2))/L_all;                        % frequency
% 
% X_OG2_FFT_out_all(:,2) = X_OG2_scale_all(1:round(L_all/2+1));                  % single sided
% X_OG2_FFT_out_all(2:end-1,2) = 2*X_OG2_FFT_out_all(2:end-1,2);      % doubled
% X_OG2_FFT_out_all(:,1)=Fs*(0:(L_all/2))/L_all;                        % frequency
% 
% 
% figure
% plot(abs(X_BP_FFT_out_all(:,1)),abs(X_BP_FFT_out_all(:,2)));
% xlabel('Frequency [Hz]')
% ylabel('Velocity [mm/s]')
% title('ALL X-Direction BP Frequency Spectrum for Event on November 29, 2023 at 13:01')
% 
% figure
% plot(abs(X_OG1_FFT_out_all(:,1)),abs(X_OG1_FFT_out_all(:,2)));
% xlabel('Frequency [Hz]')
% ylabel('Velocity [mm/s]')
% title('ALL X-Direction 1.OG Frequency Spectrum for Event on November 29, 2023 at 13:01')
% 
% figure
% plot(abs(X_OG2_FFT_out_all(:,1)),abs(X_OG2_FFT_out_all(:,2)));
% xlabel('Frequency [Hz]')
% ylabel('Velocity [mm/s]')
% title('ALL X-Direction 2.OG Frequency Spectrum for Event on November 29, 2023 at 13:01')


%%X-Direction Fourier Analysis Sample Area, Velocity and Displacement
Fs = 1/delta_t;                                         % Sampling frequency       
L = length(T_samp);                                     % Length of signal
t_vec = (0:L-1)*delta_t;                                % Time vector

X_BP_FFT_V = fft(V_X_BP_samp(:,2));                  % fft
X_OG1_FFT_V = fft(V_X_OG1_samp(:,2));
X_OG2_FFT_V = fft(V_X_OG2_samp(:,2));
X_BP_FFT_D = fft(D_X_BP_samp(:,2));                  % fft
X_OG1_FFT_D = fft(D_X_OG1_samp(:,2));
X_OG2_FFT_D = fft(D_X_OG2_samp(:,2));

X_BP_V_scale = X_BP_FFT_V/Fs;                            % scaling according to Parsival theorem
X_OG1_V_scale = X_OG1_FFT_V/Fs;
X_OG2_V_scale = X_OG2_FFT_V/Fs;
X_BP_D_scale = X_BP_FFT_D/Fs;                            % scaling according to Parsival theorem
X_OG1_D_scale = X_OG1_FFT_D/Fs;
X_OG2_D_scale = X_OG2_FFT_D/Fs;

X_BP_FFT_V_out(:,2) = X_BP_V_scale(1:round(L/2));             % single sided
X_BP_FFT_V_out(2:end-1,2) = 2*X_BP_FFT_V_out(2:end-1,2);      % doubled
X_BP_FFT_V_out(:,1)=Fs*(0:(L/2))/L;                         % frequency
X_BP_FFT_D_out(:,2) = X_BP_D_scale(1:round(L/2));             % single sided
X_BP_FFT_D_out(2:end-1,2) = 2*X_BP_FFT_D_out(2:end-1,2);      % doubled
X_BP_FFT_D_out(:,1)=Fs*(0:(L/2))/L;                         % frequency

X_OG1_FFT_V_out(:,2) = X_OG1_V_scale(1:round(L/2));           % single sided
X_OG1_FFT_V_out(2:end-1,2) = 2*X_OG1_FFT_V_out(2:end-1,2);    % doubled
X_OG1_FFT_V_out(:,1)=Fs*(0:(L/2))/L;                        % frequency
X_OG1_FFT_D_out(:,2) = X_OG1_D_scale(1:round(L/2));           % single sided
X_OG1_FFT_D_out(2:end-1,2) = 2*X_OG1_FFT_D_out(2:end-1,2);    % doubled
X_OG1_FFT_D_out(:,1)=Fs*(0:(L/2))/L;                        % frequency

X_OG2_FFT_V_out(:,2) = X_OG2_V_scale(1:round(L/2));           % single sided
X_OG2_FFT_V_out(2:end-1,2) = 2*X_OG2_FFT_V_out(2:end-1,2);    % doubled
X_OG2_FFT_V_out(:,1)=Fs*(0:(L/2))/L;                        % frequency
X_OG2_FFT_D_out(:,2) = X_OG2_D_scale(1:round(L/2));           % single sided
X_OG2_FFT_D_out(2:end-1,2) = 2*X_OG2_FFT_D_out(2:end-1,2);    % doubled
X_OG2_FFT_D_out(:,1)=Fs*(0:(L/2))/L;                        % frequency


%%Y-Direction Fourier Analysis, Velocity and Displacement
Y_BP_FFT_V = fft(V_Y_BP_samp(:,2));                  % fft
Y_OG1_FFT_V = fft(V_Y_OG1_samp(:,2));
Y_OG2_FFT_V = fft(V_Y_OG2_samp(:,2));
Y_BP_FFT_D = fft(D_Y_BP_samp(:,2));                  % fft
Y_OG1_FFT_D = fft(D_Y_OG1_samp(:,2));
Y_OG2_FFT_D = fft(D_Y_OG2_samp(:,2));

Y_BP_scale_V = Y_BP_FFT_V/Fs;                                 % scaling according to Parsival theorem
Y_OG1_scale_V = Y_OG1_FFT_V/Fs;
Y_OG2_scale_V = Y_OG2_FFT_V/Fs;
Y_BP_scale_D = Y_BP_FFT_D/Fs;                                 % scaling according to Parsival theorem
Y_OG1_scale_D = Y_OG1_FFT_D/Fs;
Y_OG2_scale_D = Y_OG2_FFT_D/Fs;

Y_BP_FFT_V_out(:,2) = Y_BP_scale_V(1:round(L/2));                  % single sided
Y_BP_FFT_V_out(2:end-1,2) = 2*Y_BP_FFT_V_out(2:end-1,2);      % doubled
Y_BP_FFT_V_out(:,1)=Fs*(0:(L/2))/L;                         % frequency
Y_BP_FFT_D_out(:,2) = Y_BP_scale_D(1:round(L/2));                  % single sided
Y_BP_FFT_D_out(2:end-1,2) = 2*Y_BP_FFT_D_out(2:end-1,2);      % doubled
Y_BP_FFT_D_out(:,1)=Fs*(0:(L/2))/L;                         % frequency

Y_OG1_FFT_V_out(:,2) = Y_OG1_scale_V(1:round(L/2));                  % single sided
Y_OG1_FFT_V_out(2:end-1,2) = 2*Y_OG1_FFT_V_out(2:end-1,2);      % doubled
Y_OG1_FFT_V_out(:,1)=Fs*(0:(L/2))/L;                        % frequency
Y_OG1_FFT_D_out(:,2) = Y_OG1_scale_D(1:round(L/2));                  % single sided
Y_OG1_FFT_D_out(2:end-1,2) = 2*Y_OG1_FFT_D_out(2:end-1,2);      % doubled
Y_OG1_FFT_D_out(:,1)=Fs*(0:(L/2))/L;                        % frequency

Y_OG2_FFT_V_out(:,2) = Y_OG2_scale_V(1:round(L/2));                  % single sided
Y_OG2_FFT_V_out(2:end-1,2) = 2*Y_OG2_FFT_V_out(2:end-1,2);      % doubled
Y_OG2_FFT_V_out(:,1)=Fs*(0:(L/2))/L;                        % frequency
Y_OG2_FFT_D_out(:,2) = Y_OG2_scale_D(1:round(L/2));                  % single sided
Y_OG2_FFT_D_out(2:end-1,2) = 2*Y_OG2_FFT_D_out(2:end-1,2);      % doubled
Y_OG2_FFT_D_out(:,1)=Fs*(0:(L/2))/L;                        % frequency

%%Z-Direction Fourier Analysis
Z_BP_FFT_V = fft(V_Z_BP_samp(:,2));                  % fft
Z_OG1_FFT_V = fft(V_Z_OG1_samp(:,2));
Z_OG2_FFT_V = fft(V_Z_OG2_samp(:,2));
Z_BP_FFT_D = fft(D_Z_BP_samp(:,2));                  % fft
Z_OG1_FFT_D = fft(D_Z_OG1_samp(:,2));
Z_OG2_FFT_D = fft(D_Z_OG2_samp(:,2));

Z_BP_scale_V = Z_BP_FFT_V/Fs;                                 % scaling according to Parsival theorem
Z_OG1_scale_V = Z_OG1_FFT_V/Fs;
Z_OG2_scale_V = Z_OG2_FFT_V/Fs;
Z_BP_scale_D = Z_BP_FFT_D/Fs;                                 % scaling according to Parsival theorem
Z_OG1_scale_D = Z_OG1_FFT_D/Fs;
Z_OG2_scale_D = Z_OG2_FFT_D/Fs;

Z_BP_FFT_V_out(:,2) = Z_BP_scale_V(1:round(L/2));             % single sided
Z_BP_FFT_V_out(2:end-1,2) = 2*Z_BP_FFT_V_out(2:end-1,2);      % doubled
Z_BP_FFT_V_out(:,1)=Fs*(0:(L/2))/L;                         % frequency
Z_BP_FFT_D_out(:,2) = Z_BP_scale_D(1:round(L/2));             % single sided
Z_BP_FFT_D_out(2:end-1,2) = 2*Z_BP_FFT_D_out(2:end-1,2);      % doubled
Z_BP_FFT_D_out(:,1)=Fs*(0:(L/2))/L;                         % frequency

Z_OG1_FFT_V_out(:,2) = Z_OG1_scale_V(1:round(L/2));           % single sided
Z_OG1_FFT_V_out(2:end-1,2) = 2*Z_OG1_FFT_V_out(2:end-1,2);    % doubled
Z_OG1_FFT_V_out(:,1)=Fs*(0:(L/2))/L;                        % frequency
Z_OG1_FFT_D_out(:,2) = Z_OG1_scale_D(1:round(L/2));           % single sided
Z_OG1_FFT_D_out(2:end-1,2) = 2*Z_OG1_FFT_D_out(2:end-1,2);    % doubled
Z_OG1_FFT_D_out(:,1)=Fs*(0:(L/2))/L;                        % frequency

Z_OG2_FFT_V_out(:,2) = Z_OG2_scale_V(1:round(L/2));           % single sided
Z_OG2_FFT_V_out(2:end-1,2) = 2*Z_OG2_FFT_V_out(2:end-1,2);    % doubled
Z_OG2_FFT_V_out(:,1)=Fs*(0:(L/2))/L;                        % frequency
Z_OG2_FFT_D_out(:,2) = Z_OG2_scale_D(1:round(L/2));           % single sided
Z_OG2_FFT_D_out(2:end-1,2) = 2*Z_OG2_FFT_D_out(2:end-1,2);    % doubled
Z_OG2_FFT_D_out(:,1)=Fs*(0:(L/2))/L;                        % frequency

figure
plot(abs(X_BP_FFT_D_out(:,1)),abs(X_BP_FFT_D_out(:,2)));
xlabel('Frequency [Hz]')
xlim([0 30])
ylabel('Displacement [mm]')
title('BP Frequency Spectrum (Disp FFT) for Event on 29/11/2023 at 13:01')
hold on
plot(abs(Y_BP_FFT_D_out(:,1)),abs(Y_BP_FFT_D_out(:,2)));
plot(abs(Z_BP_FFT_D_out(:,1)),abs(Z_BP_FFT_D_out(:,2)));
legend('X','Y','Z')

figure
plot(abs(X_BP_FFT_V_out(:,1)),abs(X_BP_FFT_V_out(:,2)));
xlabel('Frequency [Hz]')
xlim([0 30])
ylabel('Velocity [mm/s]')
title('BP Frequency Spectrum for Event on November 29, 2023 at 13:01')
hold on
plot(abs(Y_BP_FFT_V_out(:,1)),abs(Y_BP_FFT_V_out(:,2)));
plot(abs(Z_BP_FFT_V_out(:,1)),abs(Z_BP_FFT_V_out(:,2)));
legend('X','Y','Z')

figure
plot(abs(X_OG1_FFT_V_out(:,1)),abs(X_OG1_FFT_V_out(:,2)));
xlabel('Frequency [Hz]')
xlim([0 30])
ylabel('Velocity [mm/s]')
title('1.OG Frequency Spectrum for Event on November 29, 2023 at 13:01')
hold on
plot(abs(Y_OG1_FFT_V_out(:,1)),abs(Y_OG1_FFT_V_out(:,2)));
plot(abs(Z_OG1_FFT_V_out(:,1)),abs(Z_OG1_FFT_V_out(:,2)));
legend('X','Y','Z')

figure
plot(abs(X_OG2_FFT_V_out(:,1)),abs(X_OG2_FFT_V_out(:,2)));
xlabel('Frequency [Hz]')
xlim([0 30])
ylabel('Velocity [mm/s]')
title('2.OG Frequency Spectrum for Event on November 29, 2023 at 13:01')
hold on
plot(abs(Y_OG2_FFT_V_out(:,1)),abs(Y_OG2_FFT_V_out(:,2)));
plot(abs(Z_OG2_FFT_V_out(:,1)),abs(Z_OG2_FFT_V_out(:,2)));
legend('X','Y','Z')

%Transfer Functions
h_OG1_BP_X=X_OG1_FFT_V_out(:,2)./X_BP_FFT_V_out(:,2);
h_OG1_BP_Y=Y_OG1_FFT_V_out(:,2)./Y_BP_FFT_V_out(:,2);
h_OG1_BP_Z=Z_OG1_FFT_V_out(:,2)./Z_BP_FFT_V_out(:,2);

h_OG2_BP_X=X_OG2_FFT_V_out(:,2)./X_BP_FFT_V_out(:,2);
h_OG2_BP_Y=Y_OG2_FFT_V_out(:,2)./Y_BP_FFT_V_out(:,2);
h_OG2_BP_Z=Z_OG2_FFT_V_out(:,2)./Z_BP_FFT_V_out(:,2);

h_OG2_OG1_X=X_OG2_FFT_V_out(:,2)./X_OG1_FFT_V_out(:,2);
h_OG2_OG1_Y=Y_OG2_FFT_V_out(:,2)./Y_OG1_FFT_V_out(:,2);
h_OG2_OG1_Z=Z_OG2_FFT_V_out(:,2)./Z_OG1_FFT_V_out(:,2);

figure
plot(abs(X_BP_FFT_V_out(:,1)),abs(h_OG1_BP_X));
xlabel('Frequency [Hz]')
xlim([0 30])
ylabel('h(omega) [-]')
title('Transfer Function OG1 -> BP')
hold on
plot(abs(X_BP_FFT_V_out(:,1)),abs(h_OG1_BP_Y));
plot(abs(X_BP_FFT_V_out(:,1)),abs(h_OG1_BP_Z));
legend('X','Y','Z')

figure
plot(abs(X_BP_FFT_V_out(:,1)),abs(h_OG2_BP_X));
xlabel('Frequency [Hz]')
xlim([0 30])
ylabel('h(omega) [-]')
title('Transfer Function OG2 -> BP')
hold on
plot(abs(X_BP_FFT_V_out(:,1)),abs(h_OG2_BP_Y));
plot(abs(X_BP_FFT_V_out(:,1)),abs(h_OG2_BP_Z));
legend('X','Y','Z')

figure
plot(abs(X_BP_FFT_V_out(:,1)),abs(h_OG2_OG1_X));
xlabel('Frequency [Hz]')
xlim([0 30])
ylabel('h(omega) [-]')
title('Transfer Function OG2 -> OG1')
hold on
plot(abs(X_BP_FFT_V_out(:,1)),abs(h_OG2_OG1_Y));
plot(abs(X_BP_FFT_V_out(:,1)),abs(h_OG2_OG1_Z));
legend('X','Y','Z')


%Denoising function
function [ DSW, DSWBF, DSWB ] = DeNoising( flag,NS,type,w1,n1,w2,n2,dt )
    % This program will denois the noisy signal based on the Ansari method %
    %                                                                      %
    % Input:                                                               %
    %       flag: 1 for working on velocity, 2 for working (1 is suggested)%
    %       NS: Noisy Signal (Velocity or Acceleration)                    %
    %       type: 1 for velocities or 2 for accelarations                  %
    %       dt: Time step of signal                                        %
    %       w1: wavelet function (for acceleration): 'sym8' is suggested   %
    %       w2: wavelet function (for velocity): 'sym8' is suggested       %
    %       n1: Decomposition level (for acceleration), e.g.: 100          %
    %       n2: Decomposition level (for velocity),e.g.: 120               %
    %       v0: initial velocity                                           %
    % Output:                                                              %
    %       DSW: Denoised Signal (only wavelet implemented)                %
    %       DSWBF: Denoised Signal (Wavelet+BaselineCorrection+Filtering   %
    %              implemented)                                            %
    %       DSWB: Denoised Signal (Wavelet+BaselineCorrection implemented) %
    %----------------------------------------------------------------------%    
    %% Wavelet+Baseline
    [ fc,~ ] = Low_Cut_Freq( NS,dt );
    fc=max(fc,0.35);
    if flag==1
        if type==1
            wname = w2; lev = n2;
            [c,l] = wavedec(NS,lev,wname);
            THR2 = wthrmngr('dw1ddenoLVL','rigrsure',c,l,'mln'); 
            [~,CXC7,LXC7,~,~] = wdencmp('lvd',c,l,wname,lev,THR2,'s');
            velo2=0;
            for ii=1:lev
                dd = wrcoef('d',CXC7,LXC7,wname,ii);
                velo2=velo2+dd;
            end
            acce2=[0;diff(velo2)]./dt; 
        elseif type==2
            wname = w1; lev = n1;
            [c,l] = wavedec(NS,lev,wname);
            THR2 = wthrmngr('dw1ddenoLVL','rigrsure',c,l,'mln');
            [~,CXC7,LXC7,~,~] = wdencmp('lvd',c,l,wname,lev,THR2,'s');
            acce_mm1=waverec(CXC7,LXC7,wname);
            velo_mm1=Inte(acce_mm1,dt);
            wname = w2; lev = n2;
            [c,l] = wavedec(velo_mm1,lev,wname);
            THR2 = wthrmngr('dw1ddenoLVL','rigrsure',c,l,'mln'); 
            [~,CXC7,LXC7,~,~] = wdencmp('lvd',c,l,wname,lev,THR2,'s');
            velo2=0;
            for ii=1:lev
                dd = wrcoef('d',CXC7,LXC7,wname,ii);
                velo2=velo2+dd;
            end
            acce2=[0;diff(velo2)]./dt; 
        else 
            warning('Noisy signal type not recognized')
        end
    elseif flag==2
        if type==1
            wname = w2; lev = n2;
            [c,l] = wavedec(NS,lev,wname);
            THR2 = wthrmngr('dw1ddenoLVL','rigrsure',c,l,'mln'); 
            [~,CXC7,LXC7,~,~] = wdencmp('lvd',c,l,wname,lev,THR2,'s');
            velo_mm1=waverec(CXC7,LXC7,wname);
            disp_mm1=Inte(velo_mm1,dt);
            wname = w2; lev = n2;
            [c,l] = wavedec(disp_mm1,lev,wname);
            THR2 = wthrmngr('dw1ddenoLVL','rigrsure',c,l,'mln'); 
            [~,CXC7,LXC7,~,~] = wdencmp('lvd',c,l,wname,lev,THR2,'h');
            disp2=0;
            disp2=0;
            for ii=1:lev
                dd = wrcoef('d',CXC7,LXC7,wname,ii);
                disp2=disp2+dd;
            end
            velo2=[0;diff(disp2)]./dt;acce2=[0;diff(velo2)]./dt;
        elseif type==2
            wname = w1; lev = n1;
            [c,l] = wavedec(NS,lev,wname);
            THR2 = wthrmngr('dw1ddenoLVL','rigrsure',c,l,'mln'); 
            [~,CXC7,LXC7,~,~] = wdencmp('lvd',c,l,wname,lev,THR2,'h');
            acce_mm1=waverec(CXC7,LXC7,wname);
            velo_mm1=Inte(acce_mm1,dt);
            disp_mm1=Inte(velo_mm1,dt);
            wname = w2; lev = n2;
            [c,l] = wavedec(disp_mm1,lev,wname);
            THR2 = wthrmngr('dw1ddenoLVL','rigrsure',c,l,'mln'); 
            [~,CXC7,LXC7,~,~] = wdencmp('lvd',c,l,wname,lev,THR2,'h');
            disp2=0;
            for ii=1:lev
                dd = wrcoef('d',CXC7,LXC7,wname,ii);
                disp2=disp2+dd;
            end
            velo2=[0;diff(disp2)]./dt;acce2=[0;diff(velo2)]./dt;
        else 
            warning('Noisy signal type not recognized')
        end
    end
        DSW=velo2; % Wavelet correction only
        %DSB=BaselineCorr(NS,dt); %Base line correction only

        DSWB = velo2; % Baseline+Wavelet correction
        DSWb = DSWB;
        %% Filtering
        NF=1/2/dt; % Nyquest frequency
        f=fc/NF; %normalized cutoff frequency
        n=length(DSWb);
        DSWb=[(1:25000)'*0;DSWb;(1:25000)'*0];
        [DSWb]=LowcutFilt(DSWb,8,f,'high');
        DSWb_revers=zeros(length(DSWb),1);
        for i=1:length(DSWb)
            DSWb_revers(i)=DSWb(length(DSWb)-i+1);
        end
        [DSWb_revers]=LowcutFilt(DSWb_revers,8,f,'high');
        for i=1:length(DSWb)
            DSWb(i)=DSWb_revers(length(DSWb)-i+1);
        end
        temp_DSWb=DSWb(25001:25000+n);

        DSWBF=temp_DSWb;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions (cannot be called externally)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [a]=BaselineCorr(a,dt) % Base line correction Based on Iwan (1985) method
        % a: acceleration record data per m/s^2
        x=Inte(a,dt);
        t1=dt;
        for i=6:length(x)-6
            if abs(mean(a(i-5:i+5)))>0.5 && abs(mean(a(i-5:i+5)))<0.5
                t1=i*dt;
                break
            end
        end
        t2=dt*length(x);
        for i=length(x)-6:-1:6
            if abs(mean(a(i-5:i+5)))<0.5 && abs(mean(a(i-5:i+5)))>0.5
                t2=i*dt;
                break
            end
        end
        t2_tf=(round(t2/dt):length(x))'*dt;
        if length(t2_tf)==1
            t2=(length(x)-1)*dt;
            t2_tf=(round(t2/dt):length(x))'*dt;
        end
    
        VelFit = fit(t2_tf,x(round(t2/dt):length(x)),'poly1');
        af=VelFit.p1*[zeros(round(t2/dt),1);ones(length(x)-round(t2/dt),1)];
        am=x(round(t2/dt))/(t2-t1)*[zeros(round(t1/dt),1);ones(round(t2/dt)-round(t1/dt),1);zeros(length(x)-round(t2/dt),1)];
        a=a-af-am;

    end
%--------------------------------------------------------------------------    

    function [filteredSignal]=LowcutFilt(signal,n,f,type) % Butterworth filtering
        % Type: High, Low, Stop, Pass
        % n: order, Wn: corner frequency
        [z,p,k] = butter(n,f,type);
        sos = zp2sos(z,p,k);
        filteredSignal = sosfilt(sos,signal); %second order filtering
        %filteredSignal = filter(b,a,signal);
    end
%--------------------------------------------------------------------------    
    function [intx]=Inte(x,dt)
        intx=zeros(length(x),1);
        for i=2:length(x)
            intx(i)=intx(i-1)+(x(i)+x(i-1))*dt/2;
        end
    end    
%--------------------------------------------------------------------------        
    function [ fc,FA ] = Low_Cut_Freq( Signal,dt )
        % Inputs: 
        %        Signal: input signal
        %        dt: time step of input signal
        % Outputs
        %        fc: cut off frequency
        %        FA: fourier amplitude
    
        [f,amp]=fft_spec(Signal,dt);
        rr=0.2;
        len=length(amp);
        FA=zeros(len,1);
        dum=round((1:len)*rr);
        for i=1:len
            FA(i)=mean(amp(max(i-dum(i),1):min(i+dum(i),len)));
        end
        %semilogx(f,FA)
    
        [fc]=Find_Freq(FA,f);

    end
%--------------------------------------------------------------------------     
    function [f,fft_amp,fft_phase,fft_sig1]=fft_spec(sig,del_t) % calculate fourier traSignalforme
        n_fft=length(sig);
        fft_sig=fft(sig,n_fft)';
        f =((1/del_t)*(0:n_fft/2)/n_fft)';
        len=length(f);
        fft_amp=abs(fft_sig(1:len));
        fft_phase=angle(fft_sig(1:len));
        fft_sig1=fft_sig(1:len);
    end
%-------------------------------------------------------------------------- 
    function [yd]=dif(y,x) % Numerical differentiation
        % xmust be a column vector
        dx=[0;diff(x)];
        dy=[0;diff(y)];
        yd=dy./dx;
    end
%-------------------------------------------------------------------------- 
    function [fc]=Find_Freq(y,x) % Numerical differentiation
        % xmust be a column vector
        [yd]=dif(y,x);
        %figure
        %semilogx(x,yd)
        k=1;
        for i=4:length(x)-2
            if mean(yd(i-3:i-1))<0 && mean(yd(i:i+2))>0
                k=i;
                break
            end
        end
    
        if x(k)>0.3
            fc=x(2);
        else
            fc=x(k);
        end
    end       
%--------------------------------------------------------------------------