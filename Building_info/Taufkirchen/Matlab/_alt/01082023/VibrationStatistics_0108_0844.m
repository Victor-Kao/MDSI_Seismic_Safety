%Project Statistics
%Robert Gumuchian
%August 2023

clear
clc
close all

%% Data Import
load Data_Processing\Event_Statistics\Velocity\01082023\BP_010823_0844AM.mat
load Data_Processing\Event_Statistics\Velocity\01082023\OG1_010823_0844AM.mat
load Data_Processing\Event_Statistics\Velocity\01082023\OG2_010823_0844AM.mat


%Signal start=4.8s, end=6.1s
t1=16.612;          %Measurement start
t2=31.799;          %Measurement end
T1=21.412;          %Sampling start
T2=22.712;          %Sampling end
delta_t=0.001;      %Timestep of sensor
t=t1:0.001:t2;
T_samp=T1:0.001:T2;


V_X_BP_0844AM=table2array(BP_010823_0844AM(:,2));
V_X_BP_0844AM_samp(:,1)=T_samp;
V_X_BP_0844AM_samp(:,2)=V_X_BP_0844AM(4801:6101);

V_Y_BP_0844AM=table2array(BP_010823_0844AM(:,3));
V_Y_BP_0844AM_samp(:,1)=T_samp;
V_Y_BP_0844AM_samp(:,2)=V_Y_BP_0844AM(4801:6101);

V_Z_BP_0844AM=table2array(BP_010823_0844AM(:,4));
V_Z_BP_0844AM_samp(:,1)=T_samp;
V_Z_BP_0844AM_samp(:,2)=V_Z_BP_0844AM(4801:6101);

V_X_OG1_0844AM=table2array(OG1_010823_0844AM(:,2));
V_X_OG1_0844AM_samp(:,1)=T_samp;
V_X_OG1_0844AM_samp(:,2)=V_X_OG1_0844AM(4801:6101);

V_Y_OG1_0844AM=table2array(OG1_010823_0844AM(:,3));
V_Y_OG1_0844AM_samp(:,1)=T_samp;
V_Y_OG1_0844AM_samp(:,2)=V_Y_OG1_0844AM(4801:6101);

V_Z_OG1_0844AM=table2array(OG1_010823_0844AM(:,4));
V_Z_OG1_0844AM_samp(:,1)=T_samp;
V_Z_OG1_0844AM_samp(:,2)=V_Z_OG1_0844AM(4801:6101);

V_X_OG2_0844AM=table2array(OG2_010823_0844AM(:,2));
V_X_OG2_0844AM_samp(:,1)=T_samp;
V_X_OG2_0844AM_samp(:,2)=V_X_OG2_0844AM(4801:6101);

V_Y_OG2_0844AM=table2array(OG2_010823_0844AM(:,3));
V_Y_OG2_0844AM_samp(:,1)=T_samp;
V_Y_OG2_0844AM_samp(:,2)=V_Y_OG2_0844AM(4801:6101);

V_Z_OG2_0844AM=table2array(OG2_010823_0844AM(:,4));
V_Z_OG2_0844AM_samp(:,1)=T_samp;
V_Z_OG2_0844AM_samp(:,2)=V_Z_OG2_0844AM(4801:6101);

%BP measurments seem much smaller so they will be plotted in all directions
figure
plot(t,V_X_BP_0844AM);
xlabel('Time [s]')
ylabel('Velocity [mm/s]')
title('BP only')
hold on
plot(t,V_Y_BP_0844AM);
plot(t,V_Z_OG2_0844AM);
legend('X','Y','Z')

figure
plot(T_samp,V_X_BP_0844AM_samp(:,2));
xlabel('Time [s]')
ylabel('Velocity [mm/s]')
title('X-Direction Slab Velocity for Event on August 1, 2023 at 8:44AM')
hold on
plot(T_samp,V_X_OG1_0844AM_samp(:,2));
plot(T_samp,V_X_OG2_0844AM_samp(:,2));
legend('BP','1OG Floor', '2OG Floor')

figure
plot(T_samp,V_Y_BP_0844AM_samp(:,2));
xlabel('Time [s]')
ylabel('Velocity [mm/s]')
title('Y-Direction Slab Velocity for Event on August 1, 2023 at 8:44AM')
hold on
plot(T_samp,V_Y_OG1_0844AM_samp(:,2));
plot(T_samp,V_Y_OG2_0844AM_samp(:,2));
legend('BP','1OG Floor', '2OG Floor')

figure
plot(T_samp,V_Z_BP_0844AM_samp(:,2));
xlabel('Time [s]')
ylabel('Velocity [mm/s]')
title('Z-Direction Slab Velocity for Event on August 1, 2023 at 8:44AM')
hold on
plot(T_samp,V_Z_OG1_0844AM_samp(:,2));
plot(T_samp,V_Z_OG2_0844AM_samp(:,2));
legend('BP','1OG Floor', '2OG Floor')

% D_X_BP_010823_0844AM=cumtrapz(t_BP,V_X_BP_0844AM);
% D_Y_BP_010823_0844AM=cumtrapz(t_BP,V_Y_BP_0844AM);
% D_Z_BP_010823_0844AM=cumtrapz(t_BP,V_Z_BP_0844AM);
% 
% D_X_OG1_010823_0844AM=cumtrapz(t_OG1,V_X_OG1_0844AM);
% D_Y_OG1_010823_0844AM=cumtrapz(t_OG1,V_Y_OG1_0844AM);
% D_Z_OG1_010823_0844AM=cumtrapz(t_OG1,V_Z_OG1_0844AM);
% 
% D_X_OG2_010823_0844AM=cumtrapz(t_OG2,V_X_OG2_0844AM);
% D_Y_OG2_010823_0844AM=cumtrapz(t_OG2,V_Y_OG2_0844AM);
% D_Z_OG2_010823_0844AM=cumtrapz(t_OG2,V_Z_OG2_0844AM);
% 
% figure
% plot(t_BP,D_X_BP_010823_0844AM);
% xlabel('Time [s]')
% ylabel('Displacement [mm]')
% title('X-Direction Slab Displacement for Event on August 1, 2023 at 8:44AM')
% hold on
% plot(t_OG1,D_X_OG1_010823_0844AM);
% plot(t_OG2,D_X_OG2_010823_0844AM);
% legend('BP','1OG Floor', '2OG Floor')
% 
% figure
% plot(t_BP,D_Y_BP_010823_0844AM);
% xlabel('Time [s]')
% ylabel('Displacement [mm]')
% title('Y-Direction Slab Displacement for Event on August 1, 2023 at 8:44AM')
% hold on
% plot(t_OG1,D_Y_OG1_010823_0844AM);
% plot(t_OG2,D_Y_OG2_010823_0844AM);
% legend('BP','1OG Floor', '2OG Floor')
% 
% figure
% plot(t_BP,D_Z_BP_010823_0844AM);
% xlabel('Time [s]')
% ylabel('Displacement [mm]')
% title('Z-Direction Slab Displacement for Event on August 1, 2023 at 8:44AM')
% hold on
% plot(t_OG1,D_Z_OG1_010823_0844AM);
% plot(t_OG2,D_Z_OG2_010823_0844AM);
% legend('BP','1OG Floor', '2OG Floor')

% Xnoise_BP=snr(X_BP_010823_0844AM);
% Ynoise_BP=snr(Y_BP_010823_0844AM);
% Znoise_BP=snr(Z_BP_010823_0844AM);

%%X-Direction Fourier Analysis ALL
Fs = 1/delta_t;                                                 % Sampling frequency       
L_all = length(t);                                     % Length of signal
t_vec = (0:L_all-1)*delta_t;                                        % Time vector
X_BP_FFT_all = fft(V_X_BP_0844AM);                  % fft
X_OG1_FFT_all = fft(V_X_OG1_0844AM);
X_OG2_FFT_all = fft(V_X_OG2_0844AM);

X_BP_scale_all = X_BP_FFT_all/Fs;                                 % scaling according to Parsival theorem
X_OG1_scale_all = X_OG1_FFT_all/Fs;
X_OG2_scale_all = X_OG2_FFT_all/Fs;

X_BP_FFT_out_all(:,2) = X_BP_scale_all(1:round(L_all/2+1));                  % single sided
X_BP_FFT_out_all(2:end-1,2) = 2*X_BP_FFT_out_all(2:end-1,2);      % doubled
X_BP_FFT_out_all(:,1)=Fs*(0:(L_all/2))/L_all;                         % frequency

X_OG1_FFT_out_all(:,2) = X_OG1_scale_all(1:round(L_all/2+1));                  % single sided
X_OG1_FFT_out_all(2:end-1,2) = 2*X_OG1_FFT_out_all(2:end-1,2);      % doubled
X_OG1_FFT_out_all(:,1)=Fs*(0:(L_all/2))/L_all;                        % frequency

X_OG2_FFT_out_all(:,2) = X_OG2_scale_all(1:round(L_all/2+1));                  % single sided
X_OG2_FFT_out_all(2:end-1,2) = 2*X_OG2_FFT_out_all(2:end-1,2);      % doubled
X_OG2_FFT_out_all(:,1)=Fs*(0:(L_all/2))/L_all;                        % frequency

figure
plot(abs(X_BP_FFT_out_all(:,1)),abs(X_BP_FFT_out_all(:,2)));
xlabel('Frequency [Hz]')
ylabel('Velocity [mm/s]')
title('ALL X-Direction BP Frequency Spectrum for Event on August 1, 2023 at 8:44AM')

figure
plot(abs(X_OG1_FFT_out_all(:,1)),abs(X_OG1_FFT_out_all(:,2)));
xlabel('Frequency [Hz]')
ylabel('Velocity [mm/s]')
title('ALL X-Direction 1.OG Frequency Spectrum for Event on August 1, 2023 at 8:44AM')

figure
plot(abs(X_OG2_FFT_out_all(:,1)),abs(X_OG2_FFT_out_all(:,2)));
xlabel('Frequency [Hz]')
ylabel('Velocity [mm/s]')
title('ALL X-Direction 2.OG Frequency Spectrum for Event on August 1, 2023 at 8:44AM')


%%X-Direction Fourier Analysis Sample Area
Fs = 1/delta_t;                                         % Sampling frequency       
L = length(T_samp);                                     % Length of signal
t_vec = (0:L-1)*delta_t;                                % Time vector
X_BP_FFT = fft(V_X_BP_0844AM_samp(:,2));                  % fft
X_OG1_FFT = fft(V_X_OG1_0844AM_samp(:,2));
X_OG2_FFT = fft(V_X_OG2_0844AM_samp(:,2));

X_BP_scale = X_BP_FFT/Fs;                                 % scaling according to Parsival theorem
X_OG1_scale = X_OG1_FFT/Fs;
X_OG2_scale = X_OG2_FFT/Fs;

X_BP_FFT_out(:,2) = X_BP_scale(1:round(L/2));                  % single sided
X_BP_FFT_out(2:end-1,2) = 2*X_BP_FFT_out(2:end-1,2);      % doubled
X_BP_FFT_out(:,1)=Fs*(0:(L/2))/L;                         % frequency

X_OG1_FFT_out(:,2) = X_OG1_scale(1:round(L/2));                  % single sided
X_OG1_FFT_out(2:end-1,2) = 2*X_OG1_FFT_out(2:end-1,2);      % doubled
X_OG1_FFT_out(:,1)=Fs*(0:(L/2))/L;                        % frequency

X_OG2_FFT_out(:,2) = X_OG2_scale(1:round(L/2));                  % single sided
X_OG2_FFT_out(2:end-1,2) = 2*X_OG2_FFT_out(2:end-1,2);      % doubled
X_OG2_FFT_out(:,1)=Fs*(0:(L/2))/L;                        % frequency

figure
plot(abs(X_BP_FFT_out(:,1)),abs(X_BP_FFT_out(:,2)));
xlabel('Frequency [Hz]')
ylabel('Velocity [mm/s]')
title('X-Direction BP Frequency Spectrum for Event on August 1, 2023 at 8:44AM')

figure
plot(abs(X_OG1_FFT_out(:,1)),abs(X_OG1_FFT_out(:,2)));
xlabel('Frequency [Hz]')
ylabel('Velocity [mm/s]')
title('X-Direction 1.OG Frequency Spectrum for Event on August 1, 2023 at 8:44AM')

figure
plot(abs(X_OG2_FFT_out(:,1)),abs(X_OG2_FFT_out(:,2)));
xlabel('Frequency [Hz]')
ylabel('Velocity [mm/s]')
title('X-Direction 2.OG Frequency Spectrum for Event on August 1, 2023 at 8:44AM')


%%Y-Direction Fourier Analysis
Y_BP_FFT = fft(V_Y_BP_0844AM_samp(:,2));                  % fft
Y_OG1_FFT = fft(V_Y_OG1_0844AM_samp(:,2));
Y_OG2_FFT = fft(V_Y_OG2_0844AM_samp(:,2));

Y_BP_scale = Y_BP_FFT/Fs;                                 % scaling according to Parsival theorem
Y_OG1_scale = Y_OG1_FFT/Fs;
Y_OG2_scale = Y_OG2_FFT/Fs;

Y_BP_FFT_out(:,2) = Y_BP_scale(1:round(L/2));             % single sided
Y_BP_FFT_out(2:end-1,2) = 2*Y_BP_FFT_out(2:end-1,2);      % doubled
Y_BP_FFT_out(:,1)=Fs*(0:(L/2))/L;                         % frequency

Y_OG1_FFT_out(:,2) = Y_OG1_scale(1:round(L/2));           % single sided
Y_OG1_FFT_out(2:end-1,2) = 2*Y_OG1_FFT_out(2:end-1,2);    % doubled
Y_OG1_FFT_out(:,1)=Fs*(0:(L/2))/L;                        % frequency

Y_OG2_FFT_out(:,2) = Y_OG2_scale(1:round(L/2));           % single sided
Y_OG2_FFT_out(2:end-1,2) = 2*Y_OG2_FFT_out(2:end-1,2);    % doubled
Y_OG2_FFT_out(:,1)=Fs*(0:(L/2))/L;                        % frequency

figure
plot(abs(Y_BP_FFT_out(:,1)),abs(Y_BP_FFT_out(:,2)));
xlabel('Frequency [Hz]')
ylabel('Velocity [mm/s]')
title('Y-Direction BP Frequency Spectrum for Event on August 1, 2023 at 8:44AM')

figure
plot(abs(Y_OG1_FFT_out(:,1)),abs(Y_OG1_FFT_out(:,2)));
xlabel('Frequency [Hz]')
ylabel('Velocity [mm/s]')
title('Y-Direction 1.OG Frequency Spectrum for Event on August 1, 2023 at 8:44AM')

figure
plot(abs(Y_OG2_FFT_out(:,1)),abs(Y_OG2_FFT_out(:,2)));
xlabel('Frequency [Hz]')
ylabel('Velocity [mm/s]')
title('Y-Direction 2.OG Frequency Spectrum for Event on August 1, 2023 at 8:44AM')

%%Z-Direction Fourier Analysis
Z_BP_FFT = fft(V_Z_BP_0844AM_samp(:,2));                  % fft
Z_OG1_FFT = fft(V_Z_OG1_0844AM_samp(:,2));
Z_OG2_FFT = fft(V_Z_OG2_0844AM_samp(:,2));

Z_BP_scale = Z_BP_FFT/Fs;                                 % scaling according to Parsival theorem
Z_OG1_scale = Z_OG1_FFT/Fs;
Z_OG2_scale = Z_OG2_FFT/Fs;

Z_BP_FFT_out(:,2) = Z_BP_scale(1:round(L/2));                  % single sided
Z_BP_FFT_out(2:end-1,2) = 2*Z_BP_FFT_out(2:end-1,2);      % doubled
Z_BP_FFT_out(:,1)=Fs*(0:(L/2))/L;                         % frequency

Z_OG1_FFT_out(:,2) = Z_OG1_scale(1:round(L/2));                  % single sided
Z_OG1_FFT_out(2:end-1,2) = 2*Z_OG1_FFT_out(2:end-1,2);      % doubled
Z_OG1_FFT_out(:,1)=Fs*(0:(L/2))/L;                        % frequency

Z_OG2_FFT_out(:,2) = Z_OG2_scale(1:round(L/2));                  % single sided
Z_OG2_FFT_out(2:end-1,2) = 2*Z_OG2_FFT_out(2:end-1,2);      % doubled
Z_OG2_FFT_out(:,1)=Fs*(0:(L/2))/L;                        % frequency

figure
plot(abs(Z_BP_FFT_out(:,1)),abs(Z_BP_FFT_out(:,2)));
xlabel('Frequency [Hz]')
ylabel('Velocity [mm/s]')
title('Z-Direction BP Frequency Spectrum for Event on August 1, 2023 at 8:44AM')

figure
plot(abs(Z_OG1_FFT_out(:,1)),abs(Z_OG1_FFT_out(:,2)));
xlabel('Frequency [Hz]')
ylabel('Velocity [mm/s]')
title('Z-Direction 1.OG Frequency Spectrum for Event on August 1, 2023 at 8:44AM')

figure
plot(abs(Z_OG2_FFT_out(:,1)),abs(Z_OG2_FFT_out(:,2)));
xlabel('Frequency [Hz]')
ylabel('Velocity [mm/s]')
title('Z-Direction 2.OG Frequency Spectrum for Event on August 1, 2023 at 8:44AM')