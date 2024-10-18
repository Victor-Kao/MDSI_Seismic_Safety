folder_name = 'D:/MDSI_project/ANSYS_Building_model/Results_Ansys';
file_name = ['Disp_CH_15.csv'];
path = fullfile(folder_name,file_name );
FRF = readtable(path);

signal = FRF.REAL + 1i * FRF.IMAG;
omega = 2 * pi * FRF.Freq;
diff_signal = 1i * omega .* signal;

%plot(FRF.Freq,FRF.REAL);
%hold on 
%plot(FRF.Freq,FRF.IMAG);
%plot(FRF.Freq,FRF.AMPL);

plot(FRF.Freq,real(diff_signal));
hold on 
plot(FRF.Freq,imag(diff_signal));
plot(FRF.Freq,abs(diff_signal));