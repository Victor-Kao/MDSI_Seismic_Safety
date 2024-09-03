function data = Func_ImportMenhirData2Tab(path_csv)
    T = readtable(path_csv);
    X_mm_s = str2double(strrep(T.XInMm_s, ',', '.'));
    Y_mm_s = str2double(strrep(T.YInMm_s, ',', '.'));
    Z_mm_s = str2double(strrep(T.ZInMm_s, ',', '.'));
    date_time = T.ZeitInUTC_01_00;
    dt = seconds(date_time(3) - date_time(2));
    time = transpose(0:dt:(length(Z_mm_s)-1)*dt);
    data = table(date_time,time,X_mm_s,Y_mm_s,Z_mm_s);
end