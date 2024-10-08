
%Energetic Third Octave
v_ref=5*10^(-8);
if isequal(X_cf_BP_samp,X_cf_OG1_samp,X_cf_OG2_samp,Y_cf_BP_samp,Y_cf_OG1_samp,Y_cf_OG2_samp,Z_cf_BP_samp,Z_cf_OG1_samp,Z_cf_OG2_samp)
    cf=X_cf_BP_samp;
    fl=X_cf_BP_samp*1.9953^(-1/(2*3)); %Frequenzbänder
    fu=X_cf_BP_samp*1.9953^(1/(2*3));
    fl(end+1)=fu(end);
else
    disp('There are different central frequencies')
end

L_v_BP_edit=20*log10(abs(Z_BP_FFT_V_out(:,2)/v_ref));
L_v_OG1_edit=20*log10(abs(Z_OG1_FFT_V_out(:,2)/v_ref));
L_v_OG2_edit=20*log10(abs(Z_OG2_FFT_V_out(:,2)/v_ref));
for i=1:length(cf)
    L_v_BP_per_fband{i} = L_v_BP_edit(Z_BP_FFT_V_out(:, 1) >= fl(i) & Z_BP_FFT_V_out(:, 1) < fl(i+1));
    L_v_ges_BP(i)=10*log10(sum(10.^(L_v_BP_per_fband{i}*0.1)));
    L_v_OG1_per_fband{i} = L_v_OG1_edit(Z_OG1_FFT_V_out(:, 1) >= fl(i) & Z_OG1_FFT_V_out(:, 1) < fl(i+1));
    L_v_ges_OG1(i)=10*log10(sum(10.^(L_v_OG1_per_fband{i}*0.1)));
    L_v_OG2_per_fband{i} = L_v_OG2_edit(Z_OG2_FFT_V_out(:, 1) >= fl(i) & Z_OG2_FFT_V_out(:, 1) < fl(i+1));
    L_v_ges_OG2(i)=10*log10(sum(10.^(L_v_OG2_per_fband{i}*0.1)));
    f_per_fband_{i} = Z_OG2_FFT_V_out(Z_OG2_FFT_V_out(:, 1) >= fl(i) & Z_OG2_FFT_V_out(:, 1) < fl(i+1), 1);

%     L_v_BP{i}=20*log10(v_per_fband_BP{i}/v_ref);
%     L_v_OG1{i}=20*log10(v_per_fband_OG1{i}/v_ref);
%     L_v_OG2{i}=20*log10(v_per_fband_OG2{i}/v_ref);
       
%     v_current_BP=v_per_fband_BP{i};
%     v_squared_BP=(v_current_BP.^2)/(v_ref^2);
%     L_v_ges_BP(i)=10*log10(sum(v_squared_BP));
%     L_v_ges_BP_2(i)=10*log10(sum(10.^(L_v_BP{i}*0.1)));
%     if isequal(L_v_ges_BP(i),L_v_ges_BP_2(i))
%     else
%         disp('There are solutions for L_v_ges per 1/3 octave')
%     end
%     
%     v_current_OG1=v_per_fband_OG1{i};
%     v_squared_OG1=(v_current_OG1.^2)/(v_ref^2);
%     L_v_ges_OG1(i)=10*log10(sum(v_squared_OG1));
%     L_v_ges_OG1_2(i)=10*log10(sum(10.^(L_v_OG1{i}*0.1)));
%     if isequal(L_v_ges_OG1(i),L_v_ges_OG1_2(i))
%     else
%         disp('There are solutions for L_v_ges per 1/3 octave')
%     end
%     
%     v_current_OG2=v_per_fband_OG2{i};
%     v_squared_OG2=(v_current_OG2.^2)/(v_ref^2);
%     L_v_ges_OG2(i)=10*log10(sum(v_squared_OG2));
%     L_v_ges_OG2_2(i)=10*log10(sum(10.^(L_v_OG2{i}*0.1)));
%     if isequal(L_v_ges_OG2(i),L_v_ges_OG2_2(i))
%     else
%         disp('There are several solutions for L_v_ges per 1/3 octave')
%     end
    
%     plot(f_per_fband_{i},v_per_fband_BP{i},'r',f_per_fband_{i},v_per_fband_OG1{i},'b',f_per_fband_{i},v_per_fband_OG2{i},'g');
%     %legend('BP', '1.OG', '2.OG');
%     hold on
%     yyaxis right
%     plot(f_per_fband_{i},L_v_BP{i},'r--',f_per_fband_{i},L_v_OG1{i},'b--',f_per_fband_{i},L_v_OG2{i},'g--')
%     legend('BP', '1.OG', '2.OG','L_v BP', 'L_v 1.OG', 'L_v 2.OG');
%     xlim([0 40])
end

% hold on
% yyaxis right
% plot(cf,L_v_ges_BP(2:end-1),cf,L_v_ges_OG1(2:end-1),cf,L_v_ges_OG2(2:end-1));
    
% L_v_ges_BP=L_v_ges_BP(2:end-1)*0.00001;
% L_v_ges_OG1=L_v_ges_OG1(2:end-1)*0.00001;
% L_v_ges_OG2=L_v_ges_OG2(2:end-1)*0.00001;
% hold on
% plot(X_cf_BP_samp,L_v_ges_BP)
figure
%plot(Z_BP_FFT_V_out(:,1),L_v_BP_edit,'b--',Z_OG1_FFT_V_out(:,1),L_v_OG1_edit,'r--',Z_OG2_FFT_V_out(:,1),L_v_OG2_edit,'g--',cf,L_v_ges_BP,'b',cf,L_v_ges_OG1,'r',cf,L_v_ges_OG2,);
plot(Z_BP_FFT_V_out(:,1),L_v_BP_edit,"--",cf,L_v_ges_BP,'Color',[0 0.4470 0.7410])
hold on
plot(Z_OG1_FFT_V_out(:,1),L_v_OG1_edit,"--",cf,L_v_ges_OG1,'Color',[0.8500 0.3250 0.0980]);
hold on
plot(Z_OG2_FFT_V_out(:,1),L_v_OG2_edit,"--",cf,L_v_ges_OG2,'Color',[0.9290 0.6940 0.1250]);
xlim([0 50])
legend('BP', '1.OG', '2.OG','BP', '1.OG', '2.OG')
title('1/3 Octave Z-Direction (Sampled Time)')

%figure
%plot(Z_BP_FFT_V_out(:,1),L_v_BP_edit,Z_OG1_FFT_V_out(:,1),L_v_OG1_edit,Z_OG2_FFT_V_out(:,1),L_v_OG2_edit)%,X_cf_BP_samp,L_v_ges_BP,X_cf_BP_samp,L_v_ges_OG1,X_cf_BP_samp,L_v_ges_OG2)
%legend('BP', '1.OG', '2.OG','BP', '1.OG', '2.OG')
%Z_OG2_FFT_V_out(:,2)
%L_v=10*lg((Z_OG2_FFT_V_out(:,2))^2/v_ref)
