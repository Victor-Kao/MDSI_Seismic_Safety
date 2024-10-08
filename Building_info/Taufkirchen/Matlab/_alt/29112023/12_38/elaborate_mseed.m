classdef elaborate_mseed
   properties
      PropertyName
   end

   methods (Static)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Class functions (can be accessed externally by classname.functionname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------- 
function [Event_station_matrix] = matrix_station_events( st_name_table,event_number,st_number,Flag_station,event )
    % Create a station vs event matrix and export for latex
    % column of the station names
    Stations_list = st_name_table.Properties.VariableNames;
    Stations_list = Stations_list';
    
    Event_station_matrix = table(Stations_list);
    % Gives a x for a station that recorded an event
    for e =1:event_number
        for i = 1:st_number
                if any(Flag_station(1,:,i,e))== 1
                        %str = num2str(table2array(st_name_table(j,i)));
                        %Event_station_matrix(i,e)=strjoin([string(table2array((Event_station_matrix(i,e)))), ',', str]); % Append this string.                                     
                        Event_station_matrix(i,e+1)={'x'};
                else
                    Event_station_matrix(i,e+1)={''};
                end
        end 
    end
    % Givenames to the columns header
    Event_station_matrix.Properties.VariableNames(2:end)=cellstr(datestr(event.StartTime));
    % Export to latex 
    %table2latex(Event_station_matrix, 'tabs\Event_station_matrix.tex') 
end


%-------------------------------------------------------------------------- 
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

        DSWB = BaselineCorr(velo2,dt); % Baseline+Wavelet correction
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

        DSWBF=BaselineCorr(temp_DSWb,dt);

    end
    
%-------------------------------------------------------------------------- 
    function [EcumTH,EcumTH_norm,t_5_95,Td, arias ] = energy_duration_ai(dt, xgtt, time)

            % CUMULATIVE ENERGY
            % time history of cumulative energy
            EcumTH = cumsum(xgtt.^2).*dt;
            % Total cumulative energy at the end of the ground motion
            Ecum = EcumTH(end);
            % time history of the normalized cumulative energy
            EcumTH_norm = EcumTH/Ecum;

            % SIGNIFICANT DURATION
            % elements of the time vector which are within the significant
            % duration
            timed = time(EcumTH>=0.05*Ecum & EcumTH<=0.95*Ecum);
            % starting and ending points of the significant duration
            t_5_95 = [timed(1),timed(end)];
            % significant duration
            Td = timed(end)-timed(1)+dt;

            % ARIAS INTENSITY
            % time history of Arias Intensity
            ariasTH = 1/9.81*cumsum(xgtt(time<=Td).^2).*pi.*dt/2;
            % Total Arias Intensity at the end of the ground motion
            arias = ariasTH(end);
    end
%--------------------------------------------------------------------------    
    function [intx]=Inte(x,dt)
        intx=zeros(length(x),1);
        for i=2:length(x)
            intx(i)=intx(i-1)+(x(i)+x(i-1))*dt/2;
        end
    end    
%--------------------------------------------------------------------------
    function [H_V_ration]=HV_ratio(order, framelen, st_name_table, vel_g_bar_single_denoised,freq,Nsamples, channels, stations, events, event)
    % Suggestion
     %order=3;
     %framelen=13;
    figure;
    n_loc=0; % inizialise number of signal 
    for e = events % selected events
        for i = stations %1:st_number % selected stations
            for j= channels
                   % Geom Mean in H-dir  
                   H1 = vel_g_bar_single_denoised(1:Nsamples(1,channels(j),i,e)./2+1,channels(j),i,e);
                   H2 = interp1(freq(1:Nsamples(1,channels(j)+1,i,e)./2+1,channels(j)+1,i,e),...
                                 vel_g_bar_single_denoised(1:Nsamples(1,channels(j)+1,i,e)./2+1,channels(j)+1,i,e),...
                                 freq(1:Nsamples(1,channels(j),i,e)./2+1,channels(j),i,e));
                   % Geomterical mean
                   H_mean = (H1.^2+H2.^2).^(1/2);
                   % Smoothing functions
                    sgf_H_mean(1:Nsamples(1,channels(j),i,e)./2+1,channels(j),i,e) = sgolayfilt(abs(H_mean(1:Nsamples(1,channels(j),i,e)./2+1)),order,framelen);
                    sgf_V(1:Nsamples(1,channels(j)+2,i,e)./2+1,channels(j)+2,i,e) = sgolayfilt(abs(vel_g_bar_single_denoised(1:Nsamples(1,channels(j)+2,i,e)./2+1,channels(j)+2,i,e)),order,framelen);
                    sgf_V(1:Nsamples(1,channels(j),i,e)./2+1,channels(j),i,e) = interp1(freq(1:Nsamples(1,channels(j)+2,i,e)./2+1,channels(j)+2,i,e),...
                                                                                        sgf_V(1:Nsamples(1,channels(j)+2,i,e)./2+1,channels(j)+2,i,e),...
                                                                                        freq(1:Nsamples(1,channels(j),i,e)./2+1,channels(j),i,e));
    
                    H_V_ration=sgf_H_mean./sgf_V;
    
                    plot(freq(1:Nsamples(1,channels(j),i,e)./2+1,channels(j),i,e),H_V_ration(1:Nsamples(1,channels(j),i,e)./2+1,channels(j),i,e) )
                    n_loc=n_loc+1;
                      hold on
                    % legendCell(n_loc)=['Event:',datestr(event.StartTime(e))];
                    clear H1 H2 H_mean
            end
        end
    
    end
    %legend(legendCell)
    grid on;xlabel('$f$ [Hz]'); ylabel('$H/V$ [-]');
    title(['H/V ratio. Station '...
          ,strjoin([st_name_table.Properties.VariableNames(i) ,'-', table2array(st_name_table(j,i))])]);
         xlim([0.5 30])
    
    hold off
    end
%-------------------------------------------------------------------------- 
function [Mw] = Ml2MwGruenthal(Ml)
    % Convert according to Grünthal
    Mw = 0.67 + 0.56.*Ml+0.046.*Ml^2;
end
%-------------------------------------------------------------------------- 
function [Mw] = Ml2Mw(Ml)
    if Ml<2
    % Convert according to Allmann et al. (2010)
    Mw = 0.594*Ml + 0.985;
    elseif and(Ml>=2, Ml<4)
    Mw = 0.594*Ml + 0.985;
    elseif Ml>=4
    Mw = Ml - 0.3;
end
%-------------------------------------------------------------------------- 
    function [T,Spa,Spv,Sd,Sv,Sa]=SPEC(dt,Ag,zet,g,endp)
    
    %% Elastic Response Spectra, Version 2 
    % This function generates elastic response specra including Displacement
    % Spectrum, Pseudo Acceleration Spectrum, and Pseudo Velocity Spectrum which
    % are needed in a "Response Spectrum Analysis" of structures. To solve 
    % the "equation of motions" for different periods, the Newmark Linear Method 
    % was used. 
    
    %% Update Note:  It was clarified that Ag has the acceleration unit.
    
    %% (c) Mostafa Tazarv, South Dakota State University, May 2019
    
    %% SPEC Function Help:
    
    % INPUT:
    % dt:     Time Interval (Sampling rate) of the Ground Motion
    % Ag:     Ground Motion Acceleration in the unit of the "acceleration", e.g. m/s^2  
    % zet:    Damping Ratio in percent (%); e.g. 5
    % g:      Gravitational Constant; e.g. 9.81 m/s^2; g determines the output unit
    % endp:   End Period of the Spectra; e.g. 4 sec.
    
    % OUTPUT:
    % T:      Period of the Structure (sec.)
    % Spa:    Elastic Pseudo Acceleration Spectrum in g
    % Spv:    Elastic Pseudo Velocity Spectrum in the unit of velocity (e.g. m/s)
    % Sd:     Elastic Displacement Spectrum, in the unit of displacement (e.g. m)
    
    
    u=zeros(length(Ag),1);
    v=zeros(length(Ag),1);
    ac=zeros(length(Ag),1);
    Ag(end+1)=0;
    T(1,1)=0.00;
    for j=1:round(endp/dt)                          % equation of motion(Newmark linear method)
        omega(j,1)=2*pi/T(j);      % Natural Frequency
        m=1;       
        k=(omega(j))^2*m;
        c=2*m*omega(j)*zet/100;
        K=k+3*c/dt+6*m/(dt)^2;
        a=6*m/dt+3*c;
        b=3*m+dt*c/2;    
      for i=1:length(u)-1
         u(1,1)=0;                      %initial conditions
         v(1,1)=0;
         ac(1,1)=0;    
         df=-(Ag(i+1)-Ag(i))+a*v(i,1)+b*ac(i,1);  % delta Force
         du=df/K;
         dv=3*du/dt-3*v(i,1)-dt*ac(i,1)/2;
         dac=6*(du-dt*v(i,1))/(dt)^2-3*ac(i,1);
         u(i+1,1)=u(i,1)+du;
         v(i+1,1)=v(i,1)+dv;
         ac(i+1,1)=ac(i,1)+dac;     
      end
        Sd(j,1)=max(abs((u(:,1))));
        Sv(j,1)=max(abs(v));
        Sa(j,1)=max(abs(ac))/g;
        Spv(j,1)=Sd(j)*omega(j);
        Spa(j,1)=Sd(j)*(omega(j))^2/g;
        T(j+1,1)=T(j)+dt;
    end
    Ag(end)=[];
    T(end)=[];
    Sd(2,1)=0; Spv(1:2,1)=0;Spa(1:2,1)=max(abs(Ag))/g;
    Sv(1:2,1)=0; Sa(1:2,1)=0;
    
%     %% Plot Spectra
%     subplot(2,1,1)
%      %figure('Name','Spectral Displacement','NumberTitle','off')
%      plot(T,Sd,'LineWidth',2.)
%      grid on
%     xlabel('Period (sec)','FontSize',13);
%     ylabel('Sd (mm)','FontSize',13);
%     title('Displacement Spectrum','FontSize',13)
%     
%     subplot(2,1,2)
%      %figure('Name','Pseudo Acceleration Spectrum','NumberTitle','off')
%      plot(T,Spa,'LineWidth',2.)
%      grid on
%     xlabel('Period (sec)','FontSize',13);
%     ylabel('Spa (g)','FontSize',13);
%     title('Pseudo Acceleration Spectrum','FontSize',13)
    end  
%--------------------------------------------------------------------------  
    function [spec_mean]=spec_geom_mean(spec_1, spec_2)
        spec_mean=(spec_1.^2+spec_2.^2).^(1/2);
    end
%--------------------------------------------------------------------------  

   end % end method
%--------------------------------------------------------------------------    
end % end class
%--------------------------------------------------------------------------    

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