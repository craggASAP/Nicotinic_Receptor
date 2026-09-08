%%%%%%%%%%%%%%%%%%% Imaging data analysis %%%%%%%%%%%%%%%%%%%%
% !! Remember to modify parameters depending on acquisition protocol and stim protocol!!
% (exposure / time of stim / etc..)
%  Lucille Duquenoy 28/08/24 - added 2p100Hz, 3p100Hz 


% Before using the code you need to create two tables as following:

% an input table called A with:
% Column 1 = with the frame number 
% Column 2 to the end= raw fluorescence of transients

% an input table called "stims" with:
% column 1 = 0 
% column 2 onwards = the stim done for the corresponding fluo transient


total_number_of_stims=size(A,2)-1; %extracting how many stims were done during the whole experiment (the -1 is necessary because it counts the number of columns in table A and the first column is the frame number
total_frame_number=size(A(:,1)); %extracting the total number of frames
frame_number=A(:,1); %extracting the frame number
exposure=1/599.52; % exposure every 1.67ms (frequency of recording 600Hz)- to change if we chqnge the number of frame per seconds
time=(total_frame_number.*exposure)./100; % converting to seconds, to obtain total recording time in seconds 


first_spike_absolute_heigth_all=zeros(1,total_number_of_stims);
second_spike_absolute_heigth_all=zeros(1,total_number_of_stims);
third_spike_absolute_heigth_all=zeros(1,total_number_of_stims);
fourth_spike_absolute_heigth_all=zeros(1,total_number_of_stims);

first_spike_F0_all=zeros(1,total_number_of_stims);
second_spike_F0_all=zeros(1,total_number_of_stims);
third_spike_F0_all=zeros(1,total_number_of_stims);
fourth_spike_F0_all=zeros(1,total_number_of_stims);

first_spike_rebaselined_heigth_all=zeros(1,total_number_of_stims);
second_spike_rebaselined_heigth_all=zeros(1,total_number_of_stims);
third_spike_rebaselined_heigth_all=zeros(1,total_number_of_stims);
fourth_spike_rebaselined_heigth_all=zeros(1,total_number_of_stims);

first_spike_AUC_all=zeros(1,total_number_of_stims);
second_spike_AUC_all=zeros(1,total_number_of_stims);
third_spike_AUC_all=zeros(1,total_number_of_stims);
fourth_spike_AUC_all=zeros(1,total_number_of_stims);

first_spike_trace_all=zeros(9,total_number_of_stims);
second_spike_trace_all=zeros(9,total_number_of_stims);
third_spike_trace_all=zeros(9,total_number_of_stims);
fourth_spike_trace_all=zeros(9,total_number_of_stims);



for i=2:size(A,2); % size(A,2) is the amount of transients you have in your experiment
    
    raw_F=A(:,i); % extracting the raw fluorescence values of the recording number i

    % Creating a linear fit to correct the bleaching. 
    % The linear fit is taken to match 400 frames (~0.60 seconds) before stim onset
    % stim onset at 1.4 seconds ~ frame 839
    Fit=fit(A(400:800,1),raw_F(400:800,:), 'poly1'); % fit is a*x+b = rawF, adjusting the a and b coefs
    % getting the coefficients: slope=a=coefs(1) and F0=b=coefs(2) of the fit
    coefs=coeffvalues(Fit);
    slope_fit=coefs(1);
    F0_fit=coefs(2);
    % Compiling all F0 from all traces   
    F0_all(:,i-1)=F0_fit;

    % to correct for slope, need to subtract the slope from the fit (a*x) from the original trace  
    F_slope_corrected=raw_F-(slope_fit*frame_number); 
    % Corrects for baseline by calclating F-F0 and inverting transient to see them positive
    F_baseline_corrected=-(F_slope_corrected-F0_fit); 
    % Calculate F-F0/F0
    deltaF_over_F0=F_baseline_corrected/F0_fit;
    % Converting to % of F-F0/F0
    deltaF_over_F0_percent=deltaF_over_F0.*100;
    % Compiling all F-F0/F0 from all traces and converting to % of F-F0/F0  
    deltaF_over_F0_percent_all(:,i-1)=deltaF_over_F0_percent;

    % calculate the total AUC
    total_AUC_all(:,i-1)=trapz(deltaF_over_F0_percent(840:900)); %840 frame is when the first pulse happens, and 900 frames (100ms later) seems to be when there is a return to baseline for 1p outside of DHBE
    
    
    % Find spikes (spike heights and traces are in % of deltaF/F0 already)
    
    
    % 1st spike ( always at 1.4s = 840 frames)
    
    t1=840;
    first_spike_absolute_heigth=max(deltaF_over_F0_percent(t1-1:t1+3)); % taking a narrow window around stim (at 840) so it doesn't count the cholinergic bump when it exists
    % need to align the trace to stim onset at 840
    first_spike_F0=mean(deltaF_over_F0_percent(t1-2:t1-1)); % F0 is an average of the 2 points before time of stim t1 (at t1-2 and t1-1) 
    first_spike_trace=deltaF_over_F0_percent(t1-3:t1+5)- first_spike_F0;% taking a wider window and will cross later if necessary and re-baselining it taking 3 points before the time of stim t1 anf 5 points after
    first_spike_rebaselined_heigth=max(first_spike_trace(4:7)); % taking a narrow window around stim (at 840) so it doesn't count the cholinergic bump when it exists
    first_spike_AUC=trapz(first_spike_trace(1:7)); %from pulse onset until end 7 points later
    
    first_spike_absolute_heigth_all(:,i-1) = first_spike_absolute_heigth;
    first_spike_F0_all(:,i-1) = first_spike_F0;
    first_spike_trace_all(:,i-1) = first_spike_trace;
    first_spike_rebaselined_heigth_all(:,i-1) = first_spike_rebaselined_heigth;
    first_spike_AUC_all(:,i-1) = first_spike_AUC;
    
    
    
    % the other stims
     
    if stims(:,i)>2; % value in cell from the "stims" input table needs to be more than 2 because for a 1p it is just written 1

        % for Lulu:
        if stims(:,i) == 2100; % 2p100Hz
        t2= 846; % 1st pulse is at 1.4s which is 840 frames, and 100Hz=every 0.0102s (small 0.2ms delay due to protocol), so the second pulse is at 1.4102s*599.52fps=845 frames but 10.2 ms/1.67ms=845.44 frames
        second_spike_F0=mean(deltaF_over_F0_percent(t2-1)); % only using 1 frame as baseline for 100 Hz stims
        third_spike_F0=0;
        fourth_spike_F0=0;
        end 

        if stims(:,i) == 3100; % 3p100Hz
        t2= 846; % 1st pulse is at 1.4s which is 840 frames, and 100Hz=every 0.0102s (small 0.2ms delay due to protocol), so the second pulse is at 1.4102s*599.52fps=845 frames but 10.2 ms/1.67ms=845.44 frames
        t3= 852; % comes at 1.42s so 1.4202s*599.52fps=851.4 frames
        second_spike_F0=mean(deltaF_over_F0_percent(t2-1)); % only using 1 frame as baseline for 100 Hz stims
        third_spike_F0=mean(deltaF_over_F0_percent(t3-1));
        fourth_spike_F0=0;
        end 

        if stims(:,i) == 450; % 4p50Hz
        t2= 852; % 1st pulse is at 1.4s which is 839 frames, and 50Hz=every 0.0202s (small 0.2ms delay due to protocol), so the second pulse is at 1.4202s*599.52fps=851.4 frames
        t3= 864; % comes at 1.44s so 1.4402s*599.52fps= 863.4 frames
        t4= 876; % comes at 1.46s so 1.4602s*599.52fps= 875.4 frames
        second_spike_F0=mean(deltaF_over_F0_percent(t2-2:t2-1));
        third_spike_F0=mean(deltaF_over_F0_percent(t3-2:t3-1));
        fourth_spike_F0=mean(deltaF_over_F0_percent(t4-2:t4-1));
        end
        
        if stims(:,i) == 4100; % 4p100Hz
        t2= 846; % 1st pulse is at 1.4s which is 840 frames, and 100Hz=every 0.0102s (small 0.2ms delay due to protocol), so the second pulse is at 1.4102s*599.52fps=845 frames but 10.2 ms/1.67ms=845.44 frames
        t3= 852; % comes at 1.42s so 1.4202s*599.52fps=851.4 frames
        t4= 858; % comes at 1.43s so 1.4302s*599.52fps=857.4 frames
        second_spike_F0=mean(deltaF_over_F0_percent(t2-1)); % only using 1 frame as baseline for 100 Hz stims
        third_spike_F0=mean(deltaF_over_F0_percent(t3-1));
        fourth_spike_F0=mean(deltaF_over_F0_percent(t4-1));
        end 
        

        % % for Beth:
        % if stims (:,i) == 210; % 2p 10ms IPI, 2p100Hz
        % t2= 845; % 1st pulse is at 1.4s which is 840 frames, and 100Hz=every 0.01s, so the second pulse is at 1.41*599.52fps=845.52 frames
        % second_spike_F0=mean(deltaF_over_F0_percent(t2-1)); % only using 1 frame as baseline for 100 Hz stims
        % end
        % 
        % if stims (:,i) == 240; % 2p 40ms IPI, 2p25Hz
        % t2= 863; % 1st pulse is at 1.4s which is 840 frames, and 25Hz=every 0.04s, so the second pulse is at 1.44*599.52fps=863.3 frames
        % second_spike_F0=mean(deltaF_over_F0_percent(t2-2:t2-1));
        % end
        % 
        % if stims (:,i) == 2100; % 2p 100ms IPI, 2p10Hz
        % t2= 899; % 1st pulse is at 1.4s which is 840 frames, and 10Hz=every 0.1s, so the second pulse is at 1.50s*599.52fps=899.28 frames
        % second_spike_F0=mean(deltaF_over_F0_percent(t2-2:t2-1));
        % end
        % 
        % if stims(:,i) == 2200; % 2p 200ms IPI, 2p5Hz
        % t2= 959; % 1st pulse is at 1.4s which is 840 frames, and 5Hz=every 0.2s, so the second pulse is at 1.60s*599.52fps=95.23 frames
        % second_spike_F0=mean(deltaF_over_F0_percent(t2-2:t2-1));
        % end 

        
        % 2nd spike
        
        second_spike_absolute_heigth=max(deltaF_over_F0_percent(t2-1:t2+3)); % taking a narrow window around stim
        % need to align the trace to stim onset at t2
        second_spike_trace=deltaF_over_F0_percent(t2-3:t2+5)- second_spike_F0;% taking a wider window and will cross later if necessary and re-baselining it
        second_spike_rebaselined_heigth=max(second_spike_trace(4:7)); % taking a narrow window around stim (at 840) so it doesn't count the cholinergic bump when it exists
        second_spike_AUC=trapz(second_spike_trace(4:9)); %from pulse onset
        
        second_spike_absolute_heigth_all(:,i-1)=second_spike_absolute_heigth;
        second_spike_F0_all(:,i-1)=second_spike_F0;
        second_spike_trace_all(:,i-1)=second_spike_trace;
        second_spike_rebaselined_heigth_all(:,i-1)=second_spike_rebaselined_heigth;
        second_spike_AUC_all(:,i-1)=second_spike_AUC;
        
        if stims(:,i) == 2100; % 2p100Hz
        continue
        end
        
        % For Beth : please put the next section in %comment since you don't use 4p stims

        % 3rd spike
        
        third_spike_absolute_heigth=max(deltaF_over_F0_percent(t3-1:t3+3)); % taking a narrow window around stim 
        % need to align the trace to stim onset at t3
        third_spike_trace=deltaF_over_F0_percent(t3-3:t3+5)- third_spike_F0;% taking a wider window and will cross later if necessary and re-baselining it
        third_spike_rebaselined_heigth=max(third_spike_trace(4:7)); % taking a narrow window around stim (at 840) so it doesn't count the cholinergic bump when it exists
        third_spike_AUC=trapz(third_spike_trace(4:9)); %from pulse onset
        
        third_spike_absolute_heigth_all(:,i-1)=third_spike_absolute_heigth;
        third_spike_F0_all(:,i-1)=third_spike_F0;
        third_spike_trace_all(:,i-1)=third_spike_trace;
        third_spike_rebaselined_heigth_all(:,i-1)=third_spike_rebaselined_heigth;
        third_spike_AUC_all(:,i-1)=third_spike_AUC;
        
        if stims(:,i) == 3100; % 3p100Hz
        continue
        end
        
        %4th spike
        
        fourth_spike_absolute_heigth=max(deltaF_over_F0_percent(t4-1:t4+3)); % taking a narrow window around stim
        % need to align the trace to stim onset at t4
        fourth_spike_trace=deltaF_over_F0_percent(t4-3:t4+5)- fourth_spike_F0;% taking a wider window and will cross later if necessary and re-baselining it
        fourth_spike_rebaselined_heigth=max(fourth_spike_trace(4:7)); % taking a narrow window around stim (at 840) so it doesn't count the cholinergic bump when it exists
        fourth_spike_AUC=trapz(fourth_spike_trace(4:9)); %from pulse onset
        
        fourth_spike_absolute_heigth_all(:,i-1)=fourth_spike_absolute_heigth;
        fourth_spike_F0_all(:,i-1)=fourth_spike_F0;
        fourth_spike_trace_all(:,i-1)=fourth_spike_trace;
        fourth_spike_rebaselined_heigth_all(:,i-1)=fourth_spike_rebaselined_heigth;
        fourth_spike_AUC_all(:,i-1)=fourth_spike_AUC;
        
               

    end 
           
    t2=0;
    t3=0;
    t4=0;
        
    

% plot raw fluorescence, pre- and post-peak values and fitted exponential decay
    figure(i-1)
    plot(frame_number,raw_F)
    hold on
    plot(Fit)
    hold on
    plot(F_slope_corrected)
    yline(F0_fit)
    hold on
    %plot(F_baselinecorrected)
    hold on
    %plot (deltaF_over_Fo(001:end)*100)
    %yline(0)
    legend('raw transient','fitted curve','slope correction','b')
    %legend('baseline corrected')
    %legend('deltaFoverFz/Fit_F0')
    %legend('baseline corrected','deltaF/F0')

    
end




%to copy paste in Excel:

% For Lulu:
SPIKE_absolute_heigth_all=cat(1,first_spike_absolute_heigth_all,second_spike_absolute_heigth_all,third_spike_absolute_heigth_all,fourth_spike_absolute_heigth_all);
SPIKE_rebaselined_heigth_all=cat(1,first_spike_rebaselined_heigth_all,second_spike_rebaselined_heigth_all,third_spike_rebaselined_heigth_all,fourth_spike_rebaselined_heigth_all);
SPIKE_AUC_all=cat(1,first_spike_AUC_all,second_spike_AUC_all,third_spike_AUC_all,fourth_spike_AUC_all);
SPIKE_F0_all=cat(1,first_spike_F0_all,second_spike_F0_all,third_spike_F0_all,fourth_spike_F0_all);

TO_EXPORT_general_AUC_F0=cat(1,total_AUC_all,F0_all);
TO_EXPORT_full_traces=cat(2,frame_number*exposure,deltaF_over_F0_percent_all);
TO_EXPORT_spike_traces=cat(1,first_spike_trace_all,zeros(1,total_number_of_stims),second_spike_trace_all,zeros(1,total_number_of_stims),third_spike_trace_all,zeros(1,total_number_of_stims),fourth_spike_trace_all);
TO_EXPORT_spike_abs_rebas_height_AUC_F0=cat(1,SPIKE_absolute_heigth_all,zeros(1,total_number_of_stims),SPIKE_rebaselined_heigth_all,zeros(1,total_number_of_stims),SPIKE_AUC_all,zeros(1,total_number_of_stims),SPIKE_F0_all);


% % For Beth:
% % total_AUC_all
% deltaF_over_F0_percent_all
% F0_all
% SPIKE_absolute_heigth_all=cat(1,first_spike_absolute_heigth_all,second_spike_absolute_heigth_all);
% SPIKE_rebaselined_heigth_all=cat(1,first_spike_rebaselined_heigth_all,second_spike_rebaselined_heigth_all);
% SPIKE_AUC_all=cat(1,first_spike_AUC_all,second_spike_AUC_all);
% SPIKE_F0_all=cat(1,first_spike_F0_all,second_spike_F0_all);
% 
% TO_EXPORT_general_AUC_F0=cat(1,total_AUC_all,F0_all)
% TO_EXPORT_full_traces=cat(2,frame_number*exposure,deltaF_over_F0_percent_all)
% TO_EXPORT_spike_traces=cat(1,first_spike_trace_all,zeros(1,total_number_of_stims),second_spike_trace_all,zeros(1,total_number_of_stims))
% TO_EXPORT_spike_abs_rebas_height_AUC_F0=cat(1,SPIKE_absolute_heigth_all,zeros(1,total_number_of_stims),SPIKE_rebaselined_heigth_all,zeros(1,total_number_of_stims),SPIKE_AUC_all,zeros(1,total_number_of_stims),SPIKE_F0_all)
% 
% 
% 
% 
