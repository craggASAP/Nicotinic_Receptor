%%%%%%%%%%%%%%%%%%% Imaging data analysis %%%%%%%%%%%%%%%%%%%%
% !! Remember to modify parameters depending on acquisition protocol and stim protocol!!
% (exposure / time of stim / etc..)


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

first_spike_F0_all=zeros(1,total_number_of_stims);
second_spike_F0_all=zeros(1,total_number_of_stims);


first_spike_rebaselined_heigth_all=zeros(1,total_number_of_stims);
second_spike_rebaselined_heigth_all=zeros(1,total_number_of_stims);


first_spike_AUC_all=zeros(1,total_number_of_stims);
second_spike_AUC_all=zeros(1,total_number_of_stims);


first_spike_trace_all=zeros(9,total_number_of_stims);
second_spike_trace_all=zeros(9,total_number_of_stims);


first_stim_trace_all=zeros(71,total_number_of_stims);
second_stim_trace_all=zeros(71,total_number_of_stims);


first_stim_AUC_all=zeros(1,total_number_of_stims);
second_stim_AUC_all=zeros(1,total_number_of_stims);



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
    first_spike_absolute_heigth=max(deltaF_over_F0_percent(t1-1:t1+2)); % taking a narrow window around stim (at 840) so it doesn't count the cholinergic bump when it exists
    % need to align the trace to stim onset at 840
    first_spike_F0=mean(deltaF_over_F0_percent(t1-2:t1-1)); % F0 is an average of the 2 points before time of stim t1 (at t1-2 and t1-1) 
    first_spike_trace=deltaF_over_F0_percent(t1-3:t1+5)- first_spike_F0;% taking a wider window and will cross later if necessary and re-baselining it taking 3 points before the time of stim t1 anf 5 points after
    first_spike_rebaselined_heigth=max(first_spike_trace(4:6)); % taking a narrow window around stim (at 840) so it doesn't count the cholinergic bump when it exists
    first_spike_AUC=trapz(first_spike_trace(4:6)); %from pulse onset until end 7 points later
    
    first_stim_trace_all(:,i-1)= deltaF_over_F0_percent(t1-10:t1+60)-first_spike_F0;
    first_stim_AUC_all(:,i-1)= trapz(first_stim_trace_all(10:end,i-1));

    first_spike_absolute_heigth_all(:,i-1) = first_spike_absolute_heigth;
    first_spike_F0_all(:,i-1) = first_spike_F0;
    first_spike_trace_all(:,i-1) = first_spike_trace;
    first_spike_rebaselined_heigth_all(:,i-1) = first_spike_rebaselined_heigth;
    first_spike_AUC_all(:,i-1) = first_spike_AUC;
    
    
    
    % the other stims
     
    if stims(:,i)>2; % value in cell from the "stims" input table needs to be more than 2 because for a 1p it is just written 1

        if stims(:,i) == 10; % 2p 10 ms IPI
        t2= 846; % 1st pulse is at 1.4s which is 840 frames, add 10 ms , so the second pulse is at 1.41s*599.52fps=845.3232 frames
        end 
        
        if stims(:,i) == 25; % 2p 25 ms IPI
        t2= 854; % 1st pulse is at 1.4s which is 840 frames, add 25 ms , so the second pulse is at 1.425s*599.52fps=854.316 frames 
        end 

        if stims(:,i) == 40; % 2p 40 ms IPI
        t2= 864; % 1st pulse is at 1.4s which is 840 frames, add 40 ms, so the second pulse is at 1.44*599.52fps=863.3088 frames 
        end 

        if stims(:,i) == 80; % 2p 80 ms IPI
        t2= 888; % 1st pulse is at 1.4s which is 840 frames, add 80 ms, so the second pulse is at 1.48s*599.52fps=887.2896 frames 
        end 

        if stims(:,i) == 100; % 2p 100 ms IPI
        t2= 900; % 1st pulse is at 1.4s which is 840 frames, add 100 ms, so the second pulse is at 1.5s*599.52fps=899.28 frames 
        end 

        if stims(:,i) == 150; % 2p 150 ms IPI
        t2= 930; % 1st pulse is at 1.4s which is 840 frames, add 150 ms, so the second pulse is at 1.55s*599.52fps=929.256 frames
        end
        
        if stims(:,i) == 200; % 2p 200 ms IPI
        t2= 960; % 1st pulse is at 1.4s which is 840 frames, add 200 ms, so the second pulse is at 1.6s*599.52fps=959.232 frames
        end 

        if stims(:,i) == 250; % 2p 200 ms IPI
        t2= 990; % 1st pulse is at 1.4s which is 840 frames, add 250 ms, so the second pulse is at 1.65s*599.52fps=989.208 frames
        end

        if stims(:,i) == 500; % 2p 500ms IPI
        t2= 1140; % 1st pulse is at 1.4s which is 840 frames, add 0.5002s (small 0.2ms delay due to protocol), so the second pulse is at 1.9002s*599.52fps=1139.2079 frames 
        end 
        
        if stims(:,i) == 1000; % 2p 1s IPI
        t2= 1440; % 1st pulse is at 1.4s which is 840 frames, add 1000ms (small 0.2ms delay due to protocol), so the second pulse is at 2.4002s*599.52fps=1438.9679 frames 
        end 

        if stims(:,i) == 1500; % 2p 1.5s IPI
        t2= 1739; % 1st pulse is at 1.4s which is 840 frames, add 1500ms (small 0.2ms delay due to protocol), so the second pulse is at 2.9002s*599.52fps=1738.7279 frames
        end
        
        % 2nd spike
        
        second_spike_F0=mean(deltaF_over_F0_percent(t2-2:t2-1)); % only using 10 frames (1 frame) as baseline

        second_spike_absolute_heigth=max(deltaF_over_F0_percent(t2-1:t2+2)); % taking a narrow window around stim
        % need to align the trace to stim onset at t2
        second_spike_trace=deltaF_over_F0_percent(t2-3:t2+5)- second_spike_F0;% taking a wider window and will cross later if necessary and re-baselining it
        second_spike_rebaselined_heigth=max(second_spike_trace(4:6)); % taking a narrow window around stim (at 840) so it doesn't count the cholinergic bump when it exists
        second_spike_AUC=trapz(second_spike_trace(4:6)); %from pulse onset
        
        second_stim_trace_all(:,i-1)= deltaF_over_F0_percent(t2-10:t2+60)-second_spike_F0;
        second_stim_AUC_all(:,i-1)= trapz(second_stim_trace_all(10:end,i-1)); %from 2nd pulse onset to 100ms (60 frames) later
        
        second_spike_absolute_heigth_all(:,i-1)=second_spike_absolute_heigth;
        second_spike_F0_all(:,i-1)=second_spike_F0;
        second_spike_trace_all(:,i-1)=second_spike_trace;
        second_spike_rebaselined_heigth_all(:,i-1)=second_spike_rebaselined_heigth;
        second_spike_AUC_all(:,i-1)=second_spike_AUC;
   
        
       
    end 
           
    
 
        
    

% % plot raw fluorescence, pre- and post-peak values and fitted exponential decay
%     figure(i-1)
%     plot(frame_number,raw_F)
%     hold on
%     plot(Fit)
%     hold on
%     plot(F_slope_corrected)
%     yline(F0_fit)
%     hold on
%     %plot(F_baselinecorrected)
%     hold on
%     %plot (deltaF_over_Fo(001:end)*100)
%     %yline(0)
%     legend('raw transient','fitted curve','slope correction','b')
%     %legend('baseline corrected')
%     %legend('deltaFoverFz/Fit_F0')
%     %legend('baseline corrected','deltaF/F0')

    
end




%out put to copy paste in Excel:

total_AUC_all;
deltaF_over_F0_percent_all;
F0_all;
SPIKE_absolute_heigth_all=cat(1,first_spike_absolute_heigth_all,second_spike_absolute_heigth_all);
SPIKE_rebaselined_heigth_all=cat(1,first_spike_rebaselined_heigth_all,second_spike_rebaselined_heigth_all);
SPIKE_AUC_all=cat(1,first_spike_AUC_all,second_spike_AUC_all);
SPIKE_F0_all=cat(1,first_spike_F0_all,second_spike_F0_all);

TO_EXPORT_general_AUC_F0=cat(1,total_AUC_all,F0_all);
TO_EXPORT_full_traces=cat(2,frame_number*exposure,deltaF_over_F0_percent_all);
TO_EXPORT_spike_traces=cat(1,first_spike_trace_all,zeros(1,total_number_of_stims),second_spike_trace_all,zeros(1,total_number_of_stims));
TO_EXPORT_stim_traces=cat(1,first_stim_trace_all,zeros(1,total_number_of_stims),second_stim_trace_all,zeros(1,total_number_of_stims));
TO_EXPORT_stim_AUCs=cat(1,first_stim_AUC_all,second_stim_AUC_all);
TO_EXPORT_spike_abs_rebas_height_AUC_F0=cat(1,SPIKE_absolute_heigth_all,zeros(1,total_number_of_stims),SPIKE_rebaselined_heigth_all,zeros(1,total_number_of_stims),SPIKE_AUC_all,zeros(1,total_number_of_stims),SPIKE_F0_all);

 
