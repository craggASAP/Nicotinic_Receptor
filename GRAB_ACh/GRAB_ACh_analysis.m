%%%%%%%%%%%%%%%%%%%% Imaging data analysis %%%%%%%%%%%%%%%%%%%%
% Remember to modify parameters depending on acquisition protocol
% (exposure / time of stim / etc..)
%  Lucille Duquenoy 03/05/24 
% 

% input A with first a column with the frame number (Column 1), and then the raw fluorescence transients in the next columes

deltaF_all=zeros(size(A,1),size(A,2)-1);

for i=2:size(A,2)
    x=A(:,1); %frames
    exposure=1/600; %exposure every 100ms (frequency of recording 10Hz)
    time=(x.*exposure)./100; %total recording time in seconds 
    y=A(:,i); %fluorescence values

    %create fit
    Fit=fit(A(1:90,1),y(1:90),'log'); %here it is set for a stim arriving at frame 100. you can change this parameter depending on your stim protocol.
    %Fit=fit(A(5:38,1),y(5:38),'log'); %here it is set for a stim arriving at frame 50. you can change this parameter depending on your stim protocol.
    coefs=coeffvalues(Fit);

    %need to subtract the equation obtained from Fit function from the original trace to correct for slope 
    fitted_transient=y-(coefs(1)*log(x));
    %F0= mean(fitted_transient(250:398,:)); % F0 = avreage of few frames before the stim
    F0=coefs(2); % here we choose F0=b=coefs(2) in a*log(x)+b
    F0_all(:,i-1)=F0;
    baselinecorrected=(fitted_transient-F0); %actual calculation of deltaF 

    % Calculate F-F0/F0
    deltaF=baselinecorrected/F0;
    deltaF_percent_all(:,i-1)=deltaF.*100;
    
    % Find peak value
    forpeak=deltaF(99:250);
    %forpeak=deltaF(99:200); 
    peak=max(forpeak);
    peak_percentage=peak.*100;
    peak_percent_all(:,i-1)=peak_percentage;
    
    % Calculate Area under the curve AUC from where stimulation starts to end of transient
    AUC=trapz(deltaF(99:250)); % not sure when to stop counting the AUC?
    %AUC=trapz(deltaF(38:100)); % not sure when to stop counting the AUC?
    AUC_all(:,i-1)=AUC;

   %plot raw fluorescence, pre- and post-peak values and fitted exponential decay
    figure(i-1)
    plot(x,y)
    hold on
    plot(Fit)
    hold on
    plot(fitted_transient)
    yline(coefs(2))
    hold on
    %plot(baselinecorrected)
    hold on
    %plot (deltaF(001:end)*100)
    %yline(0)
    legend('raw transient','fitted curve','slope correction','b')
    %legend('baseline corrected')
    %legend('deltaF/F0')
    %legend('baseline corrected','deltaF/F0')
    
    
end
