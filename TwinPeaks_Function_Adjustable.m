function [swd_number,percent_time_swd,time_grid,final_selection] = TwinPeaks_Function(cx,SR,tw,ovl,lims,limppks,limPks)

%Written to search for spike-wave discharges by Ivaylo Iotchev, special thanks to
%Gilles van Luijtelaar, Dmitry Perevozniuk, Evgenia Sitnikova.

%script identifies the EEG trace by the name 'cx'
%script identifies the sampling rate by the name 'SR'
%script identifies time window length by the input 'tw', default = 1
%script identifies the overlap of time-windows requested for STFFT by the
%input 'ovl', default = 0.5
%script identifies the threshold for how many time-domain sharp peaks need
%be observed by the input 'lims', default = 2
%script identifies the threshold for how many power-peaks above 0.025 need
%be observed by the input 'limppks', default = 3
%script identifies the threshold for how many power-peaks above 0.05 need
%be observed by the input 'limPks', default = 2

if nargin < 3
    tw = 1;
    ovl = 0.5;
    lims = 2;
    limppks = 3;
    limPks = 2;
end

if nargin < 4
    ovl = 0.5;
    lims = 2;
    limppks = 3;
    limPks = 2;
end

if nargin < 5
    lims = 2;
    limppks = 3;
    limPks = 2;
end

if nargin < 6
    limppks = 3;
    limPks = 2;
end

if nargin < 7
    limPks = 2;
end

cx = zscore(cx - mean(cx)); %demeaning


cxraw = cx; %save raw version of the trace

Nyquist = SR/2; %Using SR (= Sampling Rate, provided with the file) calculate the Nyquist frequency

window = SR*tw; %Define the size of the signal segments that will be used for a sliding fft, as a magnitude of sampling rates (results in 0.5 seconds window)

normalized_cutoff = 100/Nyquist; %passband low-pass

normalized_cutoff2 = 150/Nyquist; %stopband low-pass

highpass_cutoff2 = 7/Nyquist; %passband high-pass

highpass_cutoff = 5/Nyquist; %stopband high-pass

[order,nfreq] = buttord(normalized_cutoff,normalized_cutoff2,0.5,30); %establish the best order and cut-off frequency for a low-pass filter that attenuates with 30 dB and allows for max 0.5 dB ripple

[order2,nfreq2] = buttord(highpass_cutoff2,highpass_cutoff,0.5,30); %establish the best order and cut-off frequency for a high-pass filter that attenuates with 30 dB and allows for max 0.5 dB ripple

[b,a,k1] = butter(order,nfreq,'low'); %obtaining filter coefficients with the order and cut-off specified for the low-pass filter
 
[c,d,k2] = butter(order2,nfreq2,'high'); %obtaining filter coefficients with the order and cut-off specified for the high-pass filter 

sos1 = zp2sos(b,a,k1);

sos2 = zp2sos(c,d,k2);

cx = sosfilt(sos2,cx); %transform the trace with the high-pass filter

clear b a c d rl rl2 rl3 rh rh2 rh3 %clear the filter coefficients

[result1, frequencymap1, timepoints1, power1] = spectrogram(cxraw,window,round(SR*ovl),SR*10,SR); %obtain a matrix of amplitudes, the corresponding frequencies, timepoints and power for a sliding fft with adjustable window-size (using predefined window value, see above), window-overlapp of x seconds (specified by input 'ovl'), and zero-padding for 0.1 Hz resolution (10*SR)

result1 = abs(result1');
frequencymap1 = frequencymap1';
timepoints1 = timepoints1';
power1 = power1'; %align all results to ease later indexing


[r,c] = size(result1); %obtain size parameters of the results for subsequent loop-construction

   for x = 1:r %loop repeats for the number of windows used in the fft
       
       szelet = (frequencymap1 >= 5); %originally at this spot only frequencies above 2 Hz were selected for further analysis, to avoid the estimation of the max frequency power be influenced by the DC-offset, however I later learned high-pass filtering deals with it, so the value is set to zero now
       szeletmap = frequencymap1(frequencymap1 >= 5); %here the frequency-map is adjusted for the selection above, which is now simply above 0 (but was above 2 Hz)
       selectresult = power1(x,szelet); %the power-values of each window are collected, originally only for the selected frequency search radius (> 2 Hz)
       a = x - 1;
       b = x + 1;
       if x == 1
           a = x + 2;
       end
       if x == r
           b = x - 2;
       end
       selectresult0 = power1(a,szelet);
       selectresult1 = power1(b,szelet);
       [amp,pos] = max(selectresult); %amplitude and position on the frequency-map of the highest power is selected
       maxfreq1(x) = pos; %for each window the frequency with the highest power is collected in a vector
       maxpower1 = (szeletmap <= 50); %limit space of frequency-domain search to 50 Hz
       powerrange0 = selectresult0(maxpower1);
       powerrange = selectresult(maxpower1);
       powerrange1 = selectresult1(maxpower1);
       powerpeaks0 = findpeaks(powerrange0);
       powerpeaks = findpeaks(powerrange);
       powerpeaks1 = findpeaks(powerrange1);
       i = sum(powerpeaks0 > 0.025);
       l = sum(powerpeaks0 > 0.05);
       v = sum(powerpeaks > 0.025);
       e = sum(powerpeaks > 0.05);
       o = sum(powerpeaks1 > 0.025);
       n = sum(powerpeaks1 > 0.05);
       ivo = [i v o]; %collect number of power peaks > 0.025 for 3 neighboring time-windows
       len = [l e n]; %collect number of power peaks > 0.05 for 3 neighboring time-windows
       ppks(x) = max(ivo); %make sure within the duration of 3 time-windows the sought for criterion is present at least once
       Pks(x) = max(len); %make sure within the duration of 3 time-windows the sought for criterion is present at least once
       clear amp pos selectresult selectresult0 selectresult1 maxpower1 powerrange powerrange0 powerrange1 powerpeaks powerpeaks0 powerpeaks1 a b i v o ivo l e n len %clear variables for a new loop
       
   end
   

   

      
      [estimatesamp,pci_amp] = mle(cxraw); % estimatesamp contains estimates of the true mean(1) and standard deviation(2) of the swd-candidate-amplitudes
     
      mps = estimatesamp(1); %label the mean estimate of amplitudes
      varps = estimatesamp(2); %label the std estimate of amplitudes
      
      varps = 3*varps;
      
      pilestone = mps + varps;

   

   
  forcrd = timepoints1*SR; %reconstruct the sample positions of the time-points :D
    
    
    for e = 1:length(timepoints1); %go through all time-points
        crd = forcrd(e); %find each corresponding sample
        crd1 = round(crd - (0.75*SR)); %select a point that is half of the fft-window (in samples) below the time-point (which is itself the center of the window)
        crd2 = round(crd + (0.75*SR)); %select a point that is half of the fft-window (in samples) above the time-point (which is itself the center of the window)
        crdx = round(crd - (0.75*SR));
        crdy = round(crd + (0.75*SR));
        
        if crd1 <= 0 %if the lower boundary of the real window is shorter than expected (the center is not symmetrical)
          crd1 = 1; %let the beginning of the whole trace define the beginning of the window
        end
         if crdx <= 0 %if the lower boundary of the real window is shorter than expected (the center is not symmetrical)
          crdx = 1; %let the beginning of the whole trace define the beginning of the window
        end
        if crd2 > length(cx); %if the upper boundary of the window is longer than expected (signal ends less than 0.75 seconds after the center of the window)
            crd2 = length(cx); %upper boundary of the window coincides with the end of the trace
        end
        if crdy > length(cx); %if the upper boundary of the window is longer than expected (signal ends less than 0.75 seconds after the center of the window)
            crdy = length(cx); %upper boundary of the window coincides with the end of the trace
        end
        
        ampest1(e) = rms(cx(crd1:crd2)); %calculate the root mean square (rms) for the entire window as an amplitude estimate
        spitze(e) = max(diff(abs(cx(crdx:crdy))));
        slice = cx(crd1:crd2);
        slice = slice - mean(slice);
        slice = slice*-1;
        speaks = findpeaks(slice);
        scount(e) = sum(speaks > pilestone);
        artcount(e) = sum(speaks > 10);
        
        clear crd crd1 crd2 crdx crdy slice speaks %clear variables for a new loop
    end     
    
  
  ampest1 = ampest1'; %align amplitude-estimates for better indexing
  spitze = spitze';
  scount = scount';
  artcount = artcount';
  ppks = ppks';
  Pks = Pks';
  spitze = zscore(spitze);

  maxfreq1 = maxfreq1'; 
 
  
  
  selector1 = (scount >= lims) & (ppks >= limppks) & (Pks >= limPks) & (artcount == 0); %search criteria, first round
  
  happening_times1 = timepoints1(selector1 == 1,:);
  
  h2_times = happening_times1;
  
fork = h2_times*SR;
fork = fork';

spikycx = cx(:,fork);


[dauan,daada] = mle(spikycx);

mz = dauan(1);
sz = dauan(2);

sz = 2*sz;

pstonex = mz + sz;
    
    for e = 1:length(timepoints1); %go through all time-points
        crd = forcrd(e); %find each corresponding sample
        crd1 = round(crd - (0.75*SR)); %select a point that is half of the fft-window (in samples) below the time-point (which is itself the center of the window)
        crd2 = round(crd + (0.75*SR)); %select a point that is half of the fft-window (in samples) above the time-point (which is itself the center of the window)
        crdx = round(crd - (0.75*SR));
        crdy = round(crd + (0.75*SR));
        
        if crd1 <= 0 %if the lower boundary of the real window is shorter than expected (the center is not symmetrical)
          crd1 = 1; %let the beginning of the whole trace define the beginning of the window
        end
         if crdx <= 0 %if the lower boundary of the real window is shorter than expected (the center is not symmetrical)
          crdx = 1; %let the beginning of the whole trace define the beginning of the window
        end
        if crd2 > length(cx); %if the upper boundary of the window is longer than expected (signal ends less than 0.75 seconds after the center of the window)
            crd2 = length(cx); %upper boundary of the window coincides with the end of the trace
        end
        if crdy > length(cx); %if the upper boundary of the window is longer than expected (signal ends less than 0.75 seconds after the center of the window)
            crdy = length(cx); %upper boundary of the window coincides with the end of the trace
        end
        
        ampest1(e) = rms(cx(crd1:crd2)); %calculate the root mean square (rms) for the entire window as an amplitude estimate
        spitze(e) = max(diff(abs(cx(crdx:crdy))));
        slice = cx(crd1:crd2);
        slice = slice - mean(slice);
        slice = slice*-1;
        speaks = findpeaks(slice);
        scounto(e) = sum(speaks > pstonex);
        
        clear crd crd1 crd2 crdx crdy slice speaks %clear variables for a new loop
    end 
    

    

  
    scounto = scounto';
    
  selector3 = (scounto >= lims) & (ppks >= limppks) & (Pks >= limPks) & (artcount == 0); % search criteria, second round
  
  happening_times3 = timepoints1(selector3 == 1,:);
  
  hx_times = happening_times3;
  
    
    
    

    
  for t = 1:length(hx_times) %plotting loop
      
      mark1 = hx_times(t); %go through all time-points corresponding to detections as defined after the second round of selection
      dvalue = mark1*SR; %find the corresponding sample in the original signal to each time-point
      svalue = dvalue + (SR*3); %extend 1.5 seconds from the time-point of detection and use the correpsonding sample for the upper boundary of the plotting window
      evalue = dvalue - (SR*3); %go 1.5 seconds down from the time-point of detection and use the correpsonding sample for the lower boundary of the plotting window
      if evalue <= 0 %if the time-point of detection is less than 1.5 seconds from the start of the signal...
          evalue = 1; %...then start plotting from the beginning of the signal
      elseif svalue >= length(cxraw) %if the signal ends less than 1.5 seconds from the time-point of detection...
          svalue = length(cxraw); %...then only plot as far as the signal reaches
      end
      xax = 0:1/SR:6; %a scale in seconds for the x-axis
      xax = xax - 3; %let the detection be time-point zero on the scale
      plotsignal = cxraw(evalue:svalue); %determine the part of the signal you want to plot

      if length(xax) < length(plotsignal)||length(xax) > length(plotsignal) %if the portion you want to plot is longer or shorter than the 3 second time-window around detection...
          subplot(2,1,1)
          plot(plotsignal) %...then just plot without scaling

      else
      subplot(2,1,1)
      plot(xax,plotsignal) %if scale and signal align always scale
      axis([-3 3 -5 5]) %limit the y-axis to diplay the signal only between 20 and -20 µV
      title('cortex')
      subplot(2,1,2)
      sval = dvalue + (SR*3); %extend 1.5 seconds from the time-point of detection and use the correpsonding sample for the upper boundary of the plotting window
      eval = dvalue - (SR*3); %go 1.5 seconds down from the time-point of detection and use the correpsonding sample for the lower boundary of the plotting window
      if eval <= 0 %if the time-point of detection is less than 1.5 seconds from the start of the signal...
          eval = 1; %...then start plotting from the beginning of the signal
      elseif sval >= length(cxraw) %if the signal ends less than 1.5 seconds from the time-point of detection...
          sval = length(cxraw); %...then only plot as far as the signal reaches
      end
      plotsignal2 = cxraw(eval:sval); %determine the part of the signal you want to plot
      [result3, frequencymap3, timepoints3, power3] = spectrogram(plotsignal2,window,round(SR*0.5),SR*10,SR);
      plotsignal2 = power3;
      plotsignal2 = plotsignal2(1:500,:);
      frequencymap3 = frequencymap3(1:500,:);
      timepoints3 = timepoints3 - 3;
      imagesc(timepoints3,frequencymap3,plotsignal2,[0 0.1]);
      hold on
      c = colorbar;
      xlabel 'Time (seconds)'
      ylabel 'Frequency (Hz)'
      ylabel(c, 'Power (W/Hz)')
      title(['time at zero ' num2str(mark1)])
      hold off

      end
      pause %pause after showing every detection-image
      clear mark1 dvalue svalue evalue %eliminate variables that could interfer with the loop
      
  end     
  
  
 HD = diff(hx_times); %calculate a vector of distances between detection times
 UjHD = HD > 2; %count 1 separate swd event for every two detections that are more than 2 seconds apart
 swd_number = sum(UjHD); %sum the previous ones to estimate the number of swd

 if (length(hx_times)) - 1 == sum(HD > 0) %if the number of detections and events estimated as separate is equal...
     swd_number = swd_number + 1 %...take the number of estimated separate events and add a one to adjust for the shorter diff-vector
 elseif sum(HD > 0) == 1 %if there is only one 1 in UjHD, this still means 2 separate swd
     swd_number = 2
 elseif length(hx_times) > 0 & sum(HD > 0) == 0 %if no event is marked as separate, but the detections are larger zero...
     swd_number = 1 %this is still one swd
 else
     swd_number = swd_number %if none of above complications occur, swd_number is a good enough estimate :)
 end

 
 percent_time_swd = (sum(selector3)/length(selector3))*100 %total duration swd as % of total duration recording 
 time_grid = timepoints1;
 final_selection = selector3;
