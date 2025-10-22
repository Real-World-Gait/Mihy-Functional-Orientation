function [strt,nd,fig] = gait_ID_fft(signal,Fs,windowsize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Functional sensor-to-segment orientation code provided as a supplement to
% JA Mihy, M Wagatsuma, SM Cain, JF Hafer, A functional sensor-to-segment 
% calibration method reduces the effects of varied sensor placement on 
% estimates of segment angular excursion, J Appl Biomech
%
% See notes below "Outputs" if using this function with the example data
% provided with Mihy et al.
% 
% This function determines windows of probable walking based on frequency
% analysis. Note that further processing is necessary to determine whether
% identified windows are level, in a straight line, of desired length, etc.
%
% As currently written, the function will output windows of data that
% include 2.5 second buffers added to the start and end of
% frequency-identified walking bouts. The function will also collapse any
% overlapping windows into single, larger windows
%
% Inputs:
% signal: sensor data to run detection algorithm from. Recommend using
% shank angular velocity
% Fs: data collection frequency
% windowsize: size of search window, in seconds. Recommend using 5
%
% Outputs:
% strt: start frames for windows of interest
% nd: end frames for windows of interest
% fig: figure depicting identified bouts
%
% JF Hafer, 02/2021
%
% If using with Mihy et al. example data, GaitTrials is a structure
% containing data collected during a standard in-lab gait analysis. Each
% sensor field includes raw, IMU-frame acceleration (a), angular velocity
% (w), and orientation data in quaternion format (q), as well as rotation
% matrices for transforming acceleration and angular velocity data to X and
% Z functional reference frames. Running this function on the GaitTrials
% data will provide bout start and end times that correspond to the
% frames in the Bout_time field.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set basic parameters
window = windowsize*Fs; %frames of data within set window size
da_len = length(signal); %length of raw data
da_len_use = da_len-window; %length of data to use as window start points (makes all windows the same length instead of having progressively shorter end windows)
num_win = floor(da_len_use/window); %number of consecutive windows in data
n=2^nextpow2(window);

%Initialize variables for holding power data for each axis
P11 = NaN(n/2+1,num_win);
P21 = NaN(n/2+1,num_win);
P31 = NaN(n/2+1,num_win);
for i = 1:num_win
   
   temp = signal((i*window-window+1):(i*window),1); %first axis of data
   Y=fft(temp,n,1);
   Pb = abs(Y/window);
   Pa = Pb(1:n/2+1,:);
   Pa(2:end-1,:) = 2*Pa(2:end-1,:);
   P11(:,i) = Pa;  
   
   temp = signal((i*window-window+1):(i*window),2); %second axis of data
   Y=fft(temp,n,1);
   Pd = abs(Y/window);
   Pc = Pd(1:n/2+1,:);
   Pc(2:end-1,:) = 2*Pc(2:end-1,:);
   P21(:,i) = Pc;  
   
   temp = signal((i*window-window+1):(i*window),3); %third axis of data
   Y=fft(temp,n,1);
   Pf = abs(Y/window);
   Pe = Pf(1:n/2+1,:);
   Pe(2:end-1,:) = 2*Pe(2:end-1,:);
   P31(:,i) = Pe;  
    
end

Pall(:,:) = P11+P21+P31; %add frequency content across axes

scale = max(Pall,[],'all'); %quantity to scale magnitude search by below

%find windows that contain a frequency power > 1/3 the maximum power of the
%data
scalemag = scale/3;
Q = NaN(n/2,num_win);
keep = NaN(1,num_win);
freq = 0:(Fs/n):(Fs/2-Fs/n);
Pfreq = Pall(1:n/2,:);
x = size(Pfreq,1);

for j = 1:num_win
   
    %find values of frequency peaks within each P column
    %keep windows that have frequency >0.5 and <2.2 Hz
    [pks,locs] = findpeaks(Pfreq(:,j),freq,'MinPeakProminence',scalemag);
    
   locsfreq = double(and(locs>=0.5,locs<=2.2)); %are there peaks between 0.5 and 2.2 Hz?
   freqlog = mean(locsfreq); %if the location of any frequencies is in desired range, average will be >0
   if and(max(Pfreq(:,j)>scalemag),freqlog>0)
       Q(:,j) = Pfreq(:,j);
       keep(:,j) = j;
   else
       Q(:,j) = NaN(x,1);
       keep(:,j) = 0;
   end
    
end

%find column indices of start and end of windows of interest
keep_ind = find(keep); %identify all frames within windows

%If windows not identified, exit function
if length(keep_ind)==1
    strt = 0;
    nd = 0;
    fig = 0;
    return
end

%Identify windows to keep
for k = 1:(length(keep_ind))

    if k == 1
        keep_win(:,k) = keep_ind(1); %beginning of first window
            ktemp = keep_ind(:,k+1) - keep_ind(:,k); %gap between current and next keep windows
          if (ktemp ~= 1)%singleton frame within windows
          sngltn(:,k) = keep_ind(k);    
          end
    elseif k == (length(keep_ind))
        keep_win(:,length(keep_ind)) = keep_ind(end); %end of last window
        ktemp2 = keep_ind(:,k) - keep_ind(:,k-1); %gap between current and previous keep windows
          if (ktemp2 ~= 1)%singleton frames within windows
          sngltn(:,k) = keep_ind(k);    
          end
    else
    ktemp = keep_ind(:,k+1) - keep_ind(:,k); %gap between current and next keep windows
    ktemp2 = keep_ind(:,k) - keep_ind(:,k-1); %gap between current and previous keep windows
    if and((ktemp == 1),(ktemp2 == 1)) %frame within window
        keep_win(:,k) = 0; %don't keep as independent window
        elseif and((ktemp ~= 1),(ktemp2 ~= 1))%singleton window
            sngltn(:,k) = keep_ind(k);
            keep_win(:,k) = keep_ind(k);
    elseif or((ktemp ~= 1),(ktemp2 ~= 1)) %end of window
        keep_win(:,k) = keep_ind(k);
    end

    end

end

%If windows to keep not identified, exit function
if ~exist('keep_win','var')
    strt = 0;
    nd = 0;
    fig = 0;
    return
end

%Find frames of window borders
[row,col,keep_win_ind] = find(keep_win);
if exist('sngltn') == 1
[row,col,sngltn_win_ind] = find(sngltn);
sngltn_win_frame = sngltn_win_ind*window-window+1;
end
keep_win_frame = keep_win_ind*window-window+1; %frames of identifying window borders


%split window frames into start frames and end frames
strt = []; nd = [];
for m = 1:length(keep_win_frame)

   temp = keep_win_frame(m);
   if exist('sngltn') == 1
    if ismember(temp,sngltn_win_frame) %if frame identifies a singleton window, frame is start, frame + window length is end
       strt(m) = temp;
       nd(m) = temp+window;
    elseif length(strt) == length(nd) %if not a singleton and length of start and end indices are already equal, this must be a new start
       strt(m) = temp;
       strt(m+1) = 0;
    elseif length(strt) > length(nd) %if not a singleton and length of start is longer than length of end indices, this must be the end of the previous start
       nd(m) = temp+window; %frame index indicates beginning of last window; add one window length
    end
   else
    if length(strt) == length(nd) %if not a singleton and length of start and end indices are already equal, this must be a new start
       strt(m) = temp;
       strt(m+1) = 0;
    elseif length(strt) > length(nd) %if not a singleton and length of start is longer than length of end indices, this must be the end of the previous start
       nd(m) = temp+window; %frame index indicates beginning of last window; add one window length   
    end
   end

end

strt = strt(strt~=0); nd = nd(nd~=0); %remove placeholder 0s

%Add 2.5 second buffers to beginning and ends of windows
strt = strt-2.5*Fs;
nd = nd+2.5*Fs;

%collapse overlapping windows
for i = 2:length(strt)

    if strt(i) <= nd(i-1)
        strt(i) = 0;
        nd(i-1) = 0;
    end

end

strt = strt(strt~=0); nd = nd(nd~=0); %remove placeholder 0s

%Correct for buffer of first window going into negative time
if strt(1) < 1
    strt(1) = 1;
end

%Plot start and end identifiers on data
fig = figure();
plot(signal,'k','Linewidth',0.02)
hold on
xline(strt,'b')
xline(nd,'r')
xlabel('Frames')
ylabel('Angular velocity')
title('Probable walking windows. Blue = start, red = end')

end