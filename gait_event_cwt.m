function [TO,HS,ftd] = gait_event_cwt(signal_a,signal_w,filter,Fs,LPfreq)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Functional sensor-to-segment orientation code provided as a supplement to
% JA Mihy, M Wagatsuma, SM Cain, JF Hafer, A functional sensor-to-segment 
% calibration method reduces the effects of varied sensor placement on 
% estimates of segment angular excursion, J Appl Biomech
%
% See notes below "Output" if using this function with the example data
% provided with Mihy et al.
%
% This function provides frames of toe off and heel strike events using a
% continuous wavelet transform technique.
%
% Input:
% signal_a: acceleration data from foot sensor in world reference frame
% signal_w: angular velocity data from foot sensor in world reference frame
% filter: 1 or []. 1 indicates use a lowpass filter on CWT output. If
% filter is 1, you must enter values for Fs and LPfreq
% Fs: sample frequency
% LPfreq: desired lowpass frequency cutoff
%
% Output:
% TO: frames of toe-offs
% HS: frames of heel strikes
% ftd: foot displacement from ZUPT procedures
%
% Required functions:
% find_zero_velocity.m
% zupt_displacement.m
%
% Code developed by JF Hafer, 2018
% 
% If using with Mihy et al. example data, BoutExample is a structure
% containing a single bout of walking extracted using gait_ID_fft.m for use
% identifying gait events and performing other example procedures.
% gait_events field includes the frames of gait events within strides
% (column 1 = heel strike 1, column 2 = toe-off, column 3 = heel strike 2).
% Spatiotemporal variables were derived from each stride using zero
% velocity procedures (see additional functions). Each
% sensor field includes raw, IMU-frame acceleration (a), angular velocity
% (w), orientation data in quaternion format (q), rotation
% matrices for transforming acceleration and angular velocity data to X and
% Z functional reference frames, acceleration and angular velocity data
% transformed to X and Z functional reference frames (_bodyX and _bodyZ
% variables), quaternion data transformed to rotation matrix format (or
% variable), and acceleration and angular velocity data transformed to
% world reference frame using the or variable (_world variables). Running
% this function on the BoutExample data will provide frames of gait events
% in the gait_event variable.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Continuous wavelet transform procedure
%Identify probable gait events based on high power in high-frequency data
cfs = cwt(signal_a(:,3));

if filter == 1
%filter signal
[b,a] = butter(2,LPfreq/(Fs/2)); %second-order low pass filter at whatever is entered for LPfreq
sig = filtfilt(b,a,abs(cfs(1,:)));
else
    sig = abs(cfs(1,:));
end

%Determine peak magnitudes to adjust gait event detection threshold
[m, t] = findpeaks(sig);
threshold = median(m)*0.5; %median(m)*0.2;

%Gait event magnitude (gem) and time (get) based on filtered cwt peaks
[gem, get] = findpeaks(sig,'MinPeakProminence',threshold); %use threshold


%% Velocity and displacement around gait events for event context
%Determine slope of foot displacement data at gait events

%first, ZUPT procedures
[zupt_index] = find_zero_velocity(signal_a,signal_w,1,1);
ndtime = length(signal_a)*1/Fs-1/Fs;
t = linspace(0,ndtime,length(signal_a))';

%remove gravity from signal_a
signal_a_gravcor = signal_a;
signal_a_gravcor(:,3) = signal_a(:,3) - 9.81;

[ftv,ftd] = zupt_displacement(signal_a_gravcor,t,zupt_index); %foot velocity (ftv) and displacement (ftd)

%slope of vertical displacement at gait event; central difference method
for d = 1:(length(get)-1)
ges(d) = (ftd(get(d)+1,3) - ftd(get(d)-1,3))/3; 
end

if sum(ges) == 0 %If no gait events, stop code
    figure() %make blank figure
    TO = nan;
    HS = nan;
else

%% Define HSs and TOs
%get indices of TOs and HSs depending on slope of vertical displacement
%also pull slope
%NOTE: This code is meant to work for relatively level gait
for e = 1:(length(get)-1)
    if ges(e)>0 %TOs generally have positive slope/vertical velocity
        TO(e)=get(e);
        gesTO(e) = ges(e);
    else %HSs generally have negative slope/vertical velocity
        HS(e)=get(e);
        gesHS(e) = ges(e);
    end
end

%above loop leaves 0s at loop iterations between gait events. Remove 0
%cells
for k = 1:length(TO)
    if gesTO(k) == 0
        TO(k)=NaN;
        gesTO(k)=NaN;
    end
end

TO=(TO(~isnan(TO)));
gesTO=(gesTO(~isnan(gesTO)));

for l = 1:length(HS)
    if gesHS(l) == 0
        HS(l)=NaN;
        gesHS(l)=NaN;
    end
end

HS=(HS(~isnan(HS)));
gesHS=(gesHS(~isnan(gesHS)));

%Remove common sources of incorrect gait events (repeated TOs, HSs)
%NOTE: These procedures will not correct all errors
for f = 1:length(HS)
    if gesHS(f)==0 %if heel strike on flat slope
            HS(f) = NaN; %delete
    elseif f==1 %if first HS
        if mean(and(TO>HS(f),TO<HS(f+1)))==0 %if two HSs in a row with no TO in between
            if or(gesHS(f)>gesHS(f+1),gesHS(f)==0)  %and on a flat
            HS(f) = NaN; %delete
            end
        end
    elseif f==length(HS) %if last HS
        if mean(and(TO<HS(f),TO>HS(f-1)))==0 %and no previous TO
            if isnan(HS(f-1)) %in the case that the previous HS was deleted
                if mean(and(TO<HS(f),TO>HS(f-2)))==0     
                HS(f) = NaN;
                end
            else
            HS(f)=NaN;
            end
        end
    elseif and(mean(and(TO<HS(f),TO>HS(f-1)))==0,mean(and(TO>HS(f),TO<HS(f+1)))==0) %if no TO between current HS and next AND previous HS
        if and(f==2,isnan(HS(f-1))) %if this is the second HS and previous HS was deleted
           if mean(TO<HS(f))~=0     
            HS(f) = NaN;
           end 
        elseif and(f>2,isnan(HS(f-1))) %in the case that the previous HS was deleted
            if mean(and(TO<HS(f),TO>HS(f-2)))==0     
            HS(f) = NaN;
            end
        else
            HS(f)=NaN;
        end
    elseif mean(and(TO<HS(f),TO>HS(f-1)))==0 %if only no TO between current HS and previous HS
        if and(gesHS(f)>gesHS(f+1),gesHS(f)>gesHS(f-1))  %and on a flat
            HS(f) = NaN; %delete current HS
        else %and not on a flat
            HS(f-1) = NaN; %delete previous HS
        end

    end
end

HS=(HS(~isnan(HS)));

t=1;
while and(t<length(TO),TO(t)<HS(end)) 
   if mean(TO(t)<HS)==1 %get rid of leading TOs
       TO(t) = [];
       t = t;
   elseif t<length(TO)
       if and(TO(t)<HS(t+1),TO(t+1)<HS(t+1)) %more than 1 TO in a row
       TO(t) = [];
       t=t;
       else
           t=t+1;
       end
   else
       t=t+1;
   end
end

for g=t:length(TO)
   if mean(TO(g)>HS)==1 %get rid of trailing TOs
           TO(g) = NaN;
   end
end

TO=(TO(~isnan(TO)));

%Plot candidate gait events on vertical foot displacement time series and
%visualize wavelet peaks used for gait events
figure()
ax = subplot(211);
plot(ftd(:,3))
hold on
    scatter(HS,ftd(HS,3),'k','*')
    scatter(TO,ftd(TO,3),'r','x')
title('vertical foot trajectory')
ylabel(ax,'vertical displacement (m)')
xlabel('frame #')

bx = subplot(212);
findpeaks(sig,'MinPeakProminence',threshold);
title('Wavelet gait event identification')
ylabel('wavelet power')
xlabel('frame #')
linkaxes([ax bx], 'x');
end

end
