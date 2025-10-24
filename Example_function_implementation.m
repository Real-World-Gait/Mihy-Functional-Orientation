%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Example implementation code provided as a supplement to
% JA Mihy, M Wagatsuma, SM Cain, JF Hafer, A functional sensor-to-segment 
% calibration method reduces the effects of varied sensor placement on 
% estimates of segment angular excursion, J Appl Biomech
% 
% See notes within functions and example_data.mat for details
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Example func_S2S_orientation implementation
%load data
load('Cal.mat')

align = 1;
a = Cal.Pelvis1.a;
w = Cal.Pelvis1.w;
gravi = Cal.gravi;
roti = Cal.rotiP;

Pelvis1.R_bodyX = func_S2S_orientation(align,a,w,gravi,roti);

align = 2;
a = Cal.RShank1.a;
w = Cal.RShank1.w;
gravi = Cal.gravi;
roti = Cal.rotiR;

Shank1.R_bodyZ = func_S2S_orientation(align,a,w,gravi,roti);

%% Example gait_ID_fft implementation
%load data
load('GaitTrials.mat')

%Note that the GaitTrials data includes only a single sensor per segment
signal = GaitTrials.RShank.w;
Fs = 128;
windowsize = 5;

[strt,nd,fig] = gait_ID_fft(signal,Fs,windowsize);

%% Example gait_event_cwt implementation
%load data
load('BoutExample.mat')

signal_a = BoutExample.RFoot.a_world;
signal_w = BoutExample.RFoot.w_world;
filter = 1;
Fs = 128;
LPfreq = 4;

[TO,HS,ftd] = gait_event_cwt(signal_a,signal_w,filter,Fs,LPfreq);
