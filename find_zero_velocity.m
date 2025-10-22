function [zupt_index] = find_zero_velocity(a,w,a_limit,w_limit)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Stephen Cain (smcain@umich.edu) 2017 
% This function was written for an IMU tutorial at the 2017 meeting of the
% American Society of Biomechanics, held in Boulder, CO.
%
% This function is used to find candidate times of zero-velocity using a
% threshold based method. When both acceleration and angular velocity
% magnitudes are less than define thresholds, it is assumed that the sensor
% is at zero velocity.
%
% Input:
% a: sensor-fixed acceleration (m/s^2)
% w: sensor-fixed angular rate (rad/s)
% a_limit: threshold for acceleration (m/s^2)
% w_limit: threshold for angular rate (rad/s)
%
% Output:
% zupt_index: index values of zero velocity
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grav = 9.80665; % m/s^2

% Calculate angular velocity magnitude and acceleration magnitude
wmag = (sum(w.^2,2)).^(1/2);

amag = (sum(a.^2,2)).^(1/2) - grav;

% Identify zero points that meet threshold conditions
if ~exist('a_limit','var')
    a_limit = 1;
else
end
if ~exist('w_limit','var')
    w_limit = 10;
else
end

zupt_index = find(wmag < w_limit & amag < a_limit & amag > -a_limit);
end