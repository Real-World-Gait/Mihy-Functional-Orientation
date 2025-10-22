function [v,d] = zupt_displacement(aworld,t,zupt_index)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Stephen Cain (smcain@umich.edu) 2017 
% This function was written for an IMU tutorial at the 2017 meeting of the
% American Society of Biomechanics, held in Boulder, CO.
%
% This function is used to calculate drift corrected velocity and position
% using a zero-velocity update approach. Between times of zero-velocity,
% linear drift is assumed. Also assumed is that zero-velocity points all
% occur at the same height, or Z position.
% 
% The basic method used in this function is documented in the following
% publications:
% Ojeda L and Borenstein J, J of Navigation 60, 391-407, 2007.
% Rebula JR, et al. Gait Posture 38, 974-980, 2013.
%
% Input:
% aworld: foot sensor acceleration in world-fixed frame with gravity
% removed (m/s^2)
% t: time (s)
% zupt_index: index values for zero-velocity points (from
% find_zero_velocity.m)
% 
% Output:
% v: velocity (m/s)
% d: displacement (m)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Force first and last sample to be zero points if they aren't already
if zupt_index(1) ~= 1
    zupt_index = [1; zupt_index];
else
end    
if zupt_index(end) ~= length(t)
    zupt_index = [zupt_index; length(t)];
else
end   

% Initialize variables
v = zeros(length(t),3);
d = zeros(length(t),3);

% Integrate between zero velocity points and apply drift corrections
for uu = 1:(length(zupt_index) - 1)
    acc = aworld(zupt_index(uu):zupt_index(uu+1),:);
    time = t(zupt_index(uu):zupt_index(uu+1));
    int_time = t(zupt_index(uu+1)) - t(zupt_index(uu));

    vx = cumtrapz(time,acc(:,1));
    vy = cumtrapz(time,acc(:,2));
    vz = cumtrapz(time,acc(:,3));

    % Assume linear drift velocities and correct
    vx_slope = (vx(end) - vx(1))/int_time;
    vx_drift = vx_slope*(time - time(1));
    vx = vx - vx_drift;

    vy_slope = (vy(end) - vy(1))/int_time;
    vy_drift = vy_slope*(time - time(1));
    vy = vy - vy_drift;

    vz_slope = (vz(end) - vz(1))/int_time;
    vz_drift = vz_slope*(time - time(1));
    vz = vz - vz_drift;
    
    % Store velocity values
    v(zupt_index(uu):zupt_index(uu+1),:) = [vx,vy,vz];
       
    % Integrate to obtain position
    dx = cumtrapz(time,vx);
    dy = cumtrapz(time,vy);
    dz = cumtrapz(time,vz);
        
    % Add initial position to integrated values
    if uu == 1       
    elseif uu > 1
        dx = dx + d(zupt_index(uu)-1,1);
        dy = dy + d(zupt_index(uu)-1,2);
        dz = dz + d(zupt_index(uu)-1,3);
    else
    end
    
    % Store position values
    d(zupt_index(uu):zupt_index(uu+1),:) = [dx,dy,dz];
        
end

end
