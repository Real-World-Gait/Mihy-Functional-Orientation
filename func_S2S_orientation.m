function R_body_IMU = func_S2S_orientation(align,a,w,gravi,roti)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Functional sensor-to-segment orientation code provided as a supplement to
% JA Mihy, M Wagatsuma, SM Cain, JF Hafer, A functional sensor-to-segment 
% calibration method reduces the effects of varied sensor placement on 
% estimates of segment angular excursion, J Appl Biomech
%
% See notes below "Output" if using this function with the example data
% provided with Mihy et al.
%
% This function is used to create a rotation matrix that defines the
% orientation of a body-fixed reference frame (e.g., functional anatomical
% axes) relative to the IMU-fixed reference frame. To perform
% sensor-to-segment orientation using the rotation matrix, pre-multiply
% IMU-frame data by the R_body matrix to get body-frame data.
%
% The approach used here is documented in Supplement 1 of:
% Cain SM, et al. Gait Posture 43, 65-69, 2016.
%
% X, Y, Z body-fixed axes estimate sagittal (medial-lateral), frontal
% (anterior-posterior), and longitudinal (cranial-caudal) anatomical axes.
% A period of static posture is used to define the body-fixed Z axis,  a
% period of approximately planar rotation is used to define the body-fixed X
% axis, and the body-fixed Y axis is defined with the cross-product of Z
% and X.
% 
% There are two options for creating an orthogonal reference frame after
% the above operations:
% X alignment:
% Hold original alignment of X axis, correct Z to X & Y as a final step
% Z alignment:
% Hold original alignment of Z axis, correct X to Y & Z as a final step
%
% Input: calibration data
% align: X or Z alignment (x=1, z=2)
% a: raw (IMU-frame) sensor acceleration
% w: raw (IMU-frame) sensor angular velocity
% t: time vector for data
% gravi: indices of static trial within data
% roti: indices of functional rotation within data
% 
% Output:
% R_body_IMU: Rotation matrix defining the rotation required to
% go from the IMU reference frame to the body-fixed reference frame. 
%
% If using with Mihy et al. example data, Cal is a structure containing
% data for functional orientation procedures. rotiP contains the indices of
% functonal rotation for the pelvis sensors; rotiR contains the indices of
% functional rotation for the thigh, shank, and foot sensors. Sensor fields
% contain R_body matrices that were created from the included gravi and
% roti indices and the Cal sensor data using these algorithms.
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define orientation of body segment relative to the IMU reference frame

% Using X alignment 
if align == 1
    
% Use gravity vector to define segment Z axis
ma = mean(a(gravi,:));
Z = ma./norm(ma);      

% Use rotation to define segment X axis
PCAcoeff = pca(w(roti,:));
X_pca = (PCAcoeff(:,1))';
X_pca = X_pca./norm(X_pca);

% Calculate Y according to Y = Z x X
Y_pca = cross(Z,X_pca);
Y_pca = Y_pca./norm(Y_pca);

X = X_pca;

% Re-calculate Z to ensure orthogonality (Z = X x Y)
Z = cross(X,Y_pca);
Z = Z./norm(Z);

%Define initial direction cosine matrix
R_body_IMU = [X; Y_pca; Z];

else %Using Z alignment

% Use gravity vector to define segment Z axis
ma = mean(a(gravi,:));
Z = ma./norm(ma);      

% Use rotation to define segment X axis
PCAcoeff = pca(w(roti,:));
X_pca = (PCAcoeff(:,1))';
X_pca = X_pca./norm(X_pca);

% Calculate Y according to Y = Z x X
Y_pca = cross(Z,X_pca);
Y_pca = Y_pca./norm(Y_pca);

% Re-calculate X to ensure orthogonality (X = Y x Z)
X = cross(Y_pca,Z);
X = X./norm(X);

%Define initial direction cosine matrix
R_body_IMU = [X; Y_pca; Z];  
    
end
