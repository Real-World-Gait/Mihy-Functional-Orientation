# Mihy-Functional-Orientation

This repository includes the functions used for functional sensor-to-segment orientation, walking identification, and gait event identification as performed in Mihy et al., "A functional sensor-to-segment calibration method reduces the effects of varied sensor placement on estimates of segment angular excursion", as well as example data files from a single subject. Functions for performing zero-velocity update procedures as originally shared by Stephen Cain at a 2017 American Society of Biomechanics workshop are also included as the gait event function depends on these functions. Code files include details on the function inputs, outputs, and dependencies. Some instructions are also included below.

## Code File Details:

**1. Example_function_implementation.m**  
   This script includes the syntax for each example function and its matching data file.

**2. func_S2S_orientation.m**  
   Use this function with the Cal.mat data file. The function creates a rotation matrix that defines the orientation of a body-fixed reference frame (e.g., functional anatomical axes) relative to the IMU-fixed reference frame. To perform sensor-to-segment orientation using the rotation matrix, pre-multiply IMU-frame data by the R_body matrix to get body-frame data.  
    The approach used here is documented in Supplement 1 of <ins>Cain SM, et al. Gait Posture 43, 65-69, 2016.</ins>

**3. gait_ID_fft.m**  
Use this function with the GaitTrials.mat data file. The function identifies windows of probable walking based on frequency analysis. Note that further processing is necessary to determine whether identified windows are level, in a straight line, of desired length, etc.

**4. gait_event_cwt.m**  
Use this function with the BoutExample.mat data file. The function identifies frames of toe off and heel strike events using a continuous wavelet transform technique. This function requires find_zero_velocity.m and zupt_displacement.m

**5. find_zero_velocity.m**  
This function identifies frames of zero velocity. This is generally used to determine when the foot is static for accurate spatiotemporal variable calculations.

**6. zupt_displacement.m**  
This function provides drift-corrected linear velocity and displacement based on the zero frames identified in find_zero_velocity.m

## Data File Details:

**1. CalibrationData.mat**  
A structure containing data for functional orientation procedures. rotiP contains the indices of functonal rotation for the pelvis sensors; rotiR contains the indices of functional rotation for the thigh, shank, and foot sensors. Sensor fields contain R_body matrices that were created from the included gravi and roti indices and the Cal sensor data using these algorithms.

**2. GaitTrials.mat**  
A structure containing data collected during a standard in-lab gait analysis. Each sensor field includes raw, IMU-frame acceleration (a), angular velocity (w), and orientation data in quaternion format (q), as well as rotation matrices for transforming acceleration and angular velocity data to X and Z functional reference frames.

**3. BoutExample.mat**  
A structure containing a single bout of walking extracted using gait_ID_fft.m for use identifying gait events and performing other example procedures. gait_events field includes the frames of gait events within strides (column 1 = heel strike 1, column 2 = toe-off, column 3 = heel strike 2). Spatiotemporal variables were derived from each stride using zero velocity procedures (see additional functions). Each sensor field includes raw, IMU-frame acceleration (a), angular velocity (w), orientation data in quaternion format (q), rotation matrices for transforming acceleration and angular velocity data to X and Z functional reference frames, acceleration and angular velocity data transformed to X and Z functional reference frames (_bodyX and _bodyZ variables), quaternion data transformed to rotation matrix format (or variable), and acceleration and angular velocity data transformed to world reference frame using the or variable (_world variables). 
