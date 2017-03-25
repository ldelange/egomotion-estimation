# Kinect egomotion estimation for MATLAB

Visual odometry method, which estimates the position of a Kinect camera based on sparse features. The external kinect_mex library is used to grab data.

+ Extracting SURF features in two frames
+ 6DoF linearized motion model estimates camera position
+ Unregisterd! depth is used

Algorithm is based on the paper:
Lui, W. L. D., Tang, T. J. J., Drummond, T., and Li, W. H. (2012). Robust egomotion estimation using icp
in inverse depth coordinates. In Robotics and Automation (ICRA), 2012 IEEE International Conference on,
pages 1671–1678. IEEE.

## Results

![Egomotion estimation](https://raw.githubusercontent.com/ldelange/egomotion-estimation/master/egomotion.png)

Combination of egomotion and plane detection.

![egmotion and plane detection](https://raw.githubusercontent.com/ldelange/egomotion-estimation/master/egmotion and plane detection.png
)
