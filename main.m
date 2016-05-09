%%  Kinect egomotion estimation
%   Author:  	    Leon de Lange 
%                   TuDelft, BMD Master Thesis
%                   Multi-view object retrieval
%   Last revision:  15 march 2016
%                   UNREGISTERED DEPTH IS USED!

clc
clear all
close all


%% constants
res = [640 480];            % [pixels] kinect resolution
bound = [0.75 4];           % [m] reliable sensor depth boundaries
mtresh = 0.05;              % [m] threshold matching distance
itresh = 1e-3;              % iteration threshold
iter = 1;                   % iteration count
optimal = true;             % optimal alignment boolean
frame = 1;                  % number of frames processed


%% variables
mot = eye(4);           % kinect (estimated) motion matrix
pos = [0 0 0]';         % kinect position


%% kinect camera intrinsic parameters
f = [572.88     572.74];            % [pixels] focal point
c = [240.16     319.65];            % [pixels] principle point
k = [-0.00307   3.33095];           % depth calibration


%% Create listener for kinect camera
fprintf('Making first call to initalize the kinect driver and listening thread\n');
kinect_mex();
pause(5)

fprintf('Making second call starts getting data\n');
kinect_mex();
pause(5)

figure(1)
hold on


%% create video player
videoPlayer  = vision.VideoPlayer('Position',[100 100 [res(1), res(2)]]);


%% Create a point tracker
pointTracker = vision.PointTracker('NumPyramidLevels',7,'MaxBidirectionalError', 8, 'MaxIterations',50,'BlockSize',[5 5]);

points.Count = 0;
while(points.Count < 1)
    
    % grab frame
    [t1_d, t1_rgb] = kinect_mex();

    % reshape
    t1_rgb = permute(reshape(t1_rgb,[3 res]),[3 2 1]);

    % extract SURF corners
    points = detectSURFFeatures(rgb2gray(t1_rgb));
    
end
    
% initialize tracker
initialize(pointTracker,points.Location,t1_rgb);


%% ICP algorithm - forever loop
 while (1 == 1)
    
    % clear screen
    clc
    
    % start stopwatch for frame rate
    tic
    
    % grab next frame
    [t1_d, t1_rgb] = kinect_mex();
          
    % reshape datastream using kinect resolution
    t1_d = permute(reshape(t1_d,res),[2 1]);
    t1_d = im2double(t1_d, 'indexed');
    t1_rgb = permute(reshape(t1_rgb,[3 res]),[3 2 1]);
        
    % calculate inverse depth coordinates u, v and q
    t1_v = (repmat((1:res(2))',[1,res(1)]) - c(1))./f(1);
    t1_u = (repmat((1:res(1)),[res(2),1]) - c(2))./f(2);
    t1_q = k(1) .* t1_d + k(2);
    
    % obtain euclidean coordinates x, y and z
    t1_z = 1 ./ t1_q;
    t1_x = t1_z .* t1_u;
    t1_y = t1_z .* t1_v;
     
    % clip depth data outside reliable range
    t1_z = t1_z.*(t1_z>bound(1) & t1_z<bound(2));
        
    % extract SURF corners
    points = detectSURFFeatures(rgb2gray(t1_rgb));
    
    % add points to tracker
    setPoints(pointTracker, points.Location);
    
    % update tracking step
    [isPoints, isFound] = step(pointTracker, t1_rgb);
    
    % matched points in both frames
    oldPoints = points.Location(isFound, :);
    newPoints = isPoints(isFound, :);
                 
    % show video output
    t1_rgb = insertMarker(t1_rgb, floor(newPoints), '+','Color', 'red');
    t1_rgb = insertMarker(t1_rgb, floor(oldPoints), '+','Color', 'white');
    step(videoPlayer, t1_rgb);
           
    % when frame t0 exists
    if(exist('t0_x','var')&&exist('t0_y','var')&&exist('t0_z','var'))
              
        % get indices matched points
        p0_ind = sub2ind([res(2) res(1)],floor(oldPoints(:,2)),floor(oldPoints(:,1)));
        p1_ind = sub2ind([res(2) res(1)],floor(newPoints(:,2)),floor(newPoints(:,1)));
        
        % obtain euclidean coordinates
        p0_xyz = [t0_x(p0_ind) t0_y(p0_ind) t0_z(p0_ind)];
        p1_xyz = [t1_x(p1_ind) t1_y(p1_ind) t1_z(p1_ind)];
        
        % check if a depth is assigned for each point pair
        mask = (p0_xyz(:,3) > 0)&(p1_xyz(:,3) > 0);
        p0_xyz = p0_xyz(mask,:);
        p1_xyz = p1_xyz(mask,:);
        
        % check if distance of point pair is below treshold
        mask2 = sqrt(sum((p0_xyz - p1_xyz).^2,2)) < mtresh;
        p0 = p0_xyz(mask2,:);
        p1 = p1_xyz(mask2,:);
        
        % add a colom of ones to pointlists
        p0(:,4) = 1;
        p1(:,4) = 1;
              
        % when there are enough points to solve the equations
        if((size(p0,1) > 6)&&(size(p1,1) > 6))

            % initialize total point distance for iteration criterea
            dist = inf;

            % iterative update pose estimation
            while(dist > itresh)

                 % initialize stacked Jacobian
                J = [];

                % initialize stacked point differences
                Y = [];

                % for each pointpair
                for i = 1:size(p0,1)

                    % construct Jacobian
                    Jp = [  1   0   0   0           p0(i,3)     -p0(i,2);
                            0   1   0   -p0(i,3)    0           p0(i,1);
                            0   0   1   p0(i,2)     -p0(i,1)    0           ];

                    % stack Jacobians
                    J = [J; Jp];

                    % calculate point difference
                    Yp = [  p1(i,1) - p0(i,1);
                            p1(i,2) - p0(i,2);
                            p1(i,3) - p0(i,3)   ];

                    % stack point differences
                    Y = [Y; Yp]; 

                end

                % error of each point to its objective
                err = sqrt(sum((p0-p1).^2,2));

                % standard deviation of error
                std = sqrt(mean(sum((err - mean(err)).^2)));

                % obtain outlier weighting matrix using M-estimator
                W = diag(repelem((1-(err.^2)./(std.^2 + err.^2)),3));

                % estimate motion parameters
                B = (J'*W*J)\J'*W*Y;

                % alternating for each iteration (optimal rotation then
                % translation)
                if(optimal)
                    B(1:3) = 0;
                    optimal = false;
                else
                    B(4:6) = 0;
                    optimal = true;
                end

                % estimated motion matrix
                M = [   1           -B(6)       B(5)        B(1);
                        B(6)        1           -B(4)       B(2);
                        -B(5)       B(4)        1           B(3);
                        0           0           0           1       ];

                % calculate new estimated point list
                pe = (M * p0')';

                % calculate distance point list update (termination criterea)
                dist = sum(sqrt(sum((p0 - pe).^2,2)));

                % set estimate point list as new point list
                p0 = pe;

                % update estimated motion parameters
                mot = M * mot;

                % increase iteration (termination criterea)
                iter = iter + 1;
                if(iter > 15)
                    break;
                end

            end

            % body to world orientation
            wmot = inv(mot);
        
            % update position estimation
            pos = [pos wmot(1:3,4)];

        else

            display('too few points')

        end
        
    end
                
    % set current frame to previous frame
    t0_x = t1_x;
    t0_y = t1_y;
    t0_z = t1_z;
                    
    % reset paramters
    iter = 1;                   % iteration
    optimal = true;             % reset optimal alignment boolean
    
    % stopwatch (fps)
    fps = 1/toc;
    display(fps)
             
    % draw trajectory
    if(exist('pos','var')&&exist('wmot','var'))
        clf;
        plot3(pos(1,:),pos(2,:),pos(3,:),'.')
        xlabel('x')
        ylabel('y')
        zlabel('z')
        hold on
        plotCamera('Location',wmot(1:3,4),'Orientation',mot(1:3,1:3),'Opacity',0,'Size',0.05);
        hold off
        view(0,0);
        axis('equal');           
        drawnow
    end
             
    frame = frame + 1;
        
end
    
%% close listener
kinect_mex('q');


%% release videoplayer
release(videoPlayer);