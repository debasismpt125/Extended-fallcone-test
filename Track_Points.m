%% Track points on a object from video (Fall cone test)
% Read the MP4 file and convert into frames
% Debasis Mohapatra, Aalto University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fall cone tests on sensitive marine clay: a comprehensive experimental study and its replication with the Generalized Interpolation Material Point Method%
%Debasis Mohapatra, Saeideh Mohammadi, Maarit Saresma, Joonas J.
%Virtasalo,Wojciech T. Sołowski
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                      

clear
close all
framerate = input('Insert frames per second: ');   % e.g. 1000 frames/second
scale = input('Insert scale value meter/pixel: ') ; % e.g. scale =4.72E-05;   % meter/pixel  % 30/09/2022

videoReader = VideoReader('path_to_video\video.mp4');
totalFrames = videoReader.NumberOfFrames;
videoPlayer = vision.VideoPlayer('Position',[0,0,1000,2000]);
%% 
% Read the first video frame, which contains the object, define the region.
%%
objectFrame = readFrame(videoReader);
figure; 
imshow(objectFrame); 
objectRegion=round(getPosition(imrect))
% 
% Show initial frame with a red bounding box.
%%
objectImage = insertShape(objectFrame,'Rectangle',objectRegion,'Color','red'); 
%figure;
imshow(objectImage);
title('Red box shows object region');
%% 
% Automatically detect interest points in the object region.
points = detectMinEigenFeatures(im2gray(objectFrame),'ROI',objectRegion,'MinQuality',0.0000001);
%% 
% Display the detected points in green.
%%
pointImage = insertMarker(objectFrame,points.Location,'+','Color','green');
figure;
imshow(pointImage);
title('Detected interest points');
%% 
% Create a tracker object.
%%
tracker = vision.PointTracker('MaxBidirectionalError', 1);
%% 
% Initialize the tracker.
%%
initialize(tracker,points.Location,objectFrame);
%% 
% Read, track, display points, and results in each video frame.
points_n_1 = points.Location;
points_n = points.Location;
np = size(points_n,1);
k =1;
t(1) = 0;
dt = 0.001;
for i=1:np
frames(1).acc(i) = 0;
frames(1).position(i, 1:2)=points.Location(i,1:2);
frames(2).position(i, 1:2)=points.Location(i,1:2);
frames(1).disp(i) = 0;
frames(2).disp(i) = 0;
frames(1).vel(i)=0;

end

while k<totalFrames 
      
      frames(k).position(1:np, 1:2)=points_n_1;
      frames(k+1).position(1:np, 1:2)=points_n;
      frame = readFrame(videoReader);
      [points,validity] = tracker(frame);
      points_n1 = points;
      frames(k+2).position(1:np, 1:2)=points_n1;
      t(k+1) = t(k)+dt;
      for i=1:np
          frames(k+2).dispx(i)=((points_n1(i,1)-points_n(i,1)));
          frames(k+2).dispy(i)=((points_n1(i,2)-points_n(i,2)));
          frames(k+2).disp(i) =  frames(k+1).disp(i)+frames(k+2).dispy(i)*scale;
          frames(k+1).vel(i) = ((points_n1(i,2)-points_n_1(i,2))/2) * framerate*scale;
          vel2 = ((points_n1(i,2)-points_n(i,2))/2) * framerate*scale;
          vel1 = ((points_n(i,2)-points_n_1(i,2))/2) * framerate*scale;
          frames(k+1).acc(i) =  (vel2 - vel1)/0.001;            
      end

      out = insertMarker(frame,points(validity, :),'+');
      videoPlayer(out);
      x = ['frame number ', num2str(k)];
      disp(x) 
      points_n_1=points_n;
      points_n=points_n1;
      k=k+1;
      for i = 1:np
      startpoint=[points_n(i,1) points_n(i,2)];endpoint=[points_n_1(i,1) points_n_1(i,2)];headsize=1.1;
      v1 = headsize*(startpoint-endpoint)/4;
      theta = 15*pi/180; 
      theta1 = -1*15*pi/180; 
      rotMatrix = [cos(theta) -sin(theta) ; sin(theta) cos(theta)]; 
      rotMatrix1 = [cos(theta1) -sin(theta1) ; sin(theta1) cos(theta1)]; 
      v2 = v1*rotMatrix; 
      v3 = v1*rotMatrix1; 
      x1 = endpoint; 
      x2 = x1 + v2; 
      x3 = x1 + v3; 
      hold on; 
      % below line fills the arrowhead (black) 
      %fill([x1(1) x2(1) x3(1)],[x1(2) x2(2) x3(2)],[0 0 1] ); 
      % below line draws line (black) 
      plot([startpoint(1) endpoint(1)],[startpoint(2) endpoint(2)],'linewidth',0.5,'color',[0 0 1]);
      hold off
      end
      
end
    frames(1).velsum = 0.;
    frames(1).dispsum = 0;
    frames(1).accsum = 0;
for i=1:np
    for k=1:totalFrames
    if (frames(k).disp(i)<1e-4)
        frames(k).disp(i) = 0;
    end
    end

end

a =1;
for k=1:totalFrames
    frames1(k).dispsum = 0;
    frames1(k).velsum = 0;
    frames1(k).accsum = 0;
end
for i=1 :np
    if(frames(totalFrames).disp(i)>1e-3)
       for k=1:totalFrames
          frames1(k).disp(a)=frames(k).disp(i);
          frames1(k).vel(a)=frames(k).vel(i);
          frames1(k).acc(a)=frames(k).acc(i);
          frames1(k).dispsum = frames1(k).dispsum+frames1(k).disp(a);
          frames1(k).velsum = frames1(k).velsum+frames1(k).vel(a);
          frames1(k).accsum = frames1(k).accsum+frames1(k).acc(a);
       end
       a=a+1;
    end
end
for k=1:totalFrames
    frames1(k).dispavg = frames1(k).dispsum/(a-1) ;
    frames1(k).velavg = frames1(k).velsum/(a-1);
    frames1(k).accavg = frames1(k).accsum/(a-1);
end




%%
release(videoPlayer);
%% Plotting the Displacement, Velocity and acceleration of cone
for i=1:totalFrames-1
    
    vel1(i)=frames1(i).vel(1);
    vel2(i)=frames1(i).vel(2);
    vel3(i)=frames1(i).vel(3);

    velavg(i)= frames1(i).velavg;

    disp1(i)= frames1(i).disp(1);
    disp2(i)= frames1(i).disp(2);
    dispavg(i) = frames1(i).dispavg;

    acc1(i)= frames1(i).acc(1);
    acc2(i)= frames1(i).acc(2);
    accavg(i)= frames1(i).accavg;

end

figure
plot(t(1:totalFrames-1), disp1,t(1:totalFrames-1), disp2, '-r','LineWidth',0.3)
hold on
plot(t(1:totalFrames-1), dispavg, '-k','LineWidth',1);
xlabel('Time (s)', 'FontSize',12,'FontWeight','bold','Color','b')
ylabel('Displacement (m)', 'FontSize',12,'FontWeight','bold','Color','b')
legend('point 1','point 2', 'Average displacement')

figure
plot(t(1:totalFrames-1), vel1,t(1:totalFrames-1), vel2,'-r','LineWidth',0.5);
hold on
plot(t(1:totalFrames-1), velavg, '-k','LineWidth',1.5);
xlabel('Time (s)', 'FontSize',12,'FontWeight','bold','Color','b')
ylabel('Velocity (m/s)', 'FontSize',12,'FontWeight','bold','Color','b')
legend('point 1','point 2', 'Average velocity')
hold off


figure
plot(t(1:totalFrames-1), acc1,t(1:totalFrames-1), acc2, '-r','LineWidth',0.3)
hold on
plot(t(1:totalFrames-1), accavg, '-k','LineWidth',1);
xlabel('Time (s)', 'FontSize',12,'FontWeight','bold','Color','b')
ylabel('Acceleration (m/s^2)', 'FontSize',12,'FontWeight','bold','Color','b')
legend('point 1','point 2', 'Average acceleration')