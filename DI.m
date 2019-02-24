%this code to use if you want to know shape index of every frames
imagefiles = dir('*.jpg');      
nfiles = length(imagefiles);    % Number of files found
fontSize = 16;
for j=1:2
for i=1:nfiles
   currentfilename = imagefiles(i).name;
   %currentimage = imread(currentfilename);
   %images{i} = currentimage;
   %figure
   %A= imshow(images{i})
   grayImage = imread(currentfilename);

imshow(grayImage, []);
axis on;
title('Original Grayscale Image', 'FontSize', fontSize);
%set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
%message = sprintf('Left click and hold to begin drawing.\nSimply lift the mouse button to finish');
%uiwait(msgbox(message));


   if i==1 
      hFH = imfreehand();

% Create a binary image ("mask") from the ROI object.
      binaryImage = hFH.createMask();
      xy= hFH.getPosition;

     labeledImage = bwlabel(binaryImage);
     measurements = regionprops(binaryImage, grayImage, ...
    'area', 'Centroid', 'WeightedCentroid', 'Perimeter','MajorAxisLength','MinorAxisLength');
    area(i,j) = measurements.Area
    centroid = measurements.Centroid;
   centerOfMass = measurements.WeightedCentroid;
   perimeter(i,j) = measurements.Perimeter
   D1(i,j)=measurements.MajorAxisLength;
   D2(i,j)=measurements.MinorAxisLength;

% Calculate the area, in pixels, that they drew.
   numberOfPixels1(i,j) = sum(binaryImage(:));
% Another way to calculate it that takes fractional pixels into account.
  numberOfPixels2(i,j) = bwarea(binaryImage);
   
else 
  w = input('type number');
  if w==1
    break;
  else    %1 to continue, other to stop
  
        
        hFH = imfreehand();

% Create a binary image ("mask") from the ROI object.
binaryImage = hFH.createMask();
xy= hFH.getPosition;

% Now make it smaller so we can show more images.
%subplot(2, 2, 1);
%imshow(grayImage, []);
%axis on;
%title('Original Grayscale Image', 'FontSize', fontSize);

% Display the freehand mask.
%subplot(2, 2, 2);
%imshow(binaryImage);
%axis on;
%title('Binary mask of the region', 'FontSize', fontSize);

% Label the binary image and computer the centroid and center of mass.
labeledImage = bwlabel(binaryImage);
measurements = regionprops(binaryImage, grayImage, ...
    'area', 'Centroid', 'WeightedCentroid', 'Perimeter','MajorAxisLength','MinorAxisLength');
area(i,j) = measurements.Area;
centroid = measurements.Centroid;
centerOfMass = measurements.WeightedCentroid;
perimeter(i,j) = measurements.Perimeter;
D1(i,j)=measurements.MajorAxisLength;
D2(i,j)=measurements.MinorAxisLength;

% Calculate the area, in pixels, that they drew.
numberOfPixels1(i,j) = sum(binaryImage(:));
% Another way to calculate it that takes fractional pixels into account.
numberOfPixels2(i,j) = bwarea(binaryImage);
    
  end
  end
end
end
