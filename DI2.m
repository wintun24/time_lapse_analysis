%this code to use if you want to get shape factor of first and last frames.

imagefiles = dir('*.jpg'); 
filename='xy_02';%%%%%%%%%%%change file name to save%%%%%%%%%
nfiles = length(imagefiles);    % Number of files found
fontSize = 16;
for j=1:35% how many number of cell you want to track in one video
for i=1:nfiles%number of frames
   currentfilename = imagefiles(i).name;
   Image = imread(currentfilename);

   imshow(Image, []);
   axis on;
%
%set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
%message = sprintf('Left click and hold to begin drawing.\nSimply lift the mouse button to finish');
%uiwait(msgbox(message));


   if i==1 || i==nfiles %only like to calculate shape factor of 1st and last frame
       grayImage=rgb2gray(Image);
       imshow(grayImage,[]);
      title('Original Grayscale Image', 'FontSize', fontSize);
      hFH = imfreehand(); %outline interested cell
      
    % Create a binary image ("mask") from the ROI object.
      binaryImage = hFH.createMask();
      xy= hFH.getPosition;

     labeledImage = bwlabel(binaryImage);
     measurements = regionprops(binaryImage, grayImage, ...
    'area', 'Centroid', 'WeightedCentroid', 'Perimeter','MajorAxisLength','MinorAxisLength');
    area(i,j) = measurements.Area;
    centroid = measurements.Centroid;
   centerOfMass = measurements.WeightedCentroid;
   perimeter(i,j) = measurements.Perimeter;
   Dmax(i,j)=measurements.MajorAxisLength;
   Dmin(i,j)=measurements.MinorAxisLength;

% Calculate the area, in pixels, that they drew.
   numberOfPixels1(i,j) = sum(binaryImage(:));
% Another way to calculate it that takes fractional pixels into account.
  numberOfPixels2(i,j) = bwarea(binaryImage);
   
else 
  w = input('type number');
  if w==1 %type one if the cell disappear before end of frames otherwise press enter
    break;
          
  end
  end
end
end

%AREA=area(1,:;last,:);
area1=area(1,:);%areas of cell in frame 1
area2=area(end,:); %areas of cells in last frame
final_area=[area1;area2];

Dmax1=Dmax(1,:);%major axis of first frame
Dmax2=Dmax(end,:);%major axis of last frame
final_Dmax=[Dmax1;Dmax2];
Dmin1=Dmin(1,:);%minor aixis of first frame
Dmin2=Dmin(end,:);%minor axis of last frame
final_Dmin=[Dmin1;Dmin2];

perimeter1=perimeter(1,:);%perimeter of first frame
perimeter2=perimeter(end,:);%perimeter of last frame
final_perimeter=[perimeter1;perimeter2];

%removing zero columns

%SF1=(4*pi*A)/p^2
zero1=any(final_area==0);
final_final_area=final_area(:,~zero1);
zero2=any(final_Dmax==0);
final_final_Dmax=final_Dmax(:,~zero2);
zero3=any(final_Dmin==0);
final_final_Dmin=final_Dmin(:,~zero3);
zero4=any(final_perimeter==0);
final_final_perimeter=final_perimeter(:,~zero4);

pi=3.14159;
%shape factor for first frame
for j=1:size(final_final_area,2)
    for i=1
        SF1(i,j)=(4*pi*final_final_area(i,j))/final_final_perimeter(i,j)^2;
    end
end
%shape facotr for last frames
for j=1:size(final_final_area,2)
    for i=1       

   SF2(i,j)=(4*pi*final_final_area(i+1,j))/final_final_perimeter(i+1,j)^2;
    end
end

%SI=dmin/dmax
%shape index for first frame
for j=1:size(final_final_area,2);
    for i=1;
        SI1(i,j)=final_final_Dmin(i,j)/final_final_Dmax(i,j);
    end
end
%shape index for last frames
for j=1:size (final_final_area,2)
    for i=1
        SI2(i,j)=final_final_Dmin(i+1,j)/final_final_Dmax(i+1,j);
    end
end

ave_SF1=mean(SF1);
ave_SF2=mean(SF2);
ave_SI1=mean(SI1);
ave_SI2=mean(SI2);

save(filename);














