clear all
close all
%%
%calculating shape index THIS IS SAME AS DI and DI_2 but much easier
load('xy_01.mat')
filename='xy_01_SF.mat';
for i=1:size(index,1)
    frame1X(i)=XX{index(i)}(1);
    frame65X(i)=XX{index(i)}(65);
    frame1Y(i)=YY{index(i)}(1);
    frame65Y(i)=YY{index(i)}(65);
end
    
   
     I=imread('first.png');
     position= [frame1X' frame1Y'];
     value=1:1:size(position,1);
     RGB = insertText(I,position,value,'FontSize',18);
     figure()
     imshow(RGB)
     grayImage1=rgb2gray(RGB);
     figure()
     imshow(grayImage1,[]);
   
     for i=1:size(position,1)
     hFH1 = imfreehand();
     
     
     binaryImage1 = hFH1.createMask();
     % xy= hFH.getPosition;

     %labeledImage1 = bwlabel(binaryImage1);
     measurements1 = regionprops(binaryImage1, grayImage1, ...
    'area', 'Perimeter','MajorAxisLength','MinorAxisLength');
     
   area1(i) = measurements1.Area;
   perimeter1(i) = measurements1.Perimeter;
   Dmax1(i)=measurements1.MajorAxisLength;
   Dmin1(i)=measurements1.MinorAxisLength;
   
     end
  
   
   
     I=imread('last.png');
     position=    [frame65X' frame65Y'];
     %value=1:1:size(position,1);
     RGB = insertText(I,position,value,'FontSize',18);
     figure()
     imshow(RGB);
     grayImage2=rgb2gray(RGB);
     figure()
     imshow(grayImage2,[]);
     
     for i=1:size(position,1)
     hFH2 = imfreehand(); 
     
     
     binaryImage2 = hFH2.createMask();
    %  xy= hFH.getPosition;

     %labeledImage2 = bwlabel(binaryImage2);
     measurements2 = regionprops(binaryImage2, grayImage2, ...
    'area', 'Perimeter','MajorAxisLength','MinorAxisLength');
     
   area2(i) = measurements2.Area;
   perimeter2(i) = measurements2.Perimeter;
   Dmax2(i)=measurements2.MajorAxisLength;
   Dmin2(i)=measurements2.MinorAxisLength;
   
     end
     
    
     prompt = 'Do you want to delete any cell? ';
     for s=1:1000
     Remove(s)=input(prompt);
       if Remove(s)==0 %type zero when you finish
       break;
          
      end
     end
  
     Remove(:,end) = [];
     Remove=unique(Remove);
     area1(Remove)=[];
     area2(Remove)=[];
     perimeter1(Remove)=[];
     perimeter2(Remove)=[];
     Dmin1(Remove)=[];
     Dmax1(Remove)=[];
     Dmin2(Remove)=[];
     Dmax2(Remove)=[];
     pi=3.14159;
%shape factor for first frame
for j=1:size(area1,2)
    
        SF1(j)=(4*pi*area1(j))/perimeter1(j)^2;
    
end
%shape facotr for last frames
for j=1:size(area2,2)
        

   SF2(j)=(4*pi*area2(j))/perimeter2(j)^2;
   
end

%SI=dmin/dmax
%shape index for first frame
for j=1:size(area1,2);
    
        SI1(j)=Dmin1(j)/Dmax1(j);
   
end
%shape index for last frames
for j=1:size (area2,2)
  
        SI2(j)=Dmin2(j)/Dmax2(j);
    
end
  

ave_SF1=mean(SF1);
ave_SF2=mean(SF2);
ave_SI1=mean(SI1);
ave_SI2=mean(SI2);

save(filename);
