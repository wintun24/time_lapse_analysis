
%this code is used to produce to the same data as Combined_RW but it is
%less complicated. 
clear all
close all

%%
%Parameters and reading original file %%%%%%%%%%%%%%%%%
filename='xy_22c.xls';%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%change file name%%%%%%%]
filename1='xy_22c_notT.mat'
A=xlsread(filename);%the whole matrix
%break down the matrix
cellno=A(:,2);%cell number it will repeate 1 1 1  2 2 2 ..etc
frame=A(:,3);% number of frame it is 1 2 3 4 5 till 65 or 97, may stop before that too. 
X=A(:,4);%X coordinate
Y=A(:,5);%Y coordinate
distance=A(:,6);%auto calculate distance
velocity=A(:,7);%auto calculate velocity
matrixsize=size(cellno,1);% total number of whole matrix
%%
%%%looking for end point of each cells by getting the diffrence col 1. 

for i=1:matrixsize-1
    DIFF=frame(i+1)-frame(i);
    if DIFF < 0
        endOfCell = frame(i);
    else endOfCell = 0;
  
    end
  RESULT(i)=endOfCell; 
end

%endCellPoint=[reshape(result,[],1)]; E=endCellPoint
E=RESULT(RESULT~=0);%location of end points
LAST=A(end,3);
E=[E,LAST];%the cutting points of individual cell
cellcount=size(E,2);%counting number of cells. 
%%

for i=1:matrixsize-1
    diff=cellno(i+1)-cellno(i);
    if diff > 0
        endCell = i;
    else endCell = 0;
   
    end
  result(i)=endCell;
end
%endCellPoint=[reshape(result,[],1)]; E=endCellPoint
%e is quite similar to E
e=result;
e=result(result~=0);%location of end points
last=A(end,1);
e=[e,last];
e=e';

%C=mat2cell(A,E);
%celldisp(C);
%fileID=('ans.txt');
%fid=fopen(fileID,'w');


 %for j=1:cellcount
    % for i=1:size(C{j})-1
    %   H(j, i)=C{j}(i)+C{j}(i+1);
    %end 
  %end

%fprintf(fid, '%10f %10f %10f %10f\n', H);



%%
%%calculating average velocity of all cells in one video based on velocity
%per frame
VELOCITY=velocity(find(velocity~=-1));
mean_velocity=mean(VELOCITY);

%%
%%calculating displacement of individual cells each clicking points
displacement1=sum(distance(2:e(1),1));
for i=1:cellcount-1
    displacement(i)=sum(distance((e(i)+2):(e(i+1)),1));
end
DP=[displacement1,displacement];
final_DP=DP';
%%
%%Calculating distance between inital point and final point of each cells%%

dist1=sqrt((X(e(1))-X(1))^2+(Y(e(1))-Y(1))^2) ; 
dist1=dist1*0.37;
for i=1:cellcount-1
dist(i)= sqrt(( X((e(i))+1)-X(e(i+1)))^2+ ( Y((e(i))+1)-Y(e(i+1)))^2); 
dist(i)=dist(i)*0.37;%0.37micrometer per pixel
end

d=[dist1,dist];
final_d=d';
%%
%calculatinng d/D
d1BYD1=dist1/displacement1;
dBYD=dist./displacement;

finaldBYD=[d1BYD1,dBYD];
final_dBYD=finaldBYD';
%%
%calculating angle of individual steps
x=mat2cell(X,E);
y=mat2cell(Y,E);

for j=1:cellcount
for i=1:size(x{j})-1
angle(j,i)=atan((y{j}(i+1)-y{j}(i))/(x{j}(i+1)-x{j}(i)));


pi=3.14159;

     if y{j}(i+1)>y{j}(i)&& x{j}(i+1)<x{j}(i);
         angleInRadian(j,i)=(pi)+angle(j,i);
         
     elseif y{j}(i+1)>y{j}(i)&& x{j}(i+1)>x{j}(i);
         angleInRadian(j,i)=angle(j,i);
         
     elseif y{j}(i+1)<y{j}(i)&& x{j}(i+1)>x{j}(i);
         angleInRadian(j,i)=(2*pi)+(angle(j,i));
         
     elseif y{j}(i+1)<y{j}(i)&& x{j}(i+1)<x{j}(i);
         angleInRadian(j,i)=(pi)+(angle(j,i));
         
     elseif y{j}(i+1)==y{j}(i)&& x{j}(i+1)<x{j}(i);
         angleInRadian(j,i)=pi;
         
     elseif y{j}(i+1)>y{j}(i)&& x{j}(i+1)==x{j}(i);
         angleInRadian(j,i)=pi/2;
         
     elseif y{j}(i+1)==y{j}(i)&& x{j}(i+1)>x{j}(i);
         angleInRadian(j,i)=2*pi;%%%%%%%%%%%%%
         
     elseif y{j}(i+1)<y{j}(i)&& x{j}(i+1)==x{j}(i);
         angleInRadian(j,i)=(3*pi)/2;
         
     elseif y{j}(i+1)==y{j}(i)&& x{j}(i+1)==x{j}(i);
          angleInRadian(j,i)=0;  
     end
    
end
end
AR=angleInRadian';


%%
%below is the same calculation for angle but it took the data as in one
%columne directly from the sheet



for i=1:size(X)-1
angle_col(i)=atan((Y(i+1)-Y(i))/(X(i+1)-X(i)));


pi=3.14159;

     if Y(i+1)>Y(i)&& X(i+1)<X(i);
         angleInRadian_col(i)=(pi)+angle_col(i);
         
     elseif Y(i+1)>Y(i)&& X(i+1)>X(i);
         angleInRadian_col(i)=angle_col(i);
         
     elseif Y(i+1)<Y(i)&& X(i+1)>X(i);
         angleInRadian_col(i)=(2*pi)+(angle_col(i));
         
     elseif Y(i+1)< Y(i)&& X(i+1)< X(i);
         angleInRadian_col(i)= 3.14159 + angle_col(i);
         
     elseif Y(i+1)==Y(i)&& X(i+1)<X(i);
         angleInRadian_col(i)=pi;
         
     elseif Y(i+1)>Y(i)&& X(i+1)==X(i);
         angleInRadian_col(i)=pi/2;
         
     elseif Y(i+1)==Y(i)&& X(i+1)>X(i);
         angleInRadian_col(i)=2*pi;%%%%%%%%%%%%%
         
     elseif Y(i+1)<Y(i)&& X(i+1)==X(i);
         angleInRadian_col(i)=(3*pi)/2;
         
     elseif Y(i+1)==Y(i)&& X(i+1)==X(i);
          angleInRadian_col(i)=0;  
     end
    
end

AR_col=angleInRadian_col';
count=0;

for i=1:size(e)-1
    AR_col(e(i)-count)=[];
    count=count+1;
end
%%
%%calculating deviation the following angle from previous angle

for j=1:size(AR,2)
    for i=1:size(AR,1)-1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
   Dev(i,j)=AR((i+1),j)-AR(i,j);
   
   if AR((i+1),j)==0 ||AR(i,j)==0%which means cell stop moving
       trueDev(i,j)=0;
       
   elseif  Dev(i,j)<-pi
      trueDev(i,j)=((2*pi)-AR(i,j)+AR(i+1,j));
   elseif Dev(i,j)>pi
       trueDev(i,j)=((2*pi)-AR(i+1,j)+AR(i,j)).*(-1);
       
        
   else trueDev(i,j)=Dev(i,j);
       

   end
   end
end

%%
%below deviatin calculation is same as above but it calcualte from col
%direct from the sheet

    for i=1:size(AR_col,1)-1
        
   Dev_col(i)=AR_col((i+1))-AR_col(i);
   
   if AR_col((i+1))==0 ||AR_col(i)==0%which means cell stop moving
       trueDev_col(i)=0;
       
   elseif  Dev_col(i)<-pi
      trueDev_col(i)=((2*pi)-AR_col(i)+AR_col(i+1));
   elseif Dev_col(i)>pi
       trueDev_col(i)=((2*pi)-AR_col(i+1)+AR_col(i)).*(-1);
       
        
   else trueDev_col(i)=Dev_col(i);
       

   end
    end
   trueDev_col=trueDev_col';
   

  
dev_cut=e-1;
Size=size(dev_cut,1);
counter=0;
for i=1:Size-1
   trueDev_col(dev_cut(i)-counter)=[];
    counter=counter+2;
end

%%

%%calculating deviation the all angles from the very frist angle
%for j=1:size(AR,2)
   
   %Dev_first(i,j)=AR((i+1),j)-AR(1,j);
   %if Dev_first(i,j)<-pi
      %trueDev_first(i,j)=((2*pi)-AR(1,j)+AR(i+1,j));
      
   %elseif Dev_first(i,j)>pi
       %trueDev_first(i,j)=((2*pi)-AR(i+1,j)+AR(1,j)).*(-1);
   
   %else trueDev_first(i,j)=Dev_first(i,j);
       

   %end
   %end
   %%
   %%%%%%%%%%%%%%%%%%%%%%new adding in version 2################
   %Correlate angle and velocity

AVpair=[AR_col,VELOCITY];
AV=sortrows(AVpair,1);



cut=0.7853975;
for i=1:8
    
    count=AV(find(AV(:,1)<=cut ));
    S(i)=size(count,1);
    
    cut=cut+0.7853975;
end


AVE1=mean(AV(1:S(1),2));
STD1=std(AV(1:S(1),2));
for i=1:7
     AVE(i)=mean(AV(S(i)+1:S(i+1),2));
     STD(i)= std(AV(S(i)+1:S(i+1),2));
end
averageVelocityEach=[AVE1,AVE];
STDeach=[STD1,STD];
VE_eachQuadrant=averageVelocityEach';%%%%%%%%%%%%TAKE THIS ###############
final_STD=STDeach';

save(filename1);
