%this code will draw the trajectory of cells which stay the whole lenght of
%durination.
clear all
close all


filename='xy_02c.xls';%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%change file name%%%%%%%
A=xlsread(filename);%the whole matrix
%break down the matrix
cellno=A(:,2);%cell number it will repeate 1 1 1  2 2 2 ..etc
frame=A(:,3);% number of frame it is 1 2 3 4 5 till 65 or 97, may stop before that too. 
X=A(:,4);%X coordinate
Y=A(:,5);%Y coordinate
distance=A(:,6);%auto calculate distance
velocity=A(:,7);%auto calculate velocity
matrixsize=size(cellno,1);% total number of whole matrix

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
x=mat2cell(X,E);
y=mat2cell(Y,E);

%select only x and y of cells whihc has complete 65 or 97 frames
x=x(cellfun(@(x)size(x,1)==65,x));%%%%%%%%%%%%%%%%may be 65 or 97
y=y(cellfun(@(y)size(y,1)==65,y));%%%%%%%%%%%%%%%%%%%%%


for j=1:size(x,1)
for i=1:size(x{j},1)
    aa(i,j)=x{j}(i)-x{j}(1);
end
end

for j=1:size(x,1)
    for i=1:size(y{j},1)
        bb(i,j)=y{j}(i)-y{j}(1);
    end
end

aa=aa';
bb=bb';
FT=reshape([aa(:) bb(:)]',2*size(aa,1), [])';

%plot 50 cells which has full trajectories

plot(FT(:,1),FT(:,2),'y-')
xlim([-600 600])%%%%%%%%%%%%%%%%
grid on
hold on

plot(FT(:,3),FT(:,4),'m-');
plot(FT(:,5),FT(:,6),'c-');
plot(FT(:,7),FT(:,8),'r-');
plot(FT(:,9),FT(:,10),'g-');
plot(FT(:,11),FT(:,12),'b-');
plot(FT(:,13),FT(:,14),'k-');
plot(FT(:,15),FT(:,16),'y-');
plot(FT(:,17),FT(:,18),'m-');
plot(FT(:,19),FT(:,20),'c-');
plot(FT(:,21),FT(:,22),'r-');
plot(FT(:,23),FT(:,24),'g-');
plot(FT(:,25),FT(:,26),'b-');
plot(FT(:,27),FT(:,28),'k-');
plot(FT(:,29),FT(:,30),'y-');
plot(FT(:,31),FT(:,32),'m-');
plot(FT(:,33),FT(:,34),'c-');
plot(FT(:,35),FT(:,36),'r-');
plot(FT(:,37),FT(:,38),'g-');
plot(FT(:,39),FT(:,40),'b-');
plot(FT(:,41),FT(:,42),'k-');
plot(FT(:,43),FT(:,44),'y-');
plot(FT(:,45),FT(:,46),'m-');
plot(FT(:,47),FT(:,48),'c-');
plot(FT(:,49),FT(:,50),'r-');
plot(FT(:,51),FT(:,52),'g-');
plot(FT(:,53),FT(:,54),'b-');
plot(FT(:,55),FT(:,56),'k-');
plot(FT(:,57),FT(:,58),'y-');
plot(FT(:,59),FT(:,60),'m-');
plot(FT(:,61),FT(:,62),'c-');
plot(FT(:,63),FT(:,64),'r-');
plot(FT(:,65),FT(:,66),'g-');
plot(FT(:,67),FT(:,68),'b-');
plot(FT(:,69),FT(:,70),'k-');
plot(FT(:,71),FT(:,72),'y-');
plot(FT(:,73),FT(:,74),'m-');
plot(FT(:,75),FT(:,76),'c-');
plot(FT(:,77),FT(:,78),'r-');
plot(FT(:,79),FT(:,80),'g-');
plot(FT(:,81),FT(:,82),'b-');
plot(FT(:,83),FT(:,84),'k-');
plot(FT(:,85),FT(:,86),'y-');
plot(FT(:,87),FT(:,88),'m-');
plot(FT(:,89),FT(:,90),'c-');
plot(FT(:,91),FT(:,92),'r-');
plot(FT(:,93),FT(:,94),'g-');
plot(FT(:,95),FT(:,96),'b-');
plot(FT(:,97),FT(:,98),'k-');
plot(FT(:,99),FT(:,100),'k-');

