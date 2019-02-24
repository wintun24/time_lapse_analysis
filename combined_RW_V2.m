%this code generates speed, angle, displacment, distance and witht time
%this code based on Rojan's code where she applied threshold
%the output from this code and that from MigrationAnalysis_V2.m is not
%significantly different
%But this code has some extra feature that Rojan analysed for her thesis
%but I don't analsye those feacture 

clear all
close all

%distance=one cell represents distance of one cells in each time step after
%Threshold (including zero movement)
%theta=one cell represnet direction of one cell in each time step
%(direction from orginal) after threhold (afer removing zero/nan)
%thetatovelocity=one cell represent direction of cell from flow direction
%after threshold (after removing nan/zero)
%velocity=velocity of each cell in each step (this calculate on
%displacment in ecah step) after T
%averagevelocity=average velocity of one cell 
%AveVelOvercells=averge velocity of all cells in one video
%AveJustDisp= collectinon of ave distance of all cell after 1 frame, 2 frame,
%3.... upto final frame;here sqrt is below;
%AveDirecSpeedPerTime = collection of ave speed of all cell after 1 frame,2
%frame, 3.. but this speed calcuate on distance (not displacment)
%AveDirectSpeed= this is average speed of all speed (after 1,2,3 and upt to
%final frames,? useful)..
%%
%#############READING FILE AND DEFINE PARAMETERS####################
filename='xy_22c.xls';%change file name
filename1='xy_22.mat';
frameinterval=15;%%%%%%%%%change interval if requried%%%%%%%%%%%%%
frame_num=(960/frameinterval)+1;
A=xlsread(filename);%the whole matrix

cellno=A(:,2);%cell number it will repeate 1 1 1  2 2 2 ..etc
frame=A(:,3);% number of frame 
X=A(:,4);%X coordinate
Y=A(:,5);%Y coordinate

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


E=RESULT(RESULT~=0);%location of end frame for each cell

LAST=A(end,3);
E=[E,LAST];%the cutting points of individual cell
cellcount=size(E,2);%counting number of cells. 


XX=mat2cell(X,E);
YY=mat2cell(Y,E);

%%
%finding distance, velocity and angle
%Initializing
  for j = 1 : cellcount    %cell numbers
      nofRows= size(XX{j},1);
      %d2(1,j)=0;
      for i=2: nofRows
         distance{j}(i-1)= sqrt((XX{j}(i)- XX{j}(i-1))^2+(YY{j}(i)- YY{j}(i-1))^2);
         time(i,j)= (i-1)*frameinterval;
         %distanceInum(i-1,j)=distance{j}(i-1).*3.08; %convert displacement from No. of traveled cells to um                          
         
      end
      AveDis{j}=mean(distance{j}(:));
      summ=0;
      count=0;
      for i=1: size(distance{j},2)
          if distance{j}(i)<= AveDis{j};
             summ= summ+distance{j}(i);
             count=count+1;
          end          
      end
      aveAveDis(j)=summ/count; %finding Threshold less than which the computed distance is not regarded as cell displacement but the user error in clicking
      Theta{j}(1)=atan2((YY{j}(2)-YY{j}(1)),(XX{j}(2)-XX{j}(1))); % angle of a cell between each movement at each time interval
      if Theta{j}(1)>=0 
        ThetatoVelocity{j}(1)=pi-Theta{j}(1);
      else
        ThetatoVelocity{j}(1)=-pi-Theta{j}(1);
      end
      Velocity{j}(1)=distance{j}(1)./frameinterval; % Time step in the source code is frameintervalmin
      Velocity2{j}(1)=distance{j}(1)./frameinterval;
      EnsembleTheta{j}(1)=atan2((YY{j}(2)-YY{j}(1)),(XX{j}(2)-XX{j}(1))) ; % Total angle of a cell position at each time with respect to initial position
      if EnsembleTheta{j}(1)>=0 
        ThetatoVelocity2{j}(1)=pi-EnsembleTheta{j}(1); % angle between the call position and blood flow
      else
        ThetatoVelocity2{j}(1)=-pi-EnsembleTheta{j}(1);
      end
      DDispTime{1,j}(1)=distance{j}(1).^2;
      XX1{j}(1)=XX{j}(1);
      YY1{j}(1)=YY{j}(1);
      XX1{j}(2)=XX{j}(2);
      YY1{j}(2)=YY{j}(2);
      ArrestTime(j)=0;
  end
MeanaveAveDis = mean(aveAveDis); %Threshold for one cell
summ=0;
count=0;
for i=1: size(aveAveDis,2)
    if aveAveDis(i)<= MeanaveAveDis;
       summ= summ+aveAveDis(i);
       count=count+1;
    end          
end
MinaveAveDis=summ/count; % Threshold for all cells, averaged over all cells
% finding cell velocity, displacement and angles after excluding the user error

countDis=1;
  for j = 1 : cellcount    %cell numbers
      for i=2: size(distance{j},2)
          if distance{j}(i)>= MinaveAveDis && (isnan((YY{j}(i+1)-YY{j}(1))/(XX{j}(i+1)-XX{j}(1)))~=1)
            XX1{j}(i+1)=XX{j}(i+1);
            YY1{j}(i+1)=YY{j}(i+1);
            distance{j}(i)=sqrt((XX1{j}(i+1)- XX1{j}(i))^2+(YY1{j}(i+1)- YY1{j}(i))^2);
            Theta{j}(i)=atan2((YY1{j}(i+1)-YY1{j}(i)),(XX1{j}(i+1)-XX1{j}(i))); %generating angle between successive displacement for timestep 
            if Theta{j}(i)>=0 
                ThetatoVelocity{j}(i)=pi-Theta{j}(i);
            else
                ThetatoVelocity{j}(i)=-pi-Theta{j}(i);
            end        
            Velocity{j}(i)=distance{j}(i)./frameinterval; 
            Velocity2{j}(i)=distance{j}(i)./frameinterval;
            EnsembleTheta{j}(i)=atan2((YY{j}(i+1)-YY{j}(1)),(XX{j}(i+1)-XX{j}(1))) ; % Total angle         
            if EnsembleTheta{j}(i)>=0 
                ThetatoVelocity2{j}(i)=pi-EnsembleTheta{j}(i);
            else
                ThetatoVelocity2{j}(i)=-pi-EnsembleTheta{j}(i);
            end
            
          else
              
               XX1{j}(i+1)=XX1{j}(i);
               YY1{j}(i+1)=YY1{j}(i);
               distance{j}(i)=sqrt((XX1{j}(i+1)- XX1{j}(i))^2+(YY1{j}(i+1)- YY1{j}(i))^2);
               Theta{j}(i)=NaN; % no angle is recorded is the cell does not move
               ThetatoVelocity{j}(i)=ThetatoVelocity{j}(i-1); 
               Velocity{j}(i)=0; % cell velocity is zero when cell does not move
               Velocity2{j}(i)=NaN;
               EnsembleTheta{j}(i)=EnsembleTheta{j}(i-1) ; % Total angle of a cell position at each time with respect to initial position         
               ThetatoVelocity2{j}(i)=ThetatoVelocity2{j}(i-1);
               ArrestTime(j)=ArrestTime(j)+1;
          end
          DDispTime{1,j}(i)=distance{j}(i).^2;
          final_DP{j}=sum(distance{j});%total displacemnt of each cell
          
          
      end
      JustDispX{1}(j) = XX{j}(1);
      JustDispY{1}(j)= YY{j}(1);
      % Finding Average
      Velocity2{j}(isnan(Velocity2{j}))=[];
      AverageVelocity(j,1) = mean((Velocity{j}(:)));
      MaxVel(j)=max(Velocity{j}(:));
      MinVel(j)=min(Velocity{j}(:));
      Theta2{j}=Theta{j};%theta2 is all the angles includein nan which represnets for non-moving cell angle
      Theta{j}(isnan(Theta{j}))=[];%theta is the angles which exclude the cells which are not moving and this theta is not corrected yet.
      FirstAngle(j)=Theta{j}(1);
      
      
      
      if size(distance{j}(:),1)==42
          sumdis(countDis)=sum(distance{j});
          countDis=countDis+1;
      end
  end
  
final_displacement=final_DP';
final_dp=cat(1,final_displacement{:});
final_dp_mm=final_dp.*0.37;%%%%%%%%final ###########TAKE THIS##########THIS IS THE DISPLACEMNT OF ONE CELL (REAL PATH LENGHT)

%%
%theta with time(along or against the flow with time)%after calculation, I
%don't use it cos it seems not very meaningful.

for j=1:size(Theta2,2)
for i=1:frame_num-1
    if i<=size(Theta2{j},2)
        if Theta2 {j}(i)> 1.570796 || Theta2{j}(i)<-1.570796
            Cell_Dir(i,j)=1;%1 is along the flow
        elseif isnan( Theta2{j}(i));
            Cell_Dir(i,j)=nan;
            
        elseif Theta2{j}(i)== 1.570796 || Theta2{j}(i)== -1.570796
            Cell_dir(i,j)=3;%perpendicular to flow
        else
            Cell_Dir(i,j)=2;%2 is against the flow
       
        end
    else
        Cell_Dir(i,j)=nan;
    end
 end
end
for j=1
for i=1:size(Cell_Dir,1)
    countone(i,j)=sum(Cell_Dir(i,:)==1);
    counttwo(i,j)=sum(Cell_Dir(i,:)==2);
    countthree(i,j)=sum(Cell_Dir(i,:)==3);
   alongF(i,j)=(countone(i)*100)/(countone(i)+ counttwo(i)+countthree(i)); 
   againstF(i,j)=(counttwo(i)*100)/(countone(i)+counttwo(i)+countthree(i));
   perpenF(i,j)=(countthree(i)*100)/(countone(i)+counttwo(i)+countthree(i));
end
end
%%
%trying to taking out theta from cell and combine
 
Tpose_theta=cellfun(@transpose,Theta,'UniformOutput',false);%transpose all contant in each cell
final_theta=cat(1,Tpose_theta{:});%combine all cell vertically 

%%
%looking for displacment with time of each cell
for j=1:size(distance,2)
for i=1:frame_num-1
    if i<=size(distance{j},2)
    
    DPwithT (j,i)=sum(distance{j}(1:i));
    
    else
        DPwithT(j,i)=nan;
    end
    
end

end
%%
%looking for average speed with time for all cell
for i=1:frame_num-1
    mean_DPwithT(i)=nanmean(DPwithT(:,i));
    mean_speedwithT(i)=mean_DPwithT(i)/(frameinterval*i);
end
final_mean_speedwithT=mean_speedwithT';
final_mean_speedwithT_mm=final_mean_speedwithT*0.37;%############TAKE THIS
%%
%trying to correct theta to make inot 2pi complete direction

for i=1:size(final_theta)
    if  final_theta(i)<0
        corrected_final_theta(i)=6.28318-(final_theta(i)*(-1));
    else
        corrected_final_theta(i)=final_theta(i);
    end
end
corrected_theta=corrected_final_theta';%%%%%%%%%% this is the corrected theta without non-moving cells%%%%%% TAKE THIS
%%
  MaxMaxVel=max(MaxVel);
  MinMinVel=max(MinVel);
  AveVelOverCells = mean(AverageVelocity);
  AVeVelOverCells_mm=AveVelOverCells*0.37;%%%################TAKE THIS
  %ArrestCoe=mean(ArrestTime)*10/max(max(time));
  countDis=0;
 
  %MeanSumDis=mean(sumdis);
  %SpeedofMeanSumDis = MeanSumDis/420;
  %%-----------------------------------------------------------------------
%Staright Distance (straight-line length ) from the firt position (start point) and the endpoint after time step: 10, 20, 30, ... min 
for k =2 : size(time,1)
for j = 1 : cellcount    %cell numbers
    if size(XX1{j},2)>=k
        JustDispX{k-1}(j) = XX1{j}(k)- XX1{j}(1);
        JustDispY{k-1}(j)= YY1{j}(k)- YY1{j}(1);
        JustDisp{k-1}(j)= JustDispX{k-1}(j)^2+JustDispY{k-1}(j)^2;
        DwithT(k-1,j)=sqrt(JustDisp{k-1}(j));
     
        
    else
       JustDispX{k-1}(j) = NaN;
       JustDispY{k-1}(j)= NaN ;
       JustDisp{k-1}(j)=NaN;%distance after each time step
       DwithT(k-1,j)=nan;
       
    end
end
%JustDisp{k-1}(isnan(JustDisp{k-1}))=[];
AveJustDisp(k-1)=nanmean(JustDisp{k-1});% average distance after k time step of all frist frame , all 2nd ...3,4,5th frame
AveDirecSpeedPerTime(k-1)=nanmean(sqrt(JustDisp{k-1}))/((k-1)*frameinterval);%average speed of cell after 1,2,3,4 frame etc 
end
AveDirectSpeed=(nanmean(AveDirecSpeedPerTime));
%%
%looking for dbyD with time
Tpose_DP=DPwithT';
dbyDPwithtime=DwithT./Tpose_DP;%dbyD with time for each cell
for j=1
for i=1:size(dbyDPwithtime)
    ave_dbyDwithtime(i,j)=nanmean(dbyDPwithtime(i,:));%dbyD with tiem for all cells in movie#################TAKE THIS
end
end

%%
%%%%%%%%%%%% I DON'T USE BELOW %%%%%%%%%%%%%%%%%%%%
% %BELOW IS TO DRAW LOG LOG PLOT
% %AveJustDisp(size(time,1))=nanmean(JustDisp{size(time,1)-1});
%  timme= frameinterval:frameinterval:960;%%%%%%%%%%%%%%%%%
% % gaussEqn = 'a*exp(-(x/b))+c*x'
% % startPoints = [AveJustDisp(1) timme(1) 19];
% % f = fit(timme.',AveJustDisp.',gaussEqn,'Start', startPoints)
% %f=fit(timme.',AveJustDisp.','poly1','Exclude', [1 2 3 4 ])
% f=fit(timme.',AveJustDisp.','poly1');
% figure
% plot(f,timme,AveJustDisp)
% title('MSD curve','FontSize', 20);
% xlabel('Time (min)','FontSize', 20);
% ylabel('Displacement^2 (pixel^2)','FontSize', 20);
% %text(300,1.7*10^4,['<d^2> = ' num2str(f.p1,'%.2f'),' t ',num2str(f.p2,'%.1f')])
% 
% 
% figure
% loglog(timme,AveJustDisp,'.')
% title('MSD log curve','FontSize', 20);
% xlabel('log Time (min)','FontSize', 20);
% ylabel('log Displacement^2 (pixel^2)','FontSize', 20);
% p = polyfit(log(timme), log(AveJustDisp), 1);
% % compute fit in linear domain
% y_hat = exp(p(1) * log(timme) + p(2));
% % make log log plot
% hold on
% loglog(timme, y_hat);
% 
% label = ['log(y) = ' num2str(p(1)) 'log(x) + ' num2str(p(2))];
% legend('data', label);
% 
% 
% %----------------finding Cm by bestfit using least squares--------------
% dt=10; %delta t is 10 min in experiment
% sumt=0;
% sumt2=0;
% sumy = 0;
% sumyt=0;
% n1 = size(time,1)-1;
% n2 = 5;
% n = n1-n2+1
% for i=n2:n1
%     t=i*dt;
%     y = AveJustDisp(i);
%     sumt = sumt+t;
%     sumt2 = sumt2 + t^2;
%     sumy = sumy + y;
%     sumyt = sumyt + y*t;
% end
% a=(n*sumyt - sumy*sumt)/(n*sumt2 - sumt*sumt);
% b=(sumy - a*sumt)/n;
% Cm = a/4 %our model is 2D for 3D it should be changed to a/6

%%%%%%%%%%%%%%%%%%%%%%

%%
%THE OUTPUT TO TAKE %%%%%%%%%%%%%%%%
%1. Distance = distance one cell, one frame (in pixel)
%2. Theta = angle one cell, one frame (not corrected, not include non- move
%3. corrected_theta = corrected in one circle direction
%3. velocity = velocity one cell one move (in pixel)
%4. final_dp_mm = diplacement (one cell, all move in mm)
%5. DPwithT = displacement one cell, with time (in pixel)
%6. DwithT = distance one cell, with time (in pixel)
%7. ave_dbyDwithtime = d/Dp (all cell, all move)- ratio
%8. AVeVelOverCells_mm = mean velocity over movie (in mm)
%9. final_mean_speedwithT_mm = speed all cell, with time (in mm)
%%
%mean AVERAGE displacement of all cells. This include only the cells thet
%present throuhg out the movies;

 index=find(E==65);
 index=index';
DDPP=final_dp_mm(index);
mean_DDPP=mean(DDPP);%%% TAKE THIS#####

%%
%to find relatinship between angle and velocity (corrected_theta Vs final V
Tpose_V=cellfun(@transpose,Velocity2,'UniformOutput',false);%transpose all contant in each cell
final_V=cat(1,Tpose_V{:});%combine all cell vertically 
final_V=final_V.*0.37;


AVpair=[corrected_theta,final_V];
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
VE_eachQuadrant=averageVelocityEach';
final_STD=STDeach';

save(filename1);



