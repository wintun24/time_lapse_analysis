%this code is to calculte the displacement for the cells that exist
%throughoout the frames
totalframenumber=65;
count=1;
filename1='xy_22_displacemnt.mat';
for i= 1:size(distance,2)
    if size(distance{i},2)==totalframenumber-1;
        Total_displacement(count)=sum(distance{i})*0.37;
        count=count+1;
    end
end
T=Total_displacement';
save(filename1,'T')