% genration of wireless sensor network
% given range of from 0 to xmax, 0 to ymax and number N of wireless sensors,
% randomly genrate WSN and plot topological map.

% HM20150324

clear;
close all;

xmax=100;                  % x limit
ymax=100;                    % y limit
n=50;                          % z limit
srange=20;                    % sensor com range

node=[];                % nodes coordinates
for i=1:n
    node(i,1)=rand(1,1)*xmax;
    node(i,2)=rand(1,1)*ymax;
end

% plot nodes and id of node
plot(node(:,1),node(:,2),'bo');
labels=cellstr(num2str([1:length(node)]'));
text(node(:,1),node(:,2), labels);
hold on;

% connect nodes, and restore neighbors in com range
nodeNeighbor=[];
for i=1:n
    cnt=1;
    for j=i+1:n
        if(sqrt((node(i,1)-node(j,1)).^2+(node(i,2)-node(j,2)).^2)<srange)
            nodeNeighbor(i,cnt)=j;
            cnt=cnt+1;
        end
    end
end

% display random genenrated wireless sensor network in grid    
for i=1:length(nodeNeighbor)
    for j=1:length(nodeNeighbor(i,:))
        if(nodeNeighbor(i,j)~=0)
            endIndex=nodeNeighbor(i,j);
            plot([node(i,1);node(endIndex,1)], [node(i,2);node(endIndex,2)],'r');
            hold on;
        end
    end
end

grid on;
