% genration of wireless sensor network
% given range of from 0 to xmax, 0 to ymax and number N of wireless sensors,
% randomly genrate WSN and plot topological map.

% HM20150324

clear;
close all;

xmax=100;                  % x limit
ymax=70;                    % y limit
n=20;                          % z limit
srange=15;                    % sensor com range

node=[];                % nodes coordinates
for i=1:n
    node(i,1)=rand(1,1)*xmax;
    node(i,2)=rand(1,1)*ymax;
end
plot(node(:,1),node(:,2),'ro');

% connect nodes, and restore neighbors in com range
nodeNeighbor=[];
for i=1:n
    for j=1:n
        if(j~=j)
            
        end
    end
end
    
% plot the whole sensor network graph