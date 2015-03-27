% genration of wireless sensor network
% given range of from 0 to xmax, 0 to ymax and number N of wireless sensors,
% randomly genrate WSN and plot topological map.

% HM20150324

% clear;
% close all;
% 
% xmax=100;                  % x limit
% ymax=100;                    % y limit
% n=50;                          % z limit
% srange=20;                    % sensor com range
function node=wsngeneration(xmax,ymax,n,srange)
    node=[];                % nodes coordinates
    for i=1:n
        node(i,1)=rand(1,1)*xmax;
        node(i,2)=rand(1,1)*ymax;
    end
end
