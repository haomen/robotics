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
% flag is to switch generation betweend different topologies
% flag=0:   random topo generation
% flag=1:   grid format

function node=wsngeneration(xmax,ymax,n,srange,flag,noiseRatio)
    node=[];                % nodes coordinates
    if (flag==0)
       disp('random wsn generation');
        for i=1:n
            node(i,1)=rand(1)*xmax;
            node(i,2)=rand(1)*ymax;
        end
    end
    
%  grid formatation comes out
%  suggest n=m^2. example 100
    if(flag==1)
        disp('grid wsn genration');
        m=cast(sqrt(n),'uint8');
%         fprintf('n=%d, m=%d\n',n,m);
        xi=single(xmax/(m+2));
        yi=single(ymax/(m+2)); 
%         fprintf('xi=%f, yi=%f\n',xi,yi);
        for i=1:m
            for j=1:m
                node(m*(i-1)+j,1)=xi*single(i)*(1+noiseRatio*rand(1));
                node(m*(i-1)+j,2)=yi*single(j)*(1+noiseRatio*rand(1));
            end
        end
    end
end