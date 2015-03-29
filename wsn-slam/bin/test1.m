% integration test as test1 for all module functons
% HM20150327

% call node generation
xmax=100;
ymax=100;
num=80;
srange=20;

nodes=wsngeneration(xmax,ymax,num,srange);

% calculate nearest neighbors
nn=findneighborNodes(nodes,srange);

% plot wsn
close all;
plotwsn(nodes,nn);

% calculate measurements of every paird sensors
nM=wsnMeasurement(nodes,nn,0);