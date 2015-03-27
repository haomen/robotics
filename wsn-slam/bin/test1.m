% integration test as test1 for all module functons
% HM20150327

% call node generation
xmax=100;
ymax=100;
num=50;
srange=15;

nodes=wsngeneration(xmax,ymax,num,srange);

% calculate nearest neighbors
nn=findneighborNodes(nodes,srange);

% plot wsn
plotwsn(nodes,nn);