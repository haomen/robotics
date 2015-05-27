% wsnSolver.m
% HM20150329

% resolve location for each element in given wireless sensor network
% iteratively.

% in which, nodes is given nodes, nn is neighbors, nm is neighbor
% measurements

function solvedNodes=wsnSolver(nodes,nn,nm)
    solvedNodes=[];
    lamda=[];
    u=[];
%     parameters definitions for solver
    delta=5;             %anglular measurement error 5 degree
    eta=0.05;           %distance error 5%
    epsilon=0.05;       %corresponding to eq.(19)
    k0=0.5;              % k value definition
    k1=0.5;
    k2=0.5;
    k3=0.5;
                                                                                                                                                                                                                                                                                                                                                                        
    
    for t=1:10
        for i=1:length(nodes)
            for j=1:length(nn)
                thetaij=nm(i,j*2);
                dij=nm(i,j*2-1);
                Aij=[-sind(thetaij+delta),cosd(theataij+delta);sind(thetaij-delta),-cosd(thetaij-delta);-cosd(thetaij),-sind(thetaij)];
                bij=[0;0;dij*(1-eta)*cosd(delta)];
                dLdxi=0;
                lamdaij=
                if(j~=0)

                end
            end
        end
    end
end