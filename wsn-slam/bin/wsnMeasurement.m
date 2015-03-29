% wsnMeasurement.m
% HM20150327

% this function is to measure nodes and nodes neighbors distance and
% orientation angle, return as longer array.

% input node is Nx2 matrix, nn uis NxN matrix, noise range defines scope of
% measurement, therefore noise is 0~noisescope
% output nm is Nx2N matrix, restoring (d, theta) for every neighbor points


function nm=wsnMeasurement(node, nn, noiserange) 
    nm=[];
    for i=1:length(node)
        for j=1:length(nn(i,:))
%             if neighbor is valid, then calculate measurement and return
            if(nn(i,j)~=0)
                nx=node(nn(i,j),1)+rand(1)*noiserange;
                ny=node(nn(i,j),2)+rand(1)*noiserange;
                dx=nx-node(i,1);
                dy=ny-node(i,2);
                r=sqrt(dx.^2+dy.^2);
                theta=atand(dy/dx);
                if(dx<0)
                    theta=theta+180;
                end
                if(dx>0 && dy<0)
                    theta=theta+360;
                end
                    
                nm(i,j*2-1)=r;
                nm(i,j*2)=theta;
            end
        end
    end
end