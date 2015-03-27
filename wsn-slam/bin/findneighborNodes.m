% findneighborNodes.m
% HM20150327

% find nearest nodes within srange and return matrix for nearest neighbors,
%  because of no hashmap functions in 2008a(7.6), i have to use matrix to
% restore neighbor nodes index

function nn=findneighborNodes(node,srange)
    nn=[];
    n=length(node);
    for i=1:n
        cnt=1;
        for j=1:n
            if(sqrt((node(i,1)-node(j,1)).^2+(node(i,2)-node(j,2)).^2)<srange)
                nn(i,cnt)=j;
                cnt=cnt+1;
            end
        end
    end
end

