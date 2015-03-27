% plotwsn.m
% HM20150327

% plot wsn with neighbors

function plotwsn(node,nn)

    
    figure;
    % plot nodes and id of node
    plot(node(:,1),node(:,2),'bo');
    labels=cellstr(num2str([1:length(node)]'));
    text(node(:,1),node(:,2), labels);
    hold on;
    
    % display random genenrated wireless sensor network in grid    
    for i=1:length(nn)
        for j=1:length(nn(i,:))
            if(nn(i,j)~=0)
                endIndex=nn(i,j);
                plot([node(i,1);node(endIndex,1)], [node(i,2);node(endIndex,2)],'r');
                hold on;
            end
        end
    end

    grid on;
    axis equal;
    
end