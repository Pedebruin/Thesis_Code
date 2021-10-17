function P = plotBeam(nodes,links,faces,plotSettings,ax)
    numNodes = length(nodes);
    numLinks = length(links);
    numFaces = length(faces);
    P = [];
    
    if plotSettings.fill == true
        cm = colormap(turbo);
        x = zeros(numNodes,1);
        y = zeros(numNodes,1);
        
        f = zeros(length(faces{1}.neighbours),numFaces);
        
        for i = 1:numNodes
            x(i) = nodes{i}.pos(1);
            y(i) = nodes{i}.pos(2);
        end

        for i = 1:numFaces
            f(:,i) = faces{i}.neighbours;
        end
        
        
        switch plotSettings.color
            case {'disp'}
                caxis([0,15])
                colorbar(ax)
                disp = zeros(numNodes,1);
                for i = 1:numNodes
                   disp(i) = sqrt(nodes{i}.disp(1)^2+nodes{i}.disp(2)^2);
                end
                NodeColor = abs(disp);
                EdgeColor = 'interp';
                FaceColor = 'interp';  
                
            case {'deform'}
                caxis([0,1.5])
                deform = zeros(numLinks,1);
                for i = 1:numLinks
                    deform(i) = links{i}.deform;
                end
                NodeColor = abs(deform)/max(deform);
                EdgeColor = 'interp';
                FaceColor = 'interp';

            case {'none'}
                NodeColor = [];
                EdgeColor = [0.3010 0.7450 0.9330];     % A nice blue
                FaceColor = [0.3010 0.7450 0.9330];
                
            otherwise                    
        end
        
        if plotSettings.links == true
            EdgeAlpha = 1;
            FaceAlpha = 0.1;
        else
            EdgeAlpha = 0;
            FaceAlpha = 1;
        end
        coloring = patch(ax,'Faces',f','Vertices',[x,y],...
                                                'FaceVertexCData',NodeColor,...
                                                'EdgeColor',EdgeColor,...
                                                'FaceColor',FaceColor,...
                                                'EdgeAlpha',EdgeAlpha,...
                                                'FaceAlpha',FaceAlpha);                                
        P = [P,coloring];
    end

    % Plot Link numbers
    for i = 1:numLinks
        % Link Numbers
        if plotSettings.linkNumbers == true
            linkNumberPlots = links{i}.plotNumber(ax,'k');
            P = [P,linkNumberPlots];
        end
    end

    % Plot Nodes
    for i = 1:numNodes
        % Nodes
        if plotSettings.nodes == true
            nodePlots = nodes{i}.plotNode(ax,'k');
            P = [P,nodePlots];
        end

        % Node Numbers
        if plotSettings.nodeNumbers == true
            nodeNumberPlots = nodes{i}.plotNumber(ax,'k');
            P = [P,nodeNumberPlots];
        end
    end
end
