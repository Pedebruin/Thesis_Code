function P = plotModel(beam,plotSettings,simulationSettings,ax)
    if isempty(ax)
        figure()
        hold on
        axis equal
        ax = gca;
    end
    numNodes = length(beam.nodes);
    numLinks = length(beam.links);
    numFaces = length(beam.faces);
    P = [];
    
    if plotSettings.fill == true
        colormap(turbo);
        x = zeros(numNodes,1);
        y = zeros(numNodes,1);
        
        f = zeros(length(beam.faces{1}.neighbours),numFaces);
        
        for i = 1:numNodes
            x(i) = beam.nodes{i}.pos(1);
            y(i) = beam.nodes{i}.pos(2);
        end

        for i = 1:numFaces
            f(:,i) = beam.faces{i}.neighbours;
        end
        
        switch plotSettings.color
            case {'disp'}
                caxis([0,7.5])
                colorbar(ax)
                disp = zeros(numNodes,1);
                for i = 1:numNodes
                   disp(i) = sqrt(beam.nodes{i}.disp(1)^2+beam.nodes{i}.disp(2)^2);
                end
                NodeColor = abs(disp);
                EdgeColor = 'interp';
                FaceColor = 'interp';  
                
            case {'deform'}
                caxis([0,5e-4])
                deform = zeros(numLinks,1);
                for i = 1:numLinks
                    deform(i) = beam.links{i}.deform;
                end
                
                d = zeros(numNodes,1);
                for i = 1:numNodes
                    d(i) = mean(deform(beam.nodes{i}.neighbours));
                end
                NodeColor = abs(d);
                EdgeColor = 'interp';
                FaceColor = 'interp';

            case {'none'}
                NodeColor = [0.3010 0.7450 0.9330];
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
            linkNumberPlots = beam.links{i}.plotNumber(ax,'k');
            P = [P,linkNumberPlots];
        end
    end

    % Plot Nodes
    for i = 1:numNodes
        % Nodes
        if plotSettings.nodes == true
            nodePlots = beam.nodes{i}.plotNode(ax,'k');
            P = [P,nodePlots];
        end

        % Node Numbers
        if plotSettings.nodeNumbers == true
            nodeNumberPlots = beam.nodes{i}.plotNumber(ax,'k');
            P = [P,nodeNumberPlots];
        end
    end
end
