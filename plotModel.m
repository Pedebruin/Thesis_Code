function P = plotModel(bodies,plotSettings,simulationSettings,ax)
    if isempty(ax)
        figure()
        hold on
        axis equal
        ax = gca;
    end
    % Get node locations!
    nodes = unique([bodies.nodes]);
    nodes = [nodes.number;nodes.pos]';
    nodes = sortrows(nodes,1);
    nodes = nodes(:,[2,3])';
    
    numNodes = length(nodes);
    numLinks = length([bodies.links]);
    P = [];
    
    if plotSettings.fill == true
        colormap(turbo);
        
        beam = [bodies(1).elements];
        beamElements = [beam.neighbours]';
        beamElements = beamElements(:,1:3);
        beamElements = [beamElements(:,end),beamElements];
        
        beamColor = [0.3010 0.7450 0.9330];
        
        if plotSettings.links == true
            EdgeAlpha = 1;
            FaceAlpha = 0.1;
        else
            EdgeAlpha = 0;
            FaceAlpha = 1;
        end
        beamPatch = patch(ax,'Faces',beamElements,'Vertices',nodes',...
                                                'FaceVertexCData',beamColor,...
                                                'EdgeColor',beamColor,...
                                                'FaceColor',beamColor,...
                                                'EdgeAlpha',EdgeAlpha,...
                                                'FaceAlpha',FaceAlpha); 
        P = [P,beamPatch];
        
        if simulationSettings.Nsensors > 0
            sensors = [bodies(2:simulationSettings.Nsensors+1).elements];
            sensorElements = [sensors.neighbours]';
            sensorElements = sensorElements(:,1:3);
            sensorElements = [sensorElements(:,end),sensorElements];

            sensorColor = [0.8500 0.3250 0.0980];   % Nice red
            sensorPatch = patch(ax,'Faces',sensorElements,'Vertices',nodes',...
                                                     'FaceVertexCData',sensorColor,...
                                                    'EdgeColor',sensorColor,...
                                                    'FaceColor',sensorColor,...
                                                    'EdgeAlpha',EdgeAlpha,...
                                                    'FaceAlpha',FaceAlpha);  
            P = [P,sensorPatch];
        end
        if simulationSettings.Nactuators > 0
            actuators = [bodies(simulationSettings.Nsensors+2:end).elements];
            actuatorElements = [actuators.neighbours]';
            actuatorElements = actuatorElements(:,1:3);
            actuatorElements = [actuatorElements(:,end),actuatorElements];

            actuatorColor = [0.4660 0.6740 0.1880]; % Nice green
            actuatorPatch = patch(ax,'Faces',actuatorElements,'Vertices',nodes',...
                                                     'FaceVertexCData',actuatorColor,...
                                                    'EdgeColor',actuatorColor,...
                                                    'FaceColor',actuatorColor,...
                                                    'EdgeAlpha',EdgeAlpha,...
                                                    'FaceAlpha',FaceAlpha);
           P = [P,actuatorPatch];
        end
    end

    % Plot Link numbers
    if plotSettings.linkNumbers == true
        for i = 1:numLinks
            % Link Numbers

            linkNumberPlots = beam.links{i}.plotNumber(ax,'k');
            P = [P,linkNumberPlots];
        end
    end

    % Plot Nodes
    if plotSettings.nodes == true
        for i = 1:numNodes
            nodePlots = beam.nodes{i}.plotNode(ax,'k');
            P = [P,nodePlots];
        
            % Node Numbers
            if plotSettings.nodeNumbers == true
                nodeNumberPlots = beam.nodes{i}.plotNumber(ax,'k');
                P = [P,nodeNumberPlots];
            end
        end
    end
end
