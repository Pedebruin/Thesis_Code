function plotBeam(nodes,links,plotSettings)
    numNodes = length(nodes);
    numLinks = length(links);
    
    figure()
    hold on
    axis equal
    xlabel 'mm'
    ylabel 'mm'
    title 'Beam plot'
    ax = gca;
    P = [];

    % Plot Links
    for i = 1:numLinks
        % Links
        if plotSettings.links == true
            linkPlots = links{i}.plotLink(ax,[0 0.4470 0.7410]);
            P = [P,linkPlots];
        end

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
