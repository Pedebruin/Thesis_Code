function [FEM,nodes,links,faces,beam] = buildBeam(L,h,N,E,mu,rho,plotSettings)
%% Set up the model using PDE toolbox!
% This approach uses the meshing tool form the pde toolbox to mesh a
% geometry also defined in this file. 
% Inputs: Parameters of the beam
% Outputs: Mesh matrices, nodal vector and also the entire beam model

    % Create container
    beam = createpde('structural','modal-planestress');

    % Define Geometry
    rects = zeros(10,N-1);
    ns = [];
    sf = [];
    for i = 1:N-1
        rects(:,i) = [3,4,-h/2,h/2,h/2,-h/2,(i-1)*L/(N-1),(i-1)*L/(N-1),i*L/(N-1),i*L/(N-1)]';
        if length(num2str(i))==1
            n = ['Seg','00',num2str(i)];
        elseif length(num2str(i)) == 2
            n = ['Seg','0',num2str(i)];
        else 
            n = ['Seg',num2str(i)];
        end
        ns = [ns;n];

        if ~isempty(sf)
            sf = [sf,'+',n];
        else
            sf = n;
        end
    end

    g = decsg(rects,sf,ns');
    geometryFromEdges(beam,g);

    % Apply material to structure
    structuralProperties(beam,'YoungsModulus', E, ...
                               'PoissonsRatio', mu,...
                               'massDensity', rho);

    structuralBC(beam,'edge',1,'Constraint','roller');
 
    % Generate Mesh
    generateMesh(beam,'Hmax',L/N,...
                        'Hmin',L/N,...
                        'GeometricOrder','linear');

    % Obtain the FEM matrices!
    FEM = assembleFEMatrices(beam,'KM');
    
%% Plot original mesh (Just for validation)
    if plotSettings.realMesh == true
         pdeplot(beam,'NodeLabels','on'); % You can  plot the mesh if you want to         
    end
    
%% Extract mesh and put it in the objects!
% But get connectivity graph first
    numNodes = size(beam.Mesh.Nodes,2);

    numEl = size(beam.Mesh.Elements,2);
    connections = zeros(numNodes);
    for i = 1:numEl
        edges = nchoosek(beam.Mesh.Elements(:,i),2);
        for j = 1:size(edges,1)
            from = edges(j,1);
            to = edges(j,2);
            connections(from,to) = 1;
            connections(to,from) = 1;
        end
    end
    
% Get the nodes for plotting and stuff
    nodes = {};
    for i = 1:numNodes
        neighbours = find(connections(i,:));  % Neighbouring nodes

        initPos = beam.Mesh.Nodes(:,i);       % Position of the node

        nod = node(i,neighbours,initPos);   % Make a node object
        nodes{i} = nod;                         % Place node with other nodes
    end 
    
% Do same for links 
    connections = tril(connections,-1);
    k = 1;
    links = {};
    for i = 1:numNodes
        for j = 1:numNodes
            if connections(i,j) == 1
                from = nodes{i}.pos;
                to = nodes{j}.pos;
                initPos = [from,to];
                
                links{k} = link(k,[i j],initPos);
                k = k+1;
            end
        end
    end 
    
% And the faces
    faces = {};
    numEl = size(beam.Mesh.Elements,2);
    for i = 1:numEl
        neighbours = beam.Mesh.Elements(:,i);
        initPos = zeros(2,length(neighbours));
        for j = 1:length(neighbours)
            n = neighbours(j);
            initPos(:,j) = nodes{n}.initPos;
        end
        faces{i} = face(i,neighbours,initPos);
    end
    
end