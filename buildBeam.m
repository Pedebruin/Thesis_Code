function [K,M,R,nodes,links,beam] = buildBeam(L,h,N,E,mu,rho,alpha,beta)
%% Set up the model using PDE toolbox!
% This approach uses the meshing tool form the pde toolbox to mesh a
% geometry also defined in this file. 
% Inputs: Parameters of the beam
% Outputs: Mesh matrices, nodal vector and also the entire beam model

    % Create container
    beam = createpde('structural','frequency-planestress');

    % Define Geometry
    rects = zeros(10,N-1);
    ns = [];
    sf = [];
    for i = 1:N-1
        rects(:,i) = [3,4,(i-1)*L/(N-1),(i-1)*L/(N-1),i*L/(N-1),i*L/(N-1),0, h, h, 0]';
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
    structuralDamping(beam,'Alpha',alpha,...
                            'Beta',beta);

    % Set up boundary conditions
    structuralBC(beam,'edge',1,'Constraint','fixed');

    % Generate Mesh
    generateMesh(beam,'Hmax',L/N,...
                        'Hmin',L/N,...
                        'GeometricOrder','linear');
        pdeplot(beam,'NodeLabels','on'); % You can  plot the mesh if you want to

    % Obtain the FEM matrices!
    FEM = assembleFEMatrices(beam,'KMR');

    K = FEM.K;
    M = FEM.M;
    R = FEM.R;
    
%% ALso get the nodes for plotting and stuff
    numNodes = size(beam.Mesh.Nodes,2);
    nodes = {};
    
    connections = K~=0;
    connections = full(connections(1:numNodes,1:numNodes));
    
    for i = 1:numNodes
        neighbours = find(connections(i,:));  % Neighbouring nodes

        initPos = beam.Mesh.Nodes(:,i);       % Position of the node

        nod = node(i,neighbours,initPos);   % Make a node object
        nodes{i} = nod;                         % Place node with other nodes
    end 
    
%% Do same for links 
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
end