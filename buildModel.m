function [beam,sensors,actuators,FEM,PDEbeam] = buildModel(beam,sensor,actuator,plotSettings,simulationSettings)
%% Set up the model using PDE toolbox!
% This approach uses the meshing tool form the pde toolbox to mesh a
% geometry also defined in this file. 
% Inputs: Parameters of the beam
% Outputs: Mesh matrices, nodal vector and also the entire beam model

%% Build beam geometry
    disp('Building Beam...')  
    % Create container
    PDEbeam = createpde('structural','modal-planestress');

    % Define beam Geometry
    Brects = zeros(10,beam.N-1);
    Bns = [];
    Bsf = [];
    for i = 1:beam.N-1
        Brects(:,i) = [3,4,-beam.h/2,beam.h/2,beam.h/2,-beam.h/2,... % x values
                        (i-1)*beam.L/(beam.N-1),(i-1)*beam.L/(beam.N-1),... % y values (1,2)
                        i*beam.L/(beam.N-1),i*beam.L/(beam.N-1)]';          % y values(3,4)
        if length(num2str(i))==1
            Bn = ['B1','00',num2str(i)];
        elseif length(num2str(i)) == 2
            Bn = ['B1','0',num2str(i)];
        else 
            Bn = ['B1',num2str(i)];
        end
        Bns = [Bns;Bn];

        if ~isempty(Bsf)
            Bsf = [Bsf,'+',Bn];
        else
            Bsf = Bn;
        end
    end
    
%% Add sensors
    Nsensors = simulationSettings.Nsensors;
    if Nsensors*sensor.L >= beam.L
        error('Too much sensors to fit on the beam')
    end
    
    Srects = [];
    Sns = [];
    Ssf = [];
    for i = 1:Nsensors
        sns = [];
        ssf = [];
        srects = zeros(10,sensor.N-1);
        dist = (beam.L-Nsensors*sensor.L)/Nsensors;
        bottom = (i-1)*(dist+sensor.L);
        for j = 1:sensor.N-1
            srects(:,j) = [3,4,-sensor.h/2,sensor.h/2,sensor.h/2,-sensor.h/2,... % x values
                            (j-1)*sensor.L/(sensor.N-1),(j-1)*sensor.L/(sensor.N-1),... % y values (1,2)
                            j*sensor.L/(sensor.N-1),j*sensor.L/(sensor.N-1)]';          % y values(3,4) 
            if length(num2str(j))==1
                sn = ['S',num2str(i),'00',num2str(j)];
            elseif length(num2str(j)) == 2
                sn = ['S',num2str(i),'0',num2str(j)];
            else 
                sn = ['S',num2str(i),num2str(j)];
            end
            Sns = [Sns;sn];
            
            if i == 1
                if ~isempty(ssf)
                    ssf = [ssf,'+',sn];
                else
                    ssf = sn;
                end
            else
                if ~isempty(ssf)
                    ssf = [ssf,'+',sn];
                else
                    ssf = ['+',sn];
                end
            end
        end
        Sns = [Sns,sns];
        Ssf = [Ssf,ssf];  
        
        srects(3:6,:) = srects(3:6,:) - (sensor.h/2+beam.h/2);  % Translate x positions
        srects(7:10,:) = srects(7:10,:) + bottom;               % Translate y positions
        Srects = [Srects,srects];                               % Add to rest of sensors
    end
    
    
%% Add actuators
    Nactuators = simulationSettings.Nactuators;
    if Nactuators*actuator.L >= beam.L
        error('Too much actuators to fit on the beam')
    end
    
    Arects = [];
    Ans = [];
    Asf = [];
    for i = 1:Nactuators
        ans = [];
        asf = [];
        arects = zeros(10,actuator.N-1);
        dist = (beam.L-Nactuators*actuator.L)/Nactuators;
        bottom = (i-1)*(dist+actuator.L);
        for j = 1:actuator.N-1
            arects(:,j) = [3,4,-actuator.h/2,actuator.h/2,actuator.h/2,-actuator.h/2,... % x values
                            (j-1)*actuator.L/(actuator.N-1),(j-1)*actuator.L/(actuator.N-1),... % y values (1,2)
                            j*actuator.L/(actuator.N-1),j*actuator.L/(actuator.N-1)]';          % y values(3,4)      
            if length(num2str(j))==1
                an = ['A',num2str(i),'00',num2str(j)];
            elseif length(num2str(j)) == 2
                an = ['A',num2str(i),'0',num2str(j)];
            else 
                an = ['A',num2str(i),num2str(j)];
            end
            Ans = [Ans;an];
            
            if i == 1
                if ~isempty(asf)
                    asf = [asf,'+',an];
                else
                    asf = an;
                end
            else
                if ~isempty(asf)
                    asf = [asf,'+',an];
                else
                    asf = ['+',an];
                end
            end
        end
        Ans = [Ans,ans];
        Asf = [Asf,asf]; 
        
        arects(3:6,:) = arects(3:6,:) + (actuator.h/2+beam.h/2);    % Translate x positions
        arects(7:10,:) = arects(7:10,:) + bottom;                   % Translate y positions
        Arects = [Arects,arects];                                   % Add to rest of sensors
    end
    
%% Create geometry
    if simulationSettings.Actuators == true && simulationSettings.Sensors ==  true
        rects = [Brects,Srects,Arects];
        sf = [Bsf,'+',Ssf,'+',Asf];
        ns = [Bns',Sns',Ans'];
    elseif siulationSettings.Actuators == true && simulationSettings.Sensors == false
        rects = [Brects,Arects];
        sf = [Bsf,'+',Asf];
        ns = [Bns',Ans'];
    elseif simulationSettings.Actuators == false && simulationSettings.Sensors == true
        rects = [Brects,Srects];
        sf = [Bsf,'+',Ssf];
        ns = [Bns',Sns'];
    end

    g = decsg(rects,sf,ns);

    geometryFromEdges(PDEbeam,g);

%% Apply structural properties
    %pdegplot(PDEbeam,'FaceLabels','on');
    % Beam
    NbeamFaces = (beam.N-1);
    beamFaces = 1:NbeamFaces;
    structuralProperties(PDEbeam,'Face',beamFaces,...
                                'YoungsModulus', beam.E, ...
                               'PoissonsRatio', beam.mu,...
                               'massDensity', beam.rho);
    % Sensors                       
    NsensorFaces = (sensor.N-1)*simulationSettings.Nsensors;
    sensorFaces = NbeamFaces+1:NbeamFaces+NsensorFaces;
    structuralProperties(PDEbeam,'Face',sensorFaces,...
                                'YoungsModulus', sensor.E,...
                                'PoissonsRatio', sensor.mu,...
                                'massDensity', sensor.rho);
    % Actuators
    NactuatorFaces = (actuator.N-1)*simulationSettings.Nactuators;  
    actuatorFaces = NbeamFaces+NsensorFaces+1:NbeamFaces+NsensorFaces+NactuatorFaces;
    structuralProperties(PDEbeam,'Face',actuatorFaces,...
                                'YoungsModulus', actuator.E,...
                                'PoissonsRatio', actuator.mu,...
                                'massDensity', actuator.rho);

    %pdegplot(PDEbeam,'EdgeLabels','on','FaceLabels','on');
    bottomEdge = nearestEdge(PDEbeam.Geometry,[0,0]);
    switch simulationSettings.Input
        case 'force'
            structuralBC(PDEbeam,'edge',bottomEdge,'Constraint','fixed');
        case 'disp'
            structuralBC(PDEbeam,'edge',bottomEdge,'Constraint','roller');
    end
    

%% Mesh 
    % Generate Mesh
    generateMesh(PDEbeam,'Hmax',beam.L/beam.N,...
                        'Hmin',beam.L/beam.N,...
                        'GeometricOrder','quadratic');

    % Obtain the FEM matrices!
    FEM = assembleFEMatrices(PDEbeam,'KM');

    if plotSettings.realMesh == true
        figure();
        pdeplot(PDEbeam,'NodeLabels','on'); % You can  plot the mesh if you want to 
%         figure();
%         pdeplot(PDEbeam,'ElementLabels','on');
    end

%% Extract mesh and put it in the objects!
% But get connectivity graph first
    numNodes = size(PDEbeam.Mesh.Nodes,2);

    numEl = size(PDEbeam.Mesh.Elements,2);
    connections = zeros(numNodes);
    for i = 1:numEl
        elNodes = PDEbeam.Mesh.Elements(:,i);
        edges = nchoosek(elNodes(1:3),2);
        for j = 1:size(edges,1)
            from = edges(j,1);
            to = edges(j,2);
            connections(from,to) = 1;
            connections(to,from) = 1;
        end
    end
    
% Get all the nodes for plotting and stuff
    nodes = node.empty;
    for i = 1:numNodes
        neighbours = find(connections(i,:));  % Neighbouring nodes
        initPos = PDEbeam.Mesh.Nodes(:,i);       % Position of the node
        nodes(i) = node(i,neighbours,initPos);   % Place node with other nodes
    end 
    
% Do same for all links 
    connections = tril(connections,-1);
    k = 1;
    links = link.empty;
    for i = 1:numNodes
        for j = 1:numNodes
            if connections(i,j) == 1
                from = nodes(i).pos;
                to = nodes(j).pos;
                initPos = [from,to];
                
                links(k) = link(k,[i j]',initPos);
                k = k+1;
            end
        end
    end 
    
% And the all elements
    Elements = element.empty;
    numEl = size(PDEbeam.Mesh.Elements,2);
    for i = 1:numEl
        neighbours = PDEbeam.Mesh.Elements(:,i);
        initPos = zeros(2,length(neighbours));
        for j = 1:length(neighbours)
            Bn = neighbours(j);
            initPos(:,j) = nodes(Bn).initPos;
        end
        Elements(i) = element(i,neighbours,initPos);
    end

% Find nodes, faces and links for all bodies
    % Beam
    beamNodes = findNodes(PDEbeam.Mesh,'region','Face',beamFaces);
    beam.nodes = nodes(beamNodes);
   
    beam.nodeLocations = [nodes.pos];
    
    beamElements = findElements(PDEbeam.Mesh,'region','Face',beamFaces);
    beam.elements = Elements(beamElements);
    
    beamLinks = [];
    for j = 1:length(links)
        if ismember(links(j).neighbours,beamNodes)
            beamLinks = [beamLinks,links(j).number];
        end
    end
    beam.links = links(beamLinks);
    
    % Sensors
    sensors = sensor.empty;
    for i = 1:simulationSettings.Nsensors
        sensors(i) = copy(sensor);
        sensors(i).name = [sensors(i).name,num2str(i)];
        sensors(i).number = i;
        
        % nodes
        sensorNodes = findNodes(PDEbeam.Mesh,'region','Face',sensorFaces(i:i+sensor.N-2));
        sensors(i).nodes = nodes(sensorNodes);

        % elements
        sensorElements = findElements(PDEbeam.Mesh,'region','Face',sensorFaces(i:i+sensor.N-2));
        sensors(i).elements = Elements(sensorElements);
        
        % links
        sensorLinks = [];
        for j = 1:length(links)
            if ismember(links(j).neighbours,sensorNodes)
                sensorLinks = [sensorLinks,links(j).number];
            end
        end
        sensors(i).links = links(sensorLinks);
    end

    % Actuators
    actuators = actuator.empty;
    for i = 1:simulationSettings.Nactuators
        actuators(i) = copy(actuator);
        actuators(i).name = [actuators(i).name,num2str(i)];
        actuators(i).number = i;

        % nodes
        actuatorNodes = findNodes(PDEbeam.Mesh,'region','Face',actuatorFaces(i:i+sensor.N-2));
        actuators(i).nodes = nodes(actuatorNodes);

        % elements
        actuatorElements = findElements(PDEbeam.Mesh,'region','Face',actuatorFaces(i:i+sensor.N-2));
        actuators(i).elements = Elements(actuatorElements);
        
        % links
        actuatorLinks = [];
        for j = 1:length(links)
            if ismember(links(j).neighbours,actuatorNodes)
                actuatorLinks = [actuatorLinks,links(j).number];
            end
        end
        actuators(i).links = links(actuatorLinks);
    end  
end