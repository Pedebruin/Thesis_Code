function bodeFigure(K,M,Beam,elements,modelSettings,modes)
    fullSys = modalDecomp(K,M,Beam,elements,modelSettings,1:5); % First element is full system
    fullSys.userData = 'Full System';
    
    systems = cell(length(modes)+2,1);
    
    for j = 1:length(modes)
        mode = modes(j);
        sys = modalDecomp(K,M,Beam,elements,modelSettings,mode);
        sys.userData = ['Mode: ',num2str(mode)];
        systems{j} = sys;
    end
    
    systems{end} = fullSys;
    
    truncated = [];
    for j = 1:length(modes)
        if j == 1
            truncated = systems{1};
        else
            truncated = truncated+systems{j};
        end
    end
    truncated.userData = 'Truncated';
    systems{end-1} = truncated;
    
    plotoptions = bodeoptions;
    plotoptions.Title.String = 'Modal decomposition of flexible structure';
    plotoptions.Title.Interpreter = 'latex';
    plotoptions.XLabel.Interpreter = 'latex';
    plotoptions.YLabel.Interpreter = 'latex';
    plotoptions.XLabel.FontSize = 9;
    plotoptions.YLabel.FontSize = 9;
    plotoptions.FreqUnits = 'Hz';
    plotOptions.grid = 'on';
    plotOptions.PhaseWrapping = 'off';
    plotOptions.XLimMode = {'manual','manual'};    
    
    figure()
    hold on;
    grid on;
    bodeAx = gca;
    
    %% Plot modes
    names = cell(length(modes),1);
    for i = 1:length(systems)-2
        linespec = [];
        bodeplot(bodeAx,systems{i},linespec,plotoptions);
        names{i} = systems{i}.userData;
    end
    % Change opacity
    Fh = gcf;                                                   % Handle To Current Figure
    Kids = Fh.Children;                                         % Children
    AxAll = findobj(Kids,'Type','Axes');                        % Handles To Axes
    Ax1 = AxAll(1);                                             % First Set Of Axes
    LinesAx1 = findobj(Ax1,'Type','Line');                      % Handle To Lines
    Ax2 = AxAll(2);                                             % Second Set Of Axes
    LinesAx2 = findobj(Ax2,'Type','Line');                      % Handle To Lines
    
    Opacity = 0.5;
    for j = 1:length(LinesAx1)
        LinesAx1(j).Color = [LinesAx1(j).Color,Opacity];                    
        LinesAx2(j).Color = [LinesAx2(j).Color,Opacity];                     
    end
   
    
    %% Plot truncated and true model
    for i = length(modes)+1:length(systems)
        if i == length(systems)-1       % For the truncated model
            linespec = '-.';
        else                            % For the true model
            linespec = '--k';
        end
        
        bodeplot(bodeAx,systems{i},linespec,plotoptions);
        names{i} = systems{i}.userData;
    end
    
    %% General changes
    % Change thickness
    Fh = gcf;                                                   % Handle To Current Figure
    Kids = Fh.Children;                                         % Children
    AxAll = findobj(Kids,'Type','Axes');                        % Handles To Axes
    Ax1 = AxAll(1);                                             % First Set Of Axes
    LinesAx1 = findobj(Ax1,'Type','Line');                      % Handle To Lines
    Ax2 = AxAll(2);                                             % Second Set Of Axes
    LinesAx2 = findobj(Ax2,'Type','Line');                      % Handle To Lines
    
    LineWidth = 1.2;
    for j = 1:length(LinesAx1)
        LinesAx1(j).LineWidth = LineWidth;                                  % Set ‘LineWidth’
        LinesAx2(j).LineWidth = LineWidth;                                  % Set ‘LineWidth’
    end
    
    names = string(names');
    legend(bodeAx,names);
    
    % Function to obtain the model for a single or multiple modes
    function sys = modalDecomp(K,M,Beam,elements,modelSettings,modes)
        nPatches = length(modelSettings.patches);           % Number of patches
        nsElementsP = modelSettings.nsElementsP;            % Number of elements


        % Get some parameters
        numEl = modelSettings.numEl;
        numNodes = modelSettings.numNodes;
        L = modelSettings.L;
        Nmodes = length(modes);

        % Find mode shapes and natural frequencies
        if modelSettings.Nmodes <= numEl
                [Phi,omega2] = eig(K,M);
        else
            error('Number of elements is smaller then the number of modes requested. Modal decomposition will not work!')
        end

        Phi = Phi(:,modes+2);
        omega2 = omega2(modes+2,modes+2);

        % Check eigenfrequencies
        f = diag(sqrt(omega2)/(2*pi));

        I = (Beam.h^3)/12;
        S = Beam.h;         % Need to check this formulation still!
        lambda = 1.87510407;
        analyticalf1 = lambda^2/(2*pi*L^2)*sqrt(Beam.E*I/(Beam.rho*S));

        delta = 1*L^3/(3*Beam.E*I);

        if abs(f(1) - analyticalf1) > 0.001 && nPatches == 0 && modelSettings.dispInput == false  % Check only when there are no patches
            warning('Eigenfrequencies do not correspond')
        end

        % Damping matrix 
        Cmodal = 2*Beam.zeta*sqrt(omega2);

        % External input B matrix (in d coordinates)
        Bext = zeros(numNodes*2,1);
        fNode = 3;                                  % force Node
        % Best(fNode*2) = 1;
        Bext(end-1) = 1;                          % Force at last node in x direction 


        % Measurement C matrix (in d coordinates) (have to interpolate between two
        % nodes)
            % Find interpolate nodes
        height = modelSettings.measurementHeight*L;
        pos = unique([elements.pos]');
        diff = pos-height;
        lowerNodes = find(diff<0);
        interpNodes = [lowerNodes(end),lowerNodes(end)+1];  

            % Corresponding elements
        for i = 1:numEl
            if interpNodes == elements(i).neighbours          
                interpEl = i;
            elseif interpNodes(1) == elements(i).neighbours(2) && interpNodes(2) > numNodes % If exactly at last node
                interpEl = i-1;
                interpNodes = interpNodes-1;
            end
        end

            % natural coordinate alpha
        alpha = height - pos(interpNodes(1));
        N = [1, alpha, alpha^2, alpha^3]*elements(interpEl).Ainv;
        Cmeasurement = zeros(1,numNodes*2);                         
        Cmeasurement(interpNodes(1)*2-1:interpNodes(1)*2+2) = N;  % Construct C from this (cubic interpolation)


        % Construct full system 
        A = [zeros(Nmodes),eye(Nmodes);
            -omega2,-Cmodal];
        B = [zeros(Nmodes,1);                          
            Phi'*Bext];                       % [Force input, Base force, Base force, Piezo inputs]
        C = [Cmeasurement*Phi, zeros(1,Nmodes)];                 % Laser measurement

        sys = ss(A,B,C,[]);
        end
end