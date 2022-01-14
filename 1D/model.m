classdef model
    % a 'model' contains all information on a model
    %{
        This class is made to bundle all information relating a system in
        an easy to move way. This allows for easy transfer and saving of
        different models. 
        
        A model object also contains some of its plotting functions like
        the beam plot in model.showBeam() or the bode plot in
        model.showBode().
    %}
    
    properties
        descr;                  % description of this system
        elements;               % Elements of the system in current orientation

        analyticalf1;           % analytical f
        numEl;                  % Number of elements
        numNodes;               % Number of nodes
        Nmodes;                 % Number of modes modelled
        Phi;                    % Eigenvectors
        omega2;                 % Eigenvalues

        plotSettings;           % Plot settings
        modelSettings;          % Model settings
        simulationSettings;     % Simulation settings

        simulationData;         % Data from a simulation

        Q;                      % Process covariance matrix
        R;                      % Measurement covariance matrix
        S;                      % Measurement and process cross covariance
        Bw;                     % Process noise influence matrix
        Dv;                     % Measurement noise influence matrix

        sys;                    % Actual state space system!
    end
    
    methods
        function obj = model(descr)
            obj.descr = descr;
        end
        
        %% Plot the full beam plot for a given state q (modal coordinates)
        function simPlots = showBeam(obj,Ax,q)
            if isempty(Ax)
                figure()
                hold on
                grid on
                xlabel 'm'
                ylabel 'm'
                axis equal
                xlim([-obj.modelSettings.L/6,obj.modelSettings.L/6]);
                ylim([0,1.2*obj.modelSettings.L])
                title 'Beam plot'
                Ax = gca;
            end

            simPlots = [];

            % update elements
            dfull = [obj.Phi, zeros(obj.numNodes*2,obj.Nmodes);         % full states in d space
                    zeros(obj.numNodes*2,obj.Nmodes),obj.Phi]*q;

            for j = 1:length(obj.elements)
                obj.elements(j).update(dfull(:));  
            end

            % Plot beam itself
            for j = 1:obj.numEl
                pel = obj.elements(j).show(Ax,obj.plotSettings); 
                simPlots = [simPlots,pel];
            end

            % Plot laser in beam plot
            if obj.plotSettings.sensor == true
                laserx = obj.sys.C(1,:)*q;
                laser = plot(Ax,[Ax.XLim(1) laserx],[1,1]*obj.modelSettings.measurementHeight*obj.modelSettings.L,'r','lineWidth',2);
                simPlots = [simPlots, laser];

                laserText = text(Ax,Ax.XLim(1),obj.modelSettings.measurementHeight*obj.modelSettings.L,'Laser',...
                    'Color','r',...
                    'HorizontalAlignment','Left',...
                    'VerticalAlignment','Bottom');
                simPlots = [simPlots,laserText];
            end

            % Plot input force in beam plot
            if obj.plotSettings.inputForce == true
                forcex = obj.sys.B(:,1)'*q;
                forcey = obj.modelSettings.forceHeight*obj.modelSettings.L;

                f1 = [forcex-obj.modelSettings.L/15,forcey];
                f2 = [forcex,forcey];
                df = f2-f1;
                
                color = [0.6350 0.0780 0.1840];
                force = quiver(Ax,f1(1),f1(2),df(1),df(2),0,'LineWidth',2,...
                    'MaxHeadSize',2,...
                    'Color',color);
                simPlots = [simPlots, force];

                forceText = text(Ax,f1(1),f1(2),'F','HorizontalAlignment',...
                    'right',...
                    'Color',color);
                simPlots = [simPlots, forceText];
            end
        end

        %% Plot the bode plot of this model
        function bod = showBode(obj,Ax)
            if isempty(Ax)
                figure();
                hold on
                Ax = gca;
            end

            plotoptions = bodeoptions;
            plotoptions.Title.String = '';
            plotoptions.Title.Interpreter = 'latex';
            plotoptions.XLabel.Interpreter = 'latex';
            plotoptions.YLabel.Interpreter = 'latex';
            plotoptions.XLabel.FontSize = 9;
            plotoptions.YLabel.FontSize = 9;
            plotoptions.FreqUnits = 'Hz';
            plotoptions.grid = 'on';
            plotoptions.PhaseWrapping = 'off';

            bod = bodeplot(Ax,obj.sys(1,1),plotoptions);

            % Plot analytical omega1
            xline(Ax,obj.analyticalf1)
        end

        %% Plot a given mode shape
        function modePlot = showMode(obj,Ax,n)
            q = zeros(obj.Nmodes*2,1);
            q(n) = obj.plotSettings.modeAmp*1;
            modePlot = obj.showBeam(Ax,q);
        end
    end
end

