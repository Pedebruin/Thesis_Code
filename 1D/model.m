classdef model < handle & dynamicprops & matlab.mixin.Copyable
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
        name;                   % Model name
        number;
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
        observers;              % Observers for this model!

        simulationSettings;     % Simulation settings

        simulationData;         % Data from a simulation

        Q;                      % Process covariance matrix
        R;                      % Measurement covariance matrix
        S;                      % Measurement and process cross covariance
        Bw;                     % Process noise influence matrix
        Dv;                     % Measurement noise influence matrix

        ny;                     % Number of outputs (Total)
        nu;                     % Number of inputs 
        nq;                     % Number of states

        MF;                     % Filters
        LO;
        KF;
        AKF;
        DKF;
        GDF;

        sys;                    % Actual state space system!
        dsys_sim;               % discrete time ss model for simulation
        dsys_obs;               % discrete tyme ss model for observers
    end
    
    methods
        %% Constructor
        function obj = model(m,n,z)
             if nargin ~= 0
                obj(m,n,z) = obj;
             end
         end

        %% Simulate the beam model
        function obj = simulate(obj,Udist,m)
           
            % Handy defenitions & Unpacking
            T = obj.simulationSettings.T;
            dt = obj.simulationSettings.dt;
            t = 0:dt:T;

            % Some initialisations!!!
            yfull = zeros(obj.ny,length(t));            % System output
            yfull_MF = zeros(obj.ny,length(t));
            yfull_LO = zeros(obj.ny,length(t)); 
            yfull_KF = zeros(obj.ny,length(t));
            yfull_AKF = zeros(obj.ny,length(t));
            yfull_DKF = zeros(obj.ny,length(t));
            yfull_GDF = zeros(obj.ny,length(t));
        
            U = zeros(obj.nu,1);                                        % Input vector
        
            qfull = zeros(obj.nq,length(t));                            % Full state vector over time
            qfull_MF = zeros(obj.nq,length(t));
            qfull_LO = zeros(obj.nq,length(t));                         % Also for the observers
            qfull_KF = zeros(obj.nq,length(t));
            qfull_AKF = zeros(obj.nq+obj.nu*(obj.AKF.nd+1),length(t));          % Augmented state for obj.AKF
            qfull_DKF = zeros(obj.nq,length(t));
            qfull_GDF = zeros(obj.nq,length(t));
        
            % ufull_AKF can be found in the last two states of qfull_AKF (due tobthe augmented nature)
            ufull_DKF = zeros(obj.nu,length(t));
            ufull_GDF = zeros(obj.nu,length(t));

            q1 = zeros(obj.nq,1);                                       % Normal time simulation    
            q1_LO = ones(obj.nq,1)*obj.simulationSettings.obsOffset;        % Initial state estimate for LO
            q1_KF = ones(obj.nq,1)*obj.simulationSettings.obsOffset;        % Initial state estimate for obj.KF
            q1_AKF = zeros(obj.nq+obj.nu*(obj.AKF.nd+1),1);                         % Initial sate estimate for obj.AKF 
                q1_AKF(1:obj.nq) = ones(obj.nq,1)*obj.simulationSettings.obsOffset; % Only  beam states offset error
            q1_DKF = ones(obj.nq,1)*obj.simulationSettings.obsOffset;       % Initial state estimate for DKF
                u1_DKF = zeros(obj.nu,1);
            q1_GDF = ones(obj.nq,1)*obj.simulationSettings.obsOffset;       % Initial state estimate for GDF
        
            P1_KF = zeros(obj.nq);                                        % Initial P matrix kalman filter
            P1_AKF = zeros(obj.nq+obj.nu*(obj.AKF.nd+1));
            P1_DKF = zeros(obj.nq);
                Pu1_DKF = zeros(obj.nu);
            Pq1_GDF  = zeros(obj.nq);
            
            W = mvnrnd(zeros(obj.nq,1),obj.Q,T/dt+1)';
            V = mvnrnd(zeros(obj.ny,1),obj.R,T/dt+1)';
            
            if m == 1
                fprintf(['\n Simulating %s ... \n'...
                        '    patchCov: %1.1e \n'...
                        '    accCov: %1.1e \n'],obj.name,obj.modelSettings.patchCov,obj.modelSettings.accCov);
            end
            startTime = tic;
            
            % Simulation loop!!!! Here, the system is propagated.%%%%%%%%%%%%%%%%%%
            for i = 1:T/dt+1
                % Shift one time step
                q = q1;             % Previous next state is the current state. (Yes, very deep indeed)
                q_LO = q1_LO;       % Also for the observers
                q_KF = q1_KF;
                q_AKF = q1_AKF;
                q_DKF = q1_DKF;
                    u_DKF = u1_DKF;
                q_GDF = q1_GDF;
        
                P_KF = P1_KF;
                P_AKF = P1_AKF;
                P_DKF = P1_DKF;
                    Pu_DKF = Pu1_DKF;
                Pq_GDF = Pq1_GDF;
        
                % Simulate system
                    % Pick input
                    U(obj.simulationSettings.distInput) = Udist(i);
        
                    % Pick noise
                    if obj.simulationSettings.noise == true
                        w = W(:,i);                         % mvnrnd to allow for different covariances of different sensors
                        v = V(:,i);
                    else
                        w = zeros(obj.nq,1);
                        v = zeros(obj.ny,1);
                    end

                    % Propagate discrete time dynamic system
                    q1 = obj.dsys_sim.A*q + obj.dsys_sim.B*U + obj.Bw*w;           % Propagate dynamics
                    y = obj.dsys_sim.C*q + obj.dsys_sim.D*U + obj.Dv*v;       % Measurement equation
        
                    % Also save full states for plotting (in q space, so modal)
                    qfull(:,i) = q;
                    yfull(:,i) = y;
        
                % Run state estimators
                    % obj.MF ----------------------------------------------------------
                    if any(ismember(obj.simulationSettings.observer,'MF'))
                        % Modal filter
                        q_MF = obj.MF.Psi*y(2:end);
        
                        % Save state and output
                        qfull_MF(:,i) = q_MF;
                        yfull_MF(:,i) = obj.dsys_sim.C*q_MF;
                    end
        
                    % LO ----------------------------------------------------------
                    if any(ismember(obj.simulationSettings.observer,'LO'))
                        q1_LO = obj.dsys_obs.A*q_LO + obj.dsys_obs.B*U + obj.LO.L*(y(2:end)-obj.dsys_obs.C*q_LO-obj.dsys_obs.D*U);
                        yfull_LO(:,i) = obj.dsys_sim.C*q_LO + obj.dsys_sim.D*U; % Estimated output
                        qfull_LO(:,i) = q_LO;       % Save estimated state
                    end 
        
                    % obj.KF ----------------------------------------------------------           
                    if any(ismember(obj.simulationSettings.observer,'KF'))
                        if obj.KF.stationary == true
                            % Steady state kalman filter (Same as LO, but with kalman gain)
                            q1_KF = obj.dsys_obs.A*q_KF + obj.dsys_obs.B*U + obj.KF.K*(y(2:end)-obj.dsys_obs.C*q_KF-obj.dsys_obs.D*U);
        
                            % Save state and output
                            qfull_KF(:,i) = q_KF;                    
                            yfull_KF(:,i) = dsys.C*q_KF + dsys.D*U;
        
                        else
                            % Measurement update
                            q_KF = q_KF + P_KF*obj.dsys_obs.C'/(obj.dsys_obs.C*P_KF*obj.dsys_obs.C'+obj.KF.R)*(y(2:end)-obj.dsys_obs.C*q_KF-obj.dsys_obs.D*U);
                            P_KF = P_KF - P_KF*obj.dsys_obs.C'/(obj.dsys_obs.C*P_KF*obj.dsys_obs.C'+obj.KF.R)*obj.dsys_obs.C*P_KF;
        
                            % Time update
                            q1_KF = obj.dsys_obs.A*q_KF + obj.dsys_obs.B*U;
                            P1_KF = obj.dsys_obs.A*P_KF*obj.dsys_obs.A'+obj.KF.Bw*obj.KF.Q*obj.KF.Bw';
        
                            % Save state and output
                            qfull_KF(:,i) = q_KF;                    
                            yfull_KF(:,i) = obj.dsys_sim.C*q_KF + obj.dsys_sim.D*U;
        
                        end
                    end
        
                    % obj.AKF----------------------------------------------------------
                    if any(ismember(obj.simulationSettings.observer,'AKF'))
                        if obj.AKF.stationary == true
                            % Steady state kalman filter (Same as LO, but with kalman gain)
                            q1_AKF = obj.AKF.A*q_AKF + obj.AKF.K*(y(2:end,i)-obj.AKF.C*q_AKF);
        
                            % Save state and output
                            qfull_AKF(:,i) = q_AKF;                    
                            yfull_AKF(:,i) = [dsys.C(1,:),zeros(1,obj.nu);
                                            obj.AKF.C]*q_AKF;                   
                        else
                            % Measurement update
                            q_AKF = q_AKF + P_AKF*obj.AKF.C'/(obj.AKF.C*P_AKF*obj.AKF.C'+obj.AKF.R)*(y(2:end)-obj.AKF.C*q_AKF);
                            P_AKF = P_AKF - P_AKF*obj.AKF.C'/(obj.AKF.C*P_AKF*obj.AKF.C'+obj.AKF.R)*obj.AKF.C*P_AKF;
        
                            % Time update
                            q1_AKF = obj.AKF.A*q_AKF;
                            P1_AKF = obj.AKF.A*P_AKF*obj.AKF.A'+obj.AKF.Bw*obj.AKF.Q*obj.AKF.Bw';
        
                            % Save state and output
                            qfull_AKF(:,i) = q_AKF;                    
                            yfull_AKF(:,i) = [obj.dsys_sim.C(1,:),zeros(1,obj.nu*(obj.AKF.nd+1));
                                            obj.AKF.C]*q_AKF;
                        end
                    end
        
                    % DKF----------------------------------------------------------
                    if any(ismember(obj.simulationSettings.observer,'DKF'))
                        % Measurement update of INPUT estimate
                        Ku_DKF = Pu_DKF*obj.DKF.D'/(obj.DKF.D*Pu_DKF*obj.DKF.D' + obj.DKF.R);
        
                        u_DKF = u_DKF + Ku_DKF*(y(2:end)-obj.DKF.C*q_DKF-obj.DKF.D*u_DKF);
                        Pu_DKF = Pu_DKF - Ku_DKF*obj.DKF.D*Pu_DKF;
        
                        % Measurement update STATE estimate
                        K_DKF = P_DKF*obj.DKF.C'/(obj.DKF.C*P_DKF*obj.DKF.C' + obj.DKF.R);
        
                        q_DKF = q_DKF + K_DKF*(y(2:end)-obj.DKF.C*q_DKF-obj.DKF.D*u_DKF);
                        P_DKF = P_DKF - K_DKF*obj.DKF.C*P_DKF;
        
                        % Time update of INPUT estimate
                        u1_DKF = u_DKF;                 % Random walk!
                        Pu1_DKF = Pu_DKF + obj.DKF.Qu;
        
                        % Time update STATE estimate
                        q1_DKF = obj.DKF.A*q_DKF + obj.DKF.B*u_DKF;  % Dynamics propagation using estimated input
                        P1_DKF = obj.DKF.A'*P_DKF*obj.DKF.A + obj.DKF.Q;
        
                        % Save state,estimated input and output 
                        qfull_DKF(:,i) = q_DKF;    
                        ufull_DKF(:,i) = u_DKF;
                        yfull_DKF(:,i) = [obj.dsys_sim.C(1,:);
                                        obj.DKF.C]*q_DKF; 
                    end
                    
                    % GDF----------------------------------------------------------
                    if any(ismember(obj.simulationSettings.observer,'GDF'))
                        % Input estimation
                        Rt_GDF = obj.GDF.C*Pq_GDF*obj.GDF.C' + obj.GDF.R;
                        M_GDF = (obj.GDF.D'/Rt_GDF*obj.GDF.D)\obj.GDF.D'/Rt_GDF;
                        u_GDF = M_GDF*(y(2:end)-obj.GDF.C*q_GDF);
                        Pu_GDF = eye(obj.nu)/(obj.GDF.D'/Rt_GDF*obj.GDF.D);
        
                        % Measurement update
                        K_GDF = Pq_GDF*obj.GDF.C'/obj.GDF.R;
                        q_GDF = q_GDF + K_GDF*(y(2:end)-obj.GDF.C*q_GDF-obj.GDF.D*u_GDF);
                        Pq_GDF = Pq_GDF - K_GDF*(Rt_GDF - obj.GDF.D*Pu_GDF*obj.GDF.D')*K_GDF';
                        Pqu_GDF = -K_GDF*obj.GDF.D*Pu_GDF;
        
                        % Time update
                        q1_GDF = obj.GDF.A*q_GDF + obj.GDF.B*u_GDF;
                        Pq1_GDF = [obj.GDF.A obj.GDF.B]*[Pq_GDF ,Pqu_GDF; Pqu_GDF' Pu_GDF]*[obj.GDF.A'; obj.GDF.B'] + obj.GDF.Q;
        
                        % Save state, estimated input and output
                        qfull_GDF(:,i) = q_GDF;    
                        ufull_GDF(:,i) = u_GDF;
                        yfull_GDF(:,i) = [obj.dsys_sim.C(1,:);
                                        obj.GDF.C]*q_GDF; 
                    end
        
                % Update waitbar and check cancel button
                if obj.simulationSettings.waitBar == true
                    waitbar(i/(T/dt+1),f,sprintf('Step %d/%d',i,T/dt+1))
                    if getappdata(f,'canceling')
                        break
                    end
                end
            end
        
            % Save simulation data in the larger model struct. 
            obj.simulationData.t = t;
            obj.simulationData.Udist = Udist;
            
            obj.simulationData.qfull(:,:,m) = qfull;
            obj.simulationData.qfull_MF(:,:,m) = qfull_MF;
            obj.simulationData.qfull_LO(:,:,m) = qfull_LO;
            obj.simulationData.qfull_KF(:,:,m) = qfull_KF;
            obj.simulationData.qfull_AKF(:,:,m) = qfull_AKF;
            obj.simulationData.qfull_DKF(:,:,m) = qfull_DKF;
            obj.simulationData.qfull_GDF(:,:,m) = qfull_GDF;
        
            obj.simulationData.yfull(:,:,m) = yfull;
            obj.simulationData.yfull_MF(:,:,m) = yfull_MF;
            obj.simulationData.yfull_LO(:,:,m) = yfull_LO;
            obj.simulationData.yfull_KF(:,:,m) = yfull_KF;
            obj.simulationData.yfull_AKF(:,:,m) = yfull_AKF;
            obj.simulationData.yfull_DKF(:,:,m) = yfull_DKF;
            obj.simulationData.yfull_GDF(:,:,m) = yfull_GDF;
        
            obj.simulationData.ufull_DKF(:,:,m) = ufull_DKF;
            obj.simulationData.ufull_GDF(:,:,m) = ufull_GDF;

            elapsed = toc(startTime);
            obj.simulationData.fit(:,1,m) = obj.evaluateSimulation(obj.simulationSettings.observer,m);

            fprintf(['    Iteration %d \n'...
                     '        Simulation time: %.2f s \n'],m,elapsed)
            fprintf(['        NRMSE: \n']);
            for z = 1:length(obj.simulationSettings.observer)
                observer = obj.simulationSettings.observer(z);
                if strcmp(observer,"AKF")
                    observer = [observer,num2str(obj.AKF.nd)];
                end
                if strcmp(observer,"DKF")
                    observer = [observer,num2str(obj.DKF.nd)];
                end

                fprintf(join(["            ",observer,'-> %.3f \n']),obj.simulationData.fit(z,1,m));
            end
            
            if obj.simulationSettings.waitBar == true
                delete(f);
            end
        end
        
        %% Evaluate the simulated response
        function [fits] = evaluateSimulation(obj,filters,m)
            startTime = 0;          % Start of the RMSE (To ignore transient errors)
            
            if startTime == 0
                startTimeStep = 1;
            else
                startTimeStep = startTime/obj.simulationSettings.dt;
            end

            y_true =  obj.simulationData.yfull(1,startTimeStep:end,m);
            fits = zeros(length(filters),1);
            i = 1;

            Method = 'NRMSE';

            if any(ismember(filters,'MF'))
                y_MF = obj.simulationData.yfull_MF(1,:,m);
                fits(i,1) = goodnessOfFit(y_MF',y_true',Method);
                i = i+1;
            end
            if any(ismember(filters,'LO'))
                y_LO = obj.simulationData.yfull_LO(1,:,m);
                fits(i,1) = goodnessOfFit(y_LO',y_true',Method);
                i = i+1;
            end
            if any(ismember(filters,'KF'))
                y_KF = obj.simulationData.yfull_KF(1,:,m);
                fits(i,1) = goodnessOfFit(y_KF',y_true',Method);                
                i = i+1;
            end
            if any(ismember(filters,'AKF'))
                y_AKF = obj.simulationData.yfull_AKF(1,startTimeStep:end,m);
                fits(i,1) = goodnessOfFit(y_AKF',y_true',Method);
                i = i+1;
            end
            if any(ismember(filters,'DKF'))
                y_DKF = obj.simulationData.yfull_DKF(1,:,m);
                fits(i,1) = goodnessOfFit(y_DKF',y_true',Method);
                i = i+1;
            end
            if any(ismember(filters,'GDF'))
                y_GDF = obj.simulationData.yfull_GDF(1,:,m);
                fits(i,1) = goodnessOfFit(y_GDF',y_true',Method);                
                i = i+1;
            end
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
                
                q = zeros(obj.nq,1);
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

            q = zeros(obj.Nmodes*2,1);
            q(n) = obj.plotSettings.modeAmp*1;
            modePlot = obj.showBeam(Ax,q);
        end
    end
end

