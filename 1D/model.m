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
        simulated = 0;          % Is this model simulated or not?
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
        dsys_obs;               % discrete tyme ss model for observers (Possibly erroeunous)
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
            switch obj.simulationSettings.data
                case {"Real"}
                    addpath('./../c2000/Experiments/')
                    load(obj.simulationSettings.dataset);
                    t = data{1}.Values.Time;
                    dt = mean(diff(t));
                    T = t(end)-dt;

                    %{ 
                    Scale strain output, as the scale used in
                    the simulink file is out of date. Changing
                    the gain there would however require all
                    experiments to be re-done, so the scale
                    here is adjusted for calibration of the
                    strain gauges. 
                    %}

                    positionMeasurement = data{5}.Values.Data';
                    strainMeasurement = data{2}.Values.Data'*obj.modelSettings.strainCorrectionGain;

                    % This is required for the final datasets, which do not
                    % have a calibrated strain measurement..
                    strainMeasurement = strainMeasurement - mean(strainMeasurement); % Because I didn't want to calibrate every time


                    % Change some settings
                    obj.simulationSettings.initialError = false;

                    if obj.simulationSettings.trueAccMeasurement == true
                        % Select, scale and lowpass the acceleration
                        % measurements!
                        accelerationMeasurement = double(data{1}.Values.Data')/obj.simulationSettings.accSensitivity;
                        accelerationMeasurement = accelerationMeasurement-mean(accelerationMeasurement); % Needs to be solved for real-time application!
                        accelerationMeasurement = lowpass(accelerationMeasurement,1e3,1/dt);
                    else
                        % Fake acceleration measurements by deriving the
                        % position measurement
                        vMeasurement = diff([positionMeasurement,positionMeasurement(end)])/dt;
                        aMeasurement = diff([vMeasurement,vMeasurement(end)])/dt;
                        accelerationMeasurement = lowpass(aMeasurement,1e3,1/dt);
                    end
                    
                    % Select outputs based on settings
                    if ~isempty(obj.modelSettings.Acc) && ~isempty(obj.modelSettings.patches)
                    yfull = [positionMeasurement; % Position
                            strainMeasurement; % Strain
                            accelerationMeasurement]; % Acceleration

                    elseif isempty(obj.modelSettings.Acc) && ~isempty(obj.modelSettings.patches)
                    yfull = [positionMeasurement; % Position
                            strainMeasurement]; % Strain

                    elseif ~isempty(obj.modelSettings.Acc) && isempty(obj.modelSettings.patches)
                    yfull = [positionMeasurement; % Position
                            accelerationMeasurement]; % Acceleration                        
                    end



                case {"Simulated"}
                    T = obj.simulationSettings.T;
                    dt = obj.simulationSettings.dt;
                    t = 0:dt:T;
                    
                    W = mvnrnd(zeros(obj.nq,1),obj.Q,T/dt+1)';
                    V = mvnrnd(zeros(obj.ny,1),obj.R,T/dt+1)';
                    
                    yfull = zeros(obj.ny,length(t));            % System output
           end
            
            % Some initialisations!!!
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
            ufull_DKF = zeros(obj.nu*(obj.DKF.nd+1),length(t));
            ufull_GDF = zeros(obj.nu,length(t));
            
            q = zeros(obj.nq,1);                                            % Normal time simulation    
            q1_LO = ones(obj.nq,1)*obj.simulationSettings.obsOffset;        % Initial state estimate for LO
            q_KF = ones(obj.nq,1)*obj.simulationSettings.obsOffset;         % Initial state estimate for obj.KF
            q_AKF = zeros(obj.nq+obj.nu*(obj.AKF.nd+1),1);                         % Initial sate estimate for obj.AKF 
                q_AKF(1:obj.nq) = eye(obj.nq,1)*obj.simulationSettings.obsOffset; % Only  beam states offset error
            q_DKF = ones(obj.nq,1)*obj.simulationSettings.obsOffset;       % Initial state estimate for DKF
                u_DKF = zeros(obj.nu*(obj.DKF.nd+1),1);
            q_GDF = zeros(obj.nq,1)*obj.simulationSettings.obsOffset;       % Initial state estimate for GDF
                u_GDF = zeros(obj.nu,1);
            P_KF = zeros(obj.nq);                                        % Initial P matrix kalman filter
            P_AKF = [eye(obj.nq)*obj.AKF.P0, zeros(obj.nq,obj.nu*(obj.AKF.nd+1));
                    zeros(obj.nu*(obj.AKF.nd+1),obj.nq),eye(obj.nu*(obj.AKF.nd+1))*obj.AKF.Pu0];
            P_DKF = eye(obj.nq)*obj.DKF.P0;
                Pu_DKF = eye(obj.nu*(obj.DKF.nd+1))*obj.DKF.Pu0;
            Pq_GDF  = eye(obj.nq)*obj.GDF.Pq0;
                Pu_GDF = eye(obj.nu)*obj.GDF.Pu0;
                Pqu_GDF = eye(obj.nq,obj.nu)*obj.GDF.Pqu0;
            
            if m == 1
                fprintf(['\n Simulating %s ... \n'],obj.name);
                if obj.simulationSettings.batch == true
                    fprintf('    QuTune: %1.1e \n',obj.AKF.QuTune)
                end
            end
            
            startTime = tic;
            
            % Simulation loop!!!! Here, the system is propagated%%%%%%%%%%%%%%%%%%
            for i = 1:length(t)
                % Simulate system
            
                    switch obj.simulationSettings.data
                        case {"Real"}
                            U_1 = U;
                            U(obj.simulationSettings.distInput) = Udist(i);
                            
                            y = yfull(:,i);
                            qfull(:,i) = zeros(obj.nq,1);

                            if obj.modelSettings.posMeasurement == true
                                y_obs = y;
                            else
                                y_obs = y(2:end);
                            end

                        case {"Simulated"}
                            % Pick input
                            U_1 = U;        % Save last input aswell
                            U(obj.simulationSettings.distInput) = Udist(i);
                
                            % Pick noise
                            if obj.simulationSettings.noise == true
                                w = W(:,i);                         % mvnrnd to allow for different covariances of different sensors
                                v = V(:,i);
                            else
                                w = zeros(obj.nq,1);
                                v = zeros(obj.ny,1);
                            end
                            
                            q_1 = q;             % Previous state is the past state. (Yes, very deep indeed)
                            
                            % Propagate discrete time dynamic system
                            q = obj.dsys_sim.A*q_1 + obj.dsys_sim.B*U_1 + obj.Bw*w;           % Propagate dynamics
                            y = obj.dsys_sim.C*q + obj.dsys_sim.D*U + obj.Dv*v;       % Measurement equation
                
                            % Also save full states for plotting (in q space, so modal)
                            qfull(:,i) = q;
                            yfull(:,i) = y;
                            
                            % 
                            if obj.modelSettings.posMeasurement == true
                                y_obs = y;
                            else
                                y_obs = y(2:end);
                            end
                    end

                % Run state estimators
                    % obj.MF ----------------------------------------------------------
                    if any(ismember(obj.simulationSettings.observer,'MF'))
                        % Modal filter
                        q_MF = obj.MF.Psi*y_obs;
        
                        % Save state and output
                        qfull_MF(:,i) = q_MF;
                        yfull_MF(:,i) = obj.dsys_sim.C*q_MF;
                    end
        
                    % LO ----------------------------------------------------------
                    if any(ismember(obj.simulationSettings.observer,'LO'))
                        q_LO = q1_LO;       % Also for the observers
                        q1_LO = obj.dsys_obs.A*q_LO + obj.dsys_obs.B*U + obj.LO.L*(y_obs-obj.dsys_obs.C*q_LO-obj.dsys_obs.D*U);
                        yfull_LO(:,i) = obj.dsys_sim.C*q_LO + obj.dsys_sim.D*U; % Estimated output
                        qfull_LO(:,i) = q_LO;       % Save estimated state
                    end 
        
                    % obj.KF ----------------------------------------------------------           
                    if any(ismember(obj.simulationSettings.observer,'KF'))
                        q_1_KF = q_KF;
                        P_1_KF = P_KF;
                        
                        if obj.KF.stationary == true
                            % Steady state kalman filter (Same as LO, but with kalman gain)
                            q_KF = obj.dsys_obs.A*q_1_KF + obj.dsys_obs.B*U + obj.KF.K*(y_obs-obj.dsys_obs.C*q_KF-obj.dsys_obs.D*U_1);
        
                            % Save state and output
                            qfull_KF(:,i) = q_KF;                    
                            yfull_KF(:,i) = obj.dsys_sim.C*q_KF + obj.dsys_sim.D*U_1;
        
                        else
                            % Time update
                            q_KF = obj.dsys_obs.A*q_1_KF + obj.dsys_obs.B*U_1;
                            P_KF = obj.dsys_obs.A*P_1_KF*obj.dsys_obs.A'+obj.KF.Bw*obj.KF.Q*obj.KF.Bw';

                            % Measurement update
                            q_KF = q_KF + P_KF*obj.dsys_obs.C'/(obj.dsys_obs.C*P_KF*obj.dsys_obs.C'+obj.KF.R)*(y_obs-obj.dsys_obs.C*q_KF-obj.dsys_obs.D*U);
                            P_KF = P_KF - P_KF*obj.dsys_obs.C'/(obj.dsys_obs.C*P_KF*obj.dsys_obs.C'+obj.KF.R)*obj.dsys_obs.C*P_KF;
              
                            % Save state and output
                            qfull_KF(:,i) = q_KF;                    
                            yfull_KF(:,i) = obj.dsys_sim.C*q_KF + obj.dsys_sim.D*U;
                        end
                    end
        
                    % AKF----------------------------------------------------------
                    if any(ismember(obj.simulationSettings.observer,'AKF'))
                        q_1_AKF = q_AKF;
                        P_1_AKF = P_AKF;

                        if obj.AKF.stationary == true
                            % Steady state kalman filter (Same as LO, but with kalman gain)
                            q_AKF = obj.AKF.A*q_1_AKF + obj.AKF.K*(y(2:end,i)-obj.AKF.C*q_1_AKF);
        
                            % Save state and output
                            qfull_AKF(:,i) = q_AKF;                    
                            yfull_AKF(:,i) = [dsys.C(1,:),zeros(1,obj.nu);
                                            obj.AKF.C]*q_AKF;                   
                        else
                            % Time update
                            q_AKF = obj.AKF.A*q_1_AKF;
                            P_AKF = obj.AKF.A*P_1_AKF*obj.AKF.A'+obj.AKF.Bw*obj.AKF.Q*obj.AKF.Bw';

                            % Measurement update
                            q_AKF = q_AKF + P_AKF*obj.AKF.C'/(obj.AKF.C*P_AKF*obj.AKF.C'+obj.AKF.R)*(y_obs-obj.AKF.C*q_AKF);
                            P_AKF = P_AKF - P_AKF*obj.AKF.C'/(obj.AKF.C*P_AKF*obj.AKF.C'+obj.AKF.R)*obj.AKF.C*P_AKF;
               
                            % Save state and output
                            qfull_AKF(:,i) = q_AKF;     

                            if obj.modelSettings.posMeasurement == true
                                yfull_AKF(:,i) = obj.AKF.C*q_AKF;
                            else
                                yfull_AKF(:,i) = [obj.dsys_sim.C(1,:),zeros(1,obj.nu*(obj.AKF.nd+1));
                                                obj.AKF.C]*q_AKF;
                            end
                        end
                    end
        
                    % DKF----------------------------------------------------------
                    if any(ismember(obj.simulationSettings.observer,'DKF'))
                        % Shift one time step
                        q_1_DKF = q_DKF;
                        u_1_DKF = u_DKF;
                        P_1_DKF = P_DKF;
                        Pu_1_DKF = Pu_DKF;

                        % Time update of INPUT estimate
                        u_DKF = u_1_DKF;               
                        Pu_DKF = Pu_1_DKF + obj.DKF.Qu';
                      
                        % Measurement update of INPUT estimate
                        Ku_DKF = Pu_DKF*obj.DKF.D'/(obj.DKF.D*Pu_DKF*obj.DKF.D' + obj.DKF.R);
        
                        u_DKF = u_DKF + Ku_DKF*(y_obs-obj.DKF.C*q_1_DKF-obj.DKF.D*u_DKF);
                        Pu_DKF = Pu_DKF - Ku_DKF*obj.DKF.D*Pu_DKF;
                        
                        % Time update STATE estimate
                        q_DKF = obj.DKF.A*q_1_DKF + obj.DKF.B*u_DKF;  % Dynamics propagation using estimated input
                        P_DKF = obj.DKF.A*P_1_DKF*obj.DKF.A' +  obj.DKF.Q + obj.DKF.B*Pu_DKF*obj.DKF.B';
        
                        % Measurement update STATE estimate
                        K_DKF = P_DKF*obj.DKF.C'/(obj.DKF.C*P_DKF*obj.DKF.C' + obj.DKF.R);
                        
                        q_DKF = q_DKF + K_DKF*(y_obs-obj.DKF.C*q_DKF-obj.DKF.D*u_DKF);
                        P_DKF = P_DKF - K_DKF*obj.DKF.C*P_DKF;
             
                        % Save state,estimated input and output 
                        qfull_DKF(:,i) = q_DKF;    
                        ufull_DKF(:,i) = u_DKF;
                        if obj.modelSettings.posMeasurement == true
                            yfull_DKF(:,i) = obj.DKF.C*q_DKF;
                        else
                            yfull_DKF(:,i) = [obj.dsys_sim.C(1,:); % Add laser measurement output
                                            obj.DKF.C]*q_DKF; 
                        end
                    end
                    
                    % GDF----------------------------------------------------------
                    if any(ismember(obj.simulationSettings.observer,'GDF'))
                        lastwarn(''); % Clear warnID to check for singularity warning

                        q_1_GDF = q_GDF;
                        Pq_1_GDF = Pq_GDF;
                        u_1_GDF = u_GDF;
                        
                        % Time update
                        q_GDF = obj.GDF.A*q_1_GDF + obj.GDF.B*u_1_GDF;
                        Pq_GDF = [obj.GDF.A obj.GDF.B]*[Pq_1_GDF ,Pqu_GDF; Pqu_GDF' Pu_GDF]*[obj.GDF.A'; obj.GDF.B'] + obj.GDF.Q;
                        
                        % Input estimation
                        F = obj.GDF.D;
                        Rt_GDF = obj.GDF.C*Pq_GDF*obj.GDF.C' + obj.GDF.R;
                        M_GDF = (F'/Rt_GDF*F)\F'/Rt_GDF;
                        u_GDF = M_GDF*(y_obs-obj.GDF.C*q_GDF);
                        Pu_GDF = eye(obj.nu)/(obj.GDF.D'/Rt_GDF*obj.GDF.D);

                        % Measurement update
                        K_GDF = Pq_GDF*obj.GDF.C'/Rt_GDF;
                        q_GDF = q_GDF + K_GDF*(y_obs-obj.GDF.C*q_GDF-obj.GDF.D*u_GDF);
                        Pq_GDF = Pq_GDF - K_GDF*(Rt_GDF - obj.GDF.D*Pu_GDF*obj.GDF.D')*K_GDF';
                        Pqu_GDF = -K_GDF*obj.GDF.D*Pu_GDF;
        
                        % Save state, estimated input and output
                        qfull_GDF(:,i) = q_GDF;    
                        ufull_GDF(:,i) = u_GDF;

                        if obj.modelSettings.posMeasurement == true
                            yfull_GDF(:,i) = obj.GDF.C*q_GDF;
                        else
                            yfull_GDF(:,i) = [obj.dsys_sim.C(1,:);
                                            obj.GDF.C]*q_GDF; 
                        end
                    end
        
                    
                    [~,warnID] = lastwarn;
                    if strcmp(warnID,'MATLAB:illConditionedMatrix')
                        obj.simulationData.broken = true;
                        break;
                    else
                        obj.simulationData.broken = false;
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
        
        %% Plot the final results after simulation (different from plotter())
        function [meanAccelerationSections,meanErrSections,stdErrSections] = plotFinalResults(obj)
        set(0,'defaultTextInterpreter','latex');  
        if strcmp(obj.simulationSettings.data,"Real")
            % sensor Dependancy plot & Results!
            plotSensorDependancy = true;
            if plotSensorDependancy == true
                if strcmp(obj.simulationSettings.dataset,"FINALCL1_1")
                    start = 4.82;
                    stop = 53;
                    startControl = 14;
                    stopControl = 25.4;
                    sections = ["Section 1","Section 2","Section 3"];
                elseif strcmp(obj.simulationSettings.dataset,"FINALCL2")
                    start = 4.43;
                    stop = 57.5;
                    startControl = 10.66;
                    stopControl = 22.45;
                    sections = ["Section 4","Section 5","Section 6"];
                else
                    error('Only for sdatasets 1 and 2!')
                end

                simData = obj.simulationData;
                T = 1e4;

                load(obj.simulationSettings.dataset);
                reference = data{3}.Values.Data';
                t = data{1}.Values.Time;

                accelerationMeasurement = double(data{1}.Values.Data')/obj.simulationSettings.accSensitivity;
                accelerationMeasurement = accelerationMeasurement-mean(accelerationMeasurement); % Needs to be solved for real-time application!
                accelerationMeasurement = lowpass(accelerationMeasurement,1e3,1/obj.simulationSettings.dt);
                strainMeasurement = data{2}.Values.Data'*obj.modelSettings.strainCorrectionGain;
                strainMeasurement = strainMeasurement - mean(strainMeasurement); % Because I didn't want to calibrate every time
       
                err = [];
                color = [];
                if ismember("AKF",obj.simulationSettings.observer)
                    AKF_err = simData.yfull_AKF(1,:) - simData.yfull(1,:);
                    err = [err;AKF_err];
                    color = [color; [0.4940 0.1840 0.5560]];
                end
                if ismember("DKF",obj.simulationSettings.observer)
                    DKF_err = simData.yfull_DKF(1,:) - simData.yfull(1,:);
                    err = [err;DKF_err];
                    color = [color; [0.3010 0.7450 0.9330]];
                    
                end
                if ismember("GDF",obj.simulationSettings.observer)
                    GDF_err = simData.yfull_GDF(1,:) - simData.yfull(1,:);
                    err = [err;GDF_err];
                    color = [color; [0.6350 0.0780 0.1840]];
                end
                

                % Devide measurements and errors in 3 sections 
                startSample = cast((start-start)*T+1,"int64");
                startControlSample = cast((startControl-start)*T+1,"int64");
                stopControlSample = cast((stopControl-start)*T+1,"int64");
                stopSample = cast((stop-start)*T+1,"int64");

                tSections = {t([startSample:startControlSample]);
                                t([startControlSample:stopControlSample]);
                                t([stopControlSample:stopSample])};
                strainSections = {strainMeasurement([startSample:startControlSample]);
                                strainMeasurement([startControlSample:stopControlSample]);
                                strainMeasurement([stopControlSample:stopSample])};
                accelerationSections = {accelerationMeasurement([startSample:startControlSample]);
                                accelerationMeasurement([startControlSample:stopControlSample]);
                                accelerationMeasurement([stopControlSample:stopSample])};
                errSections = {};
                errSections{1} = err(:,[startSample:startControlSample]);               % 3 filters for section 1
                errSections{2} = err(:,[startControlSample:stopControlSample]);         % 3 filters for section 3
                errSections{3} = err(:,[stopControlSample:stopSample]);                 % 3 filters for section 4

                meanStrainSections = [mean(abs(strainSections{1})),mean(abs(strainSections{2})),mean(abs(strainSections{3}))];
                meanAccelerationSections = [mean(abs(accelerationSections{1})),mean(abs(accelerationSections{2})),mean(abs(accelerationSections{3}))];
                meanErrSections = [mean(abs(errSections{1}')); % [sections, Filters]
                                    mean(abs(errSections{2}'));
                                    mean(abs(errSections{3}'))];

                stdErrSections = [std(abs(errSections{1}'-repmat(mean(errSections{1},2)',length(errSections{1}),1)));
                                    std(abs(errSections{2}'-repmat(mean(errSections{2},2)',length(errSections{2}),1)));
                                    std(abs(errSections{3}'-repmat(mean(errSections{3},2)',length(errSections{3}),1)));];

                AKFcorrCoeff = corrcoef(meanErrSections(1,:)',meanAccelerationSections');
                DKFcorrCoeff = corrcoef(meanErrSections(2,:)',meanAccelerationSections');
                GDFcorrCoeff = corrcoef(meanErrSections(3,:)',meanAccelerationSections');

                figure()
                subplot(2,1,1)
                    hold on
                    xlabel('time [s]')
                    ylabel('Position [m]')
                    title('True position and reference')
                    xlim([start,stop])
                    plot(simData.t,reference,'Color',[0.8500 0.3250 0.0980])
                    plot(simData.t,simData.yfull(1,:),'Color',[0 0.4470 0.7410])
                    xline(startControl,'r--')
                    xline(stopControl,'r--')
                    legend({'Reference','Position'})
                    text((start+startControl)/2,-4.5e-4,0,sections(1),'horizontalAlignment','center','VerticalAlignment','bottom','color','r');
                    text((startControl+stopControl)/2,-4.5e-4,0,sections(2),'horizontalAlignment','center','VerticalAlignment','bottom','color','r');
                    text((stopControl+stop)/2,-4.5e-4,0,sections(3),'horizontalAlignment','center','VerticalAlignment','bottom','color','r');

                subplot(2,1,2)
                    hold on
                    xlabel('time [s]')
                    ylabel('Err [m]')
                    xlim([start,stop])
                    xline(startControl,'r--')
                    xline(stopControl,'r--')  
                    title('Estimation error over time')
                    
                    for i = 1:length(obj.simulationSettings.observer)
                        plot(simData.t,abs(err(i,:)),'color',color(i,:))
                    end

                    legend(['','',obj.simulationSettings.observer])   
            else
                meanAccelerationSections = 0;
                meanErrSections = 0;
            end
        else
            % When data is simulated, also make nice plots
            obj.simulationSettings.observer = ["AKF" "DKF" "GDF"];
            simData = obj.simulationData;
            
            truePos = simData.yfull(1,:);
            trueStates = simData.qfull;
            trueInput = simData.Udist;

            strainMeasurement = simData.yfull(2,:);
            accMeasurement = simData.yfull(3,:);
            
            MF_pos = simData.yfull_KF(1,:);
            LO_pos = simData.yfull_KF(1,:);            
            KF_pos = simData.yfull_KF(1,:);
            AKF_pos = simData.yfull_AKF(1,:);
            DKF_pos = simData.yfull_DKF(1,:);
            GDF_pos = simData.yfull_GDF(1,:);

            MF_states = simData.qfull_MF;
            LO_states = simData.qfull_LO;
            KF_states = simData.qfull_KF;
            AKF_states = simData.qfull_AKF;
            DKF_states = simData.qfull_DKF;
            GDF_states = simData.qfull_GDF;

            AKF_u = simData.qfull_AKF(end,:);
            DKF_u = simData.ufull_DKF;
            GDF_u = simData.ufull_GDF;

            Method = 'NRMSE';
            MFposfit = goodnessOfFit(MF_pos',truePos',Method);
            LOposfit = goodnessOfFit(LO_pos',truePos',Method);
            KFposfit = goodnessOfFit(KF_pos',truePos',Method);
            AKFposfit = goodnessOfFit(AKF_pos',truePos',Method);
            DKFposfit = goodnessOfFit(DKF_pos',truePos',Method);
            GDFposfit = goodnessOfFit(GDF_pos',truePos',Method);

            MFstate1fit = goodnessOfFit(MF_states(1,:)',trueStates(1,:)',Method);
            LOstate1fit = goodnessOfFit(LO_states(1,:)',trueStates(1,:)',Method);
            KFstate1fit = goodnessOfFit(KF_states(1,:)',trueStates(1,:)',Method);
            AKFstate1fit = goodnessOfFit(AKF_states(1,:)',trueStates(1,:)',Method);
            DKFstate1fit = goodnessOfFit(DKF_states(1,:)',trueStates(1,:)',Method);
            GDFstate1fit = goodnessOfFit(GDF_states(1,:)',trueStates(1,:)',Method);

            MFstate2fit = goodnessOfFit(MF_states(2,:)',trueStates(2,:)',Method);
            LOstate2fit = goodnessOfFit(LO_states(2,:)',trueStates(2,:)',Method);
            KFstate2fit = goodnessOfFit(KF_states(2,:)',trueStates(2,:)',Method);           
            AKFstate2fit = goodnessOfFit(AKF_states(2,:)',trueStates(2,:)',Method);
            DKFstate2fit = goodnessOfFit(DKF_states(2,:)',trueStates(2,:)',Method);
            GDFstate2fit = goodnessOfFit(GDF_states(2,:)',trueStates(2,:)',Method);

            AKFInputFit = goodnessOfFit(AKF_u',trueInput',Method);
            DKFInputFit = goodnessOfFit(DKF_u',trueInput',Method);
            GDFInputFit = goodnessOfFit(GDF_u',trueInput',Method);

            simFits = [];
            colors = [];
            if ismember("MF",obj.simulationSettings.observer)
                simFits = [simFits;
                            MFposfit,MFstate1fit,MFstate2fit,NaN];
                colors = [colors; 0 0.4470 0.7410];
            end
            if ismember("LO",obj.simulationSettings.observer)
                simFits = [simFits;
                            LOposfit,LOstate1fit,LOstate2fit,NaN];
                colors = [colors; 0.8500 0.3250 0.0980];
            end
            if ismember("KF",obj.simulationSettings.observer)
                simFits = [simFits;
                            KFposfit,KFstate1fit,KFstate2fit,NaN];  
                colors = [colors; 0.9290 0.6940 0.1250];
            end
            if ismember("AKF",obj.simulationSettings.observer)
                simFits = [simFits;
                            AKFposfit,AKFstate1fit,AKFstate2fit,AKFInputFit];  
                colors = [colors; 0.4940 0.1840 0.5560];
            end
            if ismember("DKF",obj.simulationSettings.observer)
                simFits = [simFits;
                            DKFposfit,DKFstate1fit,DKFstate2fit,DKFInputFit]; 
                colors = [colors; 0.3010 0.7450 0.9330];
            end     
            if ismember("GDF",obj.simulationSettings.observer)
                simFits = [simFits;
                            GDFposfit,GDFstate1fit,GDFstate2fit,GDFInputFit];     
                colors = [colors; 0.6350 0.0780 0.1840];
            end            

            % Fit score plot!
            figure()
                hold on
                grid on
                b = bar(1:3,simFits(:,1:3)','FaceColor','flat');

                xtips = [];
                ytips = [];
                for j = 1:3
                    for k = 1:length(obj.simulationSettings.observer)
                        b(k).CData(j,:) = colors(k,:);
                        b(k).EdgeColor = [0,0,0];
                        xtips(k,:) = b(k).XEndPoints;
                        ytips(k,:) = b(k).YEndPoints;
                        labels(k,:) = string(round(b(k).YData,3));

                        text(xtips(k,:),ytips(k,:),labels(k,:),'HorizontalAlignment','center',...
                        'VerticalAlignment','bottom');
                    end
                end

                meanFits = mean(simFits,2);
                for j = 1:length(obj.simulationSettings.observer)
                    yline(meanFits(j),'--',obj.simulationSettings.observer(j),'Color',colors(j,:),...
                        'LabelHorizontalAlignment','left','LabelVerticalAlignment','middle','LineWidth',1.5)
                end

                xlim([0.4,3.5]);
                xticks(1:3);
                xticklabels({'Tip position','\eta_1','\eta_2','Input'})
                legend(obj.simulationSettings.observer,'location','northwest')
                ylabel('NRMSE [-]')
                title('Filter fit performance on 3 simulated signals')

            t = simData.t;
            T = 1e4;
            
            N = 3;

            startsLeft = [0.475,0.525];
            startsRight = [0.508,0.511];
            
            selectedLeft = startsLeft(1)*T:startsLeft(2)*T;
            selectedRight = startsRight(1)*T:startsRight(2)*T;

            figure()
            hold on
            sgtitle('Simulated system response')
            subplot(3,1,1)
                hold on
                ylabel('Position [m]')
                plot(t,truePos,'k')
            subplot(3,1,2)
                hold on
                ylabel('Acceleration [m/s$^2$]')
                plot(t,accMeasurement,'k')
            subplot(3,1,3)
                hold on
                ylabel('Strain [-]')
                plot(t,strainMeasurement,'k')

            figure()
            sgtitle('Simulated filter response')
            subplot(3,N,(1:N-1)+(1-1)*N) %
                hold on
                ylabel('Position [m]')
                xlim(startsLeft)
                plot(t(selectedLeft),truePos(selectedLeft),'k')
                plot(t(selectedLeft),AKF_pos(selectedLeft),'Color',[0.4940 0.1840 0.5560])
                plot(t(selectedLeft),DKF_pos(selectedLeft),'Color',[0.3010 0.7450 0.9330])
                plot(t(selectedLeft),GDF_pos(selectedLeft),'Color',[0.6350 0.0780 0.1840])
                xline(startsRight,'Color',[0.4660 0.6740 0.1880],'lineStyle','--','lineWidth',1.5)
                xline(startsRight,'Color',[0.4660 0.6740 0.1880],'lineStyle','--','lineWidth',1.5)
                legend('True','AKF','DKF','GDF','','','Location','northwest')
            subplot(3,N,(1:N-1)+(2-1)*N)
                hold on
                ylabel('$\eta_1$ [-]')
                xlim(startsLeft)
                plot(t(selectedLeft),trueStates(1,selectedLeft),'k')
                plot(t(selectedLeft),AKF_states(1,selectedLeft),'Color',[0.4940 0.1840 0.5560])
                plot(t(selectedLeft),DKF_states(1,selectedLeft),'Color',[0.3010 0.7450 0.9330])
                plot(t(selectedLeft),GDF_states(1,selectedLeft),'Color',[0.6350 0.0780 0.1840])
                xline(startsRight,'Color',[0.4660 0.6740 0.1880],'lineStyle','--','lineWidth',1.5)
                xline(startsRight,'Color',[0.4660 0.6740 0.1880],'lineStyle','--','lineWidth',1.5)
            subplot(3,N,(1:N-1)+(3-1)*N)
                hold on
                ylabel('$\eta_2$ [-]')
                xlim(startsLeft)
                plot(t(selectedLeft),trueStates(2,selectedLeft),'k')
                plot(t(selectedLeft),AKF_states(2,selectedLeft),'Color',[0.4940 0.1840 0.5560])
                plot(t(selectedLeft),DKF_states(2,selectedLeft),'Color',[0.3010 0.7450 0.9330])
                plot(t(selectedLeft),GDF_states(2,selectedLeft),'Color',[0.6350 0.0780 0.1840])
                xline(startsRight,'Color',[0.4660 0.6740 0.1880],'lineStyle','--','lineWidth',1.5)
                xline(startsRight,'Color',[0.4660 0.6740 0.1880],'lineStyle','--','lineWidth',1.5)
%             subplot(3,N,(1:N-1)+(4-1)*N)
%                 hold on
%                 ylabel('Input [N]')
%                 xlim(startsLeft)
%                 plot(t(selectedLeft),trueInput(selectedLeft),'k')
%                 plot(t(selectedLeft),AKF_u(selectedLeft),'Color',[0.4940 0.1840 0.5560])
%                 plot(t(selectedLeft),DKF_u(selectedLeft),'Color',[0.3010 0.7450 0.9330])
%                 plot(t(selectedLeft),GDF_u(selectedLeft),'Color',[0.6350 0.0780 0.1840])
%                 xline(startsRight,'Color',[0.4660 0.6740 0.1880],'lineStyle','--','lineWidth',1.5)
%                 xline(startsRight,'Color',[0.4660 0.6740 0.1880],'lineStyle','--','lineWidth',1.5)

            subplot(3,N,N)
                hold on
                xlim(startsRight)
                plot(t(selectedRight),truePos(selectedRight),'k')
                plot(t(selectedRight),AKF_pos(selectedRight),'Color',[0.4940 0.1840 0.5560])
                plot(t(selectedRight),DKF_pos(selectedRight),'Color',[0.3010 0.7450 0.9330])
                plot(t(selectedRight),GDF_pos(selectedRight),'Color',[0.6350 0.0780 0.1840])
            subplot(3,N,2*N)
                hold on
                xlim(startsRight)
                plot(t(selectedRight),trueStates(1,selectedRight),'k')
                plot(t(selectedRight),AKF_states(1,selectedRight),'Color',[0.4940 0.1840 0.5560])
                plot(t(selectedRight),DKF_states(1,selectedRight),'Color',[0.3010 0.7450 0.9330])
                plot(t(selectedRight),GDF_states(1,selectedRight),'Color',[0.6350 0.0780 0.1840])
            subplot(3,N,3*N)
                hold on
                xlim(startsRight)
                plot(t(selectedRight),trueStates(2,selectedRight),'k')
                plot(t(selectedRight),AKF_states(2,selectedRight),'Color',[0.4940 0.1840 0.5560])
                plot(t(selectedRight),DKF_states(2,selectedRight),'Color',[0.3010 0.7450 0.9330])
                plot(t(selectedRight),GDF_states(2,selectedRight),'Color',[0.6350 0.0780 0.1840])
%             subplot(3,N,4*N)
%                 hold on
%                 xlim(startsRight)
%                 plot(t(selectedRight),trueInput(selectedRight),'k')
%                 plot(t(selectedRight),AKF_u(selectedRight),'Color',[0.4940 0.1840 0.5560])
%                 plot(t(selectedRight),DKF_u(selectedRight),'Color',[0.3010 0.7450 0.9330])
%                 plot(t(selectedRight),GDF_u(selectedRight),'Color',[0.6350 0.0780 0.1840])


                
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
                ylim([0,1.05*obj.modelSettings.L])
                title 'Beam plot'
                Ax = gca;
                
                q = zeros(obj.nq,1);
            end

            xlim([-obj.modelSettings.L/6,obj.modelSettings.L/6]);
            ylim([0,1.05*obj.modelSettings.L]);

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

                laserText = text(Ax,Ax.XLim(1),obj.modelSettings.measurementHeight*obj.modelSettings.L,'Pos.',...
                    'Color','r',...
                    'HorizontalAlignment','Left',...
                    'VerticalAlignment','Bottom');
                simPlots = [simPlots,laserText];
            end

            % Plot input force in beam plot
            if obj.plotSettings.inputForce == true
                forcex = obj.sys.B(obj.Nmodes+1:end,1)'*q(1:obj.Nmodes);
                forcey = obj.modelSettings.forceHeight*obj.modelSettings.L;

                f2 = [forcex+obj.modelSettings.L/15,forcey];
                f1 = [forcex,forcey];
                df = f2-f1;
                
                color = [0.6350 0.0780 0.1840];
                force = quiver(Ax,f1(1),f1(2),df(1),df(2),0,'LineWidth',2,...
                    'MaxHeadSize',2,...
                    'Color',color);
                simPlots = [simPlots, force];

                forceText = text(Ax,f2(1),f2(2),'f','HorizontalAlignment',...
                    'left',...
                    'Color',color);
                simPlots = [simPlots, forceText];
            end
        end
        
        %% Plot the bode plot of this model
        function bod = showBode(obj,Ax,output,input,varargin)
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
            plotoptions.XLabel.FontSize = 12;
            plotoptions.YLabel.FontSize = 12;
            plotoptions.Title.FontSize = 12;
            plotoptions.FreqUnits = 'Hz';
            plotoptions.grid = 'on';
            plotoptions.PhaseWrapping = 'off';

%             plotoptions.Xlim = {[0,10]};
            plotoptions.Ylim = {[-120,20],[-190,190]};

            bod = bodeplot(Ax,obj.sys(output,input),plotoptions,'r');
            hold on

            if ~isempty(varargin)
                color = varargin{1};
                lineHandle = findobj(gcf,'Type','line','-and','Color','r');
                
                set(lineHandle,'Color',color);
            end
            
            if obj.plotSettings.plotBodePeak == true
                [peakMag,wpeak] = getPeakGain(obj.sys(output,input));
                peakFreq = wpeak/(2*pi);
                plot(peakFreq,mag2db(peakMag),'or');
                text(peakFreq*0.9,mag2db(peakMag),['[',num2str(round(peakFreq,2)),',',num2str(round(mag2db(peakMag),2)),']'],'HorizontalAlignment','right');
            end 
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
                ylim([0,1.05*obj.modelSettings.L])
                title 'Beam plot'
                Ax = gca;
            end

            q = zeros(obj.Nmodes*2,1);
            q(n) = obj.plotSettings.modeAmp*5;
            modePlot = obj.showBeam(Ax,q);
        end

        %% Plot multiple mode shapes (no other place to put this)
        function modesplot = showModes(obj,n)
            L = obj.modelSettings.L;
            obj.plotSettings.modeAmp = 0.75e-4;
            obj.plotSettings.sensor = false;
            d = figure('Name','Mode Shapes');
            sgtitle('Mode shapes')
            freqs = [52.9, 286];
            for i = 1:n
                subplot(1,n,i)
                    hold on
                    grid on
                    xlabel 'm'
                    ylabel 'm'
                    axis equal
                    xlim([-L/6,L/6]);
                    ylim([0,1.2*L])
                    title({['Mode ',num2str(i)],[num2str(freqs(i)),' Hz']})
                    simPlots = obj.showMode(gca,i);
            end
            movegui(d,"southwest");
        end
    end
end

