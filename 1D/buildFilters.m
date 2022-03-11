function [SYS] = buildFilters(SYS,LO,KF,AKF,DKF,GDF)
    nPatches = length(SYS.modelSettings.patches);           % Number of patches

    % MF setup-------------------------------------------------------------
    if any(ismember(SYS.simulationSettings.observer,"MF"))
        MF.Psi = pinv(SYS.dsys_obs.C);
        SYS.MF = MF;
    else
        MF = [];
    end
    
    % LO setup-------------------------------------------------------------
    if any(ismember(SYS.simulationSettings.observer,"LO"))
        poles = eig(SYS.dsys_obs.A)*(1-LO.poleMovement);                   % Place poles in Discrete time (TODO)
        LO.L = place(SYS.dsys_obs.A',SYS.dsys_obs.C',poles)'; 
    
        SYS.LO = LO;
        if SYS.modelSettings.l == 1
        fprintf(['    Luenberger Observer-> \n'...
                 '        Pole speed: %3.1f%% increase \n'],LO.poleMovement*100)
        end
    end
    
    % KF setup-------------------------------------------------------------
    if any(ismember(SYS.simulationSettings.observer,"KF"))
        KF.Q = SYS.Q*KF.QTune;
        if SYS.modelSettings.posMeasurement == true
            KF.R = SYS.R*KF.RTune;
            KF.S = SYS.S;
            KF.Dv = SYS.Dv;
        else
            KF.R = SYS.R(2:end,2:end)*KF.RTune;                          % Snip off laser measurement
            KF.S = SYS.S(:,2:end);
            KF.Dv = SYS.Dv(2:end,2:end); 
        end
        KF.Bw = SYS.Bw;
        
        if KF.stationary == true
            % Find stationary kalman gain by solving discrete algebraic ricatti equation
            % (p.162 M.Verhaegen
            [~,Kt,KF.poles] = idare(SYS.dsys_obs.A, SYS.dsys_obs.C', KF.Q, KF.R, KF.S, []);
            KF.K = Kt';
        else
            % No additional setup for non-stationary kalman filter.
        end
        
        SYS.KF = KF;
        if SYS.modelSettings.l == 1
            fprintf(['    Kalman Filter-> \n'...
                     '        Stationary: %d \n'],KF.stationary)
        end
    end
    
    % AKF setup------------------------------------------------------------
   
    if any(ismember(SYS.simulationSettings.observer,"AKF"))
        % q -> [q;qd;u] ->nq+nu*nd
    
        if SYS.simulationSettings.batch ==  false
            switch AKF.nd
                case 0
                    % Use zeroth derivative tuning parameter
                    AKF.QuTune = AKF.QuTuned0;
                case 1
                    % Use another tuning parameter
                    AKF.QuTune = AKF.QuTuned1;
                case 2
                    % Use yet another tuning parameter
                    AKF.QuTune = AKF.QuTuned2;
            end
        end
    
    
        % Generate temporary ss model for input dynamics (higher order
        % derivatives can then be modelled!
            uA = [zeros(SYS.nu*AKF.nd,SYS.nu),eye(SYS.nu*AKF.nd);
                    zeros(SYS.nu,SYS.nu*(AKF.nd+1))];
            uB = [zeros(AKF.nd*SYS.nu,SYS.nu);
                    eye(SYS.nu)];
            uC = [eye(SYS.nu),zeros(SYS.nu,(AKF.nd)*SYS.nu)];
            uD = [];
            Uss = ss(uA,uB,uC,uD);
            Uss = c2d(Uss,SYS.simulationSettings.dt,SYS.modelSettings.c2dMethod); % Checked!
    
        if SYS.modelSettings.posMeasurement == true % Depending on settings, snip off laser measurement. 
            AKF.R = SYS.R*AKF.RTune; 
            AKF.S = [SYS.S; zeros(SYS.nu*(AKF.nd),SYS.ny-1)];
            AKF.Dv = SYS.Dv;
        else
            AKF.R = SYS.R(2:end,2:end)*AKF.RTune;                                       % Do you trust your measurements?
            AKF.S = [SYS.S(:,2:end); zeros(SYS.nu*(AKF.nd),SYS.ny-1)];
            AKF.Dv = SYS.Dv(2:end,2:end);
        end
    
        AKF.Bw = [SYS.Bw,zeros(SYS.nq,SYS.nu);
            zeros(SYS.nu*(AKF.nd+1),SYS.nq),Uss.B];
        AKF.A = [SYS.dsys_obs.A, SYS.dsys_obs.B,zeros(SYS.nq,SYS.nu*AKF.nd);
            zeros(SYS.nu*(AKF.nd+1),SYS.nq),Uss.A];
        AKF.C = [SYS.dsys_obs.C, SYS.dsys_obs.D, zeros(SYS.ny-1,SYS.nu*AKF.nd)];
                                  
        AKF.Q = eye(size(SYS.Q,1))*AKF.QTune;                                           % Do you trust your model?
        AKF.Qu = [1,zeros(1,nPatches);
                  zeros(nPatches,1),zeros(nPatches,nPatches)]*AKF.QuTune;                  % Input covariance!!
        AKF.Qu(end-nPatches+1:end,end-nPatches+1:end) = eye(nPatches,nPatches)*1e-20;
    
        AKF.Q = [AKF.Q,zeros(SYS.nq,SYS.nu); % Assemble augmented Q matrix
            zeros(SYS.nu,SYS.nq),AKF.Qu];
    
        if AKF.stationary == true
            % Find stationary kalman gain by solving discrete algebraic ricatti equation
            % (p.162 M.Verhaegen)
            if AKF.nd > 0 || AKF.nd > 1
                error('Ricatti does not work for this augmented system.. (Already in backlog, working on it)')
            end
            [~,Kt,AKF.poles,INFO] = idare(AKF.A, AKF.C', AKF.Q, AKF.R, AKF.S, []);
            AKF.K = Kt';
        else
            % No additional setup for non-stationary augmented kalman filter.
        end
    
        if SYS.modelSettings.l == 1
            fprintf(['    Augmented Kalman Filter-> \n'...
                     '        Stationary: %d \n'...
                     '        nd: %d\n '],AKF.stationary,AKF.nd)
        end
    end
    SYS.AKF = AKF;
    
    % DKF setup------------------------------------------------------------
    
    if any(ismember(SYS.simulationSettings.observer,"DKF"))
        DKF.A = SYS.dsys_obs.A;
        DKF.C = SYS.dsys_obs.C;
        DKF.Bw = SYS.Bw;
                    
        % Pick tuning parameter
        if SYS.simulationSettings.batch == false
            switch DKF.nd
                case 0
                    % Use zeroth derivative tuning parameter
                    DKF.QuTune = DKF.QuTuned0;
                case 1
                    % Use another tuning parameter
                    DKF.QuTune = DKF.QuTuned1;
                case 2
                    % Use yet another tuning parameter
                    DKF.QuTune = DKF.QuTuned2;
            end
        end
   
            % Generate temporary ss model for input dynamics (higher order
            % derivatives can then be modelled!
            uA = [zeros(SYS.nu*DKF.nd,SYS.nu),eye(SYS.nu*DKF.nd);
                    zeros(SYS.nu,SYS.nu*(DKF.nd+1))];
            uB = [zeros(DKF.nd*SYS.nu,SYS.nu);
                    eye(SYS.nu)];
            uC = zeros(1,SYS.nu*(DKF.nd+1));
            uC((0:SYS.nu-1)*(DKF.nd+1)+1) = ones(1,SYS.nu);
    
            uD = [];
            Uss = ss(uA,uB,uC,[]);
            Uss = c2d(Uss,SYS.simulationSettings.dt,SYS.modelSettings.c2dMethod);
        
        if SYS.modelSettings.posMeasurement == true
            DKF.R = SYS.R*DKF.RTune;                             
        else
            DKF.R = SYS.R(2:end,2:end)*DKF.RTune;                             
        end
    
        DKF.uA = Uss.A;
        DKF.uB = Uss.B; 
        DKF.uC = [eye(SYS.nu),zeros(SYS.nu,(DKF.nd)*SYS.nu)];
        
        DKF.D = SYS.dsys_obs.D;
        DKF.B = SYS.dsys_obs.B;
        DKF.Q = eye(size(SYS.Q,1))*DKF.QTune;
        DKF.Qu = eye(SYS.nu)*DKF.QuTune;
    
        if SYS.modelSettings.l == 1
            fprintf(['    Dual Kalman Filter-> \n'...
             '        QuTune: %1.1e \n'],DKF.QuTune)
        end
    end
    SYS.DKF = DKF;
    
    % GDF setup------------------------------------------------------------
    if any(ismember(SYS.simulationSettings.observer,"GDF"))
        GDF.A = SYS.dsys_obs.A;
        GDF.B = SYS.dsys_obs.B;
        GDF.C = SYS.dsys_obs.C;
        GDF.D = SYS.dsys_obs.D;
    
        GDF.Q = eye(size(SYS.Q,1))*GDF.QTune;
        if SYS.modelSettings.posMeasurement == true
            GDF.R = SYS.R*GDF.RTune;
        else
            GDF.R = SYS.R(2:end,2:end)*GDF.RTune;                             % Snipp off laser measurement
        end
    
        if SYS.modelSettings.l == 1
            fprintf('    Giljins de Moor filter-> \n')
        end
    end
    SYS.GDF = GDF;
end