% Node class defenition
classdef element < handle & dynamicprops & matlab.mixin.Copyable
    properties
        name;                               % Name of the patch
        number;                             % Number of the patch
        neighbours;                         % Neighbouring nodes
                
        sBeam = false;                      % Switch to know if it is a smart element
        
        N;                                  % Number of interpolation points (for plotting)
        Ainv;                               % Inverse of interpolation function

        % Mechanical beam properties
        L;                                  % Its length
        h;                                  % Its height
        b;                                  % Its width
        E;                                  % Its E modulus
        mu;                                 % Its poission's ratio
        rho;                                % Its density
        zeta;                               % Its modal damping ratio
        
        % Mechanical piëzo properties
        ph = 0;
        pb = 0;
        pe = 0;
        pE = 0;
        pmu = 0;
        prho = 0;
        
        % Elemental matrices        
        Ke;                                 % Elemental stiffness matrix
        Me;                                 % Elemental mass matrix
        
        % State
        initPos;                            % Initial position 3x2
        pos;                                % Current position 3x2
    end
   
    methods
        % Constructor
        function obj = element(name)
            obj.name = name;
        end
        
        % Assignment function, to quickly assign parameters
        function obj = assign(obj,number,neighbours,initPos,Ke,Me)
            obj.number = number;
            obj.neighbours = neighbours;
            obj.initPos = initPos;
            obj.pos = initPos;
            obj.Ke = Ke;
            obj.Me = Me;
            
            
            % ALso already calculate the inverse of the A matrix for
            % interpolation so that that is not necessary in the loop
            % later
            
           A = [1, 0, 0, 0;                 % p87
               0, 1, 0, 0;
               1, obj.L, obj.L^2, obj.L^3;
               0, 1, 2*obj.L, 3*obj.L^2]; 
           obj.Ainv = eye(4)/A;
        end
        
        % update(), Updates the position of the node with the displacement
        function update(obj,d)
            n1 = obj.neighbours(1);
            n2 = obj.neighbours(2);
            
            xpos1 = d(n1*2-1);
            xpos2 = d(n2*2-1);
            ypos1 = obj.initPos(2,1);
            ypos2 = obj.initPos(2,2);
            tpos1 = d(n1*2);
            tpos2 = d(n2*2);
            
            obj.pos = [xpos1,xpos2;
                        ypos1,ypos2;
                        tpos1,tpos2];
        end
        
        % Interpolation between nodes
        function x = interp(obj,eta)
            % Use cubic interpolation function to find x value inbetween
            % nodes (on paga 87 of FE book)
            
           p = [obj.pos(1,1);
               obj.pos(3,1);
               obj.pos(1,2);
               obj.pos(3,2)];
           
           a = obj.Ainv*p;
           
           y = eta*obj.L;
           x = [1, y, y^2, y^3]*a;
        end
        
        % Plot the element
        function p = show(obj,ax,plotSettings)
            if isempty(ax)
                ax = gca;
            end
            
            % Plot this element
           eta = linspace(0,1,10);
           x = zeros(1,10);
           y = eta*obj.L + obj.pos(2,1);
           for i = 1:10
               x(i) = obj.interp(eta(i));
           end
           
           if obj.sBeam == true
               color = [0.8500 0.3250 0.0980];
           else
               color = [0 0.4470 0.7410];
           end
           
           bplot = plot(ax,x,y,'Color',color);
           p = bplot;
           if plotSettings.plotNodes == true
               nplot = plot(ax,x(end),y(end),'o','Color',[0 0.4470 0.7410]);
               p = [p,bplot, nplot];
               
               if plotSettings.nodeNumbers == true
                   ntext = text(ax,x(end)+0.01,y(end),num2str(obj.number),'horizontalAlignment','left');
                   p = [p,ntext];
               end
           end
           
           if plotSettings.elementNumbers == true   
               if obj.sBeam == false
                   color = 'k';
               end
               etext = text(ax,mean(x)+0.01,mean(y),num2str(obj.number),'horizontalAlignment','left','Color',color);
               p = [p,etext];
           end
           
        end
   end
end