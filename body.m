classdef body < handle & dynamicprops & matlab.mixin.Copyable
    properties
        name;                               % Name of the body object
        number;
        nodes;                              % Neighbours
        links;                              % Initial position
        elements;                              % Current position

        % Mechanical properties
        L;
        h;
        b;
        E;
        G;
        N;
        mu;
        rho;
        alphaC;
        betaC;
        zeta;
        
        % Piezo properties
        d31 = 0;
        d33 = 0;
        s11 = 0;
        eta33 = 0;
        eta0 = 8.854187817620e-12;
        initA;
        A;
        
        % Additional properties
        nodeLocations;
        plotElements;
    end
   
    methods
        % node(), Constructor (set parameters at initialisation)
        function obj = body(name)
            obj.name = name;   
        end   
        
        function update(obj,q)
            % d is a displacement vector the same size as obj.nodes

            numNodes = length([obj.nodes]);
            % Nodes
            for i = 1:length([obj.nodes])
                obj.nodes(i).update(q);
            end
            % Links
            for i = 1:length([obj.links])
                obj.links(i).update(obj.nodes);
            end
            % Elements
            for i = 1:length([obj.elements])
                obj.elements(i).update(obj.nodes);
            end
            
            % Area
            obj.initA = sum([obj.elements.initA]);
            obj.A = sum([obj.elements.A]);   
        end
        
        function [eta,e] = setupPiezo(obj)
            
            e = [obj.d33, obj.d31, 0;
                0, 0, 0]'/obj.s11;

            % From page 20 of paper!
            eta = eye(2)*(obj.eta33-obj.d31^2/obj.s11); % eye*eta33 
        end
   end
end