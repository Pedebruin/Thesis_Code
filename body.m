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
        s11 = 0;
        eta33 = 0;
        eta0 = 8.854187817620e-12;
        
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
            % Faces
            for i = 1:length([obj.elements])
                obj.elements(i).update(obj.nodes);
            end
        end
        
        function [C,G0,G1] = setupPiezo(obj,bodies)
            beam = bodies(1);
            
            % From page 20 of paper!
            e31_ = obj.d31/obj.s11;
            eta33_ = obj.eta33-obj.d31^2/obj.s11;
            c11_ = 1/obj.s11;
            
            
            G0 = obj.b*e31_;
            G1 = obj.b*e31_*(obj.h+beam.h)/2;
            C = (eta33_*obj.b*obj.L)/obj.h; 
        end
   end
end