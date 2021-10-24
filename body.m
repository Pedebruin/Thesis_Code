classdef body < handle & dynamicprops & matlab.mixin.Copyable
    properties
        name;                               % Name of the body object
        number;
        nodes;                              % Neighbours
        links;                              % Initial position
        elements;                              % Current position
        L;
        h;
        E;
        N;
        mu;
        rho;
        alphaC;
        betaC;
        zeta;
        
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
   end
end