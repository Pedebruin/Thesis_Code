classdef body < handle & dynamicprops & matlab.mixin.Copyable
    properties
        name;                               % Name of the body object
        nodes;                              % Neighbours
        links;                              % Initial position
        faces;                              % Current position
        state;                              % Displacement from initial
        L;
        h;
        E;
        N;
        mu;
        rho;
        alphaC;
        betaC;
        zeta;
    end
   
    methods
        % node(), Constructor (set parameters at initialisation)
        function obj = body(name)
            obj.name = name;            
        end   
        
        function update(obj,q)
            % d is a displacement vector the same size as obj.nodes
            
            numNodes = length(obj.nodes);
            % Nodes
            for i = 1:numNodes
                obj.nodes{i}.update(q);
            end
            % Links
            for i = 1:length(obj.links)
                obj.links{i}.update(obj.nodes);
            end
            % Faces
            for i = 1:length(obj.faces)
                obj.faces{i}.update(obj.nodes);
            end
        end  
   end
end