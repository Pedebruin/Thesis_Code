% Node class defenition
classdef element < handle & dynamicprops & matlab.mixin.Copyable
    properties
        number;                             % Name of the body object
        neighbours;                         % Neighbours
        initPos;                            % Initial position
        pos;                                % Current position
    end
   
    methods
        % node(), Constructor (set parameters at initialisation)
        function obj = element(number,neighbours,initPos)
            obj.number = number;
            obj.neighbours = neighbours;
            obj.initPos = initPos;
            obj.pos = initPos;
        end
        
        % update(), Updates the position of the node with the displacement
        function update(obj,nodes)
            for i = 1:length(obj.neighbours)
                n = obj.neighbours(i);
                n = find([nodes.number] == n);
                obj.pos(:,i) = nodes(n).pos;
            end
        end
        
        % plotNodeName(), Plot the name next to the node
        function P = plotNumber(obj,axName,color)
            if isempty(axName)
                figure()
                hold on
                axis equal
                axName = gca;
            end     
            P = text(obj.pos(1), obj.pos(2), num2str(obj.number),'color',color,...
                                                            'VerticalAlignment','bottom');
        end
   end
end