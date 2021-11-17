% Node class defenition
classdef node < handle & dynamicprops & matlab.mixin.Copyable
    properties
        number;                             % Name of the body object
        neighbours;                         % Neighbours
        initPos;                            % Initial position
        pos;                                % Current position
        disp;                                  % Displacement from initial
    end
   
    methods
        % node(), Constructor (set parameters at initialisation)
        function obj = node(number,neighbours,initPos)
            obj.number = number;
            obj.neighbours = neighbours;
            obj.initPos = initPos;
            obj.pos = initPos;
            obj.disp = obj.pos-obj.initPos;
        end
        
        % update(), Updates the position of the node with the displacement
        function update(obj,q)
            d(1) = q(obj.number);
            d(2) = q(obj.number+length(q)/2);
            obj.disp = d';
            obj.pos = obj.initPos + d';

        end
        
        % plotNode(), Plot the node in current location
        function P = plotNode(obj,axName,color)
            if isempty(axName)
                figure()
                hold on
                axis equal
                axName = gca;
            end 
            P = plot(axName,obj.pos(1),obj.pos(2),'.','color',color);
        end
        
        % plotNodeName(), Plot the name next to the node
        function P = plotNumber(obj,axName,color)
            if isempty(axName)
                figure()
                hold on
                axis equal
                axName = gca;
            end     
            P = text(axName,obj.pos(1), obj.pos(2), num2str(obj.number),'color',color,...
                                                            'VerticalAlignment','bottom');
        end
   end
end