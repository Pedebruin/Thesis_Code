% Node class defenition
classdef element < handle & dynamicprops & matlab.mixin.Copyable
    properties
        number;                             % Name of the body object
        neighbours;                         % Neighbours
        initPos;                            % Initial position
        initA;
        pos;                                % Current position
        A;
    end
   
    methods
        % node(), Constructor (set parameters at initialisation)
        function obj = element(number,neighbours,initPos)
            obj.number = number;
            obj.neighbours = neighbours;
            obj.initPos = initPos;
            obj.pos = initPos;
            
            obj.initA = obj.Area;
            obj.A = obj.Area;
        end
        
        % update(), Updates the position of the node with the displacement
        function update(obj,nodes)
            for i = 1:length(obj.neighbours)
                n = obj.neighbours(i);
                n = find([nodes.number] == n);
                obj.pos(:,i) = nodes(n).pos;
                
                obj.A = obj.Area;
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
        
        function A = Area(obj)
            x1 = obj.pos(1,1);
            x2 = obj.pos(1,2);
            x3 = obj.pos(1,3);
            y1 = obj.pos(2,1);
            y2 = obj.pos(2,2);
            y3 = obj.pos(2,3);
            
            A = 1/2*det([1, x1, y1;         % Area of the triangle
                         1, x2, y2;
                         1, x3, y3]);            
        end
   end
end