% Node class defenition
classdef link < handle & dynamicprops & matlab.mixin.Copyable
    properties
        number;                             % Name of the body object
        neighbours;                         % Neighbours
        initPos;
        pos;
    end
   
    methods
        % link(), Constructor (set parameters at initialisation)
        function obj = link(number,neighbours,initPos)
            obj.number = number;
            obj.neighbours = neighbours;
            obj.initPos = initPos;
            obj.pos = initPos;
        end
        
        % update(), Updates the position of the link
        function update(obj,nodes)
            n1 = obj.neighbours(1);
            n2 = obj.neighbours(2);
            
            from = nodes{n1}.pos;
            to = nodes{n2}.pos;         
            
            obj.pos = [from, to];
        end
        
        % plotLink()
        function P = plotLink(obj,axName,color)
            if isempty(axName)
                figure()
                hold on
                axis equal
                axName = gca;
            end
            x = [obj.pos(1,1),obj.pos(1,2)];
            y = [obj.pos(2,1),obj.pos(2,2)];
            P = plot(x,y,'-','color',color);
        end
        
        % plotLinkName(), Plot the name next to the node
        function P = plotNumber(obj,axName,color)
            if isempty(axName)
                figure()
                hold on
                axis equal
                axName = gca;
            end
            x = (obj.pos(1,1)+obj.pos(1,2))/2;
            y = (obj.pos(2,1)+obj.pos(2,2))/2;
            P = text(x, y, num2str(obj.number),'color',color,...
                                            'HorizontalAlignment','left');
        end
   end
end