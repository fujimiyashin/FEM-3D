classdef Node   
    properties
        id;
        x;
        y;
        z;
    end
    
    methods
        function obj = Node(id, x, y, z)
            obj.id = id;
            obj.x = x;
            obj.y = y;
            obj.z = z;
        end
    end
end

