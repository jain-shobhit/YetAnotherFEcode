classdef DielectricMaterial < KirchoffMaterial
    properties 
        % Along with standard properties for the Material class, we have:
        permittivity % [F/m]
    end

    methods
        function self = DielectricMaterial(varargin)
            % call Material Class constructor
            self@KirchoffMaterial(varargin{:})
        end
        
        function permittivity = get.permittivity(self)
            permittivity = self.VACUUM_PERMITTIVITY * self.RELATIVE_PERMITTIVITY;
        end
    end
end