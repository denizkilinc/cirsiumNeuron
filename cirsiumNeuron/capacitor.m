classdef capacitor < device
    % CAPACITOR Device class for a simple capacitor.
    % Characteristic equation: Q = Cv.
    
    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------
    
    properties (Constant)
        numTerms = 2;
        numCurrentVars = 0;
    end
    
    properties (SetAccess = private)
        C = 0;
    end
    
    methods
        % The constructer for this class
        function obj = capacitor(C, varargin)
            % varargin can be name and terminals
            obj = obj@device(varargin{:});
            if nargin > 0
                obj.C = C;
            end
        end
        
        function s = Q(thisCapacitor, varargin)
            x = varargin{1};
            v = x(1) - x(2);
            q = capacitor.devEq(v, thisCapacitor.C);
            s = [ q;
                 -q];
        end
        
        function s = dQ(thisCapacitor, varargin)
            G = thisCapacitor.C;
            s = [ G, -G;
                 -G,  G]; 
        end

        
        % Setter for C property
        function set.C(thisCapacitor, C)
            if C < 0
                error('circuitElement:badDeviceProperty',...
                    'A capacitor must have a positive capacitance value');
            end
            thisCapacitor.C = C;
        end
    end
    
    methods (Static)
        
        function q = devEq(v, C)
            q = C*v;
        end
        
        function s = I(varargin)
            % this function intentinally does nothing
            s = [0;
                 0];
        end
        
        function s = J(varargin)
            % this function intentinally does nothing
            s = [0;
                 0];
        end
        
        function s = dI(varargin)
            % this function intentinally does nothing
            s = [0, 0;
                 0, 0];
        end
        
        
    end
end

