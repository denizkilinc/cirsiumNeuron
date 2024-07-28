classdef resistor < device
    % RESISTOR Device class for a simple resistor.
    % Characteristic equation: i = v/R.
    
    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------
    
    properties (Constant)
        numTerms = 2;
        numCurrentVars = 0;
    end
    
    properties (SetAccess = private)
        R = 0;
    end
    
    methods
        % The constructer for this class
        function obj = resistor(R, varargin)
            % varargin can be name
            % 2 is the number of terminals
            % 0 is the number of current variables
            obj = obj@device(varargin{:});
            if nargin > 0
                obj.R = R;
            end
        end
        
        function s = I(thisResistor, x, ~)
            v = x(1) - x(2);
            i = resistor.devEq(v, thisResistor.R);
            s = [ i;
                 -i];
        end
        
        function s = dI(thisResistor, ~, ~)
            G = 1/thisResistor.R;
            s = [ 1, -1;
                 -1,  1]*G;
        end
        
                
        % Setter for R property
        function set.R(thisResistor, R)
            if R <= 0
                error('circuitElement:badDeviceProperty',...
                    'A resistor must have a positive resistance value!');
            end
            thisResistor.R = R;
        end

    end
    
    methods (Static)
        
        function i = devEq(v, R)
            i = v/R;
        end
        
        function s = Q(~, ~)
            % this function intentinally does nothing
            s = zeros(2,1);
        end
        
        function s = J(~, ~)
            % this function intentinally does nothing
            s = zeros(2,1);
        end
        
        function s = dQ(~, ~)
            % this function intentinally does nothing
            s = zeros(2);
        end
        
        
    end
        
end

