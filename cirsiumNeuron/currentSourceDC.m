classdef currentSourceDC < device
    % CURRENTSOURCEDC Device class for an ideal DC current source.
    
    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------
    
    properties (Constant)
        numTerms = 2;
        numCurrentVars = 0;
    end
    
    properties (SetAccess = private)
        A;      % Amplitude of the current
        % f must be of form [y,yp] = f(t)
        % here, y is the current value and yp is its time derivative
    end
    
    methods
        function obj = currentSourceDC(A, varargin)
            % varargin can be name
            obj = obj@device(varargin{:});
            if nargin > 0
                obj.A = A;
            end
        end 
        
        function s = J(thisCurrentSource, ~, ~)
            i = thisCurrentSource.A;
            s = [-i;
				  i];
        end
        
    end
    
    methods (Static)
        
        function s = Q(~, ~)
            % this function intentionally does nothing
            s = zeros(2,1);
        end
        
        function s = I(~, ~)
            % this function intentionally does nothing
            s = zeros(2,1);
        end
        
        function s = dQ(~, ~)
            % this function intentionally does nothing
            s = zeros(2);
        end
        
        function s = dI(~, ~)
            % this function intentionally does nothing
            s = zeros(2);
        end
        
    end
    
end

