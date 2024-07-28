classdef inductor < device
    % INDUCTOR Device class for a simple inductor.
    % Characteristic equation: v = L*di/dt.
    
    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------
    
    properties (Constant)
        numTerms = 2;
        numCurrentVars = 1;
    end
    
    properties (SetAccess = private)
        L = 0;
    end
    
    methods 
        function obj = inductor(L, varargin)
            obj = obj@device(varargin{:});
            if nargin > 0
                obj.L = L;
            end
        end
        
        function s = Q(thisInductor, x, ~)
            G = thisInductor.L;
            s = [ 0;
                  0;
                -x(3)*G];
        end
        
        function s = dQ(thisInductor, ~, ~)
            G = thisInductor.L;
            s = [0, 0,  0;
                 0, 0,  0;
                 0, 0, -G];
        end
    end
    
    methods (Static)
        function fi = devEq(i, L)
            fi = L*i;
        end
        
        function s = I(x, ~)
            s = [   x(3);
                   -x(3);
                 x(1)-x(2)];
        end
        
        function s = J(~, ~)
            s = zeros(3,1);
        end
        
        function s = dI(~, ~)
            s = [0,  0,  1;
                 0,  0, -1;
                 1, -1,  0];
        end
    end
    
end

