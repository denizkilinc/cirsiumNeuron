classdef currentSource < device
    % CURRENTSOURCE Device class for an ideal current source.
    
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
        f;      % Handle to time dependent current function
        % f must be of form [y,yp] = f(t)
        % here, y is the current value and yp is its time derivative
    end
    
    methods
        function obj = currentSource(A, f, varargin)
            % varargin can be name
            obj = obj@device(varargin{:});
            if nargin > 0
                obj.A = A;
                obj.f = f;
            end
        end 
        
        function s = J(thisCurrentSource, ~, t)
            I = thisCurrentSource.getCurrent(t);
            s = [-I;
				  I];
        end
        
        function [I, Ip] = getCurrent(thisCurrentSource, t)
            try
                [I, Ip] = thisCurrentSource.f(t);
                I = I*thisCurrentSource.A;
                Ip = Ip*thisCurrentSource.A;
            catch 
                I = thisCurrentSource.f(t)*thisCurrentSource.A;
                Ip = NaN;
            end
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

