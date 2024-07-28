classdef voltageSource < device
    % VOLTAGESOURCE Device class for an ideal voltage source.
    
    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------
    
    properties (Constant)
        numTerms = 2;
        numCurrentVars = 1;
    end
    
    properties (SetAccess = private)
        A;      % Amplitude of the voltage
        f;      % Handle to time dependent voltage function
        % f must be of form [y yp] = f(t)
        % here, y is the voltage value and yp is its time derivative
    end
    
    methods
        function obj = voltageSource(A, f, varargin)
            % varargin can be name
            obj = obj@device(varargin{:});
            if nargin > 0
                obj.A = A;
                obj.f = f;
            end
            obj.numVoltageSources = 1;
            obj.voltageSourceNodeNums = [1, 2];
            obj.voltageSources = obj;
        end 
        
        function s = J(thisVoltageSource, ~, t)
            s = zeros(3,1);
            s(3) = thisVoltageSource.getVoltage(t);
        end
        
        function [V Vp] = getVoltage(thisVoltageSource, t)
            try
                [V Vp] = thisVoltageSource.f(t);
                V = V*thisVoltageSource.A;
                Vp = Vp*thisVoltageSource.A;
            catch 
                V = thisVoltageSource.f(t)*thisVoltageSource.A;
                Vp = NaN;
            end
        end
    end
    
    methods (Static)
        
        function s = Q(~, ~)
            % this function intentionally does nothing
            s = zeros(3,1);
        end
        
        function s = I(x, ~)
            s = [x(3);
                -x(3);
                x(1) - x(2)];
        end
        
        function s = dQ(~, ~)
            % this function intentionally does nothing
            s = zeros(3);
        end
        
        function s = dI(~, ~)
            
            s = [0  0  1;
                 0  0 -1;
                 1 -1  0];      
        end
        
    end
    
end

