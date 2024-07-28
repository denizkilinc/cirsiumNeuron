classdef voltageSourceAC < voltageSource
    % VOLTAGESOURCEAC Device class for an ideal AC voltage source.
    
    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------
    
    methods
        function obj = voltageSourceAC(A, T, varargin)
            % varargin can be name
            obj = obj@voltageSource(A, @(t)(voltageSourceAC.f(t, T)), varargin{:});
        end 
    end
    
    methods (Static)
        function [y yp] = f(t, T)
            y = sin(2*pi/T*t);
            yp = cos(2*pi/T*t);
        end
    end
end