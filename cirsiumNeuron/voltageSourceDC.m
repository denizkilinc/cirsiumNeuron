classdef voltageSourceDC < voltageSource
    % VOLTAGESOURCEDC Device class for an ideal DC voltage source.
    
    %-------------------------------------------------------------------------------------
    % Copyright 2018 by Koc University and Deniz Kilinc, A. Gokcen Mahmutoglu, Alper Demir 
    % All Rights Reserved
    %-------------------------------------------------------------------------------------
    
    methods
        function obj = voltageSourceDC(A, varargin)
            % varargin can be name
            obj = obj@voltageSource(A, @voltageSourceDC.f, varargin{:});
        end 
    end
    
    methods (Static)
        function [y yp] = f(t)
            y = 1;
            yp = 0;
        end
    end
end

