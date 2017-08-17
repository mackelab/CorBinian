classdef ProbDist2D < handle
    %PROBDIST2D Represents a 2D probability distribution.
    
    properties (SetAccess = private);
        sup_x
        sup_y
        pdf_xy
    end
    
    methods
        function obj = ProbDist2D(x, y, pdf_xy)
            obj.sup_x = x(:);
            obj.sup_y = y(:);
            f = trapz(y, trapz(x, pdf_xy));
            obj.pdf_xy = pdf_xy / f;
        end
        
        function Ex = exp_x(obj)
            x = obj.sup_x;
            Ex = trapz(x, bsxfun(@times, x, obj.pdf_xy), 1);
        end
        
        function Ey = exp_y(obj)
            y = obj.sup_y';
            Ey = trapz(y, bsxfun(@times, y, obj.pdf_xy), 2);
        end
        
        function [Ex, Ey] = expectation(obj)
            Ex = obj.exp_x; Ex = Ex(:);
            Ey = obj.exp_y; Ey = Ey(:);
        end
        
        function COVxy = covariance(obj)
            x = obj.sup_x;
            y = obj.sup_y;
            
            tmp = (x * y') .* obj.pdf_xy;
            
            Exy = trapz(y, trapz(x, tmp, 1));
            [Ex, Ey] = obj.expectation();
            
            COVxy = Exy - Ex' * Ey;
        end
        
    end
    
end

