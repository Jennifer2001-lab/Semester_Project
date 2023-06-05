% clear all
% close all
% Ragu 

% - Quadratic constraint 
% con [Bp^2 + Bq^2 <= Converter_capacity^2] is replaced by 


%  [slope_B, constant_B] = Linearize_quadratic(Converter_capacity, no_linear);
% for sl = 1:numel(slope_B)/2
% 
%     con = [con Bp <= slope_B(sl)*Bq + constant_B(sl)];
%     con = [con Bp >= slope_B(sl)*Bq - constant_B(sl)];
% end

function [slope, constant] = Linearize_quadratic(Converter_capacity, no_linear)
    r = Converter_capacity;
    th = 0:pi/no_linear:2*pi;
    xunit = r * cos(th) + 0;
    yunit = r * sin(th) + 0;

    for j = 1:numel(yunit)-1
        slope(j) = (yunit(j+1)-yunit(j))/(xunit(j+1)-xunit(j));
        constant(j) = (yunit(j)*(xunit(j+1)-xunit(j))- xunit(j)*(yunit(j+1)-yunit(j)))/(xunit(j+1)-xunit(j));
    end

%     x = -r:0.05*r:r;
%     % x = -r_value:0.05*r_value:r_value;
% 
%     % figure
%     circle(0,0,r)
%     hold on
%     for jj = 1:numel(slope)/2
% 
%         yy = slope(jj)*x + constant(jj);
%         zz = slope(jj)*x - constant(jj);
%         plot(x,yy,'b--')
%         hold on
%         plot(x,zz,'b--')
%     end
%     axis([-Converter_capacity, Converter_capacity, -Converter_capacity, Converter_capacity])

function h = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit);
hold off
end

end
