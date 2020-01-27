% map range [a,b] to range [c,d]
function [x] = mapScale(x,a,b,c,d)
x = (x-a)*(d-c)/(b-a) + c;
end