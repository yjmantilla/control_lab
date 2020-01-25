function [x] = mapScale(x,low_in,hi_in,low_out,hi_out)
x = (x-low_in)*(hi_out-low_out)/(hi_in-low_in) + low_out;
end