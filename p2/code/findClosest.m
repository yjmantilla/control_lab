function [closest,ix] = findClosest(x,val)
[ d, ix ]  = min( abs( x-val ) );
closest = x(ix);
end
