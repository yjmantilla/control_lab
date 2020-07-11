function [idx] = findNearest(n,a)
dist    = abs(a - n);
minDist = min(dist);
idx     = find(dist == minDist);
end