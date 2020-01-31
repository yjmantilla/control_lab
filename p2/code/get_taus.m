function [taus] = get_taus(den)
ro = roots(den);
taus = [];
for i= 1:length(ro)
taus(i) = -1/ro(i);
end
end