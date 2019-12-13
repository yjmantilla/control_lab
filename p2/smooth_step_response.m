% Loading of the data
path = 'Y:\y\sync\current\courses\control\lab\p2\';
filename = 'to_matlab.csv';
data = csvread([path filename]);

% Partitioning
time = data(:,1)';  % time vector
p = data(:,2)'; % process value, in this case pressure
% set value is column 3, not needed since there is no control being done
vop = data(:,4)'; % valve opening percentage

% Normalization
pv_1 = p/norm(p);
vop_1 = vop/norm(vop);

% Smoothing
span = 0.1 * length(time); % span of 10% of the total number of data points
pv_s = smooth(time,pv,span,'moving'); % smoothing
pv_1s = pv_s/norm(pv_s);

% Plotting
figure(1)
plot(time,vop_1,'r-')
hold on
plot(time,pv_1,'g.')
hold on 
plot(time,pv_1s,'b-')
