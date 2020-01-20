function [stableTime,stableVal,stableIndex] = findStablePoint(time,signal,timeThresh)
    tol = 0.005*(max(signal) - min(signal));
    for i = 1:length(signal)
        for j = i:length(signal)
            delta = signal(j) - signal(i);
            if abs(delta) >= tol
                break
            end
            if time(j)-time(i)>= timeThresh
                stableVal = signal(j);
                stableTime = time(j);
                stableIndex = j;
                return
            end
        end
    end
    
    stableVal = false;
    stableTime = false;
    stableIndex = false;
    
end