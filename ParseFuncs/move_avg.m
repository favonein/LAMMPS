function [avgOut] = move_avg(inArr,numAvg,centering)
%MOVE_AVG - Take a moving average of inArray using numAvg samples.
%      Centering doesn't properly handle beginning/end averaging. 
%   INPUTS:
%   inArr       - 1D array to be filtered
%   numAvg      - number of samples to be contained in moving AVG
%   centering   - If true, moving avg will be done with values +/- numAvg/2
%   OUTPUTS:
%   avgOut      - Filtered array with moving avg.

if (centering); offset = ceil(numAvg/2); else; offset = 0; end

movingAvgArr = inArr(1).*ones(numAvg,1);
avgOut = NaN(length(inArr),1);
idx = 2;
for i = 1:length(inArr)-offset


    avgOut(i) = mean(movingAvgArr);

    movingAvgArr(idx) = inArr(i+offset);

    idx = mod(idx,numAvg)+1;
end

for i = length(inArr)-offset+1:length(inArr)
    avgOut(i) = avgOut(i-1);
end
