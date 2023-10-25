function [time,Bins] = zDens(fileName,maxSteps)
%zDens - processes a zDensity file
%   Detailed explanation goes here


    fid = fopen(fileName);
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);

    time = NaN(maxSteps ,1); %allocate some space based on how timesteps & steps per print
    Bins = cell(maxSteps,4);
    for i = 1:maxSteps
        %nextLine = transpose(split(strtrim(fgetl(fid))));
        nextLine = cell2mat(textscan(fid, '%f%f%f',1));
        time(i) = nextLine(1);
        numBins = nextLine(2);

        Bins(i,:) = (textscan(fid, '%f%f%f%f',numBins));
    end
    
    


end

