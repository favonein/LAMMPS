function [Table,time] = parse_DumpCustomOne(fileName, timeStep, dumpData)
%PARSE_DUMPCUSTOM Parses a custom dump file of per atom quantities
%   INPUTS:
%       FileName - file location
%       timeStep - timestep of dump you want to print
%
%       dumpData - take straight from the scan_DumpCustom file
%           .numAtoms - number of Atoms 
%           .time - list of timesteps for each dump
%           .tableStarts - list of dump data starting locations
%           .opts - table import options
%   OUTPUTS:
%       Table - single table at dump print closest to timeStep
%       time - actual timeStep of dump

    %extract dumpData
    time = dumpData{1}.time;
    numAtoms = dumpData{1}.numAtoms;
    tableStarts = dumpData{1}.tableStarts;
    opts = dumpData{1}.opts;
    %extract dumpData

    %Round timeStep to nearest dump timestep
    if timeStep > time(length(time))
        i = length(time);
    elseif  timeStep < 0
        i = 1;
    else
        i = find(time == round(timeStep/time(2),0)*time(2));
    end

    time = time(i);
    
    %read table @ time
    opts.DataLines = [tableStarts(i) tableStarts(i)+numAtoms-1];
    Table = readtable(fileName,opts);

end

