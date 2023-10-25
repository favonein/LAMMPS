function [Tables] = parse_DumpCustom(fileName, dumpData)
%PARSE_DUMPCUSTOM Parses a custom dump file of per atom quantities
%                   If you don't need all the data/know what you are
%                   looking for, use parse_DumpCustomOne on loop.
%   INPUTS:
%       FileName - file location
%
%       dumpData - take straight from the scan_DumpCustom file
%           .numAtoms - number of Atoms 
%           .tableStarts - list of dump data starting locations
%           .opts - table import options
%   OUTPUTS:
%       Tables - Cell array of Tables at each timestep of the dump


    opts = dumpData{1}.opts;
    tableStarts = dumpData{1}.tableStarts;
    numAtoms = dumpData{1}.numAtoms;

    Tables = cell(length(tableStarts),1);
    
    %This can take quite long for long runs with frequent prints.
    for i = 1:length(tableStarts)
        opts.DataLines = [tableStarts(i) tableStarts(i)+numAtoms-1];
        Tables{i} = readtable(fileName,opts);
    end

end

