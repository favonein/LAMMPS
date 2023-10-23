function [Table,time] = parse_DumpCustomOne(fileName, timeStep)
%PARSE_DUMPCUSTOM Parses a custom dump file of per atom quantities
%   INPUTS:
%       FileName - file location
%       timeStep - timestep of dump you want to print
%   OUTPUTS:
%       Table - single table at dump print closest to timeStep
%       time - actual timeStep of dump

    fid = fopen(fileName,"r");
    
    %% Setup Table formats
    %Parse first input for data on num atoms, box bounds
    discard1 = fgetl(fid); %ITEM: TIMESTEP
    discard2 = fgetl(fid); %timestep #
    discard2 = fgetl(fid); %ITEM: NUMBER OF ATOMS
    numAtoms = textscan(fid,'%f',1);
    numAtoms = numAtoms{1};
    discard3 = fgetl(fid); %atom #
    discard3 = fgetl(fid); %ITEM: BOX BOUNDS pp pp ff
    
    xbounds = textscan(fid,'%f',2); %xlo xhi
    ybounds = textscan(fid,'%f',2); %ylo yhi
    zbounds = textscan(fid,'%f',2); %zlo zhi
    
    %Setup table format
    nextl = fgetl(fid);
    nextl = fgetl(fid);
    nextl_parsed = textscan(nextl,'%s');
    col_names = nextl_parsed{1};
    col_names = transpose(col_names(3:end));
    %col_types = strings(length(col_names),1);
    %col_types(:) = "double"; %Should be detected automatically well
    
    %Find table start locations
    fullText = regexp(fileread(fileName),'\n','split');
    tableStarts = find(contains(fullText,'ITEM: ATOMS'))+1;
    
    %Table import options
    opts = detectImportOptions(fileName,'FileType','text'); 
    opts.Delimiter = {' '};
    opts.LeadingDelimitersRule = 'ignore';
    opts.TrailingDelimitersRule = 'ignore';
    opts.ConsecutiveDelimitersRule = 'join';
    opts.VariableNames = col_names;
    %opts.VariableTypes = col_types;
    
    %Read Tables
    timeArray = str2double(string(fullText(find(contains(fullText,'ITEM: TIMESTEP'))+1)));

    %Round timeStep to nearest dump timestep
    if timeStep > timeArray(length(timeArray))
        i = length(timeArray);
    elseif  timeStep < 0
        i = 1;
    else
        i = find(timeArray == round(timeStep/timeArray(2),0)*timeArray(2));
    end

    time = timeArray(i);
    
    %read ith table
    opts.DataLines = [tableStarts(i) tableStarts(i)+numAtoms-1];
    Table = readtable(fileName,opts);

end

