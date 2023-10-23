function [Tables,time,numAtoms,xbounds,ybounds,zbounds] = parse_DumpCustom(fileName)
%PARSE_DUMPCUSTOM Parses a custom dump file of per atom quantities
%   INPUTS:
%       FileName - file location
%   OUTPUTS:
%       Tables - Cell array of Tables at each timestep of the dump
%       time - timestep of each dump
%       numAtoms - number of atoms in each dump
%       *lo,*hi - simulation box bounds

    fid = fopen(fileName,"r");
    
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
    col_types = strings(length(col_names),1);
    col_types(:) = "double" ;
    
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
    opts.VariableTypes = col_types;
    
    %Read Tables
    time = str2double(string(fullText(find(contains(fullText,'ITEM: TIMESTEP'))+1)));
    Tables = cell(length(tableStarts),1);
    
    clear fullText; clear discard1; clear discard2; clear discard3;
    
    %This can take quite long for long runs with frequent prints. Try not
    %to do that :/
    for i = 1:length(tableStarts)
        opts.DataLines = [tableStarts(i) tableStarts(i)+numAtoms-1];
        Tables{i} = readtable(fileName,opts);
    end

end

