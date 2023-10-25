function [dumpData] = scan_DumpCustom(fileName)
%SCAN_DUMPCUSTOM Scans dump file and outputs information detailed below
%   INPUTS:
%       FileName - file location
%   OUTPUTS:
%       dumpData -- for use in parse funcs
%           .numAtoms - number of Atoms 
%           .time - list of timesteps for each dump
%           .tableStarts - list of dump data starting locations
%           .opts - table import options
%           .col_names - names of columns in the dump output
%           .*lo,*hi - simulation box bounds, assumes const volume, output is a
%           .single cell with 2x1 matrix inside each
%           .element_list - Array of elements in order of ID at start of sim

    
    try
        fid = fopen(fileName,"r");
    catch
        error('Dumpfile not found!');
    end
    
    %Parse first input for data on num atoms, box bounds
    fgetl(fid); %ITEM: TIMESTEP
    fgetl(fid); %timestep #
    fgetl(fid); %ITEM: NUMBER OF ATOMS
    numAtoms = textscan(fid,'%f',1);
    numAtoms = numAtoms{1};
    fgetl(fid); %atom #
    fgetl(fid); %ITEM: BOX BOUNDS pp pp ff
    
    xbounds = textscan(fid,'%f',2); %xlo xhi
    ybounds = textscan(fid,'%f',2); %ylo yhi
    zbounds = textscan(fid,'%f',2); %zlo zhi
    
    %Setup table format
    fgetl(fid);
    nextl = fgetl(fid);
    nextl_parsed = textscan(nextl,'%s');
    col_names = nextl_parsed{1};
    col_names = transpose(col_names(3:end));
    
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

    %Get timestep list
    time = str2double(string(fullText(find(contains(fullText,'ITEM: TIMESTEP'))+1)));


    %varargout
    dumpData{1}.time = time;
    dumpData{1}.numAtoms = numAtoms;
    dumpData{1}.tableStarts = tableStarts;
    dumpData{1}.opts = opts;
    dumpData{1}.xbounds = xbounds{1};
    dumpData{1}.ybounds = ybounds{1};
    dumpData{1}.zbounds = zbounds{1};
    dumpData{1}.col_names = col_names;
    %varargout


    %Scan first table for element list
    [Table] = parse_DumpCustomOne(fileName, time(1), dumpData);

    Table = sortrows(Table,'id');

    element_list = string(Table.element);
    
    %Add element list from the first table to dump data
    dumpData{1}.element_list = element_list;

    fclose(fid);

end

