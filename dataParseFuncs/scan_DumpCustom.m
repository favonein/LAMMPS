function [time,col_names,numAtoms,xbounds,ybounds,zbounds] = scan_DumpCustom(fileName)
%PARSE_DUMPCUSTOM Scans dump file and outputs information detailed below
%   INPUTS:
%       FileName - file location
%   OUTPUTS:
%       time - timestep of each dump
%       col_names - names of columns in the dump output
%       numAtoms - number of atoms in each dump
%       *lo,*hi - simulation box bounds, assumes const volume, output is a
%       single cell with 2x1 matrix inside each

    
    try
        fid = fopen(fileName,"r");
    catch
        error('Dumpfile not found!');
    end
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
    
    %Find table start locations
    fullText = regexp(fileread(fileName),'\n','split');
    
    %Read Tables
    time = str2double(string(fullText(find(contains(fullText,'ITEM: TIMESTEP'))+1)));

    
    fclose(fid);

end

