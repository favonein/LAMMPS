function [T, numAtoms] = parse_Out(fileName)
%PARSE_OUT parses a lammps output file for thermo output, assumes exact
%format as shown below in variable names.
%
%   T is the thermo data as a table, energy converted to eV/atom from
%   kcal/mol
%
%   Future work: parse the regexp output for table column names for
%   generalization
%   Take an input that specifies units

%Find where the thermo output starts and ends
fullText = regexp(fileread(fileName),'\n','split');
tableLine = find(contains(fullText,'KinEng'))+1;
tableLineEnd = find(contains(fullText,'Loop time of'))-1;

opts = detectImportOptions(fileName,'FileType','text'); 
opts.Delimiter = {' '};
opts.LeadingDelimitersRule = 'ignore';
opts.TrailingDelimitersRule = 'ignore';
opts.ConsecutiveDelimitersRule = 'join';
opts.DataLines = [tableLine tableLineEnd];
opts.VariableNames = {'Step','PotEng',...
                      'KinEng','TotEng','Temp', 'Press', 'Volume'};

opts.VariableTypes = {'double','double',...
                      'double','double','double','double','double'};

%Find num atoms in group "move"
numAtoms = cell2mat(textscan(string(fullText(contains(fullText,'atoms in group move'))),'%f'));

%Unit Adjustments
T = readtable(fileName,opts);
T.PotEng = T.PotEng.*0.0433634./numAtoms;
T.KinEng = T.KinEng.*0.0433634./numAtoms;
T.TotEng = T.TotEng.*0.0433634./numAtoms;
end

