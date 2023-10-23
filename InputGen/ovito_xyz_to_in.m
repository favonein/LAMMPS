fileName = "300k_ramp_2.xyz";
fileNameOut = "300k_ramp_2.in";


fid = fopen(fileName);
fid_out = fopen(fileNameOut,'W');

numAtoms = str2double(fgetl(fid));
Discard = fgetl(fid);

fprintf(fid_out,'Ice\n\n');

AtomList = cell(numAtoms,4);

for i = 1:numAtoms
    AtomList(i,:) = transpose(split(strtrim(fgetl(fid))));
end

AtomTable = cell2table(AtomList, "VariableNames",["Atom" "x" "y" "z"]);

xlo = 0;
ylo = 0;
zlo = 0;

xhi = max(str2double(AtomTable.x));
yhi = max(str2double(AtomTable.y));
zhi = max(str2double(AtomTable.z));
AtomTableRem = AtomTable;
%AtomTableRem = AtomTable((str2double(AtomTable.x) ~= xhi) & (str2double(AtomTable.y) ~= yhi),:);
%AtomTableRem = AtomTable((str2double(AtomTable.x) ~= xhi) & (str2double(AtomTable.y) ~= yhi) & (str2double(AtomTable.z) ~= zhi),:);
%AtomTableRem = AtomTable((str2double(AtomTable.x) ~= xhi) & (str2double(AtomTable.y) ~= yhi) & (str2double(AtomTable.z) > 0),:);
%AtomTableRem = AtomTableRem((str2double(AtomTableRem.x) >= 0) & (str2double(AtomTableRem.y) >= 0) & (str2double(AtomTableRem.z) >= 0),:);
numAtomsRem = height(AtomTableRem);

atomTypes = unique(AtomTableRem.Atom);
numAtomTypes = length(atomTypes);

fprintf(fid_out,'%u atoms\n\n',numAtomsRem);
fprintf(fid_out,'%u atom types\n\n',numAtomTypes);

fprintf(fid_out,'%f %f xlo xhi\n', xlo, xhi+xlo);
fprintf(fid_out,'%f %f ylo yhi\n', ylo, yhi+ylo);
fprintf(fid_out,'%f %f zlo zhi\n\n', zlo, zhi+zlo);

fprintf(fid_out,'Atoms\n\n');

for i = 1:numAtomsRem
    fprintf(fid_out,'%u   %u   %u   %f   %f   %f\n',i,find(contains(atomTypes,string(AtomTableRem.Atom(i)))),0,str2double(AtomTableRem.x(i)), ...
                str2double(AtomTableRem.y(i)),str2double(AtomTableRem.z(i)));
end



fclose('all');