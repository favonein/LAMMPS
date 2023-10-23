%Creates molecules of water based off of O/H distance
%Generates .in file -- the boundaries need manual editing, as well as any
%xy/xz/yz terms

clear variables;
fileName = "Data/H2O_weirdsquish.xyz";
fileNameOut = "Data/H2O_weirdsquish.in";


fid = fopen(fileName);
fid_out = fopen(fileNameOut,'W');

numAtoms = str2number(fgetl(fid));
Discard = fgetl(fid);

fprintf(fid_out,'Ice\n\n');

AtomList = cell(numAtoms,4);
numWater = numAtoms/3;
for i = 1:numAtoms
    AtomList(i,:) = transpose(split(strtrim(fgetl(fid))));
end

AtomTable = convertvars(cell2table(AtomList, "VariableNames",["Atom" "x" "y" "z"]),{'Atom'},'string');

%Offsets
xlo = 0;
ylo = 0;
zlo = 0;
%Not accurate usually if your cell is larger than the edge here.
xhi = max(str2double(AtomTable.x));
yhi = max(str2double(AtomTable.y));
zhi = max(str2double(AtomTable.z));

atomTypes = unique(AtomTable.Atom);
numAtomTypes = length(atomTypes);

idx = AtomTable.Atom == 'O';
OList = AtomTable(AtomTable.Atom == 'O',:);
HList = AtomTable(AtomTable.Atom == 'H',:);
Hidx = nonzeros((AtomTable.Atom == 'H').*transpose(1:1:numAtoms));
water = NaN(numWater,3);    %[O H1 H2]
water(:,1) = nonzeros(idx.*transpose(1:1:numAtoms)); %Get indexes of the O's
for i = 1:numWater
    O_x = str2number(cell2mat(OList.x(i))); %Got this str2number from matlab file exchange https://www.mathworks.com/matlabcentral/fileexchange/106265-str2number-fast-replacement-for-str2num-and-str2double
    O_y = str2number(cell2mat(OList.y(i))); %str2double works, but that command is horrendously slow
    O_z = str2number(cell2mat(OList.z(i)));
    O_xyz = [O_x;O_y;O_z]; 
    cnt = 0;
    for j = 1:2*numWater-2*(i-1)
        H_x = str2number(cell2mat(HList.x(j-cnt)));
        H_y = str2number(cell2mat(HList.y(j-cnt)));
        H_z = str2number(cell2mat(HList.z(j-cnt)));

        H_xyz = [H_x;H_y;H_z;];
        dist = norm(H_xyz-O_xyz);
        disp(strcat("H ", string(Hidx(j-cnt))," O ",string(water(i))," Dist: ",string(dist))); 
        if dist < 1.1   %Bond length ~1 angstrom
            water(i,cnt+2) = Hidx(j-cnt);

            Hidx(j-cnt) = [];       %Once assigned, remove from search list
            HList(j-cnt,:) = [];
            cnt = cnt+1;
        end
        if cnt > 1
            break
        end
    end
end


%% Print
fprintf(fid_out,'%u atoms\n',numAtoms);
fprintf(fid_out,'%u atom types\n',numAtomTypes);

fprintf(fid_out,['%u bonds\n' ...
                 '2 bond types\n'],numWater*2);

fprintf(fid_out,['%u angles\n' ...
                 '1 angle types\n\n'],numWater);


fprintf(fid_out,'%f %f xlo xhi\n', xlo, xhi+xlo);
fprintf(fid_out,'%f %f ylo yhi\n', ylo, yhi+ylo);
fprintf(fid_out,'%f %f zlo zhi\n\n', zlo, zhi+zlo);

fprintf(fid_out,'Atoms\n\n');

%generate charge style atom list
for i = 1:numAtoms
    fprintf(fid_out,'%u   %u    %u   %u   %f   %f   %f\n',i,mod(find(water == i),length(water)),find(contains(atomTypes,string(AtomTable.Atom(i)))),0,str2double(AtomTable.x(i))+xlo, ...
                str2double(AtomTable.y(i))+ylo,str2double(AtomTable.z(i))+zlo);
end

%generate Bond list
fprintf(fid_out,'\nBonds\n\n');
for i = 1:numWater
    fprintf(fid_out, '%u   1   %u   %u\n', i, water(i,1), water(i,2));
    fprintf(fid_out, '%u   1   %u   %u\n', i, water(i,1), water(i,3));
end

%generate Angle list
fprintf(fid_out,'\nAngles\n\n');
for i = 1:numWater
    fprintf(fid_out, '%u   1   %u   %u   %u\n', i, water(i,2), water(i,1), water(i,3));
end

fclose('all');