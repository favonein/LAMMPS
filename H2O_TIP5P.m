%Creates molecules of water based off of O/H distance
%Generates .in file -- the boundaries need manual editing, as well as any
%xy/xz/yz terms

clear variables;
fileName = "Data/H2O_3x2x3.xyz";
fileNameOut = "Data/H2O_3x2x3.in";


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
%Not accurate ... just replace.
xhi = max(str2double(AtomTable.x));
yhi = max(str2double(AtomTable.y));
zhi = max(str2double(AtomTable.z));

atomTypes = unique(AtomTable.Atom);
atomTypes = [atomTypes;'L'];
numAtomTypes = length(atomTypes);

idx = AtomTable.Atom == 'O';
OList = AtomTable(AtomTable.Atom == 'O',:);
HList = AtomTable(AtomTable.Atom == 'H',:);
Hidx = nonzeros((AtomTable.Atom == 'H').*transpose(1:1:numAtoms));
water = NaN(numWater,5);    %[O H1 H2 L1 L2]
water(:,1) = nonzeros(idx.*transpose(1:1:numAtoms)); %Get indexes of the O's
OH = NaN(3,2);
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
        %disp(strcat("H ", string(Hidx(j-cnt))," O ",string(water(i))," Dist: ",string(dist))); 
        if dist < 1.1   %Bond length ~1 angstrom
            water(i,cnt+2) = Hidx(j-cnt);

            %Scale bond length if not equal to 0.9572
            if dist ~= 0.9572

                scale = 0.9572/dist;
                %Write to xyz table
                AtomTable.x(Hidx(j-cnt)) = {sprintf('%.6f',O_x+(H_x-O_x)*scale)};
                AtomTable.y(Hidx(j-cnt)) = {sprintf('%.6f',O_y+(H_y-O_y)*scale)};
                AtomTable.z(Hidx(j-cnt)) = {sprintf('%.6f',O_z+(H_z-O_z)*scale)};

                %Edit H_xyz
                H_xyz = [O_x+(H_x-O_x)*scale;O_y+(H_y-O_y)*scale;O_z+(H_z-O_z)*scale];
            end
            Hidx(j-cnt) = [];       %Once assigned, remove from search list
            HList(j-cnt,:) = [];
            OH(:,cnt+1) = H_xyz-O_xyz;
            cnt = cnt+1;
        end
        if cnt > 1
            
            cosTheta = max(min(dot(OH(:,1),OH(:,2))/(norm(OH(:,1))*norm(OH(:,2))),1),-1);
            theta = real(acosd(cosTheta));
            axis_vec = transpose(cross(OH(:,1),OH(:,2)));
            %disp(strcat("theta: ",string(theta)," OH1 ",string(OH(1,1))," ",string(OH(2,1))," ",string(OH(3,1))," OH2 ",string(OH(1,2))," ",string(OH(2,2))," ",string(OH(3,2))));
            
            %Moves a hydrogen atom to match the 104.52 angle set for TIP5P
            if theta ~= 104.52
                %Generate rotation matrix using axes-angle form.
                diff = theta - 104.52;
                rot = axang2rotm([axis_vec deg2rad(diff)]);

                
                OH(:,1) = rot*OH(:,1);
                cosTheta_prime = dot(OH(:,1),OH(:,2))/(norm(OH(:,1))*norm(OH(:,2)));
                theta_prime = real(acosd(cosTheta_prime));
                
                %Changed only H1
                AtomTable.x(water(i,2)) = {sprintf('%.6f',OH(1,1)+O_x)};
                AtomTable.y(water(i,2)) = {sprintf('%.6f',OH(2,1)+O_y)};
                AtomTable.z(water(i,2)) = {sprintf('%.6f',OH(3,1)+O_z)};
                %disp(strcat("theta': ",string(theta_prime)," OH1 ",string(OH(1,1))," ",string(OH(2,1))," ",string(OH(3,1))," OH2 ",string(OH(1,2))," ",string(OH(2,2))," ",string(OH(3,2))));
            end

            OH_mid = transpose((OH(:,1) + OH(:,2))./2);
            %LOL plane is 90 degrees rotated from HOH plane
            Lplane_rot = axang2rotm([OH_mid pi/2]);
            Lplane_axis = transpose(Lplane_rot*transpose(axis_vec));
            %Angle of LOL = 109.47, so we will rotate 180-109.47/2 = 125.265 deg away from the OH_mid vector
            L_rot1 = axang2rotm([Lplane_axis 125.265*pi/180]);
            L_rot2 = axang2rotm([Lplane_axis 234.735*pi/180]);
            %O-L distance is 0.7 angstroms
            L1 = L_rot1*transpose(OH_mid)./norm(OH_mid).*0.7+O_xyz;
            L2 = L_rot2*transpose(OH_mid)./norm(OH_mid).*0.7+O_xyz;
            %Add to AtomTable
            newL = {'L',sprintf('%.6f',L1(1)),sprintf('%.6f',L1(2)),sprintf('%.6f',L1(3));
                    'L',sprintf('%.6f',L2(1)),sprintf('%.6f',L2(2)),sprintf('%.6f',L2(3))};
            AtomTable = [AtomTable;newL];
            numAtoms = numAtoms + 2;
            %add to water
            water(i,4) = height(AtomTable)-1;
            water(i,5) = height(AtomTable);
            break
        end
    end
end


%% Print
fprintf(fid_out,'%u atoms\n',numAtoms);
fprintf(fid_out,'%u atom types\n',numAtomTypes);

fprintf(fid_out,'%f %f xlo xhi\n', xlo, xhi+xlo);
fprintf(fid_out,'%f %f ylo yhi\n', ylo, yhi+ylo);
fprintf(fid_out,'%f %f zlo zhi\n\n', zlo, zhi+zlo);

fprintf(fid_out,'Atoms\n\n');

%generate charge style atom list
for i = 1:numAtoms
    if mod(find(water == i),length(water)) == 0
        waterNum = length(water);
    else
        waterNum = mod(find(water == i),length(water));
    end
    fprintf(fid_out,'%u   %u    %u   %u   %f   %f   %f\n',i,waterNum,find(contains(atomTypes,string(AtomTable.Atom(i)))),0,str2double(AtomTable.x(i))+xlo, ...
                str2double(AtomTable.y(i))+ylo,str2double(AtomTable.z(i))+zlo);
end

fclose('all');