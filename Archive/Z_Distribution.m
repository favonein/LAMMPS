%% Setup
clear variables; 
fileName = "Data/settled.xyz";

%% Read Data, get Z

fid = fopen(fileName,"r");

numAtoms = cell2mat(textscan(fid,"%f",1));

discard = fgetl(fid); %skip next line
discard = fgetl(fid); %skip next line

atomList = textscan(fid, '%s%f%f%f',numAtoms(1));

atomType = atomList{1};
z = atomList{4};

idx = atomType == "O";
zList = z(idx);
idx = atomType == "C";
cList = z(idx);

[countO,edgesO] = histcounts(zList,20);
[countC,edgesC] = histcounts(cList,5);

midO = ones(length(edgesO)-1,1);

midC = ones(length(edgesC)-1,1);
for i = 1:length(edgesO)-1
    midO(i) = (edgesO(i) + edgesO(i+1))/2;
end
for i = 1:length(edgesC)-1
    midC(i) = (edgesC(i) + edgesC(i+1))/2;
end

%% Plot
figure(1);
hold on;
histogram(zList,10);

histogram(cList,3);
hold off;

figure(2);
center = sum(transpose(midC).*countC)/sum(countC);

hold on;
plot(midO-center,countO);
bar(midC-center,countC);
hold off;
fclose(fid);
