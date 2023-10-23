% rdf_plot.m
% Plot Radial Distribution function from LAMMPS output file
% Necessary inputs in setup section :)
% Favian Sun
% 9/13/2023
%% Setup
close all; clear variables;

fileName = "Logs/Z-density/Graphene/spc_atom.dens";
fileName2 = "Logs/Z-density/Graphene/C_atom.dens";

fspc = fopen(fileName);
fC = fopen(fileName2);
%From cont.in
timestep = 0.5/1000;    %timestep in cont.in (ps)
numSteps = 500000;      %numSteps in cont.in
perStep = 10000;          %how many steps between each rdf print

numPrint = numSteps/perStep;

%Zbounds
zhi = 44;
zlo = 0;

scaling = zhi-zlo;

%x-y bounds
xlo = 0;
xhi = 22.14;
ylo = 0;
yhi = 22.14;    %MAG

area = (xhi-xlo)*(yhi-ylo);

%First 3 lines are autogenerated comments
fgetl(fspc);
fgetl(fspc);
fgetl(fspc);
fgetl(fC);
fgetl(fC);
fgetl(fC);
%Filter strength
mAvg_num = 10;

plotTime = 39;

%% Loop, Read, and Basic data analysis

time = NaN(numPrint ,1); %allocate some space based on how timesteps & steps per print
BinsSpc = cell(numPrint,4);
BinsC = cell(numPrint,4);
for i = 1:numPrint
    %nextLine = transpose(split(strtrim(fgetl(fid))));
    nextLine = cell2mat(textscan(fspc, '%f%f%f',1));
    time(i) = nextLine(1);
    numBinsSpc = nextLine(2);

    nextLine = cell2mat(textscan(fC, '%f%f%f',1));
    numBinsC = nextLine(2);

    BinsSpc(i,:) = (textscan(fspc, '%f%f%f%f',numBinsSpc));
    BinsC(i,:) = (textscan(fC, '%f%f%f%f',numBinsC));
end

firstStep = cell2mat(BinsSpc(1,:));

dz = (firstStep(2,2)-firstStep(1,2))*scaling;

%Moving Avg -- not necessary if avgs are taken within the LAMMPS simulation
BinsPrintSpc = cell2mat(BinsSpc(plotTime,:));
binsFinSpc = move_avg(BinsPrintSpc(:,4),ceil(numBinsSpc/1000)+1,true);

BinsPrintC = cell2mat(BinsC(plotTime,:));
binsFinC = move_avg(BinsPrintC(:,4),ceil(numBinsC/1000)+1,true);


%% Plot
figure(1);
hold on;
plot(BinsPrintSpc(:,2).*scaling-5,BinsPrintSpc(:,4),'r-');
plot(BinsPrintC(:,2).*scaling-5,BinsPrintC(:,4),'b');
xlabel("Z Distance (Angstrom)");
ylabel("Unit Density");
title("Density (g/cm^3)");
legend(["Water" "C"]);
axis([-2 10 0 7]);
hold off;
figure(2);
hold on;
plot(BinsPrintSpc(:,2).*scaling-5,binsFinSpc,'r-');
plot(BinsPrintC(:,2).*scaling-5,binsFinC,'b');
xlabel("Z Distance (Angstrom)");
ylabel("Density (g/cm^3)");
title("Density of C and H2O");
legend(["Water" "C"]);
axis([-2 10 0 7]);
hold off;
fclose('all');