% rdf_plot.m
% Plot Radial Distribution function from LAMMPS output file
% Necessary inputs in setup section :)
% Favian Sun
% 9/13/2023
%% Setup
close all; clear variables;

fileName = "Logs/Z-density/tip3p_atom.dens";
fileName2 = "Logs/Z-density/Pt_atom.dens";
fileName3 = "Logs/Z-density/spc_1.dens";
fileName4 = "Logs/Z-density/reax_atom.dens";
fileName5 = "Logs/Z-density/spc.dens";

ftip3p = fopen(fileName);
fpt = fopen(fileName2);
fspc = fopen(fileName3);
freax = fopen(fileName4);
fspc2 = fopen(fileName5);

%From cont.in
timestep = 0.5/1000;    %timestep in cont.in (ps)
numSteps = 500000;      %numSteps in cont.in
perStep = 10000;          %how many steps between each rdf print

numPrint = numSteps/perStep;

%Zbounds
zhi = 50;
zlo = -0.5;

scaling = zhi-zlo;

scaling_Reax = 22.9+0.5;

%x-y bounds
xlo = -0.5;
xhi = 13.44115;
ylo = -0.5;
yhi = 11.57339;

area = (xhi-xlo)^2; %parallelogram area

%First 3 lines are autogenerated comments
fgetl(ftip3p);
fgetl(ftip3p);
fgetl(ftip3p);
fgetl(fspc);
fgetl(fspc);
fgetl(fspc);
fgetl(freax);
fgetl(freax);
fgetl(freax);
fgetl(fspc2);
fgetl(fspc2);
fgetl(fspc2);
fgetl(fpt);
fgetl(fpt);
fgetl(fpt);
%Filter strength
mAvg_num = 10;

plotTime = 40;

%% Loop, Read, and Basic data analysis

time = NaN(numPrint ,1); %allocate some space based on how timesteps & steps per print
BinsTip3p = cell(numPrint,4);
BinsSpc = cell(numPrint,4);
BinsSpc2 = cell(numPrint,4);
BinsReax = cell(numPrint,4);
BinsPt = cell(numPrint,4);
for i = 1:numPrint
    nextLine = cell2mat(textscan(ftip3p, '%f%f%f',1));
    time(i) = nextLine(1);
    numBinsTip3p = nextLine(2);
    nextLine = cell2mat(textscan(fpt, '%f%f%f',1));
    numBinsPt = nextLine(2);
    nextLine = cell2mat(textscan(fspc, '%f%f%f',1));
    numBinsSpc = nextLine(2);
    nextLine = cell2mat(textscan(fspc2, '%f%f%f',1));
    numBinsSpc2 = nextLine(2);
    nextLine = cell2mat(textscan(freax, '%f%f%f',1));
    numBinsReax = nextLine(2);

    BinsTip3p(i,:) = (textscan(ftip3p, '%f%f%f%f',numBinsTip3p));
    BinsPt(i,:) = (textscan(fpt, '%f%f%f%f',numBinsPt));
    BinsSpc(i,:) = (textscan(fspc, '%f%f%f%f',numBinsSpc));
    BinsSpc2(i,:) = (textscan(fspc2, '%f%f%f%f',numBinsSpc2));
    BinsReax(i,:) = (textscan(freax, '%f%f%f%f',numBinsReax));

end

%Moving Avgs - Dont need if time LAMMPS Avg
BinsPrintTip3p = cell2mat(BinsTip3p(plotTime,:));
binsFinTip3p = move_avg(BinsPrintTip3p(:,4),ceil(numBinsTip3p/100)+1,true);


BinsPrintSpc = cell2mat(BinsSpc(plotTime,:));
binsFinSpc = move_avg(BinsPrintSpc(:,4),ceil(numBinsSpc/100)+1,true);

BinsPrintSpc2 = cell2mat(BinsSpc2(plotTime,:));
binsFinSpc2 = move_avg(BinsPrintSpc2(:,4),ceil(numBinsSpc2/100)+1,true);

BinsPrintReax = cell2mat(BinsReax(plotTime,:));
binsFinReax = move_avg(BinsPrintReax(:,4),ceil(numBinsReax/100)+1,true);

BinsPrintPt = cell2mat(BinsPt(plotTime,:));
binsFinPt = move_avg(BinsPrintPt(:,4),ceil(numBinsPt/100)+1,true);

%% Plot



figure(1);
hold on;
plot(BinsPrintTip3p(:,2).*scaling-7.348,BinsPrintTip3p(:,4),'g-');
plot(BinsPrintSpc(:,2).*scaling-7.348,BinsPrintSpc(:,4),'r-');
plot(BinsPrintSpc2(:,2).*scaling-7.348,BinsPrintSpc2(:,4),'k-.');
plot(BinsPrintReax(:,2).*scaling_Reax-7.348,BinsPrintReax(:,4),'m-');
plot(BinsPrintPt(:,2).*scaling-7.348,BinsPrintPt(:,4),'b');
xlabel("Z Distance (Angstrom)");
ylabel("Density (g/cm^3)");
title("Density of Pt and H2O");
legend(["Tip3p" "Spc" "Spc weaker mixing" "ReaxFF" "Pt"]);
axis([-2 15 0 10]);
hold off;
figure(2);
hold on;
plot(BinsPrintTip3p(:,2).*scaling-7.348,binsFinTip3p,'g-');
plot(BinsPrintSpc(:,2).*scaling-7.348,binsFinSpc,'r-');
plot(BinsPrintSpc2(:,2).*scaling-7.348,binsFinSpc2,'k-.');
plot(BinsPrintReax(:,2).*scaling_Reax-7.348,binsFinReax,'m-');
plot(BinsPrintPt(:,2).*scaling-7.348,binsFinPt,'b');
xlabel("Z Distance (Angstrom)");
ylabel("Unit Density (g/cm^3)");
title("Density of Pt and H2O, filtered");
legend(["Tip3p" "Spc" "Spc weaker mixing" "ReaxFF" "Pt"]);
axis([-2 15 0 10]);
hold off;
fclose('all');