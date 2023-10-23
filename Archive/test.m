%% Setup
clear variables; close all;
fileName = "Logs/spc_cs_2v.out";

%numAtoms = 183; %183 atoms not including the bottom layer of Pt ... 75 Pt + 36 H2O
numAtoms = 291; %291 with 2 layers of water

timestep = 0.5/1000;    %timestep in cont.in (ps)
numSteps = 250000;      %numSteps in cont.in
perStep = 500;          %how many steps between each thermo dump

mAvg_num = 50;

moving_ave_tot = NaN(numSteps/perStep +1,1);
moving_ave_pot = NaN(numSteps/perStep +1,1);
%moving_ave_tot = NaN(1655,1);
%moving_ave_pot = NaN(1655,1);


%% Read Data, Moving Average
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

T = readtable(fileName,opts);

[Steps,idx] = sort(T.Step);
T = T(idx,:);

Time = timestep*T.Step;

movingAvgList_tot = T.TotEng(1).*ones(mAvg_num,1); %moving average of mAvg_num samples
movingAvgList_pot = T.PotEng(1).*ones(mAvg_num,1); %moving average of mAvg_num samples
idx = 2;
for i = 1:length(moving_ave_tot)
    moving_ave_tot(i) = mean(movingAvgList_tot);
    moving_ave_pot(i) = mean(movingAvgList_pot);
    movingAvgList_tot(idx) = T.TotEng(i);
    movingAvgList_pot(idx) = T.PotEng(i);
    idx = mod(idx,mAvg_num)+1;
end


%% Plot
figure(1)
plot(Time,moving_ave_tot.*0.0433634./numAtoms,Time,moving_ave_pot.*0.0433634./numAtoms);
title("Energy over time");
xlabel("Time (ps)");
ylabel("Energy (eV/atom)");
legend(["Total Energy" "Potential Energy"]);

%% Plot Diff
% figure(2)
% diff_tot = diff(moving_ave_tot)./diff(Time);
% diff_tot = diff_tot(10/(perStep*timestep):end);
% diff_pot = diff(moving_ave_pot)./diff(Time);
% diff_pot = diff_pot(10/(perStep*timestep):end);
% plot(Time(10/(perStep*timestep)+1:end),diff_tot.*0.0433634./numAtoms,Time(10/(perStep*timestep)+1:end),diff_pot.*0.0433634./numAtoms);
% title("Energy over time");
% xlabel("Time (ps)");
% ylabel("dE/dt (eV/atom*s)");
%legend(["Total Energy" "Potential Energy"]);

