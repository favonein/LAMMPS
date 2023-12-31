function [time,Bins,rdf_ave_oh, rdf_ave_hh, rdf_ave_oo] = h2o_rdf(rdf,maxSteps,bins)
%H2O_RDF Parses rdf file of water and returns full data, and avgs
%   INPUTS:
%       rdf = file path of rdf
%       maxSteps = number of rdf prints expected
%       bins = number of bins in each rdf print
%       cutOH = bin number where O-H bond distance lies
%   OUTPUTS:
%       time = array of timesteps of the rdf prints
%       Bins = full data of each rdf print at each timestep
%       rdf_ave_xx = array of avg of respective radial distance

fid = fopen(rdf);

if fid == -1
 error('Could not read file: %s', rdf);
end

%First 3 lines are autogenerated comments
fgetl(fid);
fgetl(fid);
fgetl(fid);
%allocate some space based on how timesteps & steps per print
%
time = NaN(maxSteps ,1); 

rdf_ave_oo = NaN(maxSteps ,1);

rdf_ave_hh = NaN(maxSteps ,1);

rdf_ave_oh = NaN(maxSteps ,1);

Bins = NaN(bins,8,maxSteps);

for i = 1:maxSteps
    nextLine = cell2mat(textscan(fid, '%f%f',1));
    time(i) = nextLine(1);
    numBins = nextLine(2);

    %textscan is fast for reading many lines of same formatting. Hard coded
    %formatting for this
    BinsTemp = cell2mat(textscan(fid, '%f%f%f%f%f%f%f%f',numBins));
    Bins(:,:,i) = BinsTemp;

    BinsCutOH = BinsTemp((BinsTemp(:,4)) > 2,:);  %Cutoff set to when the coordination number increases above 2. 2 O-H per molec
    rdf_ave_oh(i) = sum((BinsCutOH(:,3)).*(BinsCutOH(:,2)))/sum((BinsCutOH(:,3)));
    BinsCutHH = BinsTemp((BinsTemp(:,6)) > 1.01,:);  %Cutoff set to when the coordination number increases above 1. 1 H/H interaction per molec
    rdf_ave_hh(i) = sum((BinsCutHH(:,5)).*(BinsCutHH(:,2)))/sum((BinsCutHH(:,5)));
    rdf_ave_oo(i) = sum((BinsTemp(:,7)).*(BinsTemp(:,2)))/sum((BinsTemp(:,7)));

    idxOHcut = find(diff(BinsTemp(:,4))==2,1);
    Bins(idxOHcut+1,3,i) = 0;
    idxHHcut = find(diff(BinsTemp(:,5))>=10,1);
    Bins(idxHHcut+1,5,i) = 0;
end

rdf_ave_oh = rdf_ave_oh(~isnan(rdf_ave_oh),:);

end

