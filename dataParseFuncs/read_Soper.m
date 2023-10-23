function [HH,OH,OO] = read_Soper(fileNameHH,fileNameOH,fileNameOO)
%READ_SOPER Summary of this function goes here
%   Detailed explanation goes here

    HH = readtable(fileNameHH);
    OH = readtable(fileNameOH);
    OO = readtable(fileNameOO);
end

