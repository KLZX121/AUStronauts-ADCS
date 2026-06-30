%% main attitude control script

clear;
close all;

format longG

if isempty(which('Q2Mat'))
    addpath(genpath('../SCT/SCTAcademic'));
end
