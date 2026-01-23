tic

close all
clear
clc

% addpath('/home/bieito/Documents/SCIENCE/EPFL/instruments/MicroCTD/odas_matlab_v4.3.07/odas/')
% addpath(genpath('/home/bieito/Documents/SCIENCE/EPFL/instruments/MicroCTD/myfunctions/'))
% addpath(genpath('..\..\MicroCTD\'))
addpath('../odas_v4.4/')
addpath('C:\Users\tdoda\OneDrive - Université de Lausanne\External drive\PhD\Field Work\Instruments\Microstructure profilers\Profiler MicroCTD\Scripts\Bieito_private_scripts\Functions_2020\myfunctions\')
datafolder='../../data/LakeGeneva/20190826/Level0/';



%%%%
% file0 = 'VMP_005';
% file0 = 'VMP_009';
%file0 = 'VMP_012';
%file0 = 'VMP_017';
% file0 = 'VMP_018';
%file0 = 'VMP_020';
file0 = 'VMP_023';

DATA=odas_p2mat([datafolder,file0,'.P']);
info.minvel_detect = 0.25;
info.mindur_detect = 20;
info.pmin = 0;
info.pmax = 30; %%%CHANGE: maximum depth
info.dp = 0.5;
info.dpD = 1;
info.prof_dir = 'up'; %%CHANGE: up or down
info.fmaxT = 110;
info.Tmethod = 'B'; % Bieto's method (with Tspec=Kraichnan)
info.time_res = 0.;
info.system = 'Lem';

%%from 089 on T1228 was installed in T1 and a peak on spectra has to be
%%removed
%info.peak_rem_T1 = [40,52]; %for T1228 (downward)
%info.peak_rem_T1 = [45,60]; %for T1228 (upward)


PLOT = 0;

iPs0 = get_profile(DATA.P_slow,DATA.W_slow,info.pmin,info.minvel_detect,info.prof_dir,info.mindur_detect,DATA.fs_slow);
Nprf = length(iPs0);
%inPall = 6:9

inPall = 1:Nprf; %VMP012, VMP017, VMP024
%inPall = [1,2,4,5,6]; %VMP005
%inPall = [1,2,3,4]; %VMP009
%inPall = [5,7,8,9,10,11,12]; %VMP018 
%inPall = [1,2,3,5,6];%VMP020

i=1;
for inp = inPall
    fprintf('Profile number %d of %d, ',inp, Nprf)
    tic
    maxP=max(DATA.P_slow(iPs0(1,inp):iPs0(2,inp)),[],"omitnan");
    if maxP>30
        warning('Not upward profile')
        continue
    end
    [BINNED{i},SLOW{i},FAST{i}] = resolve_microCTD_profile(DATA,inp,info,file0,PLOT);
    toc
    
    i= i+1;
    
end

save([datafolder,'../Level2/microCTD_Bieito_',file0,'_',info.prof_dir],'BINNED','SLOW','FAST')