% Reza Rashidi
% email: rashidyreza76@gmail.com
% SPP
clc;clear;close all;
% ==============================================================
% [rinexobs] = read_rinex_obs('bcov0010.21o');
load('matlab.mat')
% ==============================================================
algorithm_type = 1;
Filetype = 'navigation';% Filetype = 'sp3';
readnewfile = 'no';
usethisfile = rinexobs;
r0rcv = 'r0rinex';
[rinexobs] = T_emission_collins_klobuchar_XYZsat(algorithm_type,Filetype,readnewfile,usethisfile,r0rcv);
% ===============================================================
Find = find( rinexobs.data(:,rinexobs.col.C1)== 0 );
rinexobs.data(Find,:) =  [];
%====================cut of angle================================
Find = find( rinexobs.data(:,rinexobs.col.Elevation)< 15 );
rinexobs.data(Find,:) =  [];
%===========================================================
TOW = rinexobs.data(:,rinexobs.col.TOW); 
epoch = unique(TOW);
epoch = epoch(2870+1);
m = find(TOW ==epoch ); m = m(1)-1;
c = 299792458;
Code_obs = rinexobs.data(1:m,rinexobs.col.C1); 
Tr  = rinexobs.data(1:m,rinexobs.col.Trcollins);
Io  = rinexobs.data(1:m,rinexobs.col.I1_klobuchar);
TGD = rinexobs.data(1:m,rinexobs.col.TGD);
dtsat = rinexobs.data(1:m,rinexobs.col.dtsat);
Xsat = rinexobs.data(1:m,rinexobs.col.Xsat);
Ysat = rinexobs.data(1:m,rinexobs.col.Ysat);
Zsat = rinexobs.data(1:m,rinexobs.col.Zsat);
TOW = rinexobs.data(1:m,rinexobs.col.TOW); 
PRN = rinexobs.data(1:m,rinexobs.col.PRN);
TOW = rinexobs.data(1:m,rinexobs.col.TOW); 
delta_ro_rel = rinexobs.data(1:m,rinexobs.col.delta_ro_rel);
y = Code_obs-(-c*dtsat + Tr+Io+TGD);


rcv0 = rinexobs.r0;
cdtrcv = 1;
X0 = [rcv0;cdtrcv];
%===============================================================================
dxhat = 1;
while norm(dxhat) > 10^-5
dX = -Xsat+rcv0(1);
dY = -Ysat+rcv0(2);
dZ = -Zsat+rcv0(3);
y0 = sqrt( dX.^2 +dY.^2 + dZ.^2 ) +cdtrcv;
dy = y-y0;
A = [dX./y0 , dY./y0 , dZ./y0 , ones(m,1)];
dxhat = inv(A'*A)*A'*dy;
X0 = X0 + dxhat;
rcv0 = X0(1:3);cdtrcv = X0(4);
end
Qx = inv(A'*A);
rinexobs.r0 - rcv0