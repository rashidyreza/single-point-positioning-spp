% Reza Rashidi
% email: rashidyreza76@gmail.com

%{
inputs
                    rinex observation File
                    NAVIGATION FILE or sp3 FILE
                    Filetype: 'sp3' or 'navigation'
                    algorithm_type for calculate T_emission
                           1: for A Pseudorange-Based Algorithm
                           2: for A Purely Geometric Algorithm
                    r0rcv for second Algorithm
                          r0rcv = [X;Y;Z]
                           or
                 r0rcv = 'r0rinex' to use coordinates in rinexobs header;
                 readnewfile = 'yes' or 'no'
                 if no I will use old file of rinex observation File
%}

%{
outinputs
                 Time of emission
                 dt
                 [X,Y,Z] sat in Time of emission
                 dtsat

%}
%-------------------------read files------------------------------
function [rinexobs,satnav] = T_emission_collins_klobuchar_XYZsat(algorithm_type,Filetype,readnewfile,usethisfile,r0rcv)
readnewfile = sum(double(readnewfile));
switch Filetype
    case 'sp3'
        sp3fname = uigetfile('*.sp3*');     % sp3 FILE
        [SP3_XYZ] = ReadSP3(sp3fname);
        Filename = uigetfile('*.21n*');     % NAVIGATION FILE
        [navvv,Dateee,alfa,beta]=ReadNavigation(Filename,38);
    case 'navigation'
        Filename = uigetfile('*.21n*');     % NAVIGATION FILE
        [nav,Date,alfa,beta]=ReadNavigation(Filename,38);
end
if readnewfile == sum(double('yes'))
    filename = uigetfile('*.21o*');             % OBSERVATION FILE
    [rinexobs] = read_rinex_obs(filename);
    rec_xyz = rinexobs.r0;
else
    rinexobs = usethisfile;
    rec_xyz = rinexobs.r0;
end
%============================================================
if r0rcv== 'r0rinex'
    r0rcv_xyz = rec_xyz;
else
    r0rcv_xyz = r0rcv;
end
% GNSS Type G for GPS
Find = find( rinexobs.data(:,rinexobs.col.GNSStype)== double('G') );
rinexobs.data =  rinexobs.data(Find,:);
global GM we c
GM = 3.986004418e14;
we = 7.2921151467e-5;
c = 2.99792458e8;
dt=0;
i=4;
while find(dt==0)>0
    dt = rinexobs.data(:,i)/c;
    Find = find( dt==0 );
    i=i+1;
    if i==10
        rinexobs.data(Find,:)=[];
        dt(Find,:)=[];
        break
    end
end
switch algorithm_type
    % A Pseudorange-Based Algorithm
    case 1
        switch Filetype
            case 'sp3'
                PRN = unique(SP3_XYZ.data(:,SP3_XYZ.col.PRN))';
                for i = PRN
                    dtsat = 0;
                    while 1
                        Find = find( SP3_XYZ.data(:,SP3_XYZ.col.PRN) == i );
                        pointY = SP3_XYZ.data(Find,SP3_XYZ.col.dtsat)*1e-6;
                        pointX = SP3_XYZ.data(Find,SP3_XYZ.col.gps_seconds);
                        Find = find( rinexobs.data(:,rinexobs.col.PRN) == i );
                        Temission = rinexobs.data(Find,rinexobs.col.TOW) - dt(Find) - dtsat;
                        dtsat0 = lagrange_interpolation(Temission,pointX,pointY,10);
                        Norm = norm ( dtsat - dtsat0);
                        dtsat = dtsat0;
                        if Norm<1e-10
                            rinexobs.data(Find,rinexobs.col.dtsat) = dtsat;
                            rinexobs.data(Find,rinexobs.col.dt_sat_rcv) = dt(Find);
                            % Temission
                            Temission = rinexobs.data(Find,rinexobs.col.TOW) - dt(Find) - dtsat;
                            rinexobs.data(Find,rinexobs.col.Temission) = Temission;
                            break;
                        end
                    end
                    Findsp3 = find( SP3_XYZ.data(:,SP3_XYZ.col.PRN) == i );
                    Findobs = find( rinexobs.data(:,rinexobs.col.PRN) == i );
                    pointX = SP3_XYZ.data(Findsp3,SP3_XYZ.col.gps_seconds);
                    pointYx = SP3_XYZ.data(Findsp3,SP3_XYZ.col.X);
                    pointYy = SP3_XYZ.data(Findsp3,SP3_XYZ.col.Y);
                    pointYz = SP3_XYZ.data(Findsp3,SP3_XYZ.col.Z);
                    XX = lagrange_interpolation(Temission,pointX,pointYx,10);
                    YY = lagrange_interpolation(Temission,pointX,pointYy,10);
                    ZZ = lagrange_interpolation(Temission,pointX,pointYz,10);
                    rsat = [XX , YY , ZZ ];
                    rinexobs.data(Findobs,rinexobs.col.Xsat:rinexobs.col.Zsat) = rsat;
                end
            case 'navigation'
                PRN = unique(nav.data(:,nav.col.PRN))';
                Prn = nav.data(:,nav.col.PRN);
                for i = PRN
                    Find = find(Prn == i); Find = Find(1);
                    TGD = nav.data(Find,nav.col.TGD)*299792458;
                    rinexobs.data( find(rinexobs.data(:,rinexobs.col.PRN)==i),rinexobs.col.TGD) = TGD;
                end
                PRN = unique(nav.data(:,nav.col.PRN))';
                for i = PRN
                    dtsat = 0;
                    while 1
                        Findnav = find( nav.data(:,nav.col.PRN) == i );
                        Findobs = find( rinexobs.data(:,rinexobs.col.PRN) == i );
                        Temission = rinexobs.data(Findobs,rinexobs.col.TOW) - dt(Findobs) - dtsat;
                        TOE = nav.data(Findnav,nav.col.TOE);
                        Nearest = [];
                        for j=1:length(Temission)
                            MIN = abs(TOE-Temission(j));
                            nnn = find(MIN==min(MIN));
                            Nearest(j,1) = Findnav(nnn);
                        end
                        a0 = nav.data(Nearest,nav.col.SV_Clock_Bias);
                        a1 = nav.data(Nearest,nav.col.SV_Clock_Drift);
                        a2 = nav.data(Nearest,nav.col.SV_Clock_Drift_Rate);
                        toe = nav.data(Nearest,nav.col.TOE);
                        [AAAA,drel] =  NAV2ECEF(nav,rinexobs,i);
                        dtsat0 = a0 + a1.*(Temission - toe) + a2.*(Temission - toe).^2 + drel;
                        Norm = norm ( dtsat - dtsat0);
                        dtsat = dtsat0;
                        if Norm<1e-10
                            rinexobs.data(Findobs,rinexobs.col.dtsat) = dtsat;
                            rinexobs.data(Findobs,rinexobs.col.dt_sat_rcv) = dt(Findobs);
                            % Temission
                            Temission = rinexobs.data(Findobs,rinexobs.col.TOW) - dt(Findobs) - dtsat;
                            rinexobs.data(Findobs,rinexobs.col.Temission) = Temission;
                            break;
                        end
                    end
                    [rinexobs] = NAV2ECEF(nav,rinexobs,i);
                end
        end
    case 2
        switch Filetype
            case 'sp3'
                rinexobs.data(:,rinexobs.col.Temission) = rinexobs.data(:,rinexobs.col.TOW);
                PRN = unique(SP3_XYZ.data(:,SP3_XYZ.col.PRN))';
                for i = PRN
                    Findsp3 = find( SP3_XYZ.data(:,SP3_XYZ.col.PRN) == i );
                    Findobs = find( rinexobs.data(:,rinexobs.col.PRN) == i );
                    pointX = SP3_XYZ.data(Findsp3,SP3_XYZ.col.gps_seconds);
                    pointYx = SP3_XYZ.data(Findsp3,SP3_XYZ.col.X);
                    pointYy = SP3_XYZ.data(Findsp3,SP3_XYZ.col.Y);
                    pointYz = SP3_XYZ.data(Findsp3,SP3_XYZ.col.Z);
                    pointY = SP3_XYZ.data(Findsp3,SP3_XYZ.col.dtsat)*1e-6;
                    while 1
                        Temission = rinexobs.data(Findobs,rinexobs.col.Temission);
                        XX = lagrange_interpolation(Temission,pointX,pointYx,10);
                        YY = lagrange_interpolation(Temission,pointX,pointYy,10);
                        ZZ = lagrange_interpolation(Temission,pointX,pointYz,10);
                        rsat = [XX , YY , ZZ ];
                        rsat_r0rcv = (rsat - ones(length(Findobs),1)*r0rcv_xyz');
                        dt = sqrt(sum(rsat_r0rcv.^2,2))/c;
                        rinexobs.data(Findobs,rinexobs.col.Temission) = rinexobs.data(Findobs,rinexobs.col.TOW) - dt;
                        Temission = rinexobs.data(Findobs,rinexobs.col.Temission);
                        XX = lagrange_interpolation(Temission,pointX,pointYx,10);
                        YY = lagrange_interpolation(Temission,pointX,pointYy,10);
                        ZZ = lagrange_interpolation(Temission,pointX,pointYz,10);
                        rsat0 = [XX , YY , ZZ ];
                        Norm = norm ( rsat0 - rsat);
                        if Norm<1e-3
                            rinexobs.data(Findobs,rinexobs.col.dt_sat_rcv) = dt;
                            rinexobs.data(Findobs,rinexobs.col.Temission) = (rinexobs.data(Findobs,rinexobs.col.TOW) - dt -...
                                rinexobs.data(Findobs,rinexobs.col.dtrcv));
                            Temission = rinexobs.data(Findobs,rinexobs.col.Temission);
                            XX = lagrange_interpolation(Temission,pointX,pointYx,10);
                            YY = lagrange_interpolation(Temission,pointX,pointYy,10);
                            ZZ = lagrange_interpolation(Temission,pointX,pointYz,10);
                            rsat = [XX , YY , ZZ ];
                            rinexobs.data(Findobs,rinexobs.col.Xsat:rinexobs.col.Zsat) = rsat;
                            break;
                        end
                    end
                    dtsat = lagrange_interpolation(Temission,pointX,pointY,10);
                    rinexobs.data(Findobs,rinexobs.col.dtsat) = dtsat;
                end
            case 'navigation'
                
                PRN = unique(nav.data(:,nav.col.PRN))';
                Prn = nav.data(:,nav.col.PRN);
                for i = PRN
                    Find = find(Prn == i); Find = Find(1);
                    TGD = nav.data(Find,nav.col.TGD)*299792458;
                    rinexobs.data( find(rinexobs.data(:,rinexobs.col.PRN)==i),rinexobs.col.TGD) = TGD;
                end
                
                rinexobs.data(:,rinexobs.col.Temission) = rinexobs.data(:,rinexobs.col.TOW);
                PRN = unique(nav.data(:,nav.col.PRN))';
                for i = PRN
                    Findnav = find( nav.data(:,nav.col.PRN) == i );
                    Findobs = find( rinexobs.data(:,rinexobs.col.PRN) == i );
                    while 1
                        [rinexobs] = NAV2ECEF(nav,rinexobs,i);
                        rsat = rinexobs.data(Findobs,rinexobs.col.Xsat:rinexobs.col.Zsat);
                        rsat_r0rcv = (rsat - ones(length(Findobs),1)*r0rcv_xyz');
                        dt = sqrt(sum(rsat_r0rcv.^2,2))/c;
                        rinexobs.data(Findobs,rinexobs.col.Temission) = rinexobs.data(Findobs,rinexobs.col.TOW) - dt;
                        [rinexobs] = NAV2ECEF(nav,rinexobs,i);
                        rsat0 = rinexobs.data(Findobs,rinexobs.col.Xsat:rinexobs.col.Zsat);
                        Norm = norm ( rsat0 - rsat);
                        if Norm<1e-3
                            rinexobs.data(Findobs,rinexobs.col.dt_sat_rcv) = dt;
                            rinexobs.data(Findobs,rinexobs.col.Temission) = (rinexobs.data(Findobs,rinexobs.col.TOW) - dt -...
                                rinexobs.data(Findobs,rinexobs.col.dtrcv));
                            [rinexobs] = NAV2ECEF(nav,rinexobs,i);
                            break;
                        end
                    end
                    TOE = nav.data(Findnav,nav.col.TOE);
                    Temission = rinexobs.data(Findobs,rinexobs.col.Temission);
                    Nearest = [];
                    for j=1:length(Temission)
                        MIN = abs(TOE-Temission(j));
                        nnn = find(MIN==min(MIN));
                        Nearest(j,1) = Findnav(nnn);
                    end
                    a0 = nav.data(Nearest,nav.col.SV_Clock_Bias);
                    a1 = nav.data(Nearest,nav.col.SV_Clock_Drift);
                    a2 = nav.data(Nearest,nav.col.SV_Clock_Drift_Rate);
                    toe = nav.data(Nearest,nav.col.TOE);
                    [AAAA,drel] =  NAV2ECEF(nav,rinexobs,i);
                    dtsat = a0 + a1.*(Temission - toe) + a2.*(Temission - toe).^2 + drel;
                    rinexobs.data(Findobs,rinexobs.col.dtsat) = dtsat;
                end
        end
end
dt = rinexobs.data(:,rinexobs.col.dt_sat_rcv);
dw = we*dt;
for i = 1:length(dw)
    rsat = rinexobs.data(i,rinexobs.col.Xsat:rinexobs.col.Zsat)';
    rsat = Rotation(3,dw(i))*rsat;
    rinexobs.data(i,rinexobs.col.Xsat:rinexobs.col.Zsat) = rsat';
end
% compute Elevation & Azimuth
Res = cart_Geo(r0rcv_xyz,2);
rinexobs.Geodetic = Res;
phi = Res(1);
lambda = Res(2);
h = Res(3);
r0rcv_xyz = r0rcv_xyz(:);
K = size(rinexobs.data,1);
r_ct =  rinexobs.data(:,rinexobs.col.Xsat:rinexobs.col.Zsat)';
dr = r_ct - r0rcv_xyz*ones(1,K);
clear CC;
for k=1:K;
    CC(k,1)=norm(dr(:,k));
end
ro = dr./((CC*ones(1,3))');
eEe = [-sin(lambda),cos(lambda),0]';
nNn=[-cos(lambda)*sin(phi) ,-sin(lambda)*sin(phi) , cos(phi)]';
uUu= [cos(lambda)*cos(phi), sin(lambda)*cos(phi), sin(phi)]';
dU=[eEe nNn uUu]*dr;
Elevationa = asin( sum(ro.*(uUu*ones(1,k))) )*180/pi;
Az = atan2( sum(ro.*(eEe*ones(1,k))),sum(ro.*(nNn*ones(1,k))) )*180/pi;
for ii = 1:K
    d_ion(ii,1) = KlobucharAlgorithm(Az(ii)*pi/180,Elevationa(ii)*pi/180,phi,lambda,rinexobs.data(ii,rinexobs.col.TOW),alfa,beta);
end
rinexobs.data(:,rinexobs.col.I1_klobuchar) = d_ion;
rinexobs.data(:,rinexobs.col.Azimuth:rinexobs.col.Elevation) = [Az',Elevationa'];
Data = GPSt2Date(rinexobs.data(100,1),rinexobs.data(100,2));
Data = Data(1:3);
DOYy = DOY(Data);
for k = 1:K
    rinexobs.data(k,rinexobs.col.Trcollins) = collins( phi,h,DOYy,Elevationa(k)*pi/180 );
end
% delta ro_rel
rrcv0 = r0rcv_xyz*ones(1,K);
rsat0 = r_ct;
for k=1:K;
    rsat_rcv(k,1)=norm(dr(:,k));
    rrcv(k,1)=norm(rrcv0(:,k));
    rsat(k,1)=norm(rsat0(:,k));
end
delta_ro_rel = 2*GM*log( (rsat + rrcv + rsat_rcv)./(rsat + rrcv - rsat_rcv) )./(c^2);
delta_ro_rel = delta_ro_rel(:);
rinexobs.data(:,rinexobs.col.delta_ro_rel) = delta_ro_rel;

switch Filetype
    case 'sp3'
        satnav = SP3_XYZ;
    case 'navigation'
        satnav = nav;
end
[rinexobs]=Hatch_filter_GNSS(rinexobs);

end
function [R]=Rotation(i,j)
%Rotation matrix
if i==3
    R = [cos(j) sin(j) 0 ; -sin(j) cos(j) 0 ; 0 0 1];
else if i==2
        R = [cos(j) 0 -sin(j) ; 0 1 0 ; sin(j) 0 cos(j)];
    else if i==1
            R = [1 0 0 ; 0 cos(j) sin(j) ; 0 -sin(j) cos(j)];
        end
    end
end
end
function result = cart_Geo(XXX,type)
a = 6378137;
b = 6356752.3142;
e2=(a^2-b^2)/(a^2);
switch type
    case 1
        %%
        %phi lambda h to kartezian
        phi = XXX(1);
        Landa = XXX(2);
        h = XXX(3);
        N=a^2/sqrt(a^2*cos(phi)^2+b^2*sin(phi)^2);
        X=(N+h)*cos(phi)*cos(Landa);
        Y=(N+h)*cos(phi)*sin(Landa);
        Z=(N*(b^2/a^2)+h)*sin(phi);
        result = [X Y Z];
    case 2
        %--------------------------------------
        %kartezian to phi lambda h
        X = XXX(1);
        Y = XXX(2);
        Z = XXX(3);
        p=sqrt(X^2+Y^2);
        phi0=atan((Z/p)/(1-e2));
        N0=a/(1-e2*sin(phi0)^2)^.5;
        h0=p/cos(phi0)-N0;
        Landa=wrapTo2Pi(atan2(Y,X));
        repitition=0;delta_phi=1;delta_h=1;
        while abs(delta_phi)>10^-12 && abs(delta_h)>10^-5
            N=a/(1-e2*sin(phi0)^2)^.5;
            h=(p/cos(phi0))-N;
            phi=atan((Z/p)*((1-((e2*N)/(N+h)))^-1));
            delta_phi=phi-phi0;delta_h=h-h0;phi0=phi;h0=h;
            %             repitition=repitition+1
        end
        result = [phi Landa h];
end
end
function [ Tr ] = collins( phir,H,D,EL )

% inputs: latitude ang logitude of reciver-Day of year-Elevation
% outputs: tropospheric delay

phi = [15,30,45,60,75];
phi = deg2rad(phi);
% phir = rad2deg(phir);

%row: phi- columns:P0(pressure)-T0(Temperature)-e0(Water
%Vapor)-beta0(Temperature lapse rate)-landa0(Water Vapor lapse rate)
epsillon0 = [1013.25,299.65,26.31,6.30*10^(-3),2.77
    1017.25,294.15,21.79,6.05*10^(-3),3.15
    1015.75,283.15,11.66,5.58*10^(-3),2.57
    1011.75,272.15,6.78,5.39*10^(-3),1.81
    1013.00,263.65,4.11,4.53*10^(-3),1.55];


%columns: deltaP-deltaT-delta e-delta BETA-delta Landa
delta_epsillon = [0,0,0,0*10^-3,0
    -3.75,7.00,8.85,0.25*10^(-3),0.33
    -2.25,11.00,7.24,0.32*10^(-3),0.46
    -1.75,15.00,5.36,0.81*10^(-3),0.74
    -0.50,14.50,3.39,0.62*10^(-3),0.30];

if phir >= 0
    D_min=28;
elseif phir<0
    D_min=211;
end

% step1: interploate
% columns: P-T-e-beta-landa
for i=1:5
    epsilon0(i) = interp1(phi,epsillon0(:,i)',phir,'spline');
    delta_epsilon(i) = interp1(phi,delta_epsillon(:,i)',phir,'spline');
end

%step2 => calculation of parameters (P,T,e,beta,landa)
for i=1:5
    dD = D - D_min;
    epsilon(i) = epsilon0(i)-(delta_epsilon(i)*cos((2*pi*dD)/365.25));
end
P = epsilon(1); T = epsilon(2); e = epsilon(3); beta = epsilon(4); landa = epsilon(5);

%step3 => calculation of initial value  (Tr_z0,d  &  Tr_z0,w)
%input coefficients
k1=77.604 ; k2=382000 ; Rd=287.054 ; gm=9.784 ; g=9.80665;
Tr_d0 = ((10^-6)*k1*Rd*P) / gm;                                  %Tr_z0,d
Tr_w0 = ((10^-6)*k2*Rd*e) / (((landa+1)*gm) - (beta*Rd*T));    %Tr_z0,w

%step4 => calculation of tropospheric delay (Tr_z,d  &  Tr_z,w)
v1 = 1-((beta*H)/T);
v2 = g/(Rd*beta);
v3 = (((landa+1)*g)/(Rd*beta))-1;
Tr_d =  (v1^v2)*Tr_d0; %Tr_z,d
Tr_w = (v1^v3)*Tr_w0; %Te_z,w

%step5 => input the mapping function (M(E))
s = sin(EL);
M = 1.001/sqrt(0.002001 + s^2);

%step6 => final tropospheric delay
Tr=(Tr_d+Tr_w)*M;


end
function out=DOY(r2)
y=r2(1);m=r2(2);d=r2(3);
k=[31 28 31 30 31 30 31 31 30 31 30 31];
y=2000+y;
R1=abs(2012-y);
R2=(R1/4)-floor(R1/4);
if R2==0
    k(2)=29;
end
if m>1
    doy=sum(k(1:m-1))+d;
elseif m==1
    doy=d;
end
out=doy;
end
function Data = GPSt2Date(gps_week,gps_sec)
JD = 7*( gps_week + (gps_sec/(86400*7)) ) + 2444244.5;
[ Data ] = JD2Date(JD);
end
function [ Data ] = JD2Date(JD)
a = floor(JD + 0.5);
b = a + 1537;
c = floor( (b - 122.1)/365.25 );
d = floor(365.25*c);
e = floor( (b-d)/30.6001 );
D = b - d - floor(30.6001*e) +  (JD + 0.5) - floor(JD + 0.5);
M = e - 1 - 12*floor(e/14);
Y = c - 4715 - floor( (7 + M)/10 );
h = degrees2dms((D - floor(D))*24 );
D = floor(D);
H = h(1);
Min = h(2);
S = ceil(h(3));
S = h(3);
Data = [Y M D H Min S];
end
function d_ion = KlobucharAlgorithm(Az,El,phi0,lambda0,gps_sec,alfa,beta)
global c Re hgps phiP lambdaP
c=299792458;
Re=6378000;
hgps=350000;
phiP=deg2rad(78.3);
lambdaP=deg2rad(291);
AI=0;
PI=0;
%step1 % Earth-centred angle
psi=( (pi/2)- El - asin(Re*cos(El)/(Re+hgps)) );
%step2 % latitude of the IPP
phiIP=asin( sin(phi0)*cos(psi) + cos(phi0)*sin(psi)*cos(Az) );
%step3 % longitude of the IPP
lambdaIP=lambda0 + psi*sin(Az)/cos(phiIP);
% step4: Find the geomagnetic latitude of the IPP
phim = asin( sin(phiIP)*sin(phiP) + cos(phiIP)*cos(phiP)*cos(lambdaIP-lambdaP) );
% step5: Find the local time at the IPP
t1=43200*(lambdaIP/pi)+gps_sec;
t1 = mod(t1,86400);
% step6: Compute the amplitude of ionospheric delay
for i=1:4
    AI=AI+alfa(i,1)*(phim/pi)^(i-1);
end
if AI<0
    AI=0;
end
% step7: Compute the period of ionospheric delay
for i=1:4
    PI=PI+beta(i,1)*(phim/pi)^(i-1);
end
if PI<72000
    PI=72000;
end
%step 8: Compute the phase of ionospheric delay
XI=2*pi*(t1-50400)/PI;
% step9: Compute the slant factor (ionospheric mapping function)
F=(1-(Re*cos(El)/(Re+hgps))^2)^(-0.5);
% step10: Compute the ionospheric time delay
if abs(XI)<pi/2
    I1 =( 5e-9 + AI*cos(XI) )*F;
elseif abs(XI)>=pi/2
    I1 = 5e-9*F;
end
% Ik = ((f1/fk)^2)*I1;
d_ion=I1*c;
% freqs1=1575.42e6;   
% freqs2=1227.60e6;
% alpha1=1/((freqs1/freqs2)^2-1);
% d_ion = d_ion*alpha1;
end
function [rinexx]=Hatch_filter_GNSS(rinexx)
%--------------------------------------
N = 360;
freqs1=1575.42e6;   
freqs2=1227.60e6;
alpha1=1/((freqs1/freqs2)^2-1);
 TOW=rinexx.data(:,rinexx(1).col.TOW);
 interval=unique(TOW);
 interval=interval(2)-interval(1);
PRN=rinexx.data(:,rinexx.col.PRN);
rinex = rinexx;
k=1;
for i=[unique(PRN)]'
    rinex(k).PRN=rinexx(1).data(find(PRN==i),:);
    k=k+1;
end
%------------------------------------------------
for i=1:k-1
Arc=rinex(i).PRN(:,rinex(1).col.TOW)./interval;
Arc=[Arc;0]-[0;Arc];
Arc(1)=[];
Arc=find(abs(Arc)>1);
Arc=[0;Arc];
Arc=[Arc+1,Arc];
Arc=[Arc(1:end-1,1),Arc(2:end,2)];
rinex(i).Arc=Arc;
end
for Sat = [unique(PRN)]'
sat = find(unique(PRN)==Sat);
Arc = rinex(sat).Arc;
C1=rinex(sat).PRN(:,rinex(1).col.C1);
P2=rinex(sat).PRN(:,rinex(1).col.P2);
L1=rinex(sat).PRN(:,rinex(1).col.L1);
L2=rinex(sat).PRN(:,rinex(1).col.L2);
Time=rinex(sat).PRN(:,rinex(1).col.TOW);
L1_D=L1+2*alpha1*(L1-L2);
L1_LC=((freqs1^2)*L1-(freqs2^2)*L2)/((freqs1^2)-(freqs2^2));
C1_PC=((freqs1^2)*C1-(freqs2^2)*P2)/((freqs1^2)-(freqs2^2));
for iii = 1:size(Arc,1)
%single frequency smoothed code
 Rs_hat(Arc(iii,1):Arc(iii,2)) = HatchFilter( C1(Arc(iii,1):Arc(iii,2)),L1(Arc(iii,1):Arc(iii,2)),N);
% divergence-free carrier smoothed code
 Rd_hat(Arc(iii,1):Arc(iii,2)) = HatchFilter( C1(Arc(iii,1):Arc(iii,2)),L1_D(Arc(iii,1):Arc(iii,2)),N);
%  Ionosphere Free smoother
RI_hat(Arc(iii,1):Arc(iii,2)) = HatchFilter( C1_PC(Arc(iii,1):Arc(iii,2)),L1_LC(Arc(iii,1):Arc(iii,2)),N);
end
Rs_hat = Rs_hat';
Rd_hat = Rd_hat';
RI_hat = RI_hat';
rinexx.data(find(PRN==Sat),rinexx.col.Hatch_DF) = Rd_hat;
rinexx.data(find(PRN==Sat),rinexx.col.Hatch_IO_free) = RI_hat;
Rs_hat = [];
Rd_hat = [];
RI_hat = [];
end
end
function [ RH ] = HatchFilter( C1,L1,N)
RH = C1(1,:);
    for k=2:length(C1)
        if k<N
            n = k;
        else
            n = N;
        end
        if (C1(k-1)==0)
            RH(k) = C1(k);
        else
            RH(k) = (1/n)*C1(k) + ((n-1)/n)*(RH(k-1)+L1(k)-L1(k-1));
        end
    end
end