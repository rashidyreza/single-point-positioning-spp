function [rinex,Date,alfa,beta]=ReadNavigation(filename,Nk)
field = fopen(filename);
if field == -1
    disp(['The file ''' filename ''' does not exist.']);
    return;
end
n = 0;
frewind(field);
while 1
    n=n+1;
    line = fgetl(field);
    if strfind(line,'ION ALPHA');
        alfa = [str2num(line(4:51))]';
    end
    if strfind(line,'ION BETA');
        beta = [str2num(line(4:51))]';
    end
        
    if strfind(line,'END OF HEADER');
        break
    end
end
% END OF HEADER line n
frewind(field);
Data = textscan (field ,'%f','headerline',n);
Data=Data{1,1};
Data=reshape(Data,Nk,[]);
Data=Data';
rinex.data=Data;
rinex.data(:,2) = rinex.data(:,2) + 2000;
%---------------File Data-----------------
rinex.col.PRN = 1;
rinex.col.Year = 2;
rinex.col.Month = 3;
rinex.col.Day = 4;
rinex.col.Hour = 5;
rinex.col.Minute = 6;
rinex.col.second = 7;
rinex.col.SV_Clock_Bias = 8;
rinex.col.SV_Clock_Drift = 9;
rinex.col.SV_Clock_Drift_Rate = 10;
rinex.col.IODE = 11;
rinex.col.Crs = 12;
rinex.col.deltan = 13;
rinex.col.Mean_Anom = 14;
rinex.col.Cuc = 15;
rinex.col.e = 16;
rinex.col.Cus = 17;
rinex.col.a = 18;
rinex.data(:,18)=rinex.data(:,18).^2;
rinex.col.TOE = 19;
rinex.col.Cic = 20;
rinex.col.Right_Ascen_at_reference_time = 21;
rinex.col.Cis = 22;
rinex.col.i0 = 23;
rinex.col.Crc = 24;
rinex.col.Argument_of_Perigee = 25;
rinex.col.Rate_of_Right_Ascen = 26;
rinex.col.idot = 27;
rinex.col.L2_Codes_Channel = 28;
rinex.col.Week = 29;
rinex.col.L2_P_Data_Flag = 30;
rinex.col.SV_Accuraccy = 31;
rinex.col.SV_Health = 32;
rinex.col.TGD = 33;
rinex.col.IODC= 34;
rinex.col.Transmission_Time = 35;
rinex.col.Fit_Interval = 36;
Date = rinex.data(1,2:7);
fclose(field);
end