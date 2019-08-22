%%
% data download for RF and CCP

% Add path

clear all
%%


addpath('/Users/likelymax/Documents/Research/SAC_DATA')

%

%%  Set time database

% set fetch time

year1 = 2017;
year2 = 2018;

month_days = cell(16,3);
month = 9;
month1 = 1;
for i = 1:16
    if i <= 4
        disp(['2017' month])
        month_days{i,1,1} = num2str(year1);
        if month<10
        month_days{i,2,1} = ['0' num2str(month)];
        else
            month_days{i,2,1} = num2str(month);
        end
        month = month+1;
    elseif i>4
        %disp(month)
        disp(['2018' month1])
        month_days{i,1,1} = num2str(year2);
        if month<10
        month_days{i,2,1} = ['0' num2str(month1)];
        else
            month_days{i,2,1} = num2str(month1);
        end
        month1= month1+1;
    end 
end

month_days{1,3} = 30;
month_days{2,3} = 31;
month_days{3,3} = 30;
month_days{4,3} = 31;
month_days{5,3} = 31;
month_days{6,3} = 29;
month_days{7,3} = 31;
month_days{8,3} = 30;
month_days{9,3} = 31;
month_days{10,3} = 30;
month_days{11,3} = 31;
month_days{12,3} = 31;
month_days{13,3} = 30;
month_days{14,3} = 31;
month_days{15,3} = 30;
month_days{16,3} = 31;

%% Fetch data

tracewhole = cell(16,2);
for ii = 1:length(month_days)
    if ii == 1
        a = 15;
    else
        a = 1;
    end
    n = month_days{ii,3};
    mytrace = cell(n,2);
        for i = a:n
            
            Mon = month_days{ii,2};
            year = month_days{ii,1};
            if i <10
                m = ['0' num2str(i)];
            elseif i >=10
                m = num2str(i);
            end
            mytrace{i,1} = [year '-' Mon '-' m  ];
            %disp(['0' m]);
            mytrace1 = irisFetch.Traces('YG','OR13','LOG','EHE,EHN,EHZ',[ year '-' Mon '-' m ' 00:00:00'],[year '-' Mon '-' m ' 23:59:59']);
            disp([year '-' Mon '-' m ' 00:00:00'])
            disp(mytrace1)
            mytrace{i,2} = mytrace1;
        end
    
    tracewhole{ii,1} = [year '-' Mon];
    tracewhole{ii,2} = mytrace;
    disp(tracewhole{ii,1})
    
end

%% save data 

save ('tracewhole.mat','-v7.3')
save month_days

%sampletimes1 = linspace (mytrace.startTime, mytrace.endTime, mytrace.sampleCount);
%plot(sampletimes1, mytrace1.data)

%% fetch the data

clear all
clc

mytrace = irisFetch.Traces('XT','D20','**','BHE,BHN,BHZ', '2012-10-28 02:30:00', '2012-10-28 05:59:59');

colors=brighten(lines(numel(mytrace)),-0.33); % define line colors
for n=1:numel(mytrace)
  %figure(n);
  tr = mytrace(n);
  data=double(tr.data) ./ tr.sensitivity;    % scale the data
  sampletimes = linspace(tr.startTime,tr.endTime,tr.sampleCount);
  plot(sampletimes, data,'color', colors(n,:));
  hold on;
  %datetick;
  %ylabel(tr.sensitivityUnits); % assumes all units are the same 
%title(['UI-ANMO traces, starting ', datestr(mytrace(1).startTime)]); 
%legend(strcat({mytrace.channel},'-',{mytrace.location}),'location','northwest');
end
hold off;
datetick;
ylabel(tr.sensitivityUnits); % assumes all units are the same 
%title(['UI-ANMO traces, starting ', datestr(mytrace(1).startTime)]); 
legend(strcat({mytrace.channel},'-',{mytrace.location}),'location','northwest');

save('mytrace.mat','-v7.3')

disp('end')

%% bandpass filter

load('mytrace.mat')
samplerate = 2000;

%figure(3)
clf
count = 0;
sample_rate = samplerate;

A_stop1 = 80;		% Attenuation in the first stopband = 60 dB
F_stop1 = .5;		% Edge of the stopband  Hz
F_pass1 = 1.;	% Edge of the passband  Hz
F_pass2 =45;	% Closing edge of the passband  Hz
F_stop2 =49;	% Edge of the second stopband  Hz
A_stop2 = 50;		% Attenuation in the second stopband = 60 dB
A_pass = 400;		% Amount of ripple allowed in the passband = 5 dB

%if count==1
%     F_stop1=2.;
%     F_pass1=2.5;
% elseif count==2
%     F_stop1=4;
%     F_pass1=5;
% elseif count==3
%     F_stop1=6;
%     F_pass1=7;
% elseif count==4
%     F_stop1=20;
%     F_pass1=21;
% end


filter_object = fdesign.bandpass;

   BandPassSpecObj = ...
    fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
 		F_stop1, F_pass1, F_pass2, F_stop2, A_stop1, A_pass, ...
 		A_stop2, 100) ;
    
 %BandPassSpecObj = ...   
      % fdesign.highpass('Fst,Fp,Ast,Ap', ...
		%F_stop1, F_pass1,  A_stop1, A_pass, ...
		% 100) 
     
%       BandPassSpecObj = ...   
%        fdesign.lowpass('Fp,Fst,Ap,Ast', ...
% 		 F_pass2,F_stop2,  A_stop2, A_pass, ...
% 		 100) 



BandPassFilt = design(BandPassSpecObj,'equiripple');


if BandPassSpecObj.Response== 'Bandpass'
 disp([BandPassSpecObj.Response ' Filter: ' num2str(F_pass1) ' to ' num2str(F_pass2)  ' Hz  '])
elseif BandPassSpecObj.Response == 'Highpass'
 disp([BandPassSpecObj.Response ' Filter: ' num2str(F_pass1)  ' Hz cutoff '])
else
    disp([BandPassSpecObj.Response ' Filter: ' num2str(F_pass2)  ' Hz cutoff '])
end

filtered = filtfilt(BandPassFilt.Numerator,samplerate,mytrace(6).data);

tr = mytrace(6);
data=double(tr.data) ./ tr.sensitivity;    % scale the data
sampletimes = linspace(tr.startTime,tr.endTime,tr.sampleCount);

%subplot(2,1,1)

figure

subplot(2,1,1)

%set(gca,'fontsize',18)
hold on

plot(sampletimes,filtered/std(filtered),'r-');
grid on
%figure(2)
subplot(2,1,2)

plot(sampletimes,double(tr.data)/1e3,'b-');
grid on


%% bandpass filter

load('mytrace.mat')
%samplerate = 1000;

%figure(3)
%clf
%count = 0;
%sample_rate = samplerate;

for i = 1:2
    for j = 1:2
        count = count + 1;

A_stop1 = 80;		% Attenuation in the first stopband = 60 dB
F_stop1 = 0.5;		% Edge of the stopband  Hz
F_pass1 = 1.;	% Edge of the passband  Hz
F_pass2 =40;	% Closing edge of the passband  Hz
F_stop2 =50;	% Edge of the second stopband  Hz
A_stop2 = 80;		% Attenuation in the second stopband = 60 dB
A_pass = 5;		% Amount of ripple allowed in the passband = 5 dB

%if count==1
 %   F_stop1=2.;
  %  F_pass1=2.5;
%elseif count==2
 %   F_stop1=4;
  %  F_pass1=5;
%elseif count==3
 %   F_stop1=6;
  %  F_pass1=7;
%elseif count==4
 %   F_stop1=20;
  %  F_pass1=21;
%end


filter_object = fdesign.bandpass

   BandPassSpecObj = ...
    fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
 		F_stop1, F_pass1, F_pass2, F_stop2, A_stop1, A_pass, ...
 		A_stop2, 100) 
    
 %BandPassSpecObj = ...   
      % fdesign.highpass('Fst,Fp,Ast,Ap', ...
		%F_stop1, F_pass1,  A_stop1, A_pass, ...
		% 100) 
     
%       BandPassSpecObj = ...   
%        fdesign.lowpass('Fp,Fst,Ap,Ast', ...
% 		 F_pass2,F_stop2,  A_stop2, A_pass, ...
% 		 100) 



BandPassFilt = design(BandPassSpecObj,'equiripple')


if BandPassSpecObj.Response== 'Bandpass'
 disp([BandPassSpecObj.Response ' Filter: ' num2str(F_pass1) ' to ' num2str(F_pass2)  ' Hz  '])
elseif BandPassSpecObj.Response == 'Highpass'
 disp([BandPassSpecObj.Response ' Filter: ' num2str(F_pass1)  ' Hz cutoff '])
else
    disp([BandPassSpecObj.Response ' Filter: ' num2str(F_pass2)  ' Hz cutoff '])
end

%ANTf=cell(4,4);
% for row=1:4
%     for col=1:4
 %       ANTf{row,col}=filtfilt(BandPassFilt.Numerator,samplerate,ANT{row,col});  % Note the use of filtfilt
 %    end
% end

subplot(2,2,count)

set(gca,'fontsize',18)
hold on

plot(Td,x/1e3,'b-','linewidth',2)
plot(Tr,x/1e3,'k-','linewidth',2)
plot(Th,x/1e3,'g-','linewidth',2)

plot(-Td,x/1e3,'b-','linewidth',2)
plot(-Tr,x/1e3,'k-','linewidth',2)
plot(-Th,x/1e3,'g-','linewidth',2)



tempdist = 0;
for row=1:4
    for col=1:4
     tempdist = tempdist ;   
     plot(lags/sample_rate,0.25*ANTf{row,col}/std(ANTf{row,col}) ...
         + station_distance(row,col)/1e3+tempdist,'r-')   
%        plot(lags/sample_rate,0.25*ANTf{row,col} ...
%          + station_distance(row,col)/1e3+tempdist,'r-')   
    % plot(lags/samplerate,.25*ANTf/std(ANTf)+tempdist,'r-')
      text(35,tempdist+0.5,string2{row,col},'fontsize',14)   
    end
end
xlim([-45 45])
xlabel('time lag (s)')
ylabel('distance (km)')

ylim([0 50])

% c_i=2600;
% c_rock=4850;
% c_sw=1450;
% D=900;
% D2=100;
% 
% theta=asin(c_sw/c_rock);
% disp(['critical angle is ' num2str(theta*180/pi) ' degrees.'])
% 
% x=[0:100:50000]';
% 
% Td=x/c_i;
% 
% Tr=2*sqrt( x.^2/4 + D^2)/c_sw;
% 
% Tr2=2*sqrt( x.^2/4 + D2^2)/c_sw;
% 
% 
% Th= (x-2*D*tan(theta))/c_rock + 2*D/(c_sw*cos(theta));
% 
% Th2= (x-2*D2*tan(theta))/c_rock + 2*D2/(c_sw*cos(theta));

% c_i=2600;; % p wave phase speed in solid ice m/s
% 
% x=[0:1:50]; % spatial locations in km, source receiver separation
% T_d = x/(c_i/1e3);
% 
% % Eurika #2 Gosh, there's a slower phase than the direct wave, could it be
% %  hydro acoustic?
% 
% c_sw = 1450 ;  % sound speed in sea water at 0C
% D = 1000; % water depth
% 
% T_reflected = 2* sqrt(x.^2/4 + D^2/1e6)/( c_sw/1e3 );
% 
% % Eurika #3 Gosh, what about the "head wave" (wave in rock)?
% c_rock=4850;
% theta=asin(c_sw/c_rock);
% disp(['critical angle is ' num2str(theta*180/pi) ' degrees.'])
% 
% Th= (x-2*D/1e3*tan(theta))/(c_rock/1e3) + 2*D/(c_sw*cos(theta));
    end
end
legend('direct wave','reflected wave', 'refracted wave',...
['c ice = ' num2str(c_i) ' m/s'],['c water = ' num2str(c_sw) ' m/s'], ['c rock = ' num2str(c_rock) ' m/s'],'location','southeast','fontsize',10)

%% fetch the data (test)

clc

% set time period

mm = 1;



for lin1 = 319: len
%lin1 = 39;

%for lin1 = 1:14
disp(['working on station ', datensta{lin1,1}, datensta{lin1, 2}]);


lenevent = length(datensta{lin1,5});

%lin2 = 25;

for lin2 = 1: lenevent

tp1 = datensta{lin1,5}{lin2,2}(1,1);

NT = datensta{lin1,1};

sta_name = datensta{lin1,2};

%tp2start = datenum(datestr(tp1 + minutes(5)));
 tp2start = tp1;
% 
 tp2end = datenum(datestr(tp1 + minutes(60)));
 
 mytrace2=struct('network',[],'station',[],'location',[],...
             'channel',[],'quality',[],'latitude',[],'longitude',[],...
             'elevation',[],'depth',[],'azimuth',[],'dip',[],...
             'eventlatitude',[],'eventlongitude',[],'eventdepth',[],'eventmag',[],...
             'sensitivity',[],'sensitivityFrequency',[],'instrument','',...
             'sensitivityUnits',[],'data',[],'sampleCount',[],'sampleRate',[],...
             'startTime',[],'endTime',[],'sacpz',[]);
%mytrace = irisFetch.Traces( NT,sta_name,'*','BHE,BHN,BHZ', datestr(tp2start), datestr(tp2end),'WRITESAC:/Users/yitanwang/Documents/Research/DATA/OREGON_CCP/SAC_DATA');
%ev = irisFetch.Events('starttime',datestr(tp2start),'endtime',datestr(tp2end),'minmag', 6.0, 'limit', 1 ,'includeallorigins',true,'includeallmagnitudes',true);
mytrace = irisFetch.Traces( NT,sta_name,'*','BHE,BHN,BHZ', datestr(tp2start), datestr(tp2end));
if isempty(mytrace) == 1
 irisFetch.Trace2SAC(mytrace,'/Users/likelymax/Documents/Research/SAC_DATA/hold');   
else
 ev = irisFetch.Events('starttime',datestr(tp2start),'endtime',datestr(tp2end),'minmag', 6.0, 'limit', 1 ,'includeallorigins',true,'includeallmagnitudes',true);
for i = 1:3
    mytrace2(i).network = mytrace(i).network;
    mytrace2(i).station = mytrace(i).station;
    mytrace2(i).location = mytrace(i).location;
    mytrace2(i).channel = mytrace(i).channel;
    mytrace2(i).quality = mytrace(i).quality;
    mytrace2(i).latitude = mytrace(i).latitude;
    mytrace2(i).longitude = mytrace(i).longitude;
    mytrace2(i).elevation = mytrace(i).elevation;
    mytrace2(i).depth = mytrace(i).depth;
    mytrace2(i).azimuth = mytrace(i).azimuth;
    mytrace2(i).dip = mytrace(i).dip;
    mytrace2(i).eventlatitude = ev.PreferredLatitude;
    mytrace2(i).eventlongitude = ev.PreferredLongitude;
    mytrace2(i).eventdepth = ev.PreferredDepth;
    mytrace2(i).eventmag = ev.Magnitudes.Value;
    mytrace2(i).eventname = ev.PreferredTime;
    mytrace2(i).sensitivity = mytrace(i).sensitivity;
    mytrace2(i).sensitivityFrequency = mytrace(i).sensitivityFrequency;
    mytrace2(i).instrument = mytrace(i).instrument;
    mytrace2(i).sensitivityUnits = mytrace(i).sensitivityUnits;
    mytrace2(i).data = mytrace(i).data;
    mytrace2(i).sampleCount = mytrace(i).sampleCount;
    mytrace2(i).sampleRate = mytrace(i).sampleRate;
    mytrace2(i).startTime = mytrace(i).startTime;
    mytrace2(i).endTime = mytrace(i).endTime;
    mytrace2(i).sacpz = mytrace(i).sacpz;
    
end

irisFetch.Trace2SAC(mytrace2,'/Users/likelymax/Documents/Research/SAC_DATA/hold');
end
% if lin2 == 1 || lin2 == 2 || lin2 ==3 || lin2 == 4 || lin2 ==5 || lin2 == 7 || lin2==12 || lin2 == 13 || lin2 ==9 || lin2 == 11 || lin2 == 14 || lin2 == 6
%      tp2start = datenum(datestr(tp1 + minutes(7)));
% 
% %tp2start = tp1;
% 
%     tp2end = datenum(datestr(tp1 + minutes(37)));
% else 
%   
%     
%     tp2start = tp1;
% 
%     tp2end = datenum(datestr(tp1 + minutes(30)));
% end

%mytrace = irisFetch.Traces( NT,sta_name,'**','BHE,BHN,BHZ', datestr(tp2start), datestr(tp2end));

% colors=brighten(lines(numel(mytrace)),-0.33); % define line colors
% 
% figure(mm);

% for n=1:numel(mytrace)
%   %figure(n);
%   tr = mytrace(n);
%   data=double(tr.data) ./ tr.sensitivity;    % scale the data
%   sampletimes = linspace(tr.startTime,tr.endTime,tr.sampleCount);
%   plot(sampletimes, data,'color', colors(n,:));
%   hold on;
%   %datetick;
%   %ylabel(tr.sensitivityUnits); % assumes all units are the same 
% %title(['UI-ANMO traces, starting ', datestr(mytrace(1).startTime)]); 
% %legend(strcat({mytrace.channel},'-',{mytrace.location}),'location','northwest');
% end
% hold off;
% datetick;
% ylabel(tr.sensitivityUnits); % assumes all units are the same 
% %title(['UI-ANMO traces, starting ', datestr(mytrace(1).startTime)]); 
% %legend(strcat({mytrace.channel},'-',{mytrace.location}),'location','northwest');
% title(['station ' num2str(lin1) 'event ' num2str(lin2) ]);
% 
% %save('mytrace.mat','-v7.3')

%datensta{lin1,5}{lin2,3} = ev;

disp(['working on event ', num2str(lin2), ' out of ', num2str(lenevent)]);

%mm = mm + 10;

end


end

disp('end')
%end

% 1284

%% fetch the data (test)

clc

% set time period


%lin1 = 84;

lin1 = 321;

disp(['working on station ', datensta{lin1,1}, datensta{lin1, 2}]);

lenevent = length(datensta{lin1,5});

for lin2 = 1160 : lenevent
min1 = 0;

min2 = 60;

tp1 = datensta{lin1,5}{lin2,2}(1,1);

NT = datensta{lin1,1};

sta_name = datensta{lin1,2};

tp2start = datenum(datestr(tp1 + minutes(min1)));

%tp2start = tp1;

tp2end = datenum(datestr(tp1 + minutes(min2)));
mytrace2=struct('network',[],'station',[],'location',[],...
             'channel',[],'quality',[],'latitude',[],'longitude',[],...
             'elevation',[],'depth',[],'azimuth',[],'dip',[],'eventname',[],...
             'eventlatitude',[],'eventlongitude',[],'eventdepth',[],'eventmag',[],...
             'sensitivity',[],'sensitivityFrequency',[],'instrument','',...
             'sensitivityUnits',[],'data',[],'sampleCount',[],'sampleRate',[],...
             'startTime',[],'endTime',[],'sacpz',[]);
%mytrace = irisFetch.Traces( NT,sta_name,'*','BHE,BHN,BHZ', datestr(tp2start), datestr(tp2end),'WRITESAC:/Users/yitanwang/Documents/Research/DATA/OREGON_CCP/SAC_DATA');
mytrace = irisFetch.Traces( NT,sta_name,'*','BHE,BHN,BHZ', datestr(tp2start), datestr(tp2end));
if isempty(mytrace) == 1
 irisFetch.Trace2SAC(mytrace,'/Users/likelymax/Documents/Research/SAC_DATA/hold');   
else
 ev = irisFetch.Events('starttime',datestr(tp2start),'endtime',datestr(tp2end),'minmag', 6.0, 'limit', 1 ,'includeallorigins',true,'includeallmagnitudes',true);
for i = 1:3
    mytrace2(i).network = mytrace(i).network;
    mytrace2(i).station = mytrace(i).station;
    mytrace2(i).location = mytrace(i).location;
    mytrace2(i).channel = mytrace(i).channel;
    mytrace2(i).quality = mytrace(i).quality;
    mytrace2(i).latitude = mytrace(i).latitude;
    mytrace2(i).longitude = mytrace(i).longitude;
    mytrace2(i).elevation = mytrace(i).elevation;
    mytrace2(i).depth = mytrace(i).depth;
    mytrace2(i).azimuth = mytrace(i).azimuth;
    mytrace2(i).dip = mytrace(i).dip;
    mytrace2(i).eventlatitude = ev.PreferredLatitude;
    mytrace2(i).eventlongitude = ev.PreferredLongitude;
    mytrace2(i).eventdepth = ev.PreferredDepth;
    mytrace2(i).eventmag = ev.Magnitudes.Value;
    mytrace2(i).eventname = ev.PreferredTime;
    mytrace2(i).sensitivity = mytrace(i).sensitivity;
    mytrace2(i).sensitivityFrequency = mytrace(i).sensitivityFrequency;
    mytrace2(i).instrument = mytrace(i).instrument;
    mytrace2(i).sensitivityUnits = mytrace(i).sensitivityUnits;
    mytrace2(i).data = mytrace(i).data;
    mytrace2(i).sampleCount = mytrace(i).sampleCount;
    mytrace2(i).sampleRate = mytrace(i).sampleRate;
    mytrace2(i).startTime = mytrace(i).startTime;
    mytrace2(i).endTime = mytrace(i).endTime;
    mytrace2(i).sacpz = mytrace(i).sacpz;
    
end

irisFetch.Trace2SAC(mytrace2,'/Users/likelymax/Documents/Research/SAC_DATA/hold');
end
%irisFetch.Trace2SAC(ev,'/Users/yitanwang/Documents/Research/DATA/OREGON_CCP/SAC_DATA/hold');

% colors=brighten(lines(numel(mytrace)),-0.33); % define line colors
% 
% figure(500)
% 
% for n=1:numel(mytrace)
%   %figure(n);
%   tr = mytrace(n);
%   data=double(tr.data) ./ tr.sensitivity;    % scale the data
%   sampletimes = linspace(tr.startTime,tr.endTime,tr.sampleCount);
%   plot(sampletimes, data,'color', colors(n,:));
%   hold on;
%   %datetick;
%   %ylabel(tr.sensitivityUnits); % assumes all units are the same 
% %title(['UI-ANMO traces, starting ', datestr(mytrace(1).startTime)]); 
% %legend(strcat({mytrace.channel},'-',{mytrace.location}),'location','northwest');
% end
% hold off;
% datetick;
% ylabel(tr.sensitivityUnits); % assumes all units are the same 
% %title(['UI-ANMO traces, starting ', datestr(mytrace(1).startTime)]); 
% legend(strcat({mytrace.channel},'-',{mytrace.location}),'location','northwest');
% 
% title(['station ' num2str(lin1) 'event ' num2str(lin2) ]);
% 
% %save('mytrace.mat','-v7.3')

%datensta{lin1,5}{lin2,3} = mytrace2;

disp(['working on event ', num2str(lin2), ' out of ', num2str(lenevent)]);
end

disp('end')

%%

save('datastation.mat', '-v7.3');
disp('saved');

%end

%%

load 'datastation.mat';

