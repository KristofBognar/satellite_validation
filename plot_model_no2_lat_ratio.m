% plot ratios of evening twilight no2 for different latitudes

% load data
load('/home/kristof/work/models/NO2_box_model/no2_cols_12-40_km_2017_79lat.mat');
model79=model_no2;

load('/home/kristof/work/models/NO2_box_model/no2_cols_12-40_km_2017_80lat.mat');
model80=model_no2;

load('/home/kristof/work/models/NO2_box_model/no2_cols_12-40_km_2017_81lat.mat');
model81=model_no2;

clear model_no2

% select evening twilght
model79(model79.ampm==0,:)=[];
model80(model80.ampm==0,:)=[];
model81(model81.ampm==0,:)=[];

% d76 was recorded as 75 for some reason
model79.day(511:527)=model79.day(511:527)+1;
model79.fractional_time(511:527)=model79.fractional_time(511:527)+1;

model80.day(511:527)=model80.day(511:527)+1;
model80.fractional_time(511:527)=model80.fractional_time(511:527)+1;

model81.day(511:527)=model81.day(511:527)+1;
model81.fractional_time(511:527)=model81.fractional_time(511:527)+1;

% % time grid slightly different for each lat: interpolate to sparsest grid
% no2_79_tmp=interp1(model79.fractional_time,model79.no2,model81.fractional_time);
% sza_79_tmp=interp1(model79.fractional_time,model79.sza,model81.fractional_time);
% no2_80_tmp=interp1(model80.fractional_time,model80.no2,model81.fractional_time);
% sza_80_tmp=interp1(model80.fractional_time,model80.sza,model81.fractional_time);
% 
% no2_81_tmp=model81.no2;
% sza_81_tmp=model81.sza;

days=46:300;
no2_79=[];
no2_80=[];
no2_81=[];

% get 90 deg or closest

for i=days
    
    [~,tmp]=min(abs(model79.sza(model79.day==i)-90));
    tmp2=find(model79.day==i);
    ind=tmp+tmp2(1)-1;
    
    no2_79=[no2_79,model79.no2(ind)];

    [~,tmp]=min(abs(model80.sza(model80.day==i)-90));
    tmp2=find(model80.day==i);
    ind=tmp+tmp2(1)-1;
    
    no2_80=[no2_80,model80.no2(ind)];

    [~,tmp]=min(abs(model81.sza(model81.day==i)-90));
    tmp2=find(model81.day==i);
    ind=tmp+tmp2(1)-1;
    
    no2_81=[no2_81,model81.no2(ind)];
    
end


% plot ratios
figure()

plot(days,no2_79./no2_80,'r--'), hold on
plot(days,no2_79./no2_81,'b--'), hold on
plot(days,no2_80./no2_81,'g--'), hold on

legend('79/80','79/81','80/81','location','best')
xlim([40,306])



