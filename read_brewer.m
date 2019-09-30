% function [ output_args ] = read_brewer( input_args )
%READ_BREWER read and reformat brewer O3 data
%   


load('/home/kristof/work/brewer/Brewer69_2004_2016_all_modes.mat');
load('/home/kristof/work/brewer/o3data2017.mat');

% merge all data to one table
data=[combined_raw_2004;combined_raw_2005;combined_raw_2006;combined_raw_2007;...
      combined_raw_2008;combined_raw_2009;combined_raw_2010;combined_raw_2011;...
      combined_raw_2012;combined_raw_2013;combined_raw_2014;combined_raw_2015;...
      combined_raw_2016];

data.ColumnO3=data.ColumnO3;
data.StdDevO3=data.StdDevO3;
  
dates=mjd2k_to_date(data.UTC-yeartime(2000));

[ft,years]=fracdate(dates);

months=month(dates);
days=day(dates);


% get direct sun data
ind=find(strcmp(data.ObsCode,'DS'));

brewer_ds=table(data.UTC(ind)-yeartime(2000),years(ind),months(ind),days(ind),ft(ind),...
                data.ColumnO3(ind),data.StdDevO3(ind),data.ZA(ind));

brewer_ds.Properties.VariableNames={'mjd2k','year','month','day','fractional_time',...
                                    'tot_col','tot_col_err','sza'};            
            
% get zenith-sky data
ind=find(strcmp(data.ObsCode,'ZS'));

brewer_zs=table(data.UTC(ind)-yeartime(2000),years(ind),months(ind),days(ind),ft(ind),...
                data.ColumnO3(ind),data.StdDevO3(ind),data.ZA(ind));

brewer_zs.Properties.VariableNames={'mjd2k','year','month','day','fractional_time',...
                                    'tot_col','tot_col_err','sza'};            

                                
% add 2017 DS data
ind=find(strcmp(o3data2017.cType,'DS'));

brewer_ds2017=table(o3data2017.mjd2k(ind),o3data2017.year(ind),o3data2017.month(ind),...
                o3data2017.day(ind),o3data2017.fractional_time(ind),...
                o3data2017.o3(ind),o3data2017.err_o3(ind),o3data2017.sza(ind));

brewer_ds2017.Properties.VariableNames={'mjd2k','year','month','day','fractional_time',...
                                    'tot_col','tot_col_err','sza'};            

                                
% add 2017 ZS data
ind=find(strcmp(o3data2017.cType,'ZS'));

brewer_zs2017=table(o3data2017.mjd2k(ind),o3data2017.year(ind),o3data2017.month(ind),...
                o3data2017.day(ind),o3data2017.fractional_time(ind),...
                o3data2017.o3(ind),o3data2017.err_o3(ind),o3data2017.sza(ind));

brewer_zs2017.Properties.VariableNames={'mjd2k','year','month','day','fractional_time',...
                                    'tot_col','tot_col_err','sza'};            

      
brewer_ds=[brewer_ds; brewer_ds2017];
brewer_zs=[brewer_zs; brewer_zs2017];
                                
save /home/kristof/work/brewer/brewer_tables.mat brewer_ds brewer_zs                                
                                
% end

