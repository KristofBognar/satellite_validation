function merge_gbs_cristen()
%MERGE_GBS_CRISTEN 

load('/home/kristof/work/GBS/VCD_results/June2011_Reanalysis_All.mat')

% merge datasets
no2_uv=merge(vcd_p0_no2_ndacc_filt,vcd_u0_no2_ndacc_filt);

no2_tmp=merge(vcd_u1_no2_ndacc_filt,vcd_u2_no2_ndacc_filt);
no2=merge(no2_tmp,vcd_p1_no2_ndacc_filt);

o3_tmp=merge(vcd_u1_o3_ndacc_filt,vcd_u2_o3_ndacc_filt);
o3=merge(o3_tmp,vcd_p1_o3_ndacc_filt);

% convert to table
cristen_no2uv=create_table(no2_uv);
cristen_no2=create_table(no2);
cristen_o3=create_table(o3);

% save
save('/home/kristof/work/GBS/VCD_results/Cristen_merged.mat',...
     'cristen_no2uv','cristen_no2','cristen_o3')


end

function out = create_table(arr_in)

    tmp=array2table(arr_in);
    
    tmp.Properties.VariableNames={'year' 'day' 'ampm' 'fd_min' 'fd_max' 'fd'...
                                  'sza_min' 'sza_max' 'sza' 'saa_min' 'saa_max' 'saa'...
                                  'mean_vcd' 'sigma_w_vcd' 'std_vcd' 'sigma_mean_vcd'};
                              
    tmp.mjd2k=ft_to_mjd2k(tmp.fd-1,tmp.year);
    tmp.fractional_time=tmp.fd-1;
    
    out=tmp;

    
end

function out = merge(t1,t2)

    t1=t1(:,1:16);
    t2=t2(:,1:16);

    % find matching twilights
    [~,ind1,ind2]=intersect(t1(:,1:3),t2(:,1:3),'rows');
    
    % average matching values 
    out=(t1(ind1,:)+t2(ind2,:))/2;
    
    % replace error average with quadrature
    out(:,15)=sqrt( t1(ind1,15).^2 +t2(ind2,15).^2 )/2; % sigma_mean_vcd
    out(:,16)=sqrt( t1(ind1,16).^2 +t2(ind2,16).^2 )/2; % std_vcd
%     out(:,18)=sqrt( t1(ind1,18).^2 +t2(ind2,18).^2 )/2; % langley_vcd_err
    
    
    % add rest of data
    ind12 = setdiff(1:length(t1), ind1);
    ind22 = setdiff(1:length(t2), ind2);
    
    out=[out; t1(ind12,:); t2(ind22,:)];
    
    % sort by time
    out=sortrows(out,[1,2,3]);
    
end