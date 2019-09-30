function smooth_nosmooth_diff(  )
%SMOOTH_NOSMOOTH_DIFF Summary of this function goes here
%   Detailed explanation goes here


load('/home/kristof/work/satellite_validation/all_data_nosmooth.mat');

a_ns_o3 = find_coincidences_twilight( ace_fts_o3, ut_o3, false );
a_ns_no2 = find_coincidences_twilight( ace_fts_no2, ut_no2, false );
a_ns_no2uv = find_coincidences_twilight( ace_fts_no2uv, gbs_no2uv, false );

o_ns_o3 = find_coincidences_time( osiris_o3, ut_o3, 12,false );
o_ns_no2 = find_coincidences_time( osiris_no2, ut_no2, 12,false );
o_ns_no2uv = find_coincidences_time( osiris_no2uv, gbs_no2uv, 12,false );

[a_ns_o3bk,~,~,~,times_ns_a] = find_coincidences_time( ace_fts_o3_bk, bruker_o3, 12,false );
[o_ns_o3bk,~,~,~,times_ns_o] = find_coincidences_time( osiris_o3_bk, bruker_o3, 12,false );

[bk_ns_o3,~,~,~,times_ns_bk] = find_coincidences_time( bruker_o3_doas, ut_o3, 12,false );

load('/home/kristof/work/satellite_validation/all_data_smooth.mat');

a_s_o3 = find_coincidences_twilight( ace_fts_o3, ut_o3, false );
a_s_no2 = find_coincidences_twilight( ace_fts_no2, ut_no2, false );
a_s_no2uv = find_coincidences_twilight( ace_fts_no2uv, gbs_no2uv, false );

o_s_o3 = find_coincidences_time( osiris_o3, ut_o3, 12,false );
o_s_no2 = find_coincidences_time( osiris_no2, ut_no2, 12,false );
o_s_no2uv = find_coincidences_time( osiris_no2uv, gbs_no2uv, 12,false );

[a_s_o3bk,~,~,~,times_s_a] = find_coincidences_time( ace_fts_o3_bk, bruker_o3, 12,false );
[o_s_o3bk,~,~,~,times_s_o] = find_coincidences_time( osiris_o3_bk, bruker_o3, 12,false );

[bk_s_o3,~,~,~,times_s_bk] = find_coincidences_time( bruker_o3_doas, ut_o3, 12,false );

fprintf('\n')

disp('DOAS O3')

mm=mean((a_s_o3./a_ns_o3)-1)*100;
ss=std((a_s_o3./a_ns_o3)-1)*100;
fprintf('ACEFTS: %.1f +/- %.1f\n', [mm,ss])

mm=mean((o_s_o3./o_ns_o3)-1)*100;
ss=std((o_s_o3./o_ns_o3)-1)*100;
fprintf('OSIRIS: %.1f +/- %.1f\n', [mm,ss])

fprintf('\n')

disp('DOAS NO2')

mm=mean((a_s_no2./a_ns_no2)-1)*100;
ss=std((a_s_no2./a_ns_no2)-1)*100;
fprintf('ACEFTS: %.1f +/- %.1f\n', [mm,ss])

mm=mean((o_s_no2./o_ns_no2)-1)*100;
ss=std((o_s_no2./o_ns_no2)-1)*100;
fprintf('OSIRIS: %.1f +/- %.1f\n', [mm,ss])

fprintf('\n')

disp('DOAS NO2 UV')

mm=mean((a_s_no2uv./a_ns_no2uv)-1)*100;
ss=std((a_s_no2uv./a_ns_no2uv)-1)*100;
fprintf('ACEFTS: %.1f +/- %.1f\n', [mm,ss])

mm=mean((o_s_no2uv./o_ns_no2uv)-1)*100;
ss=std((o_s_no2uv./o_ns_no2uv)-1)*100;
fprintf('OSIRIS: %.1f +/- %.1f\n', [mm,ss])

fprintf('\n')

disp('BRUKER O3')

mm=mean((a_s_o3bk./a_ns_o3bk)-1)*100;
ss=std((a_s_o3bk./a_ns_o3bk)-1)*100;
fprintf('ACEFTS: %.1f +/- %.1f\n', [mm,ss])

[~,ia,ib]=intersect(times_ns_o,times_s_o);
mm=mean((o_s_o3bk(ib)./o_ns_o3bk(ia))-1)*100;
ss=std((o_s_o3bk(ib)./o_ns_o3bk(ia))-1)*100;
fprintf('OSIRIS: %.1f +/- %.1f\n', [mm,ss])

fprintf('\n')

disp('BRUKER O3 smoothed with DOAS AVK')

[~,ia,ib]=intersect(times_ns_bk,times_s_bk);
mm=mean((bk_s_o3(ib)./bk_ns_o3(ia))-1)*100;
ss=std((bk_s_o3(ib)./bk_ns_o3(ia))-1)*100;
fprintf('BRUKER: %.1f +/- %.1f\n', [mm,ss])

end

