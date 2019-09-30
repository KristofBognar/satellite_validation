
% % % % altitude grid
% % % z_des=[0.795, [1.5:1:59.5]];
% % % if z_des(1) < z_des(length(z_des))
% % %     z_des = fliplr(z_des);
% % % end
% % % 
% % % % list of sza
% % % % sza_des=[74:2:84];
% % % sza_des=ones(5,1)*74;
% % % 
% % % % list of scattering heights
% % % % z_scat=ones(size(sza_des))*9.8;
% % % z_scat=[5.8:9.8];
% % % 
% % % % radius of the Earth in km
% % % R_e = 6378.1;
% % % 
% % % for i = 1:length(sza_des)
% % %     for j = 1:length(z_des)
% % %         if z_des(j) < z_scat(i)
% % %             range(i,j) = 0;
% % %         else
% % %             range(i,j) = sza_des(i) - asind( (R_e+z_scat(i))./(R_e+z_des(j)) .* sind(sza_des(i)));
% % %         end
% % %     end
% % % end
% % % 
% % % range_km = range * (pi/180) * R_e;
% % % 
% % % label=[];
% % % for i=1:length(sza_des)
% % %     plot(range_km(i,:), z_des, 'x-'), hold on
% % % %     label=[label;num2str(sza_des(i))];
% % %     label=[label;num2str(z_scat(i))];
% % % end
% % % 
% % % legend(label)




