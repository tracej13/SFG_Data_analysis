global vis_freq vis_angle ir_angle sf_angle n_1 vis_n ir_n sf_n n_i ir_freq sf_freq r 

r = 2.5;
beta_ratio_asym = 3.84;

vis_angle = 60;
ir_angle = 55;
sf_angle = 59.3;

vis_freq = 18796.99; 
ir_freq = 2880;
sf_freq = vis_freq + ir_freq; 

% Air
n_1 = 1.0003;

        %Gold From online            
% vis_n = 0.59111 + 2.5204i; 
% ir_n = 3.9946 + 24.033i; 
% sf_n = 1.3802 + 2.0345i; 

        %Gold From Baldelli            
% vis_n = 0.402 +2.54i; 
% ir_n = 2.046 +21.3i; 
% sf_n = 1.426 + 1.846i; 



        %copper Babar Weaver 2015
% vis_n = 1.0087 + 2.4408i;
% ir_n =  0.94206 + 24.975i;
% sf_n = 1.1451 + 2.3807i;

        %copper ordal et al. 1985
% vis_n = 1.0854 + 2.7704i;
% ir_n =  1.8663 + 23.172i;
% sf_n = 1.16 + 2.64i;

        %copper oxide
vis_n = 2.65 + 0.046i;
ir_n =  2.65 + 0.046i;
sf_n = 2.65 + 0.046i;



      %PP
% vis_n = 1.49;
% ir_n =  1.49;
% sf_n = 1.49;

%Interface
n_i = 1.65;



B = @angle;
cos_gamma = @refracted_cos;
refract = @refractive_cu;
xx = @L_xx;
yy = @L_yy;
zz = @L_zz;

% A_ppp_sym = 45.13;
% width_ppp_sym = 5.136;
% 
% A_ssp_sym = 8.877;
% width_ssp_sym = 2.755;
% 
% A_48_ppp = (44.8/5.126)+(29.88/5.179);
% A_48_ssp = (14.09/2.857)+(12.29/8.056);
% 
% A_48_2_ppp = (45.53/5.432)+(29.17/5.248);
% A_48_2_ssp = (14.53/2.891)+(11.42/7.673);
% 
% A_48_ratio = (A_48_ppp/A_48_ssp)^2;
% 
% A_48_2_ratio = (A_48_2_ppp/A_48_2_ssp)^2;
% 
% A_49_ppp = (55.85/5.732)+(31.63/6.142);
% A_49_ssp = (7.222/3.887)+(20.31/8.013);
% 
% A_49_ratio = (A_49_ppp/A_49_ssp)^2;
% 
% % A_61_ppp = (41.56/4.5)+(16.36/2.506);
% % A_61_ssp = (12.35/2.875)+(11.39/8.321);
% % 
% % A_61_ratio = (A_61_ppp/A_61_ssp)^2;
% 
% A_61_ppp = (37.52/3.209)+(14.01/2.531);
% A_61_ssp = (12.35/2.875)+(11.39/8.321);
% 
% A_61_ratio = (A_61_ppp/A_61_ssp)^2;
% 
% A_53_ppp = (32.72/4.366)+(12.93/3.495);
% A_53_ssp = (9.078/3.644)+(10.11/6.721);
% A_53_ratio = (A_53_ppp/A_53_ssp)^2;
% 
% A_58_ppp = (14.77/1.009)+(3.75/1);
% A_58_ssp = (2.91/1.001)+(9.067/20.84);
% A_58_ratio = (A_58_ppp/A_58_ssp)^2;
% 
% A_59_ppp = (23.08/2.439)+(11.78/6.112);
% A_59_ssp = (2.996/1)+(2.687/4.756);
% A_59_ratio = (A_59_ppp/A_59_ssp)^2;
% 
% A_60_ppp = (13.14/2.331)+(2.69/2.384);
% A_60_ssp = (2.695/2.49)+(9.529/23.97);
% A_60_ratio = (A_60_ppp/A_60_ssp)^2;
% 
% A_63_ppp = (114.8/10.49)+(172.2/16.94);
% A_63_ssp = (26.9/9.211)+(32.63/14.62);
% A_63_ratio = (A_63_ppp/A_63_ssp)^2;
% 

%A_ratio_avg = (A_49_ratio+A_53_ratio+A_61_ratio)/3;

A_35_1 = coeff_values(33,18);
A_35 = coeff_values(35,18);
A_37_1 = coeff_values(34,18);
A_37 = coeff_values(37,18);
A_39 = coeff_values(39,18);
A_41 = coeff_values(41,18);
A_43 = coeff_values(43,18);
A_48 = coeff_values(48,18);
A_49 = coeff_values(49,18);
A_53 = coeff_values(53,18);
A_54 = coeff_values(54,18);
A_55 = coeff_values(55,18);
A_58 = coeff_values(58,18);
A_59 = coeff_values(59,18);
A_60 = coeff_values(60,18);
A_62 = coeff_values(62,18);
A_69 = coeff_values(69,18);
A_70 = coeff_values(70,18);
A_71 = coeff_values(71,18);
A_72 = coeff_values(72,18);
A_73 = coeff_values(73,18);
A_74 = coeff_values(74,18);
A_75 = coeff_values(75,18);
A_76 = coeff_values(76,18);

amplitudes = zeros(80,length(fit_coeff));





% raman = {[-0.219544 -0.0413986 0.000376351; -0.0413986 -0.203998 -0.00130766; 0.000376351 -0.00130766 -0.195995],[0.000392402 -0.00602168 0.0987599; -0.00602168 -0.00324209 -0.0449072; 0.0987599 -0.0449072 0.00230837],[-0.0576117 -0.150374 -0.00590667; -0.150374 -0.169684 -0.00173303; -0.00590667 -0.00173303 0.0665775]};
% dipole = {[-0.00419403 -0.000140365 -0.000008552],[-0.0000178477 -0.000108388 0.00369652],[-0.0000816781 -0.00209643 -0.000197564]};
%r = alpha(1,1,3,1,2880,raman,dipole)/alpha(3,3,3,1,2880,raman,dipole);

for i = 1:89
     [ssp_sym(i)] = L_yy(sf_freq)*L_yy(vis_freq)*L_zz(ir_freq)*sind(ir_angle)*X_yyz_sym(i);
%     [ssp_asym(i)] = abs(L_yy(sf_freq)*L_yy(vis_freq)*L_zz(ir_freq)*sind(ir_angle)*X_yyz_asym(i,raman,dipole));
%     
      [test_xxz(i)] = X_xxz_sym(i);
      [test_yyz(i)] = X_yyz_sym(i);
      [test_xzx(i)] = X_xzx_sym(i);
      [test_zxx(i)] = X_zxx_sym(i);
      [test_zzz(i)] = X_zzz_sym(i);
 
    [ppp_sym(i)] = (-L_xx(sf_freq)*L_xx(vis_freq)*L_zz(ir_freq)*cosd(sf_angle)* ...
         cosd(vis_angle)*sind(ir_angle)*X_xxz_sym(i)) ...
         +(L_zz(sf_freq)*L_zz(vis_freq)*L_zz(ir_freq)*sind(sf_angle) ...
         *sind(vis_angle)*sind(ir_angle)*X_zzz_sym(i)) ...
         -(L_xx(sf_freq)*L_zz(vis_freq)*L_xx(ir_freq)*cosd(sf_angle) ...
         *sind(vis_angle)*cosd(ir_angle)*X_xzx_sym(i)) ...
         +(L_zz(sf_freq)*L_xx(vis_freq)*L_xx(ir_freq)*sind(sf_angle) ...
         *cosd(vis_angle)*cosd(ir_angle)*X_zxx_sym(i));

     [ppp_asym(i)] = abs((-L_xx(sf_freq)*L_xx(vis_freq)*L_zz(ir_freq)*cosd(sf_angle)* ...
        cosd(vis_angle)*sind(ir_angle)*X_xxz_asym(i)) ...
        -(L_xx(sf_freq)*L_zz(vis_freq)*L_xx(ir_freq)*cosd(sf_angle) ...
        *sind(vis_angle)*cosd(ir_angle)*X_xzx_asym(i)) ...
        +(L_zz(sf_freq)*L_xx(vis_freq)*L_xx(ir_freq)*sind(sf_angle) ...
        *cosd(vis_angle)*cosd(ir_angle)*X_zxx_asym(i)) ...
        +(L_zz(sf_freq)*L_zz(vis_freq)*L_zz(ir_freq)*sind(sf_angle) ...
        *sind(vis_angle)*sind(ir_angle)*X_zzz_asym(i)));
% 
%      ssp_ratio(i) = ssp_sym(i)./ssp_asym(i);
%      ppp_ssp_sym(i) = (ppp_sym(i)./ssp_sym(i));
%      ppp_ssp_asym(i) = ppp_asym(i)./ssp_asym(i);
end

%I=(A/width)^2
ppp_ratio = abs((1/beta_ratio_asym)*(ppp_sym./ppp_asym)).^2;
ppp_ssp_ratio = (abs(ppp_sym)./abs(ssp_sym)).^2;

figure(1)
plot(ppp_ratio)
ylim([0 20])

coverage_samples = [34 33];
sample_angles = [28 24]; 
ppp_ref_measured = coeff_values(43,20);
ppp_ref_calculated = abs(ppp_sym(26))^2;

for i = 1:length(coverage_samples)
    ppp_sample_measured = coeff_values(coverage_samples(i),20);
    ppp_sample_calculated = abs(ppp_sym(sample_angles(i)))^2;
    ratio(i) = sqrt((ppp_sample_measured*ppp_ref_calculated)/(ppp_ref_measured*ppp_sample_calculated));
end




function [Beta] = angle(omega) %correct
    global vis_angle ir_angle sf_angle vis_freq
    if omega==vis_freq
        Beta = vis_angle;
    

    elseif omega>= 2800 && omega<=3000
        Beta = ir_angle;
    else
        Beta = sf_angle;
    end
end

function [n_cu] = refractive_cu(omega) %correct
    global vis_n ir_n sf_n vis_freq
    if omega==vis_freq
        n_cu = vis_n;
    

    elseif omega>= 2800 && omega<=3000
        n_cu = ir_n;
    else
        n_cu = sf_n;
    end
end

function [n_interface] = refractive_interface(omega) %correct
    global vis_n ir_n sf_n vis_freq n_1
    if omega==vis_freq
        n_interface = abs(vis_n+n_1)/2;
    

    elseif omega>= 2800 && omega<=3000
        n_interface = abs(ir_n+n_1)/2;
    else
        n_interface = abs(sf_n+n_1)/2;
    end
end


function [cos_gamma] = refracted_cos(omega) %correct
    global n_1
    
    [n_cu] = refractive_cu(omega);
    
    [Beta] = angle(omega);

    sin_sq = (sind(Beta))^2;

    n_ratio = (n_1/(n_cu))^2;

    cos_gamma = sqrt(1-(n_ratio*sin_sq));
end

function [L_xx] = L_xx(omega) %correct
    global n_1 
    [n_cu] = refractive_cu(omega);
    [Beta] = angle(omega);
    [cos_gamma] = refracted_cos(omega);

    L_xx = (2*n_1*(cos_gamma))/(n_1*(cos_gamma)+(n_cu)*cosd(Beta));
end

function [L_yy] = L_yy(omega) %correct
    global n_1 
    [n_cu] = refractive_cu(omega);
    [Beta] = angle(omega);
    [cos_gamma] = refracted_cos(omega);

    L_yy = (2*n_1*cosd(Beta))/(n_1*cosd(Beta)+((n_cu)*(cos_gamma)));
end

function [L_zz] = L_zz(omega) %correct
    global n_1 n_i
    [n_cu] = refractive_cu(omega);
    [n_interface] = refractive_interface(omega);
    [Beta] = angle(omega);
    [cos_gamma] = refracted_cos(omega);

    L_zz = (2*(n_cu)*cosd(Beta))/(n_1*(cos_gamma)+(n_cu)*cosd(Beta)) *(n_1/n_i)^2;
end
% function [ratio] = ratio_coverage(reference_signal,unknown_signal,calculated_reference,calculated_signal)
%     signal_ratio = unknow
% end

% function [alpha] = alpha(alpha_a,alpha_b,alpha_c,mode,mode_freq,raman,dipole)
%     global dielectric
%     
%     coeff = -1/(2*dielectric*mode_freq);
%     draman_dq = raman{mode}(alpha_a,alpha_b);
%     ddipole_dq = dipole{mode}(alpha_c);
%     alpha = coeff*draman_dq*ddipole_dq;
%     alpha = 1;
% end

function [X_yyz_sym] = X_yyz_sym(theta)
    global r
    %r = alpha(1,1,3,1,2880,raman,dipole)/alpha(3,3,3,1,2880,raman,dipole);
    %r = 1.8;
    X_yyz_sym = 1/2*1*1*((cosd(theta)*(1+r))-((cosd(theta)^3)*(1-r)));
end

function [X_yyz_asym] = X_yyz_asym(theta)
    X_yyz_asym = (-1/2)*1*(cosd(theta)-(cosd(theta)^3));
end

function [X_xxz_sym] = X_xxz_sym(theta)
    global r
    %r = alpha(1,1,3,1,2880,raman,dipole)/alpha(3,3,3,1,2880,raman,dipole);
    %r = 1.8;
    X_xxz_sym = 1/2*1*1*((cosd(theta)*(1+r))-((cosd(theta)^3)*(1-r)));
end

function [X_xzx_sym] = X_xzx_sym(theta)
    global r
    %r = alpha(1,1,3,1,2880,raman,dipole)/alpha(3,3,3,1,2880,raman,dipole);
    %r = 1.8;
    X_xzx_sym = 1/2*1*1*((cosd(theta)-(cosd(theta)^3))*(1-r));
end

function [X_zxx_sym] = X_zxx_sym(theta)
    global r
    %r = alpha(1,1,3,1,2880,raman,dipole)/alpha(3,3,3,1,2880,raman,dipole);
    %r = 1.8;
    X_zxx_sym = 1/2*1*1*((cosd(theta)-(cosd(theta)^3))*(1-r));
end

function [X_zzz_sym] = X_zzz_sym(theta)
    global r
    %r = alpha(1,1,3,1,2880,raman,dipole)/alpha(3,3,3,1,2880,raman,dipole);
    %r = 1.8;
    X_zzz_sym = 1*1*((r*cosd(theta))+((cosd(theta)^3)*(1-r)));
end

function [X_xxz_asym] = X_xxz_asym(theta)
    X_xxz_asym = (-1/2)*1*(cosd(theta)-(cosd(theta)^3));
end

function [X_xzx_asym] = X_xzx_asym(theta)
    X_xzx_asym = (1/2)*(cosd(theta)^3);
end

function [X_zxx_asym] = X_zxx_asym(theta)
    X_zxx_asym = (1/2)*(cosd(theta)^3);
end

function [X_zzz_asym] = X_zzz_asym(theta)
    X_zzz_asym = (cosd(theta)-(cosd(theta)^3));
end


