% density development

a_lst = zeros(length(run1background(:,1)),1);
for idx = 1:length(run1background(:,1))
    a_lst(idx) = 1/(1+run1background(idx,1));
end





rho_CDM = run1background1(:,11);
rho_run1 = run1background(:,11)+run1background(:,17)+run1background(:,22)+run1background(:,27);
rho_run2 = run2background(:,11)+run2background(:,17)+run2background(:,22)+run2background(:,27);
rho_sf = sfbackground(:,11)+sfbackground(:,17)+sfbackground(:,22);
rho_CDM1 = run1background(:,11);
rho_CDM2 = sfbackground(:,11)


figure(1)
a1=semilogx(a_lst,rho_run1./rho_CDM,"b",LineWidth=3);
a1name = "3 component SFDM";
hold on
a2=semilogx(a_lst,rho_sf./rho_CDM,"k",LineWidth=3);
a2name = "2 component SFDM";
a3=semilogx(a_lst,rho_CDM./rho_CDM,"r",LineWidth=3);
a3name = "CDM";
hold off
xlabel("a",FontSize=15);
ylabel("\Omega_{\phi+CDM}/\Omega_{\Lambda CDM}",FontSize=15);
xlim([1e-7 1e-4])
ylim([0.5 1.5])
legend([a1,a2,a3],[a1name,a2name,a3name],FontSize=15);
title("Density Evolution");