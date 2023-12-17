figure(1)
a1=loglog(run1pk(:,1),run1pk(:,2),"b",LineWidth=3);
a1name = "R=0";
hold on
a2=loglog(run2pk(:,1),run2pk(:,2),"k",LineWidth=3);
a2name = "R=0.5";
a3=loglog(run3pk(:,1),run3pk(:,2),"r",LineWidth=3);
a3name = "R=1";
hold off
xlabel("k(hMpc^{-1})",FontSize=15);
ylabel("P",FontSize=15);
legend([a1,a2,a3],[a1name,a2name,a3name],FontSize=15);
title("Linear Power Spectrum(m~10^{-22}eV, \lambda=0)");




%{
hold on
loglog(Pk2(:,1),Pk2(:,2))
loglog(Pk3(:,1),Pk3(:,2))
loglog(Pk4(:,1),Pk4(:,2))
loglog(Pk5(:,1),Pk5(:,2))
hold off




figure(2)
plot(CMB(:,1),CMB(:,2))
hold on
plot(CMB2(:,1),CMB2(:,2))
%}








