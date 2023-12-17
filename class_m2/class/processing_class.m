

figure(1)
a1=loglog(run9pk(:,1),run9pk(:,2),"b_",LineWidth=3);
a1name = "R=0";
hold on
a2=loglog(run9pk(:,1),run9pk(:,2),"k.",LineWidth=3);
%a2=loglog(sfpk(:,1),sfpk(:,2),"k",LineWidth=3);
a2name = "R=0.5";
a3=loglog(run9pk1(:,1),run9pk1(:,2),"r",LineWidth=3);
a3name = "R=1";
hold off
xlabel("k(hMpc^{-1})",FontSize=15);
ylabel("P",FontSize=15);
xlim([1e-5 1e3])
ylim([1e-20 1e6])
legend([a1,a2,a3],[a1name,a2name,a3name],FontSize=15);
title("Linear Power Spectrum(m~10^{-18}eV, \lambda=-1e2)");

%{
figure(2)
a1=plot(run7cl(:,1),run7cl(:,2),"b_",LineWidth=1);
a1name = "R=0";
hold on
a2=plot(run8cl(:,1),run8cl(:,2),"k.",LineWidth=2);
%a2=plot(sfcl(:,1),sfcl(:,2),"k.",LineWidth=2);
a2name = "R=0.5";
a3=plot(run9cl(:,1),run9cl(:,2),"r",LineWidth=1);
a3name = "R=1";
hold off
xlabel("l",FontSize=15);
ylabel("l(l+1)C_l/2\pi[\muK^2]",FontSize=15);
xlim([0 3000])
ylim([0 1e-9])
legend([a1,a2,a3],[a1name,a2name,a3name],FontSize=15);
title("CMB Power Spectrum(m~10^{-18}eV, \lambda=-1e2)");
%}


%{
figure(3)
a1=loglog(run1pk(:,1),run1pk(:,2)./run1pk1(:,2),"b",LineWidth=3);
a1name = "3 component SFDM";
hold on
%a2=loglog(run2pk(:,1),run2pk(:,2)./run1pk1(1:628,2),"k",LineWidth=3);
a2=loglog(sfpk(:,1),sfpk(:,2)./run1pk1(:,2),"k",LineWidth=3);
a2name = "2 component SFDM";
a3=loglog(run1pk1(:,1),run1pk1(:,2)./run1pk1(:,2),"r",LineWidth=3);
a3name = "CDM";
hold off
xlabel("k(hMpc^{-1})",FontSize=15);
ylabel("P/P_{CDM}",FontSize=15);
legend([a1,a2,a3],[a1name,a2name,a3name],FontSize=15);
title("Linear Power Spectrum(\lambda=0)");


figure(4)
a1=plot(run1cl(:,1),run1cl(:,2)./run1cl1(:,2),"b_",LineWidth=1);
a1name = "3 component SFDM";
hold on
%a2=plot(run2cl(:,1),run2cl(:,2)./run1cl1(:,2),"k.",LineWidth=2);
a2=plot(sfcl(:,1),sfcl(:,2)./run1cl1(:,2),"k.",LineWidth=2);
a2name = "2 component SFDM";
a3=plot(run1cl1(:,1),run1cl1(:,2)./run1cl1(:,2),"r",LineWidth=1);
a3name = "CDM";
hold off
xlabel("l",FontSize=15);
ylabel("C_l/C_{lCDM}",FontSize=15);
legend([a1,a2,a3],[a1name,a2name,a3name],FontSize=15);
title("CMB Power Spectrum(\lambda=0)");


%{
% Compute sigma(M)
% compute the case for z=0 first
deltasq = sfpk(:,1).^3/2/pi^2.*sfpk(:,2);
M_lst = 10.^(12:0.5:18);
Mkmatrix = zeros(length(M_lst), length(sfpk(:,1)));
rhom = 1.9537*10^10;
for idx = 1:length(M_lst)
    for jdx = 1:length(sfpk(:,1))
        k = sfpk(jdx,1);
        M = M_lst(idx);
        P = sfpk(jdx,2);
        Mkmatrix(idx, jdx) = k^2/2/pi^2*P*4*pi*rhom/k^3/M*((sin(k*(3*M/4/pi/rhom)^(1/3))-k*(3*M/4/pi/rhom)^(1/3)*cos(k*(3*M/4/pi/rhom)^(1/3))))^2;
    end
end
sigmasq = zeros(length(M_lst),1);
for idx = 1:length(M_lst)
    temp = 0;
    integrand_lst = Mkmatrix(idx,:);
    for jdx = 1:length(sfpk(:,1))
        p2 = integrand_lst(jdx);
        if jdx == 1
            temp = temp + p2*sfpk(jdx,1);
        else
            p1 = integrand_lst(jdx-1);
            temp = temp + (p1+p2)/2*(sfpk(jdx,1)-sfpk(jdx-1,1));
        end
    end
    sigmasq(idx) = temp;
    temp
end
sigmaM = sigmasq.^(1/2);
%}
% compute PS mass function
% compute z=0

deltac = 1.69;
n_ps_0 = zeros(length(M_lst),1);
for idx = 1:length(M_lst)
    sigmavalue = sigmaM(idx);
    M = M_lst(idx);
    if idx == 1
        derivative = (sigmaM(idx+1)-sigmaM(idx))/(M_lst(idx+1)-M);
    elseif idx == length(M_lst)
        derivative = (sigmaM(idx)-sigmaM(idx-1))/(M-M_lst(idx-1));
    else
        derivative = (sigmaM(idx+1)-sigmaM(idx-1))/(M_lst(idx+1)-M_lst(idx-1));
    end
    derivative
    if idx == 1
        n_ps_0(idx) = -(2/pi)^(1/2)*deltac/sigmavalue^2*rhom/M*exp(-deltac^2/2/sigmavalue^2)*derivative*10^9*(M_lst(idx+1)-M_lst(idx));
    elseif idx == length(M_lst)
        n_ps_0(idx) = -(2/pi)^(1/2)*deltac/sigmavalue^2*rhom/M*exp(-deltac^2/2/sigmavalue^2)*derivative*10^9*(M_lst(idx)-M_lst(idx-1));
    else
        n_ps_0(idx) = -(2/pi)^(1/2)*deltac/sigmavalue^2*rhom/M*exp(-deltac^2/2/sigmavalue^2)*derivative*10^9*((M_lst(idx+1)-M_lst(idx-1))/2);
    end
end


figure(5)
mf1 = loglog(M_lst',n_ps_0CDM./n_ps_0CDM,"b",LineWidth=3);
mf1name = "CDM";
hold on
mf2 = loglog(M_lst',n_ps_03sf./n_ps_0CDM,"r",LineWidth=3);
mf2name = "3-component SFDM";
mf3 = loglog(M_lst',n_ps_02sf./n_ps_0CDM,"k",LineWidth=3);
mf3name = "2-component SFDM";
%mf4 = loglog(M_lst',n_ps_0m18./n_ps_0CDM,"g",LineWidth=3);
%mf4name = "m~10^{-18}eV";
hold off
legend([mf1,mf2,mf3],[mf1name,mf2name,mf3name],FontSize=15);
xlabel("M(M_{sun}/h)",FontSize=15);
ylabel("n(M)/n_{CDM}(M)",FontSize=15);
title("PS Mass Function normalized by CDM MF(z=0)")

%}