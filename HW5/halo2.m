N = 1024^3;

s1 = sum(haloz0);
s2 = sum(haloz1);
%{
figure(1)
[N, edges] = histcounts(haloz0,100)
edges = edges(2:end) - (edges(2)-edges(1))/2;
plot(edges, N)
set(gca,'xscale','log')
set(gca,'yscale','log')

figure(2)
bins = 10.^(13:0.01:18);
histogram(haloz0,bins)
set(gca,'xscale','log')
set(gca,'yscale','log')
%}
m0 = s1/1024^3;
m1 = s2/1024^3;

figure(1)
bins1 = 10.^(12:0.01:18);
[N1,edges1] = histcounts(haloz0,bins1);
edges1 = edges1(2:end) - (edges1(2)-edges1(1))/2;
n1 = zeros(length(N1),1);
for idx=1:length(N1)
    if idx == 1
        n1(idx) = N1(idx)/(edges1(2)-edges1(1))
    elseif idx == length(N1)
        n1(idx) = N1(idx)/(edges1(idx)-edges1(idx-1))
    else
        n1(idx) = N1(idx)/((edges1(idx+1)-edges1(idx-1))/2)
    end
end
plot(edges1, n1,"b",LineWidth=3);
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("M (Msun/h)",FontSize=15);
ylabel("n(M)",FontSize=15);
title("Mass Function(z=0)")

figure(2)
bins1 = 10.^(12:0.01:18);
[N2,edges2] = histcounts(haloz1,bins1);
edges2 = edges2(2:end) - (edges2(2)-edges2(1))/2;
n2 = zeros(length(N2),1);
for idx=1:length(N2)
    if idx == 1
        n2(idx) = N2(idx)/(edges2(2)-edges2(1))
    elseif idx == length(N2)
        n2(idx) = N2(idx)/(edges2(idx)-edges2(idx-1))
    else
        n2(idx) = N2(idx)/((edges2(idx+1)-edges2(idx-1))/2)
    end
end
plot(edges2, n2,"b",LineWidth=3);
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("M (Msun/h)",FontSize=15);
ylabel("n(M)",FontSize=15);
title("Mass Function(z=1)")

% find number of particles below m_crit for z=0
[M1, I1] = max(N1);
numpart1 = 0;
for i=1:I1
    numpart1 = numpart1 + N1(i)*edges1(i);
end
numpart1 = numpart1/m0;

% find number of particles below m_crit for z=1
[M2, I2] = max(N2);
numpart2 = 0;
for i=1:I2
    numpart2 = numpart2 + N2(i)*edges2(i);
end
numpart2 = numpart2/m1;
1024^3


% Compute sigma(M)
% compute the case for z=0 first
deltasq = linearpk(:,1).^3/2/pi^2.*linearpk(:,2);
M_lst = 10.^(12:0.01:18);
Mkmatrix = zeros(length(M_lst), length(linearpk(:,1)));
rhom = 1.9537*10^10;
for idx = 1:length(M_lst)
    for jdx = 1:length(linearpk(:,1))
        k = linearpk(jdx,1);
        M = M_lst(idx);
        P = linearpk(jdx,2);
        Mkmatrix(idx, jdx) = k^2/2/pi^2*P*4*pi*rhom/k^3/M*((sin(k*(3*M/4/pi/rhom)^(1/3))-k*(3*M/4/pi/rhom)^(1/3)*cos(k*(3*M/4/pi/rhom)^(1/3))))^2;
    end
end
sigmasq = zeros(length(M_lst),1);
for idx = 1:length(M_lst)
    temp = 0;
    integrand_lst = Mkmatrix(idx,:);
    for jdx = 1:length(linearpk(:,1))
        p2 = integrand_lst(jdx);
        if jdx == 1
            temp = temp + p2*linearpk(jdx,1);
        else
            p1 = integrand_lst(jdx-1);
            temp = temp + (p1+p2)/2*(linearpk(jdx,1)-linearpk(jdx-1,1));
        end
    end
    sigmasq(idx) = temp;
    temp
end

sigmaM = sigmasq.^(1/2);
figure(3)
semilogx(M_lst, sigmaM,"b",Linewidth=3);
xlabel("M",FontSize=15);
ylabel("sigma(M)",FontSize=15);
title("Sigma(M) z=0")
   
%=========================================================================
% compute the z=1 case
factor2 = Dz_lst(3);
figure(4)
semilogx(M_lst, sigmaM*factor2,"b",Linewidth=3);
xlabel("M",FontSize=15);
ylabel("sigma(M)",FontSize=15);
title("Sigma(M) z=1")

%====
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
    n_ps_0(idx) = -(2/pi)^(1/2)*deltac/sigmavalue^2*rhom/M*exp(-deltac^2/2/sigmavalue^2)*derivative*10^9;

end


figure(5)
mf1 = loglog(M_lst',n_ps_0,"b",LineWidth=3);
mf1name = "PS mass function";
hold on
mf2 = loglog(edges1, N1, "g", LineWidth=3);
mf2name = "Measured mass function";
hold off
legend([mf1,mf2],[mf1name,mf2name]);
xlabel("M",FontSize=15);
ylabel("n(M)",FontSize=15);
title("Mass Function(z=0)")


% compute z=1
n_ps_1 = zeros(length(M_lst),1);
for idx = 1:length(M_lst)
    sigmavalue = sigmaM(idx)*factor2;
    M = M_lst(idx);
    if idx == 1
        derivative = (sigmaM(idx+1)-sigmaM(idx))/(M_lst(idx+1)-M);
    elseif idx == length(M_lst)
        derivative = (sigmaM(idx)-sigmaM(idx-1))/(M-M_lst(idx-1));
    else
        derivative = (sigmaM(idx+1)-sigmaM(idx-1))/(M_lst(idx+1)-M_lst(idx-1));
    end
    derivative
    n_ps_1(idx) = -(2/pi)^(1/2)*deltac/sigmavalue^2*rhom/M*exp(-deltac^2/2/sigmavalue^2)*derivative*10^9;

end


figure(6)
mf3 = loglog(M_lst',n_ps_1,"b",LineWidth=3);
mf3name = "PS mass function";
hold on
mf4 = loglog(edges2, N2, "g", LineWidth=3);
mf4name = "Measured mass function";
hold off
legend([mf3,mf4],[mf3name,mf4name]);
xlabel("M",FontSize=15);
ylabel("n(M)",FontSize=15);
title("Mass Function(z=1)")





