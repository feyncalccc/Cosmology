Omega_m = 0.3175;
Omega_b = 0.049;
h = 0.6711;
n_s = 0.9624;
Sigma8 = 0.834;


%Plot the linear power spectrum
figure(1)
a=loglog(linearpk(:,1),linearpk(:,2),"b",LineWidth=3);
aname="linear power spectrum";
xlabel("k (h/Mpc)",FontSize=15);
ylabel("P(k) (Mpc/h)^3",FontSize=15);
title("Linear Power Spectrum");


%Verify the large scale behavior
largescale = linearpk(1:400,1).^n_s;
largescale2 = linearpk(1:400,1).^1;
factor=linearpk(1,2)/largescale(1);
factor2=linearpk(1,2)/largescale2(1);
largescale = factor*largescale;
largescale2 = factor2*largescale2;
figure(2)
a=loglog(linearpk(:,1),linearpk(:,2),"b",LineWidth=3);
aname="linear power spectrum";
hold on;
b=loglog(linearpk(1:400,1),largescale,"r",LineWidth=3);
bname="large scale power law";
c=loglog(linearpk(1:400,1),largescale2,"g",LineWidth=3);
cname="Harrison-Zel'dovich";
hold off;
xlabel("k (h/Mpc)",FontSize=15);
ylabel("P(k) (Mpc/h)^3",FontSize=15);
legend([a,b,c],[aname,bname,cname]);
title("Linear Power Spectrum");

[M,I] = max(linearpk(:,2));
k_eq = linearpk(369,1);

%Compute the transfer function
P_i = linearpk(:,1).^n_s;
T_k = (linearpk(:,2)./P_i).^(1/2);
T_k = T_k/T_k(1);
figure(3)
loglog(linearpk(:,1),T_k,"b",LineWidth=3);
xlabel("k (h/Mpc)",FontSize=15);
ylabel("T(k)",FontSize=15);
title("Transfer Function");






%Compute sigma8 using tophat-filter
new_sigma_8 = 0;
for i=1:884
    if i == 1
        step = linearpk(1,1);
    else
        step = linearpk(i,1)-linearpk(i-1,1);
    end
    k = linearpk(i,1);
    temp = (3*besselj(1,k*8)/8/k)^2*k^2*linearpk(i,2)*step;
    new_sigma_8 = new_sigma_8+temp;
end

new_sigma_8 = sqrt(new_sigma_8/2/pi^2);