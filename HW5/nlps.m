Omega_m = 0.3175;
Omega_b = 0.049;
h = 0.6711;
n_s = 0.9624;
Sigma8 = 0.834;



xq = Pkmz0(:,1);
pkz0 = interp1(xq,Pkmz0(:,2),xq,'Spline');
pkz0_5 = interp1(xq,Pkmz0_5(:,2),xq,'Spline');
pkz1 = interp1(xq,Pkmz1(:,2),xq,'Spline');
pkz2 = interp1(xq,Pkmz2(:,2),xq,'Spline');
pklinear = interp1(linearpk(:,1),linearpk(:,2),xq,'Spline');







% power spectrum 
figure(1)
a=loglog(xq,pkz0,"b",LineWidth=3);
aname="z=0";
hold on;
b=loglog(xq,pkz0_5,"r",LineWidth=3);
bname="z=0.5";
c=loglog(xq,pkz1,"y",LineWidth=3);
cname="z=1";
d=loglog(xq,pkz2,"g",LineWidth=3);
dname="z=2";
f=loglog(xq,pklinear,"k",LineWidth=3);
fname="linear power spectrum";
hold off;
xlabel("k (h/Mpc)",FontSize=15);
ylabel("P(k) (Mpc/h)^3",FontSize=15);
legend([a,b,c,d,f],[aname,bname,cname,dname,fname]);
title("Power Spectrum");




%power specturm to linear ps ratio
figure(2)
a1=loglog(xq,pkz0./pklinear,"b",LineWidth=3);
a1name="z=0";
hold on;
b1=loglog(xq,pkz0_5./pklinear,"r",LineWidth=3);
b1name="z=0.5";
c1=loglog(xq,pkz1./pklinear,"y",LineWidth=3);
c1name="z=1";
d1=loglog(xq,pkz2./pklinear,"g",LineWidth=3);
d1name="z=2";
hold off;
xlabel("k (h/Mpc)",FontSize=15);
ylabel("P(k)/P_linear(k)",FontSize=15);
legend([a1,b1,c1,d1],[a1name,b1name,c1name,d1name]);
title("Nonlinear Power Spectrum to Linear Power Spectrum ratio vs z");


%D(z)
figure(3)

%use method of part3 to solve
syms y(t)
[V] = odeToVectorField(diff(y, 2) == -(2.5+1.5*(7/3*exp(3*t)/(7/3*exp(3*t)+1)))*diff(y) - 3*7/3*exp(3*t)/(7/3*exp(3*t)+1)*y)
M = matlabFunction(V,'vars', {'t','Y'})
%sol = ode45(M,[-5 0],[1 0]);
[t,y] = ode45(M,[-5 0],[1 0]);
gf = zeros(length(t),1);
a = zeros(length(t),1);
for i=1:length(t)
    gf(i) = y(i,1)*exp(t(i));
    a(i) = exp(t(i));
    ti=t(i)
end

%extract growth function from given data
pp = xq.^n_s;
Dz0 = pkz0(1)/pp(1);
Dz0_5 = pkz0_5(1)/pp(1);
Dz1 = pkz1(1)/pp(1);
Dz2 = pkz2(1)/pp(1);
%convertion in terms of a for ease of comparison
a_lst = [1 1/(1+0.5) 1/2 1/3];
Dz_lst = [1 Dz0_5/Dz0 Dz1/Dz0 Dz2/Dz0];
a_lst_interp = 0.007:0.001:1;
Da = interp1(a_lst,Dz_lst, a_lst_interp, 'Spline');
Da = Da*gf(53);
ob1 = semilogx(a_lst_interp, Da,"b",LineWidth=3);
ob1name = "extracted growth function result";
hold on
ob2 = semilogx(a, gf,"r",LineWidth=3);
ob2name = "analytic growth function result";
hold off
xlabel("a",FontSize=15);
ylabel("D(a)",FontSize=15);
legend([ob1,ob2],[ob1name,ob2name]);
title("growth function D(a)");






