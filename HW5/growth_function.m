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
plot(t, gf,"b",LineWidth=3);
xlabel("ln(a)",FontSize=15);
ylabel("g(a)",FontSize=15);
title("g(a)")
%fplot(@(x)deval(sol,x,1)*x, [-5, 0])