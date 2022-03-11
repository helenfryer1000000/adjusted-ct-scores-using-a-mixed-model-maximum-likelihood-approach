function dx=dx(t,x)

global gamma betavec t1vec t2vec

beta=sum((t>=t1vec).*(t<t2vec).*betavec);


S=x(1);
I=x(2);

dS=-beta*S*I;
dI=beta*S*I-gamma*I;
dR=gamma*I;
INC=beta*S*I;

dx=[dS;dI;dR;INC];

end