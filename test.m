clear all
clc
x=35:-0.01:20;

y=(1-0.014).^exp(x-30);

y=exp((x-20));
y(y>0)=1;

y2=1./(1+exp(-1*(x-27)));
hold on; plot(x,y,'r')
hold on; plot(x,y2,'b')
