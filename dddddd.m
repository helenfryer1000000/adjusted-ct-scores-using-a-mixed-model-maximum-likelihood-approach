clear all
clc

global log10xbarmatrix

log10xbarmatrix=zeros(4,4,4);
log10xbarmatrix(1,:,:)=-0.01;
log10xbarmatrix(2,:,:)=0;
log10xbarmatrix(3,:,:)=1;
log10xbarmatrix(4,:,:)=2;

% loglikofdataiFUNtry2([0,1,2],[-Inf,0,1],0.7,0.014)
loglikofdataiFUNtry2([0],[-Inf],0.7,0.014)