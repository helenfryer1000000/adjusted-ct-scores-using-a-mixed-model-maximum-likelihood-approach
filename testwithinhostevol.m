clear all
clc

readindivwith2timepoints = readtable('C:\Users\helenf\Dropbox\ONS_sequencing_restricted_access_shared\Paper_2_Within_Host_Evolution\IndivsWith2Timepoints.csv');

load('d.mat')

totalsnp=sum(d(:,3));

sitetime=d(:,1).*d(:,2);

d=[d,sitetime];

totalrate=totalsnp./d(:,4);

d=[d,totalrate];

% d(1:30,:)

% d(d(:,3)==0,:)=[];

d(d(:,2)<20,:)=[];


d(d(:,3)>0.1,:)=0.12;

figure;
scatter(d(:,5),d(:,2))
xlim([-0.001,0.121])
xlabel('SNP rate')
ylabel('genome coverage')

length(dddd)
is0=sum(dddd==0);
% is1=sum(dddd>0)

jj=histcounts(dddd)';

nums=(0:(length(jj)-1))';

ss=[nums,jj,cumsum(jj)];

plot(nums,cumsum(jj))

over9=sum(jj(10:end));

shortversion=[jj(1:9);over9];
% figure;
% bar(0:9,shortversion)

n=length(dddd);
% p=35/561=0.06

wwww=binopdf(0:100,561,0.000119)

newwww=[wwww(1:9),sum(wwww(10:end))];

figure
hold on; bar(0:9',shortversion','r')
hold on; bar(0:9',561*newwww','k')

figure
hold on; bar(1:9',shortversion(2:10)','r')
hold on; bar(1:9',561*newwww(2:10)','k')