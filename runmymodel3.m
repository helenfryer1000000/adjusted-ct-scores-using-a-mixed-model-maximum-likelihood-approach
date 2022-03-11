clear all
% clc

%% are we deleting too many time points
%%check issues with minusinfinity

global whatlinktofulldata maxgap
whatlinktofulldata='C:\Users\helenf\Dropbox\ONS_sequencing_restricted_access_shared\ONS_Household_Data\Archived\ONS_2021_07_27\ONS_positives_all_results_household_20210726.csv';


load('themetadata')
% fulldataXXX = readtable(whatlinktofulldata,'VariableNamingRule','modify');

maxgap=40;
WWWsd=6.7; % calistri
Hsd=1.2; %calistri
Aspace=0.25;
%Aspace=0.1
% WWWspace=1;
  WWWspace=1;
HHspace=0.1;
% HHspace=0.01;
minofD=1;
% maxofD=17;
maxofD=100
theta=1.56; % relationship between b117 and not
Hspaceend=7.5; %6.5 %IMPORTANT THIS CHANGES THE ESTIMATES DRAMATICALLY!!!!

[whattvec,whatxvec,whatlogxvec,whatDvec,whatshiftvec,whatzvec,whatageindexvec,incdaysinfz, P_a_t_s_D,verylargest,whatIDvec]=incdaysinfzFUN_patsD(theta,minofD,maxofD,Aspace);
avec=Aspace:Aspace:verylargest;
WWWvec=WWWspace:WWWspace:100; %if this is too small then can get -INf problems!!
Hvec=HHspace:HHspace:Hspaceend;

possiblelog10xb=[-1,0:HHspace:Hspaceend];
possiblelog10xbALT=[-HHspace:HHspace:Hspaceend];
possiblelogxijdataALT=[-HHspace:HHspace:Hspaceend]';

BIGpossiblelog10xb=repmat(possiblelog10xbALT,length(possiblelogxijdataALT),1);
BIGpossiblelogxijdata=repmat(possiblelogxijdataALT,1,length(possiblelog10xbALT));
rrrrr=round(10.^possiblelog10xb);




rrT=0;
rrF=0;
rrD=0;
for k=1:length(whatDvec)
    if whatzvec(k)==1
        rrT=rrT+1;
        whattvecT(rrT)=whattvec(k);
        whatxvecT{rrT}=whatxvec{k};
        whatlogxvecT{rrT}=whatlogxvec{k};
        whatDvecT(rrT)=whatDvec(k);
        whatshiftvecT{rrT}=whatshiftvec{k};
        whatzvecT(rrT)=whatzvec(k);
        whatageindexvecT(rrT)=whatageindexvec(k);
        
    elseif whatzvec(k)==2
        rrF=rrF+1;
        whattvecF(rrF)=whattvec(k);
        whatxvecF{rrF}=whatxvec{k};
        whatlogxvecF{rrF}=whatlogxvec{k};
        whatDvecF(rrF)=whatDvec(k);
        whatshiftvecF{rrF}=whatshiftvec{k};
        whatzvecF(rrF)=whatzvec(k);
        whatageindexvecF(rrF)=whatageindexvec(k);

    elseif whatzvec(k)==3
        rrD=rrD+1;
        whattvecD(rrD)=whattvec(k);
        whatxvecD{rrD}=whatxvec{k};
        whatlogxvecD{rrD}=whatlogxvec{k};
        whatDvecD(rrD)=whatDvec(k);
        whatshiftvecD{rrD}=whatshiftvec{k};
        whatzvecD(rrD)=whatzvec(k);
        whatageindexvecD(rrD)=whatageindexvec(k);
    end
end

dataT={whatxvecT,whatlogxvecT,whatDvecT,whatshiftvecT,whatzvecT,whattvecT,whatageindexvecT};
dataF={whatxvecF,whatlogxvecF,whatDvecF,whatshiftvecF,whatzvecF,whattvecF,whatageindexvecF};
dataD={whatxvecD,whatlogxvecD,whatDvecD,whatshiftvecD,whatzvecD,whattvecD,whatageindexvecD};


% for i=1:length(dataF{1})
%     firstlogx2(i)=whatlogxvecF{i}(1);
% end

% dataF=cell(1,7);
% dataF{1}{1}=[10.^[2,4,6,5.14,4.29,3.43,2.57,1.71,0.86],0.01,0.01,0.01,0.01,0.01,0.01];
% dataF{2}{1}=[2,4,6,5.14,4.29,3.43,2.57,1.71,0.86,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf];
% dataF{3}(1)=1;
% dataF{4}{1}=(2:2:30)-2;
% dataF{5}(1)=1;
% dataF{6}(1)=300;
% dataF{7}(1)=6;


%%fake (20, 6,0.3)
% dataT=cell(1,7);
% dataT{1}{1}=[10.^[2,4,6,5.14,4.29,3.43,2.57,1.71,0.86],0.01,0.01,0.01,0.01,0.01,0.01];
% dataT{2}{1}=[2,4,6,5.14,4.29,3.43,2.57,1.71,0.86,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf];
% dataT{3}(1)=1;
% dataT{4}{1}=(2:2:30)-2;
% dataT{5}(1)=2; % B117 true
% dataT{6}(1)=300; %time
% dataT{7}(1)=6;  %age
% 
% 
% 
% epsilonvec=[0.0001];
% dofTvec=[10]; %10 has a variance of 1.25. options from 3: infnnity (6 has a variance of 1.5)
% agegradvec=0;
% WsplitWWvec=0.3;

% %FALSE (not b117)
% bigwvecF=[20];
% HvecF=[6];%
% % 
%%%%TRUE (=B117)
% bigwvecT=[20];
% HvecT=[5.9,6,6.1];


%%%%SEARCH PARAMETERS

%epsilonvec=[0.004:0.001:0.006]; %0.014
epsilonvec=[0.005];
dofTvec=[10]; %10 has a variance of 1.25. options from 3: infnnity (6 has a variance of 1.5)
%agegradvec=[0.0022:0.0001:0.0024]; %0.06
agegradvec=0.0023;
% WsplitWWvec=[0.15:0.05:0.55];
WsplitWWvec=0.3;

%FALSE (not b117)
%bigwvecF=[21.2:0.1:21.4]; %21
bigwvecF=[21.3];
%HvecF=[4.32:0.01:4.34];
HvecF=[4.33];%

%%%%TRUE (=B117)
%bigwvecT=[20.7:0.1:20.9];
bigwvecT=[20.8];
% HvecT=[4.51:0.01:4.53];
HvecT=[4.52];

%%%%DELTA
%bigwvecD=[20.7:0.1:20.9];
bigwvecD=[20.8];
%HvecD=[4.51:0.01:4.53];
HvecD=[4.52];

for i00=1:length(WsplitWWvec)
    log10xbarmatrix=log10mappedFUN2(avec,WWWvec,Hvec,HHspace,WsplitWWvec(i00));
    for i0=1:length(agegradvec)
        agegradi0=agegradvec(i0);
        for i1=1:length(dofTvec)
            tof0READFROM=tcdf(possiblelog10xb-HHspace,dofTvec(i1)); %check this
            tof0READFROM(1,1)=1;
            for i2=1:length(epsilonvec)
                drdr2=[i0,length(agegradvec),i0/length(agegradvec),i2,length(epsilonvec),i2/length(epsilonvec)]
                pfalsenegREADFROM=binopdf(0,rrrrr,epsilonvec(i2));
                qqqqqREADME=tcdf(BIGpossiblelogxijdata-BIGpossiblelog10xb,dofTvec(i1))-tcdf(BIGpossiblelogxijdata-HHspace-BIGpossiblelog10xb,dofTvec(i1)); %check this
                qqqqqREADME(:,1)=0;
                for i3=1:length(bigwvecF)
                    ii333=[i3,length(bigwvecF)]
                    bigwvecFi3=bigwvecF(i3);%
                    %                     parfor i4=1:length(HvecF)
                    for i4=1:length(HvecF)
                        vvvv1=-loglikofalldata(agegradi0,dataF,bigwvecFi3,WWWsd,HvecF(i4),Hsd,P_a_t_s_D, Aspace,WWWspace,log10xbarmatrix, pfalsenegREADFROM, tof0READFROM, qqqqqREADME, HHspace,avec, WWWvec, Hvec);
                        minusLLF(i00,i0,i1,i2,i3,i4)=sum(vvvv1);
                    end
                end
                for i5=1:length(bigwvecT)
                    ii555=[i5,length(bigwvecT)]
                    bigwvecTi5=bigwvecT(i5);
                    %                    parfor i6=1:length(HvecT)
                    for i6=1:length(HvecT)
                        vvvv2=-loglikofalldata(agegradi0,dataT,bigwvecTi5,WWWsd,HvecT(i6),Hsd,P_a_t_s_D,Aspace,WWWspace,log10xbarmatrix, pfalsenegREADFROM, tof0READFROM, qqqqqREADME, HHspace,avec, WWWvec, Hvec);
                        minusLLT(i00,i0,i1,i2,i5,i6)=sum(vvvv2);
                    end
                end
                for j5=1:length(bigwvecD)
                    jj555=[j5,length(bigwvecD)]
                    bigwvecDj5=bigwvecD(j5);
                    %                    parfor j6=1:length(HvecD)
                    for j6=1:length(HvecD)
                        vvvv3=-loglikofalldata(agegradi0,dataD,bigwvecDj5,WWWsd,HvecD(j6),Hsd,P_a_t_s_D,Aspace,WWWspace,log10xbarmatrix, pfalsenegREADFROM, tof0READFROM, qqqqqREADME, HHspace,avec, WWWvec, Hvec);
                        minusLLD(i00,i0,i1,i2,j5,j6)=sum(vvvv3);
                    end
                end
                
                for i3=1:length(bigwvecF)
                    for i4=1:length(HvecF)
                        for i5=1:length(bigwvecT)
                            for i6=1:length(HvecT)
                                for j5=1:length(bigwvecD)
                                    for j6=1:length(HvecD)
                                        minusLLbig(i00,i0,i1,i2,i3,i4,i5,i6,j5,j6)=minusLLF(i00,i0,i1,i2,i3,i4)+minusLLT(i00,i0,i1,i2,i5,i6)+minusLLD(i00,i0,i1,i2,j5,j6)
                                        
                                        bigWsplitWW(i00,i0,i1,i2,i3,i4,i5,i6,j5,j6)=WsplitWWvec(i00);
                                        bigagegrad(i00,i0,i1,i2,i3,i4,i5,i6,j5,j6)=agegradvec(i0);
                                        bigdofT(i00,i0,i1,i2,i3,i4,i5,i6,j5,j6)=dofTvec(i1);
                                        bigep(i00,i0,i1,i2,i3,i4,i5,i6,j5,j6)=epsilonvec(i2);
                                        bigwF(i00,i0,i1,i2,i3,i4,i5,i6,j5,j6)=bigwvecF(i3);
                                        bigHF(i00,i0,i1,i2,i3,i4,i5,i6,j5,j6)=HvecF(i4);
                                        bigwT(i00,i0,i1,i2,i3,i4,i5,i6,j5,j6)=bigwvecT(i5);
                                        bigHT(i00,i0,i1,i2,i3,i4,i5,i6,j5,j6)=HvecT(i6);
                                        bigwD(i00,i0,i1,i2,i3,i4,i5,i6,j5,j6)=bigwvecD(j5);
                                        bigHD(i00,i0,i1,i2,i3,i4,i5,i6,j5,j6)=HvecD(j6);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

bigWsplitWW=bigWsplitWW(:);
bigagegrad=bigagegrad(:);
bigdofT=bigdofT(:);
bigep=bigep(:);
bigwF=bigwF(:);
bigHF=bigHF(:);
bigwT=bigwT(:);
bigHT=bigHT(:);
bigwD=bigwD(:);
bigHD=bigHD(:);
[whichmin,wheremin]=min(minusLLbig(:));
minWsplitWW=bigWsplitWW(wheremin)
minagegrad=bigagegrad(wheremin)
minep=bigep(wheremin)
mindofT=bigdofT(wheremin)
minwF=bigwF(wheremin)
minHF=bigHF(wheremin)
minwT=bigwT(wheremin)
minHT=bigHT(wheremin)
minwD=bigwD(wheremin)
minHD=bigHD(wheremin)


save('allenviromnent')
%%surface of wT and HT
[X,Y] = meshgrid(bigwvecT,HvecT);
XX=X;
for i=1:width(XX)
    for j=1:height(XX)
        XX(j,i)=minusLLT(1,1,1,1,i,j);
    end
end
figure;
surf(X,Y,-XX)


% %%surface of wT and HT
% [X,Y] = meshgrid(bigwvecT,HvecT);
% XX=X;
% for i=1:width(XX)
%     for j=1:height(XX)
%         XX(j,i)=minusLLT(1,1,i,j);
%     end
% end
% figure;
% surf(X,Y,-XX)

% % figure;
% % surf(X,Y,minusLLF22(:,:,5)')
% xlabel('maximum log10(mapped)')
% zlabel('log likelihood')
% ylabel('infected period (days)')
% 
% 
% figure;
% surf(X,Y,minusLLT22(:,:,3)')
% % figure;
% % surf(X,Y,minusLLT22(:,:,5)')
% 
% 
% [a,b]=find(minusLLT22(:,:,1)==min(minusLLT22(:,:,1),[],'all'));
% VLout=VLmeanvectorF(a) % 7.4
% bigMeanout=bigMMeanvecF(b) %10 %area=37
% 
% [a,b]=find(minusLLF22(:,:,1)==min(minusLLF22(:,:,1),[],'all'));
% VLout=VLmeanvectorF(a) % 6.8
% bigMeanout=bigMMeanvecF(b) %12 %area=40
% 
% 
% % [X,Y] = meshgrid(VLmeanvectorF,sigmavecF);
% % figure;
% % surf(X,Y,minusLL)
% % 
% % minLLatsigmamean=minusLL(:,b);
% % plot(VLmeanvectorF,minLLatsigmamean)
% % 
% % diffinliksminusLLatsigmamean-minusLLatsigmamean(a)
% % log(L(xML))-log(L(x))<1.92
% % 
% % save('LLetc.mat','VLmeanvectorF','sigmavecF','minusLL')
% 
% % shortmeta=metadatatable(categorical(metadatatable.ForVLb117)=='True',:);
% % data=shortmeta.log10mapped;
% % [phatT,pciT]= mle(data,'pdf',@(data,a,sig)mypdf(data,theVLmean,sig),'start',searchstart)
% % 
% % shortmeta=metadatatable(categorical(metadatatable.ForVLb117)=='False',:);
% % data=shortmeta.log10mapped;
% % [phatF,pciF]= mle(data,'pdf',@(data,a,sig)mypdf(data,theVLmean,sig),'start',searchstart)
% 
% 






