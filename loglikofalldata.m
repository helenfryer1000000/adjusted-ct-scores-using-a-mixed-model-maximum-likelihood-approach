function loglikofalldata=loglikofalldata(agegrad,dataX,WWWmean,WWWsd,Hmean,Hsd,P_a_t_s_D, Aspace, WWWspace,log10xbarmatrix, pfalsenegREADFROM, normof0READFROM, qqqqqREADME, HHspace,avec, WWWvec, Hvec)

whatlogxvec=dataX{2};
whatDvec=dataX{3};
whatshiftvec=dataX{4};
whatzvec=dataX{5};
whattvec=dataX{6};
whatageindexvec=dataX{7};

normWWW=normcdf(WWWvec,WWWmean,WWWsd)-normcdf(WWWvec-WWWspace,WWWmean,WWWsd);
thebignormWWW=zeros(length(avec),length(WWWvec),length(Hvec));

e1age=[3,7,15,25,35,45,55,65,75,85];
e2age=e1age-45;

normH=cell(length(e1age),1);
thebignormH=cell(length(e1age),1);

for e1=1:length(e1age)
    HmeanNEW=Hmean+agegrad*e2age(e1);    
    normH{e1}=normcdf(Hvec,HmeanNEW,Hsd)-normcdf(Hvec-HHspace,HmeanNEW,Hsd);
    thebignormH{e1}=zeros(length(avec),length(WWWvec),length(Hvec));
end

for i1=1:length(avec)
    for i2=1:length(WWWvec)
        for i3=1:length(Hvec)
            thebignormWWW(i1,i2,i3)=normWWW(i2);
            for e1=1:length(e1age)               
                thebignormH{e1}(i1,i2,i3)=normH{e1}(i3);
            end
        end
    end
end

loglikofdata=zeros(length(whatDvec),1);

for i=1:length(whatDvec)
%     i
    D=whatDvec(i);
    z=whatzvec(i);
    t=whattvec(i);
    ageindex=whatageindexvec(i);
    shift=whatshiftvec{i};
    logxijdata=whatlogxvec{i};
        
    loglikofdata(i)=loglikofdataiFUNtry2(D,shift,logxijdata,thebignormWWW,thebignormH{ageindex},P_a_t_s_D{D,z}(t/Aspace,:),Aspace,log10xbarmatrix, pfalsenegREADFROM,normof0READFROM,qqqqqREADME, HHspace);
end
loglikofalldata=loglikofdata;




    
    
            

