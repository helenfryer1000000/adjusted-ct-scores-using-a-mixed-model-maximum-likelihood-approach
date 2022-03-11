function [whattvec,whatxvec,whatlogxvec,whatDvec,whatshiftvec,whatzvec,whatageindexvec,incdaysinfzFUN, P_a_t_s_D,verylargest,whatIDvec]=incdaysinfzFUN_patsD(theta,minofD,maxofD,Aspace)

global whatlinktofulldata maxgap 

load('themetadata')

temptable=readtable('newdata030821.xlsx');

% metadatatable(isnan(metadatatable.mapped),:)=[];
% metadatatable(isnan(metadatatable.log10mapped),:)=[];
metadatatable(isnan(metadatatable.ct_mean),:)=[];
metadatatable(~isfinite(metadatatable.ct_mean),:)=[];


metadatatable(datenum(metadatatable.collection_date)<datenum(min(temptable.dte)),:)=[]; 
metadatatable(datenum(metadatatable.collection_date)>datenum(max(temptable.dte)),:)=[]; 


timeFROMlastNEG=datenum(metadatatable.collection_date)-datenum(metadatatable.last_negative_date);


figure;
histogram(metadatatable.collection_date)

metadatatable(~isfinite(metadatatable.collection_date),:)=[];
metadatatable(~isfinite(metadatatable.timebetweensamples),:)=[];
metadatatable=[metadatatable, array2table(datenum(metadatatable.collection_date)-datenum(min(temptable.dte)),'VariableNames',{'timefromfirst'})];
metadatatable(isnan(metadatatable.timefromfirst),:)=[];
% metadatatable(metadatatable.protocol~='veSeq',:)=[];

metadatatable.VaccDate=datetime(metadatatable.covid_vaccine_date_1(:),'InputFormat','ddMMMyyyy');
vaccbeforpositive=metadatatable.VaccDate<=metadatatable.collection_date;
metadatatable=[metadatatable, array2table(vaccbeforpositive,'VariableNames',{'vaccinebeforepositive'})];

serobeforepositive=metadatatable.first_seroPos_date<=metadatatable.collection_date;
metadatatable=[metadatatable, array2table(serobeforepositive,'VariableNames',{'serobeforepositive'})];

verylargest=max(datenum(temptable.dte))-min(datenum(metadatatable.last_negative_date))+20; %big of a fudge

uend=max(metadatatable.timefromfirst)/Aspace+1;

%%%%%%
ForVLb117=metadatatable.B_1_1_7;
ForVLb117(:)='_';
ForVLb117(categorical(metadatatable.B_1_1_7)=='True')='True';
ForVLb117(categorical(metadatatable.B_1_1_7)=='False')='False';

ForVLb117(categorical(metadatatable.B_1_1_318)=='True')='otherVOC';
ForVLb117(categorical(metadatatable.B_1_351)=='True')=	'otherVOC';
ForVLb117(categorical(metadatatable.A_23_1)=='True')='otherVOC';
ForVLb117(categorical(metadatatable.B_1_324_1)=='True')='otherVOC';
ForVLb117(categorical(metadatatable.AV_1)=='True')='otherVOC';
ForVLb117(categorical(metadatatable.B_1_525)=='True')='otherVOC';
ForVLb117(categorical(metadatatable.B_1_617_1)=='True')='otherVOC';
ForVLb117(categorical(metadatatable.P_2)=='True')='otherVOC';
ForVLb117(categorical(metadatatable.P_3)=='True')='otherVOC';
ForVLb117(categorical(metadatatable.C_36_3)=='True')='otherVOC';
ForVLb117(categorical(metadatatable.P_1)=='True')='otherVOC';
ForVLb117(categorical(metadatatable.B_1_617_3)=='True')='otherVOC';

ForVLb117(categorical(metadatatable.B_1_617_2)=='True')='DELTA';


metadatatable=[metadatatable,array2table(ForVLb117)];

% ssssss=metadatatable(categorical(metadatatable.ForVLb117)=='DELTA',:)
% histogram(ssssss.timebetweensamples)


metadatatable(find((categorical(metadatatable.ForVLb117)~='True').*(categorical(metadatatable.ForVLb117)~='False').*(categorical(metadatatable.ForVLb117)~='DELTA')),:)=[];
% head(metadatatable)

%%%use to convert between ct and log10mapped
metaforlinefit=metadatatable;
metaforlinefit(isnan(metadatatable.ct_mean),:)=[];
metaforlinefit(isnan(metadatatable.log10mapped),:)=[];

ccc = polyfit(metaforlinefit.ct_mean,metaforlinefit.log10mapped,1)

% hold on; scatter(metaforlinefit.ct_mean,metaforlinefit.log10mapped)

CTtolog10mapped=polyval(ccc,metadatatable.ct_mean); %%ORIGINAL
CTtolog10mapped<0==0;

CTtomapped=ceil(10.^CTtolog10mapped);

metadatatable=[metadatatable,array2table(CTtolog10mapped),array2table(CTtomapped)];

min(CTtomapped)
head(metadatatable)

min(metadatatable.log10mapped)
min(metadatatable.mapped)

y_est = polyval(ccc,22.3);
y_est2 = polyval(ccc,22.3+4.2);
y_est3 = polyval(ccc,22.3-4.2);
y_est4 = polyval(ccc,40);
y_est5 = polyval(ccc,19);
y_est6 = polyval(ccc,20.2);
y_est7 = polyval(ccc,19.6);
y_est8 = polyval(ccc,9.6);
y_est9 = polyval(ccc,19.9);
y_est9 = polyval(ccc,10.4);
y_est10 = polyval(ccc,21.4)
y_est11 =polyval(ccc,19.0)
y_est12 = polyval(ccc,22.0)
y_est13 =polyval(ccc,15.8)
yest_14 = polyval(ccc,10)

y_estXXXX = polyval(ccc,38)


head(temptable)

temptable(end,:)=[];
DATES=temptable.dte;
PCPOS=temptable.pcpositive;
PCPOSB117=PCPOS.*(temptable.fracB117);
PCPOSDELTA=PCPOS.*(temptable.fracDELTA);
PCPOSnotB117=PCPOS.*(1-temptable.fracB117-temptable.fracDELTA);

RTfulldates=DATES(1):DATES(end);
PCPOSB117long=NaN(size(RTfulldates));
PCPOSDELTAlong=NaN(size(RTfulldates));
PCPOSnotB117long=NaN(size(RTfulldates));
for i=1:length(RTfulldates)
     for j=2:length(DATES)
        if (RTfulldates(i)>=DATES(j-1))&&(RTfulldates(i)<=DATES(j))
%             i            
               PCPOSB117long(i)=interp1(DATES([j-1,j]),PCPOSB117([j-1,j]),   RTfulldates(i));
            PCPOSnotB117long(i)=interp1(DATES([j-1,j]),PCPOSnotB117([j-1,j]),RTfulldates(i));
              PCPOSDELTAlong(i)=interp1(DATES([j-1,j]),PCPOSDELTA([j-1,j]),  RTfulldates(i));
            
        end
    end
end

BIGRTDATAtab=array2table([RTfulldates'],'VariableNames',{'date'});
BIGRTDATAtab=[BIGRTDATAtab,array2table([PCPOSB117long',PCPOSnotB117long',PCPOSDELTAlong'],'VariableNames',{'prevb117','prevnotB117','prevDELTA'})];

BIGRTDATAtab(isnan(BIGRTDATAtab.prevb117),:)=[];

figure(22); 
hold on; plot(DATES,PCPOSB117,'ro')
hold on; plot(DATES,PCPOSnotB117,'bo')
hold on;  plot(DATES,PCPOSDELTA,'ko')

global t1vec t2vec startx gamma

S0=10000;
gamma=1/15;
startx=[S0,20,0,0];

plotthing1={'r-','b-','k-'};
plotthing2={'r:','b:','k:'};
plotthing3={'rs--','bs--','ks--'};

whatkey={'r-','b-','k-'};

for z=1:3
    RTtimes=datenum(BIGRTDATAtab.date)-datenum(BIGRTDATAtab.date(1));
    RTdates=BIGRTDATAtab.date;  
    global betavec   
    
    t1vec=RTtimes;
    t1vec(length(t1vec))=[];
    t2vec=RTtimes;
    t2vec(1)=[];
    t2vec(length(t2vec))=t2vec(length(t2vec))+1;
    ThinRTtimes=RTtimes(1):Aspace:RTtimes(end);
    ThinRTdates=RTdates(1)+ThinRTtimes;
    
    if z==1 %%b117=alpha
        Istart=BIGRTDATAtab.prevb117;
        Istart(end)=[];
        Iend=BIGRTDATAtab.prevb117;
        Iend(1)=[];
        betavec=(log(Iend./Istart)+gamma)/S0;        
        betavec(isnan(betavec))=0;
        betavec(betavec==Inf)=0;
        betavec(betavec<0)=0.00000001;
%         betavec(betavec<0)=0.00001;
        startx=[S0,0,0,0];
%         TTTTTTTTTTTTTTTTTTTTT=size(ThinRTtimes)
        [TT1,ZZ1]=ode45(@dx,ThinRTtimes(1:700),startx);
        startx=[S0,0.49,0,0];
        [TT2,ZZ2]=ode45(@dx,ThinRTtimes(701:end),startx);
        TT=[TT1;TT2];
        ZZ=[ZZ1;ZZ2];
       
    elseif z==2 %NOT VOC
        Istart=BIGRTDATAtab.prevnotB117;
        Istart(end)=[];
        Iend=BIGRTDATAtab.prevnotB117;
        Iend(1)=[];
        betavec=(log(Iend./Istart)+gamma)/S0;
        betavec(isnan(betavec))=0;
        betavec(betavec==Inf)=0;
        betavec(betavec<0)=0.00000001;    
        startx=[S0,14.25,0,0];
        [TT,ZZ]=ode45(@dx,ThinRTtimes,startx);
        
    elseif z==3 % DELTA
        Istart=BIGRTDATAtab.prevDELTA;
        Istart(end)=[];
        Iend=BIGRTDATAtab.prevDELTA;
        Iend(1)=[];
        betavec=(log(Iend./Istart)+gamma)/S0;        
        betavec(isnan(betavec))=0;
        betavec(betavec==Inf)=0;
        betavec(betavec<0)=0.00000001;
        startx=[S0,0,0,0];
        [TT1,ZZ1]=ode45(@dx,ThinRTtimes(1:1250),startx);
        startx=[S0,0.0285,0,0];
        [TT2,ZZ2]=ode45(@dx,ThinRTtimes(1251:end),startx);
        TT=[TT1;TT2];
        ZZ=[ZZ1;ZZ2];
    end
    
    
    
    SOUT=ZZ(:,1);
    IOUT=ZZ(:,2);
    ROUT=ZZ(:,3);
    INC=ZZ(:,4);
    
     figure(22);
     hold on; plot(ThinRTdates,IOUT/S0,whatkey{z})
%      
%      figure(23);
%      hold on; plot(ThinRTdates,IOUT)
%      
%      figure(24);
%      hold on; plot(ThinRTdates,IOUT)
%      
%      figure(25);
%      hold on; plot(ThinRTdates,INC)
    
    INCstart=INC;
    INCstart(end)=[];
    INCend=INC;
    INCend(1)=[];
    INCOUT=INCend-INCstart;
    INCOUT(INCOUT<0)=0.000000001;
    
%     [INCstart,INCend,INCOUT]
    
    Iprev=IOUT./(SOUT+IOUT+ROUT);
    TTshort=TT;
    TTshort(1)=[];
    
    uvector=1:uend;
    INCDAYSINF=zeros(length(INCOUT),length(uvector));
    for v=1:length(INCOUT)
        for u=uvector
            if v>=u
                INCDAYSINF(v,u)=INCOUT(v-u+1);
            end
        end
    end   
    incdaysinfzFUN{z}=INCDAYSINF;
end

Dvec=minofD:maxofD; %5:3
P_a_t_s_D=cell(length(Dvec),2);

for D=Dvec
    AA=1:(D/Aspace+1);
    for s=1:3
         step1=incdaysinfzFUN{s};
         step2=incdaysinfzFUN{s}./sum(incdaysinfzFUN{s}(:,AA),2);
         step2(step1==0)=0;
          P_a_t_s_D{D,s}=step2;
    end
end

timefromfirstcoll=datenum(metadatatable.collection_date)-datenum(RTdates(1));

metadatatable=[metadatatable,array2table(timefromfirstcoll,'VariableNames',{'timefromfirstcoll'})];

[N,X]=hist(metadatatable.person_id);
tableofpersonidcounts2=[array2table(N','VariableNames',{'counts2'}),cell2table(X','VariableNames',{'person_id2'})];

fulldataXXX = readtable(whatlinktofulldata,'VariableNamingRule','modify');

% head(fulldataXXX)
fulldataXXX.pseudo_person_id=categorical(fulldataXXX.pseudo_person_id);
fulldataXXX.result_mk=categorical(fulldataXXX.result_mk);

whatDvec=[];
whatzvec=[];
whattvec=[];

head(metadatatable)

head(fulldataXXX)

summary(categorical(metadatatable.ForVLb117))

yy=1;
for k=2:height(tableofpersonidcounts2)
    ssss=tableofpersonidcounts2(k,:);
    
    if ssss.counts2>0
        persontable=metadatatable(metadatatable.person_id==tableofpersonidcounts2.person_id2(k),:);
        persontable=sortrows(persontable,{'collection_date'});
        
        persontable(datenum(persontable.collection_date)>datenum(max(temptable.dte)),:)=[];
        
        
        if height(persontable)>1
            for r=[height(persontable):-1:2]
                if persontable.collection_date(r)==persontable.collection_date(r-1)
                    persontable(r,:)=[];
                end
            end
        end
        
        if ~isempty(persontable)
            if persontable.timebetweensamples(1)<=maxofD
                if sum(persontable.vaccinebeforepositive(1))==0
                    if sum(persontable.serobeforepositive(1))==0
                        if (persontable.source_age(1)>=3)&&(persontable.source_age(1)<=85)
                            
                            whattvec(yy)=persontable.timefromfirstcoll(1);  %need to
                            whatDvec(yy)=persontable.timebetweensamples(1); %%something like this
                            
                            whatageindexvec(yy)=find([3,7,15,25,35,45,55,65,75,85]==persontable.source_age(1));
                            
                            
                            if categorical(persontable.ForVLb117(1))=='True'
                                whatzvec(yy)= 1;
                            elseif categorical(persontable.ForVLb117(1))=='False'
                                whatzvec(yy)=2;
                            elseif categorical(persontable.ForVLb117(1))=='DELTA'
                                whatzvec(yy)=3;
                            else
                                whatzvec(yy)=0;
                            end
                            
                            WPtable=fulldataXXX(fulldataXXX.pseudo_person_id==tableofpersonidcounts2.person_id2(k),:);
                            WPtable=sortrows(WPtable,{'collection_date'});
                            for hy=height(WPtable):-1:1
                                if (WPtable.result_mk(hy)~="Positive")&&(WPtable.result_mk(hy)~="Negative")
                                    WPtable(hy,:)=[];
                                end
                            end
                            
                            
                            WPpositive=WPtable(WPtable.result_mk=="Positive",:);
                            firstpositive=min(WPpositive.collection_date);
                            zzz=WPtable(WPtable.collection_date>=firstpositive,:);
                      
                            if height(zzz)>2
                                for j=height(zzz):-1:3
                                    if (zzz.result_mk(j)=="Positive")&&(zzz.result_mk(j-1)=="Negative")&&(zzz.collection_date(j)-zzz.collection_date(j-1)>maxgap)
                                        zzz(j:end,:)=[];
                                    end
                                end
                            end

                            if height(zzz)>2
                                for j=height(zzz):-1:3
                                    if (zzz.result_mk(j)=="Positive")&&(zzz.result_mk(j-1)=="Negative")&&(zzz.collection_date(j)-zzz.collection_date(j-2)>maxgap)
                                        zzz(j:end,:)=[];
                                    end
                                end
                            end    
                            if height(zzz)>1
                                for j=height(zzz):-1:2
                                    if (zzz.result_mk(j)=="Positive")&&(zzz.result_mk(j-1)=="Negative")&&(zzz.collection_date(j)-zzz.collection_date(j-1)>maxgap)
                                        zzz(j:end,:)=[];
                                    end
                                end
                            end
 
                            
                            whatdatesXXX=zzz.collection_date;
                            timebetweensamplesXXX=datenum(zzz.collection_date)-datenum(persontable.last_negative_date(1));
                            
                            whatxvecXXX=zeros(height(zzz),1);
                            for mm=height(zzz):-1:1
                                if zzz.result_mk(mm)=='Negative'
                                    whatxvecXXX(mm)=0;
                                elseif zzz.result_mk(mm)=='Positive'
                                    correctdate=persontable.collection_date==whatdatesXXX(mm);
                                    
                                    if  sum(correctdate)>0
                                        %  whatxvecXXX(mm)=persontable.mapped(correctdate);
                                        whatxvecXXX(mm)=persontable.CTtomapped(correctdate);
                                    else
                                        whatxvecXXX(mm)=[];
                                    end
                                else
                                    timebetweensamplesXXX(mm)=[];
                                    whatxvecXXX(mm)=[];
                                end
                            end
                            
                            whatIDvec{yy}=persontable.person_id;
                            whatshiftvec{yy}=timebetweensamplesXXX-timebetweensamplesXXX(1); %%something like this
                            whatxvec{yy}=whatxvecXXX+0.01*(whatxvecXXX==0);
                            whatlogxvec{yy}=-1*(whatxvecXXX==0.01)+~(whatxvecXXX==0.01).*log10(whatxvecXXX);
                            if length(whatshiftvec{yy})>length(whatxvec{yy})
                                whatshiftvec{yy}((length(whatxvec{yy})+1):end)=[];
                            end
                            yy=yy+1;
                        end
                    end
                end
            end
        end
    end
end