clear all
clc

figcount=1;

load('themetadata')

%%%Plot how much data in terms of time since last negative
figure(1)
histogram(metadatatable.timebetweensamples,0:7:140)
xlabel('time between sample and last negative (days)')
ylabel('Number of samples')

metadatatable(isnan(metadatatable.ct_mean),:)=[];
metadatatable(isnan(metadatatable.log10mapped),:)=[];

metadatatable.VaccDate=datetime(metadatatable.covid_vaccine_date_1(:),'InputFormat','ddMMMyyyy');
vaccbeforpositive=metadatatable.VaccDate<=metadatatable.collection_date;
metadatatable=[metadatatable, array2table(vaccbeforpositive,'VariableNames',{'vaccinebeforepositive'})];

serobeforepositive=metadatatable.first_seroPos_date<=metadatatable.collection_date;
metadatatable=[metadatatable, array2table(serobeforepositive,'VariableNames',{'serobeforepositive'})];

eitherbeforepositive=categorical(metadatatable.vaccinebeforepositive);
eitherbeforepositive(categorical(metadatatable.serobeforepositive)=='true')='true';

metadatatable=[metadatatable, array2table(eitherbeforepositive)];

head(metadatatable)

%%%%%possible without defined other
ForVLb117=metadatatable.B_1_1_7;
ForVLb117(:)='_';
ForVLb117(categorical(metadatatable.B_1_1_7)=='True')='Alpha';
ForVLb117(categorical(metadatatable.B_1_1_7)=='False')='Not VOC';

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

ForVLb117(categorical(metadatatable.B_1_617_2)=='True')='Delta';

for i=1:max(metadatatable.collectionweek_HF)
     howmanytrue(i)=sum((metadatatable.collectionweek_HF==i).*(ForVLb117=='Alpha' ));
    howmanyfalse(i)=sum((metadatatable.collectionweek_HF==i).*(ForVLb117=='Not VOC'));
end
proportionstrue=howmanytrue./(howmanytrue+howmanyfalse);

figure(77); plot(1:max(metadatatable.collectionweek_HF),proportionstrue)


metadatatable=[metadatatable,array2table(ForVLb117)];


%%%%%%%%LOG10MAPPED%%%%%%%%%
metadatatable(isnan(metadatatable.ct_mean),:)=[];
trueALL=find((metadatatable.ForVLb117=='Alpha').*isfinite(metadatatable.log10mapped));
falseALL=find((metadatatable.ForVLb117=='Not VOC').*isfinite(metadatatable.log10mapped));
OtherVOCALL=find((metadatatable.ForVLb117=='otherVOC').*isfinite(metadatatable.log10mapped));
DELTAALL=find((metadatatable.ForVLb117=='Delta').*isfinite(metadatatable.log10mapped));
DELTAnovaccine=find((metadatatable.ForVLb117=='Delta').*(categorical(metadatatable.eitherbeforepositive)=='false').*isfinite(metadatatable.log10mapped));
DELTAvaccine  =find((metadatatable.ForVLb117=='Delta').*(categorical(metadatatable.eitherbeforepositive)=='true').*isfinite(metadatatable.log10mapped));
newTab2=metadatatable([trueALL',falseALL',DELTAALL'],:);%,OtherVOCALL'],:);
% newTab2=metadatatable([trueALL',falseALL',DELTAnovaccine',DELTAvaccine'],:);%,OtherVOCALL'],:);
newTab2.ForVLb117=removecats(newTab2.ForVLb117);
newTab2.ForVLb117 = reordercats(newTab2.ForVLb117,{'Not VOC','Alpha','Delta'});
figure(306)
metadatatable.vaccinebeforepositive=categorical(metadatatable.vaccinebeforepositive);
% BB=boxchart(newTab2.ForVLb117,newTab2.ct_mean)%,'linewidth',1.5,'Fontsize','12')
hold on; boxchart(newTab2.ForVLb117,newTab2.log10mapped,'GroupByColor',newTab2.vaccinebeforepositive,'boxwidth',0.9,'JitterOutliers','on')
ylabel('log_1_0(mapped)')
gtext({'vaccinated';'unvaccinated'})



metadatatable(isnan(metadatatable.ct_mean),:)=[];
trueALL=find((metadatatable.ForVLb117=='Alpha').*isfinite(metadatatable.ct_mean));
falseALL=find((metadatatable.ForVLb117=='Not VOC').*isfinite(metadatatable.ct_mean));
OtherVOCALL=find((metadatatable.ForVLb117=='otherVOC').*isfinite(metadatatable.ct_mean));
DELTAALL=find((metadatatable.ForVLb117=='Delta').*isfinite(metadatatable.ct_mean));
DELTAnovaccine=find((metadatatable.ForVLb117=='Delta').*(categorical(metadatatable.eitherbeforepositive)=='false').*isfinite(metadatatable.ct_mean));
DELTAvaccine=find((metadatatable.ForVLb117=='Delta').*(categorical(metadatatable.eitherbeforepositive)=='true').*isfinite(metadatatable.ct_mean));
newTab2=metadatatable([trueALL',falseALL',DELTAALL'],:);%,OtherVOCALL'],:);
% newTab2=metadatatable([trueALL',falseALL',DELTAnovaccine',DELTAvaccine'],:);%,OtherVOCALL'],:);
newTab2.ForVLb117=removecats(newTab2.ForVLb117);
categories(newTab2.ForVLb117)
newTab2.ForVLb117 = reordercats(newTab2.ForVLb117,{'Not VOC','Alpha','Delta','otherVOC'});
figure(307)
metadatatable.vaccinebeforepositive=categorical(metadatatable.vaccinebeforepositive);
% BB=boxchart(newTab2.ForVLb117,newTab2.ct_mean)%,'linewidth',1.5,'Fontsize','12')
hold on; boxchart(newTab2.ForVLb117,newTab2.ct_mean,'GroupByColor',newTab2.vaccinebeforepositive,'boxwidth',0.9,'JitterOutliers','on')
ylabel('ct mean')
gtext({'vaccinated';'unvaccinated'})



newtabAlpha=metadatatable([trueALL'],:);
newtabDelta=metadatatable([DELTAALL'],:);
newtabNotVOC=metadatatable([falseALL'],:);
figure;
% hold on; scatter(newtabAlpha.ct_mean, newtabAlpha.log10mapped,'rx')
hold on; scatter(newtabDelta.ct_mean, newtabDelta.log10mapped,'bx')
% hold on; scatter(newtabNotVOC.ct_mean,newtabNotVOC.log10mapped,'gx')
ylabel('ct mean')


ccc1 = polyfit(newtabAlpha.ct_mean,newtabAlpha.log10mapped,1);
ccc2 = polyfit(newtabDelta.ct_mean,newtabDelta.log10mapped, 1);
ccc3 = polyfit(newtabNotVOC.ct_mean,newtabNotVOC.log10mapped,1);

% figure;
ffff=7:35;
% hold on; plot(ffff,polyval(ccc1,ffff),'r','linewidth',2);
hold on; plot(ffff,polyval(ccc2,ffff),'b','linewidth',2);
% hold on; plot(ffff,polyval(ccc3,ffff),'b');
xlabel('ct mean')
ylabel('log_1_0(mapped)')
legend({'Alpha','Delta'})
% legend({'Alpha','Delta','Not VOC'})


figure
% hold on; scatter(newTab2.ForVLb117, newTab2.ct_mean,'x','jitter','on')
hold on; gscatter(newTab2.ForVLb117, newTab2.ct_mean,newTab2.vaccinebeforepositive)
title('All samples')
ylabel('mean(ct)')


head(metadatatable)

figcount=figcount+1;
figure(figcount)
truex=metadatatable.timebetweensamples(find((metadatatable.ForVLb117=='Alpha').*isfinite(metadatatable.ct_mean)));
falsex=metadatatable.timebetweensamples(find((metadatatable.ForVLb117=='Not VOC').*isfinite(metadatatable.ct_mean)));
deltax=metadatatable.timebetweensamples(find((metadatatable.ForVLb117=='Delta').*isfinite(metadatatable.ct_mean)));
truey=metadatatable.ct_mean(find((metadatatable.ForVLb117=='Alpha').*isfinite(metadatatable.ct_mean)));
falsey=metadatatable.ct_mean(find((metadatatable.ForVLb117=='Not VOC').*isfinite(metadatatable.ct_mean)));
deltay=metadatatable.ct_mean(find((metadatatable.ForVLb117=='Delta').*isfinite(metadatatable.ct_mean)));
truemat=[truex,truey];
falsemat=[falsex,falsey];
deltamat=[deltax,deltay];
sortedtrue=sortrows(truemat,'descend');
sortedfalse=sortrows(falsemat,'descend');
sorteddelta=sortrows(deltamat,'descend');
sortedtrue((isnan(sortedtrue(:,1))),:)=[];
sortedfalse((isnan(sortedfalse(:,1))),:)=[];
sorteddelta((isnan(sorteddelta(:,1))),:)=[];
fittrue=smooth(sortedtrue(:,1)/2,sortedtrue(:,2),0.85,'loess');
fitfalse=smooth(sortedfalse(:,1)/2,sortedfalse(:,2),0.85,'loess');
fitdelta=smooth(sorteddelta(:,1)/2,sorteddelta(:,2),0.85,'loess');
hold on; plot(sortedfalse(:,1)/2,fitfalse,'k-','linewidth',1.5)
hold on; plot(sortedtrue(:,1)/2,fittrue,'r-','linewidth',1.5)
hold on; plot(sorteddelta(:,1)/2,fitdelta,'c-','linewidth',1.5)
hold on; scatter(truex/2,truey,'kx')
hold on; scatter(falsex/2,falsey,'rx')
hold on; scatter(deltax/2,deltay,'cx')
xlim([0,50])
ylabel('mean(ct)')
xlabel('midpoint of time since last negative')
legend('Not VOC','Alpha','Delta')



%%%%PLOT how medians change with 
% cutoffvector=[8:40];
% for i=1:length(cutoffvector)    
%     truex=metadatatable.log10mapped(find((metadatatable.timebetweensamples<cutoffvector(i)).*(metadatatable.ForVLb117=='True').*isfinite(metadatatable.log10mapped)));
%     falsex=metadatatable.log10mapped(find((metadatatable.timebetweensamples<cutoffvector(i)).*(metadatatable.ForVLb117=='False').*isfinite(metadatatable.log10mapped)));    
%     mediantrue(i)=median(truex);
%     medianfalse(i)=median(falsex);
%     iqrtrue(i)=iqr(truex);
%     iqrfalse(i)=iqr(falsex);  
% end
% 
% figcount=figcount+1;
% figure(figcount)
% hold on;plot(cutoffvector,mediantrue,'rx-')
% hold on;plot(cutoffvector,medianfalse,'bx-')
% hold on;plot(cutoffvector,mediantrue+iqrtrue,'r:')
% hold on;plot(cutoffvector,mediantrue-iqrtrue,'r:')
% hold on;plot(cutoffvector,medianfalse+iqrfalse,'b:')
% hold on;plot(cutoffvector,medianfalse-iqrfalse,'b:')
% xlabel('days since last negative limit')
% ylabel('median log10mapped')
% legend('b117true','b117false')

%%%%PLOT boxplots of vl in first 20 days
 dayscutofftrue=20;
dayscutofffalse=20;
truelessthan8=find((metadatatable.timebetweensamples<dayscutofftrue).*(metadatatable.ForVLb117=='Alpha').*isfinite(metadatatable.ct_mean));
falselessthan8=find((metadatatable.timebetweensamples<dayscutofffalse).*(metadatatable.ForVLb117=='Not VOC').*isfinite(metadatatable.ct_mean));
otherVOClessthan8=find((metadatatable.timebetweensamples<dayscutofffalse).*(metadatatable.ForVLb117=='otherVOC').*isfinite(metadatatable.ct_mean));
DELTAlessthan8=find((metadatatable.timebetweensamples<dayscutofffalse).*(metadatatable.ForVLb117=='Delta').*isfinite(metadatatable.ct_mean));
newTab=metadatatable([truelessthan8',falselessthan8',DELTAlessthan8'],:);%,otherVOClessthan8',],:);
newTab.ForVLb117=removecats(newTab.ForVLb117);
% writetable(newTab,'newTab.csv')
mediantruelessx=median(metadatatable.ct_mean(truelessthan8))
medianfalselessx=median(metadatatable.ct_mean(falselessthan8))
prighttaillessthanx = ranksum(metadatatable.log10mapped(truelessthan8),metadatatable.ct_mean(falselessthan8),'tail','right')
  ptwotaillessthanx = ranksum(metadatatable.log10mapped(truelessthan8),metadatatable.ct_mean(falselessthan8),'tail','right')
figcount=figcount+1;
figure(figcount)
BB=boxchart(newTab.ForVLb117,newTab.ct_mean)%'linewidth',1.5)
hold on; scatter(newTab.ForVLb117, newTab.ct_mean,'x','jitter','on')
ylabel('mean(ct)')
title('No more than 10 days since a negative test')
xticklabels({'Not VOC','Alpha','Delta'})


mediantrueALL=median(metadatatable.log10mapped(trueALL))
medianfalseALL=median(metadatatable.log10mapped(falseALL))
prighttailALL = ranksum(metadatatable.log10mapped(trueALL),metadatatable.log10mapped(falseALL),'tail','right')
ptwotailedALL = ranksum(metadatatable.log10mapped(trueALL),metadatatable.log10mapped(falseALL),'tail','right')



%%%PLOT vl over time since last sample (in multiple samples)
allgooddata=metadatatable(find(isfinite(metadatatable.log10mapped).*(metadatatable.counts_person_id>1).*(metadatatable.counts_person_id<100).*(~isnat(metadatatable.last_negative_date))),:);

%%%%%%%%%%%%%%%%%%%%%%

for k=1:height(tableofpersonidcounts)
        findpersonsamples=find(metadatatable.person_id==tableofpersonidcounts.person_id(k));
        persontable{k}={tableofpersonidcounts(k,:),metadatatable(findpersonsamples,:)};
end

% allgooddatatrue=allgooddata(find(allgooddatatrue.ForVLb117=='True'),:);
% allgooddatafalse=allgooddata(find(allgooddatatrue.ForVLb117=='False'),:);

figcount=figcount+1;
figure(figcount)
BIGgradtrue=[];
BIGgradfalse=[];
for i=1:length(persontable)
    personcell=persontable{i};
    Personinfo=personcell{1};
    PersonTable=personcell{2};
%     if (Personinfo.counts(1)>1)&&(PersonTable.timebetweensamples(1)<17)
    if (Personinfo.counts(1)>1)  %%dont restrict who is included
        x=PersonTable.timebetweensamples;
        y=PersonTable.log10mapped;
        if (PersonTable.ForVLb117(1)=='True')
            hold on; plot(x,y,'ko-')
           [grad,S]=polyfit(x,y,1);
           if isnan(grad)
%                grad
%                x
%                y
           end
            BIGgradtrue=[BIGgradtrue,grad(1)];
        elseif (PersonTable.ForVLb117(1)=='False')
            hold on; plot(x,y,'ro-')
            grad=polyfit(x,y,1);
            [grad,S]=polyfit(x,y,1); 
            BIGgradfalse=[BIGgradfalse,grad(1)];
        end
    end
end
BIGgradtrue=rmmissing(BIGgradtrue);
BIGgradfalse=rmmissing(BIGgradfalse);
xlabel('time since last negative test')
ylabel('log10mapped')
% gtext({'black=B117';'red=not B117'})


% figcount=figcount+1;
% figure(figcount)
% togethervl=[BIGgradtrue,BIGgradfalse];
% togetherlog=[0*BIGgradtrue,1+0*BIGgradfalse];
% hold on; scatter(togetherlog,togethervl,'x','jitter','on')
% hold on; boxchart(togetherlog,togethervl)
% xticks([0,1])
% xticklabels({'probably B117','probably B1177'})
% xlim([-0.5,1.5])
% title('gradients')
% mediantrueALL2=median(BIGgradtrue)
% medianfalseALL2=median(BIGgradfalse)
% prighttailALL2 = ranksum(BIGgradtrue,BIGgradfalse,'tail','right')
% ptwotailedALL2 = ranksum(BIGgradtrue,BIGgradfalse,'tail','right')
% 
% 
% 

B117 = categorical(metadatatable.ForVLb117);
SOURCESEX = categorical(metadatatable.source_sex);
CTPATTERN = categorical(metadatatable.ct_pattern);

DATEASNUM = datenum(metadatatable.collection_date)-min(datenum(metadatatable.collection_date));
metadatatable=[metadatatable,array2table(DATEASNUM)];
metadatatable.Properties.VariableNames(width(metadatatable)) = {'dateasnum'};

whereOKDATE=isfinite(DATEASNUM);
whereOKLOG10MAPPED=isfinite(metadatatable.log10mapped);
whereOKSOURCEAGE=isfinite(metadatatable.source_age);
whereOKB117=(B117=='Alpha')+(B117=='Not VOC');
whereOKSEX=(SOURCESEX=='M')+(SOURCESEX=='F');

findwhereok=find(whereOKLOG10MAPPED.*whereOKB117);
% findwhereok=find(whereOKDATE.*whereOKLOG10MAPPED.*whereOKSOURCEAGE.*whereOKB117.*whereOKSEX);

minitabtable=metadatatable(findwhereok,:);
% writetable(minitabtable,'mintiabtable.csv')


% figure(1);
% hold on; boxchart(minitabtable.log10mapped,minitabtable.ForVLb117)
a=minitabtable.log10mapped(find(categorical(minitabtable.ForVLb117)=='Alpha'));
b=minitabtable.log10mapped(find(categorical(minitabtable.ForVLb117)=='Not VOC'));

% tcutoff=[-1,4:4:20];
tcutoff=[-1,4,12,20];
tlower=tcutoff(1:(length(tcutoff)-1));
tupper=tcutoff(2:(length(tcutoff)));
mediana=0*tlower;
medianb=0*tlower;
for i=1:length(tlower)
   wherelooktruenumber=(minitabtable.dateasnum>tlower(i)).*(minitabtable.dateasnum<=tupper(i));
   aa=minitabtable.log10mapped(find(wherelooktruenumber.*(categorical(minitabtable.ForVLb117)=='Alpha')));
   bb=minitabtable.log10mapped(find(wherelooktruenumber.*(categorical(minitabtable.ForVLb117)=='NOt VOC')));
   mediana(i)=median(aa);
   medianb(i)=median(bb);
end

% figcount=figcount+1;
% figure(figcount)
% % hold on; plot(tupper,mediana,'kx-')
% % hold on; plot([4,20],[median(b),median(b)],'r-')
% hold on; plot(tupper,medianb,'gx-')
% legend('B117','not B117')
% xlabel('time (days)')
% ylabel('median log10(mapped reads)')
% axis([4,20,2.5,5])
% title('test2')


tcutoff=[-1,4:4:20];
%tcutoff=[-1,4,12,20];
tlower=tcutoff(1:(length(tcutoff)-1));
tupper=tcutoff(2:(length(tcutoff)));
mediana=0*tlower;
medianb=0*tlower;
for i=1:length(tlower)
   wherelooktruenumber=(minitabtable.dateasnum>tlower(i)).*(minitabtable.dateasnum<=tupper(i));
   aa=minitabtable.ct_mean(find(wherelooktruenumber.*(categorical(minitabtable.ct_pattern)=='OR+N')));
   bb=minitabtable.ct_mean(find(wherelooktruenumber.*(categorical(minitabtable.ct_pattern)=='OR+N+S')));
   mediana(i)=median(aa);
   medianb(i)=median(bb);
end

% figcount=figcount+1;
% figure(figcount)
% a=minitabtable.ct_mean(find(categorical(minitabtable.ct_pattern)=='OR+N'));
% b=minitabtable.ct_mean(find(categorical(minitabtable.ct_pattern)=='OR+N+S'));
% hold on; plot(tupper,mediana,'kx-')
% hold on; plot([4,20],[median(b),median(b)],'r-')
% % hold on; plot(tupper,medianb,'gx-')
% legend('B117','not B117')
% xlabel('time (days)')
% ylabel('median log10(mapped reads)')


figcount=figcount+1;
figure(figcount)
AGE=minitabtable.source_age;
AGEasindex=AGE;
AGEasindex(AGE==3)=2;
AGEasindex(AGE==7)=4;
AGEasindex(AGE==15)=6;
AGEasindex(AGE==25)=8;
AGEasindex(AGE==35)=10;
AGEasindex(AGE==45)=12;
AGEasindex(AGE==55)=14;
AGEasindex(AGE==65)=16;
AGEasindex(AGE==75)=18;
AGEasindex(AGE==85)=20;
tableforplotting=[minitabtable,array2table(AGEasindex)];
hold on; boxchart(tableforplotting.AGEasindex,tableforplotting.log10mapped,'GroupByColor',tableforplotting.ForVLb117,'boxwidth',0.9,'JitterOutliers','on')

figure; boxchart(tableforplotting.AGEasindex,tableforplotting.log10mapped,'boxwidth',0.9,'JitterOutliers','on')


B117S = categorical(tableforplotting.ForVLb117);

for n=1:max(unique(tableforplotting.AGEasindex))
    whereB117trueAGE=find((tableforplotting.AGEasindex==n).*(B117S=='Alpha'));
    whereB117falsAGE=find((tableforplotting.AGEasindex==n).*(B117S=='Not VOC'));
    X1=tableforplotting.AGEasindex(whereB117trueAGE);
    X2=tableforplotting.AGEasindex(whereB117falsAGE);
    Y1=tableforplotting.log10mapped(whereB117trueAGE);
    Y2=tableforplotting.log10mapped(whereB117falsAGE);
%     hold on; scatter(ones(sum(X1'),1)+n-1.3,Y1,'jitter','on','JitterAmount',0.1);
%     hold on; scatter(ones(sum(X2'),1)+n-0.7,Y2,'jitter','on','JitterAmount',0.1);
end
% set(gca,"XTick", unique(AGE),"XTickLabel",categories(L))

xticks(2:2:20)
xticklabels({'3','7','15','25','35','45','55','65','75','85'})
xlabel('age')
ylabel('log_1_0(mapped)')
legend('not B117','B117')


