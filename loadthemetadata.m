clear all
clc

metadatatable = readtable('C:\Users\helenf\Dropbox\ONS_sequencing_restricted_access_shared\ONS_Summary_Spreadsheets\ONS_summary_newFormat_210727.csv');

metadatatable.log10mapped=str2double(metadatatable.log10mapped);
metadatatable.mapped=str2double(metadatatable.mapped);

metadatatable.sample_name=categorical(metadatatable.sample_name);
metadatatable.B_1_1_7=categorical(metadatatable.B_1_1_7);
metadatatable.person_id=categorical(metadatatable.person_id);

%%SORT OUT THE WRONG YEAR FORMATTING
for i=1:height(metadatatable)
    
    ggg=metadatatable.collection_date(i);
    if isnat(ggg)
    else
        newStr = strrep(datestr(ggg),'2121','2021');
        metadatatable.collection_date(i)=datetime(newStr,'Format','yyyy-MM-dd');
    end
end
for i=1:height(metadatatable)
    ggg=metadatatable.last_negative_date(i);
    if isnat(ggg)
    else
        newStr = strrep(datestr(ggg),'2121','2021');
        metadatatable.last_negative_date(i)=datetime(newStr,'Format','yyyy-MM-dd');
    end
end

%%ADD the time between samples
timebetweensamples = datenum(metadatatable.collection_date)-datenum(metadatatable.last_negative_date);
metadatatable=[metadatatable,array2table(timebetweensamples)];

%%Add counts_person_id
[N,X]=hist(metadatatable.person_id);
tableofpersonidcounts=[array2table(N','VariableNames',{'counts'}),cell2table(X','VariableNames',{'person_id'})];
countsnew=zeros(height(metadatatable),1);
for i=1:height(metadatatable)
    if isundefined(metadatatable.person_id(i))
    else
        countsnew(i)=tableofpersonidcounts.counts(find(metadatatable.person_id(i)==tableofpersonidcounts.person_id));
    end
end
metadatatable=[metadatatable,array2table(countsnew,'VariableNames',{'counts_person_id'})];

 %%Add my own collection week
 daysfrom20201214=datenum(metadatatable.collection_date)-datenum(datetime(2020,12,14));
 weeksfrom20201214=ceil(daysfrom20201214/7)+1;
 metadatatable=[metadatatable,array2table(weeksfrom20201214,'Variablenames',{'collectionweek_HF'})];

metadatatable.sequence_name=categorical(metadatatable.sequence_name);


% ff2=datetime(metadatatable.covid_vaccine_date_1,'format','ddMMMyyyy')
% %sort out the vaccine date issues
% gg2=metadatatable.covid_vaccine_date_1;
% ff2=categorical(gg2);
% hhh=ff2;
% for i=length(gg2):-1:1
%     if isnat(gg2(i))
%         hhh(i)='NoData';
%     end
% end
% metadatatable=[metadatatable,array2table(hhh,'VariableNames',{'VaccDate'})];

head(metadatatable)
save('themetadata','metadatatable','tableofpersonidcounts')