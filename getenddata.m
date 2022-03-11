clear all
clc

maxgap=20;

load C:\Users\helenf\Desktop\ONS_metadata2\themetadata.mat

fulldataXXX = readtable('C:\Users\helenf\Desktop\ONS_metadata2\metadata20210428\ONS_positives_all_results_household_20210426.csv','VariableNamingRule','modify');
fulldataXXX.pseudo_person_id=categorical(fulldataXXX.pseudo_person_id);
fulldataXXX.result_mk=categorical(fulldataXXX.result_mk);
WP=1;
% for i=1:height(fulldataXXX)
for i=1:9
    if (i==1)||(fulldataXXX.pseudo_person_id(i)~=fulldataXXX.pseudo_person_id(i-1))
        xxx=fulldataXXX(fulldataXXX.pseudo_person_id==fulldataXXX.pseudo_person_id(i),:);
        WPtable{WP}=xxx;
        yyy=xxx(xxx.result_mk=="Positive",:);
        WPtablepositive{WP}=yyy;
        firstpositive=min(yyy.collection_date);
        zzz=xxx(xxx.collection_date>=firstpositive,:);
        
        if height(zzz)>1
            for j=height(zzz):-1:2
                if (zzz.result_mk(j)=="Positive")&&(zzz.result_mk(j-1)=="Negative")&&(zzz.collection_date(j)-zzz.collection_date(j-1)>maxgap)
                    
                    zzz(j:end,:)=[];
                end
            end
        end
        allafterfirstpositive{WP}=zzz;
    end
end
        