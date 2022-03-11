% for i=1:length(whatshiftvec)
%     wwww= size(whatshiftvec{i});
%     if (wwww(1)==10)&&(whatshiftvec{i}(end)==286)
%         i
%         whatshiftvec{i}
%     end
% end

maxshift=0;
for i=1:length(whatshiftvec)
    maxshiftNEW=max(whatshiftvec{i});
    if maxshiftNEW>maxshift
        maxshift=maxshiftNEW;
    end   
end
maxshift


% whatshiftvec{289}
% whatIDvec{289}
% eeee=metadatatable(metadatatable.person_id=='ONS_000044',:)
% 
% i=289

%      0
%     23
%     51
%     79
%    134
%    162
%    193
%    228
%    258
%    286