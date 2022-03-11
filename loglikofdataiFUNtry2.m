function loglikofdataiFUNtry2=loglikofdataiFUNtry2(D,shift,logxijdata,thebiggamWWW,thebiggamH,P_a_t_s_DofDst,Aspace,log10xbarmatrix, pfalsenegREADFROM, normof0READFROM, qqqqqREADME, HHspace,age)

AA=1:(D/Aspace);

pxijgaijHWWW=1+0*log10xbarmatrix(AA,:,:);
% pxijgaij000=1+0*log10xbarmatrix(AA,:,:);

logxijdataALT=logxijdata;
logxijdataALT(logxijdataALT==-Inf)=-HHspace;

for j=1:length(logxijdata)
%     
%     AA
%     shift(j)
%     AA+shift(j)/Aspace
   

%     size(log10xbarmatrix)
%     AA
%     shift(j)
%     AA+shift(j)/Aspace
    
    whatlogxbarINDEX=round(log10xbarmatrix(AA+shift(j)/Aspace,:,:)/HHspace+2);
    pfalseneg=pfalsenegREADFROM(whatlogxbarINDEX);
           
    psample0gxijhat=normof0READFROM(whatlogxbarINDEX);  
    
    psamplexijgxijhat=qqqqqREADME(round(logxijdataALT(j)/HHspace+2),whatlogxbarINDEX);
    
    psamplexijgxijhatNEW=whatlogxbarINDEX;
    psamplexijgxijhatNEW(:)=psamplexijgxijhat;
    psamplexijgxijhat=psamplexijgxijhatNEW;

%     pxijgaij000=pxijgaij000.*(pfalseneg.*(1-psample0gxijhat)+psample0gxijhat);
    if j==1
    pxijgaij000=pfalseneg.*(1-psample0gxijhat)+psample0gxijhat;
    end

    if logxijdataALT(j)<0
        pxijgaijHWWW=pxijgaijHWWW.*(pfalseneg.*(1-psample0gxijhat)+psample0gxijhat); %%if xij is 0
    else
        pxijgaijHWWW=pxijgaijHWWW.*(1-pfalseneg).*psamplexijgxijhat;
    end
end

sumover0000=sum(sum( pxijgaij000.*thebiggamWWW(AA,:,:).*thebiggamH(AA,:,:),3),2);
sumoverHWWW=sum(sum(pxijgaijHWWW.*thebiggamWWW(AA,:,:).*thebiggamH(AA,:,:),3),2);

% figure
% hold on;
% plot(1:length(P_a_t_s_DofDst),P_a_t_s_DofDst,'r-')
% P_a_t_s_DofDst(AA)

sumovera000=sum(sumover0000'.*P_a_t_s_DofDst(AA),2); % if xij is     end
sumoveraxij=sum(sumoverHWWW'.*P_a_t_s_DofDst(AA),2); % if xij is     end

loglikofdataiFUNtry2=log(sumoveraxij/(1-sumovera000));



