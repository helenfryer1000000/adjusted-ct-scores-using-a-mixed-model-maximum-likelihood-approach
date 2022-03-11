function log10mappedFUN2=log10mappedFUN2(avec,WWWvec,Hvec,HHspace,WsplitWW)

  A=zeros(length(avec),length(WWWvec),length(Hvec));
WWW=zeros(length(avec),length(WWWvec),length(Hvec));
  H=zeros(length(avec),length(WWWvec),length(Hvec));

for i1=1:length(avec)
    for i2=1:length(WWWvec)
        for i3=1:length(Hvec)
            A(i1,i2,i3)=avec(i1);
            WWW(i1,i2,i3)=WWWvec(i2);
            H(i1,i2,i3)=Hvec(i3);
        end
    end
end

W=WsplitWW*WWW;
WW=(1-WsplitWW)*WWW;

cutoffvec=[0*W(:),W(:),W(:)+WW(:)];

q=(A(:)>cutoffvec);
q(:,1:2)=q(:,1:2)-q(:,2:3);

log10mappedFUN22=sum(q.*[A(:).*H(:)./W(:),H(:)-(A(:)-W(:)).*H(:)./WW(:),-HHspace+0*A(:)],2);

A(:)=log10mappedFUN22;
log10mappedFUN2=round(A,-log10(HHspace));
