function si=calcSI(sigt,oi,binX,binY,Yedges,Xedges)
    sigmean=mean(sigt);
    zsmap=zeros(size(oi));
    simap=zeros(size(oi));
    for y=1:size(zsmap,1)
        for x=1:size(zsmap,2)
            indt=find(binX==x & binY==y);
            if ~isempty(indt)               
                ai=mean(sigt(indt),'omitnan');
                simap(y,x)=oi(y,x)*ai*log2(ai/sigmean);
            else
                zsmap(y,x)=0;
                simap(y,x)=0;
            end
        end
    end
    simap(isnan(simap))=0;
    si=sum(simap(:));
end