function binresponse=avgbinresponse(dur,binsize,response)
    bins=discretize(1:dur,0:binsize:dur);
    bincount=dur/binsize;
    binresponse=zeros(size(response,1),bincount);
    for i = 1:bincount
        binresponse(:,i)=mean(response(:,bins==i),2);
    end
end