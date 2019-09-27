function FR=calcTF(sigtt)
    FR=zeros(size(sigtt,1),1);
    for i=1:size(sigtt,1)
        sigt=sigtt(i,:);
        bsline=median(sigt);
        stdsigt=mad(sigt);
        threshold=bsline+1.65*stdsigt;
        pks=findpeaks(sigt,'MinPeakDistance',10,'MinPeakdistance',3*stdsigt,'MinPeakHeight',threshold,'MinPeakProminence',3);
        FR(i)=numel(pks);
    end
end