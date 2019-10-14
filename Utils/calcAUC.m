function FR=calcAUC(sigtt)
    sigtt(sigtt<0)=0;
    FR=zeros(size(sigtt,1),1);
    for n=1:size(sigtt,1)
        FR(n)=trapz(sigtt(n,:));
    end
    FR=FR./size(sigtt,2);
end