function FR=calcFR(sigtt)
    sigtt(sigtt>0)=1;
    FR=sum(sigtt,2)/size(sigtt,2);
end
