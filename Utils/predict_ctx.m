function [accuracyA,accuracyB]=predict_ctx(sig,session_start,protocol,ms_ts)
    sigz=zscore(sig,[],2);
    ctx_dur=280;
    sig_ctxBpre=sigz(:,session_start((strcmp(protocol,'preB')==1)):session_start((strcmp(protocol,'preB')==1))+ctx_dur*10-1);
    sig_ctxApre=sigz(:,session_start((strcmp(protocol,'preA')==1)):session_start((strcmp(protocol,'preA')==1))+ctx_dur*10-1);
    sig_ctxBpost=sigz(:,session_start((strcmp(protocol,'postB')==1)):session_start((strcmp(protocol,'postB')==1))+ctx_dur*10-1);
    sig_ctxApost=sigz(:,session_start((strcmp(protocol,'postA')==1)):session_start((strcmp(protocol,'postA')==1))+ctx_dur*10-1);
    
    %%
    siginterp_ctxApre=interpsig(sig_ctxApre,'preA',protocol,ms_ts,ctx_dur);
    siginterp_ctxBpre=interpsig(sig_ctxBpre,'preB',protocol,ms_ts,ctx_dur);
    siginterp_ctxApost=interpsig(sig_ctxApost,'postA',protocol,ms_ts,ctx_dur);
    siginterp_ctxBpost=interpsig(sig_ctxBpost,'postB',protocol,ms_ts,ctx_dur);

    %%
    pvd=[siginterp_ctxBpre siginterp_ctxApre]';
    ctx=cell(size(pvd,1),1);
    for i=1:ctx_dur
        ctx{i}='B';
    end
    for i=ctx_dur+1:size(pvd,1)
        ctx{i}='A';
    end

    MdLinear=fitcdiscr(pvd,ctx);
    %%
    predict_postB=zeros(size(siginterp_ctxBpost,2),1);
    for n=1:size(siginterp_ctxBpost,2)
        prediction=predict(MdLinear,siginterp_ctxBpost(:,n)');
        if strcmp(prediction,'B')==1
            predict_postB(n)=1;
        else
            predict_postB(n)=0;
        end
    end
    accuracyB=length(find(predict_postB==1))/size(siginterp_ctxBpost,2);
    %%
    predict_postA=zeros(size(siginterp_ctxApost,2),1);
    for n=1:size(siginterp_ctxApost,2)
        prediction=predict(MdLinear,siginterp_ctxApost(:,n)');
        if strcmp(prediction,'A')==1
            predict_postA(n)=1;
        else
            predict_postA(n)=0;
        end
    end
    accuracyA=length(find(predict_postA==1))/size(siginterp_ctxApost,2);


    %%
    function siginterp=interpsig(sigt,trialname,protocol,ms_ts,ctx_dur)

        ms_ctx=double(ms_ts{(strcmp(protocol,trialname)==1)}(1:ctx_dur*10))/1000;
        siginterp=zeros(size(sigt,1),ctx_dur);
        ts_new=1:1:ctx_dur;
        for n=1:size(sigt,1)
            siginterp(n,:)=interp1(ms_ctx,sigt(n,:),ts_new);
        end
    end
end