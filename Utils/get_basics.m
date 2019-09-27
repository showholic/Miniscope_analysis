function newdata=get_basics(animaldata)
    newdata=animaldata;
    newdata.sig=animaldata.ms.sigraw';
    shockts=animaldata.shockts;
    ms_start=animaldata.ms_start;
    shock_start=animaldata.shock_start;
    ms_ts=animaldata.ms.ms_ts;
    session_dur=cellfun(@length,ms_ts);
    session_start=1;
    for i=1:length(session_dur)-1
        starttemp=session_start(i)+session_dur(i);
        session_start=[session_start; starttemp];
    end
    session_end=session_start+session_dur'-1;
    
    % Get shock time
    shock_ts2=shockts.time-shockts.time(ms_start(4));
    shock_bvt=shock_ts2(shock_start);
    shock_mst=ms_ts{4}/1000;
    temp1=abs(double(repmat(shock_mst,length(shock_bvt))')-shock_bvt');
    [minval,frame_shock]=min(temp1,[],1);
    newdata.session_start=session_start;
    newdata.session_end=session_end;
    newdata.frame_shock=frame_shock;
end