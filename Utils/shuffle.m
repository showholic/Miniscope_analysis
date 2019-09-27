function sigshuffle=shuffle(sigtemp,maxshift,splitsize)
    randshift=randi(maxshift,1);
    sigshift=circshift(sigtemp,randshift);
    randsplit=randperm(splitsize);
    sigsplit=reshape(sigshift,[],splitsize);
    sigshuffle=reshape(sigsplit(:,randsplit),1,[]);
end