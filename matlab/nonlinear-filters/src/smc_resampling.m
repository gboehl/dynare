function indx = smc_resampling(weights,noise,number)
    indx = zeros(number,1);
    cumweights = cumsum(weights);
    randvec = (transpose(1:number)-1+noise(:))/number;
    j = 1;
    for i=1:number
        while (randvec(i)>cumweights(j))
            j = j+1;
        end
        indx(i) = j;
    end
