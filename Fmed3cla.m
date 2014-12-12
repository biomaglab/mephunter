function [fmed, rms, fmean] = Fmed3cla(x,fsamp,epoch)
%[P,f]=psd(x,fsamp*epoch_len,fsamp,ones(fsamp*epoch_len,1),0);

if mod(length(x),epoch), error('vector length should be multiple of the epoch chosen'), return, end
if all(isnan(x)),
    fmed = nan(round(size(x,1)/epoch),1);
    rms = fmed;
    fmean = rms;
    return
end
    
ReshapedX = detrend(reshape(x,epoch,numel(x)/epoch),'linear');
rms = zeros(size(ReshapedX,2), 1);
fmean = zeros(size(ReshapedX,2), 1);
fmed = zeros(size(ReshapedX,2), 1);
for column = 1:size(ReshapedX,2)
    x = ReshapedX(:,column);
    
    rms(column,1)=norm(x)/sqrt(length(x));
    x=x-mean(x);
%     [P,f] = psd(x,fsamp,fsamp,boxcar(length(x)),0);
    % I changed psd to pwelch to supress matlab warning
    % I made tests and pwelch and psd return very similar P and f outputs
    % Doing hamming windowing besides boxcar the difference increases
    [P, f] = pwelch(x, boxcar(length(x)), 0, fsamp, fsamp);

    if sum(P)~=0
        num=sum(f.*P);
        den=sum(P);
        fmean(column,1)=num/den;
        den=sum(P);
        k=1;
        while (sum(P(1:k)))<=den/2
            k=k+1;
        end
        i=k;

        if i<length(P)
            alfa1=sum(P(1:i-1))/den;
            if P(i)<=P(i+1),
                radi=roots([abs(P(i+1)-P(i)) 2*P(i) -2*(0.5-alfa1)*den]);
            else
                radi=roots([abs(P(i+1)-P(i)) P(i)+P(i+1) -2*(0.5-alfa1)*den]);
            end
            if radi(1)>0 && radi(1)<1,
                x=radi(1);
            else
                x=radi(2);
            end
            df=f(2)-f(1);
            fmed(column,1)=f(k)+df*x;
        else
            fmean(column,1)=NaN;
            fmed(column,1)=NaN;
        end
    else
        fmean(column,1)=NaN;
        fmed(column,1)=NaN;
    end
end

