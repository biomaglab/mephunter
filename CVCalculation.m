%% CV SET OF FUNCTIONS

function [cv cc] = CVCalculation(Segna,dint,fsamp)

% Function for the estimation of CV (global or of single MUs) using Maximum Likelihood
% multiple-channel methods [1]. The computational cost is reduced due to the application of the iterative Newton
% method for efficient minimum detection (in a similar way as it is done for the 
% case of two channels [2]).
%
% Input parameters:
%
% Segna		matrix with the EMG signals in raws (no innervation zones; all signals propagate in the same direction)
% dint		interelectrode distance (in meters)
% fsamp		sampling frequency (in Hz)
%
% Output parameters:
%
% cv			the estimated conduction velocity 
% cc			the average correlation coefficient between the channels (after optimal
%				alignment) 
% References:
%
% [1] D. Farina, W. Muhammad, E. Fortunato, O. Meste, R. Merletti, H. Rix, "Estimation of 
% single motor unit conduction velocity from the surface EMG signal detected 
% with linear electrode arrays", Med. Biol. Eng. Comput., vol. 39, pp. 225-236, 2001
%
% [2] Mc Gill K,  Dorfman L. High resolution alignment of sampled waveforms. IEEE Trans Biomed Eng 1984;31:462-70.
%
% Method developed by D. Farina, E. Fortunato, R. Merletti, O. Meste, W. Muhammad, H. Rix
%
% Copyright (2000) LISiN Politecnico di Torino (lisin@polito.it)
% All right reserved
% Author(s): D. Farina and W. Muhammad
% Modified by Marco Gazzoni to return the correlation coefficient value
% Latest update: 20060226


% Set the limits for CV (the coarse delay estimation correspond to 
% CVs within these two limits)
CVmin= 2;
CVmax= 7;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This part of the function determines the propagation direction.
% To do so, the crosscorrelation between the two channels of the matrix 
% is computed and, depending on the sign of the lag corresponding to 
% the maximum crosscorrelation, the signal matrix is flipped or not. 
% In this way, Segna has always signals propagating in a specific direction.

% Compute the crosscorrelation function
corr=xcorr(Segna(2,:),Segna(1,:));
% Determine the peak of the crosscorrelation function
[cc i_cc]=max(corr);
% Flip the signal matrix if the detected lag is negative (i.e., smaller than 
% half of the crosscorrelation vector length)
if i_cc < ceil(length(corr)/2)
   Segna=flipud(Segna);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This part of the function provides a coarse estimation of the delay. 
% This is used as starting point of the iterative Newton procedure. 
% The selection of a starting point close to the true value is 
% important, especially for a large number of channels (refer to [1]). 

% ch is a vector containing two numbers which correspond to the central 
% channels of the array
ch = [ceil(size(Segna,1)/2) ceil(size(Segna,1)/2)+1];
% The coarse delay estimation is provided by a local crosscorrelation 
% method (see function localac.m for details)
start = localac(Segna(ch(1),:), Segna(ch(2),:), dint, fsamp, CVmin, CVmax);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This part of the function recall the main routine which applied the 
% iterative Newton method for the multi-channel delay estimation

% Application of the iterative procedure to estimate CV
[cv, cc] = mle3(Segna,start,dint,fsamp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function teta=localac(sig1,sig2,dint,fsamp, CVmin, CVmax)

% This function is called by mle_CV_est.m.
% The function estimates a coarse delay to start the iterative Newton procedure in
% case of multi-channel CV estimation. The method is based on the crosscorrelation 
% function between two of the channels of the array (inputs) computed for a subset 
% of delays corresponding to "physiological" CV values (2-7 m/s). Note that the CV 
% estimated by the function mle_CV_est is not limited in the range for which 
% the coarse delay is estimated
%
% Input parameters:
%
% sig1		the first signal (row vector)
% sig2		the second signal (row vector)
% dint		interelectrode distance (in meters)
% fsamp		sampling frequency (in Hz)
% CVmin     lower bound for CV value
% CVmax     upper bound for CV value
%
% Output parameters:
%
% teta		the estimated starting point for the iterative Newton method (see [1]) 
%
% References:
%
% [1] D. Farina, W. Muhammad, E. Fortunato, O. Meste, R. Merletti, H. Rix, "Estimation of 
% single motor unit conduction velocity from the surface EMG signal detected 
% with linear electrode arrays", Med. Biol. Eng. Comput., vol. 39, pp. 225-236, 2001
%
%
% Copyright (2000) LISiN Politecnico di Torino (lisin@polito.it)
% All right reserved
% Author(s): D. Farina and W. Muhammad


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This part of the function sets the limits within which the coarse 
% delay estimation is searched.

% Set the limits for CV
% Corresponding limits for the delay teta (they correspond to the limits of CV)
TETAmin=floor(dint/CVmax*fsamp);
TETAmax=ceil(dint/CVmin*fsamp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This part of the function computes the crosscorrelation function 
% only for the lags corresponding to the limits imposed on TETA.

% Compute the crosscorrelation function for the range of time delays selected
for i = TETAmin : TETAmax
   corrloc(i-TETAmin+1)=sum(sig1(1:length(sig1)-i).*sig2(i+1:length(sig2)));
end;

% Vector of the selected delays
delay=[TETAmin:TETAmax];

% Location of the maximum of the crosscorrelation function
[a b]=max(corrloc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This part of the function interpolates around the peak of the 
% crosscorrelation function with a second order polynomial function 
% (parabola) in order to improve the coarse delay estimation.

% Interpolation with a parabola using three points, with the central 
% point being the peak 
% In case the detected peak is at the first or last sample, the 
% interpolation is skipped
if b > 1 & b < length(delay),
[P S]=polyfit(delay([b-1 b b+1]),corrloc([b-1 b b+1]),2);
teta=-P(2)/(2*P(1));
else
   teta=b;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   function [cv, cc] = mle3(Segna,start,dint,fsamp)

% This function is called by mle_CV_est.m.
% The function estimates the delay between propagating signals 
% starting from a coarse estimation of the delay provided by 
% the function "localac.m" (called by mle_CV_est.m). 
% The method is based on the Newton iterative technique which is 
% based on the analytical calculation of the first and second derivative of the
% MLE mean square error function (see [1]).
%
% Input parameters:
%
% Segna		matrix with the EMG signals in raws (no innervation zones; 
%           all signals propagate in the same direction)
% start     the coarse estimation of the delay provided by the function 
%           "localac.m" (which is called by mle_CV_est.m)
% dint		interelectrode distance (in meters)
% fsamp		sampling frequency (in Hz)
%
% Output parameters:
%
% cv		the estimated conduction velocity  (in m/s)
% cc		the average correlation coefficient between the channels (after optimal
%			alignment) 
%
% References:
%
% [1] D. Farina, W. Muhammad, E. Fortunato, O. Meste, R. Merletti, H. Rix, "Estimation of 
% single motor unit conduction velocity from the surface EMG signal detected 
% with linear electrode arrays", Med. Biol. Eng. Comput., vol. 39, pp. 225-236, 2001
%
%
% Copyright (2000) LISiN Politecnico di Torino (lisin@polito.it)
% All right reserved
% Author(s): D. Farina and W. Muhammad

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inizialize the starting point
t=start;
teta=10;

% Inizialize the number of trials (the method stops after a 30 iterations)
trial=0; 

% Number of signals
num_sig=size(Segna,1);

% While cycle which constitutes the Newton method 
% The method stops when the new estimate of delay (in samples) is closer than 5e-5
% to the previous estimation or when the number of iterations is 30
while (abs(teta-t)>=5e-5) & trial < 30
    % Up-date the iterations
    trial=trial+1;
    % Up-date the estimate of teta    
    teta=t;     
    % Inizialization of the first (de1) and second (de2) derivative of the mean square error    
    de1=0;
    de2=0;
    % This cycle compute the two derivatives calling the function "derivbeam.m" (see also [1])
    for i=1:num_sig,
        [de1t, de2t]=derivbeam(Segna,i,teta);
        de1=de1+de1t;
        de2=de2+de2t;
    end;
    % Newton's criteria
    if (de2>0) 
        u=-de1/de2;
        if (abs(u)>0.5)
            u=-0.5*abs(de1)/de1;
        end
    else
        u=-0.5*abs(de1)/de1;
    end
    % Up-date the delay estimation    
    t=teta+u;
end;

% CV computed from the delay (CV is in m/s while teta is in samples)
cv=abs(dint/(teta/fsamp));

%calc CC
for i = 1 : num_sig
   sig_shift(i,:) = freshift(Segna(i,:),(i-1)*teta);   
end;

sig_aver = sum(sig_shift);
for i = 1 : num_sig
   A = corrcoef(sig_aver,sig_shift(i,:));
   cc_shift(i) = A(1,2);   
end;
cc = real(mean(cc_shift));

function [de1, de2]=derivbeam(Segna,posref,teta)

% This function is called by mle3.m, which is called by mle_CV_est.m.
% The function provides the derivative of the multi-channel mean square error function
% using the expression provided in [1]. Note that the function provides the mean square error 
% of the beamforming selecting a reference signal. The MLE mean square error is the summation of the 
% beamforming mean square errors selecting all the channels of the array as reference signals. Thus, the 
% derivative is the summation of the derivatives of the MSE of the beamformings (this summation is 
% made in the function mle3.m).
%
% Input parameters:
%
% Segna		matrix with the EMG signals in raws (no innervation zones; all signals propagate in the same direction)
% posref	the number corresponding to the reference signal for the beamforming mean square error
% teta		the delay (in samples) for which the mean square error and its derivatives are computed
%
% Output parameters:
%
% de1		first derivative of the mean square error 
% de2		second derivative of the mean square error 
%
% References:
%
% [1] D. Farina, W. Muhammad, E. Fortunato, O. Meste, R. Merletti, H. Rix, "Estimation of 
% single motor unit conduction velocity from the surface EMG signal detected 
% with linear electrode arrays", Med. Biol. Eng. Comput., vol. 39, pp. 225-236, 2001
%
%
% Copyright (2000) LISiN Politecnico di Torino (lisin@polito.it)
% All right reserved
% Author(s): D. Farina and W. Muhammad


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of signals
num_sig=size(Segna,1);
M=num_sig-1;

% Number of samples in the signals
N=size(Segna,2);

% Define a vector which contains the numbers corresponding to all the channels except the reference one
pos=find([1:M+1]~=posref);

% Matrix with the reference signal as first row and all the other signals as the following rows
Htemp=[Segna(posref,:);Segna(pos,:)];
Segna=Htemp;

% Define the positions of the signals in the array relative to the reference channel position
position(1)=posref;
for i=2:posref
   position(i)=i-(posref+1);
end;
for i=1:num_sig-posref
   position(i+posref)=i;
end;

position=position(2:num_sig);

% Define the frequency vector (bins)
k=[1:N/2];

% Inizialize the terms used for the expression of the derivatives
terminede1=zeros(1,length(k));
terminede2=zeros(1,length(k));
terminede12=zeros(1,length(k));
terminede22=zeros(1,length(k));

% Compute the Fourier transform of the signals of the array
for i=1:num_sig
   SEGNAfft(i,:)=fft(Segna(i,:));
end;

% Compute the first derivative as the summation of two terms, defined in [1]
for i=1:M
   for u=i+1:M
         terminede1=terminede1-imag(SEGNAfft(i+1,k+1).*conj(SEGNAfft(u+1,k+1)).*exp(j*2*pi*k*(position(i)-position(u))*teta/N).*2*pi.*k.*(position(i)-position(u))./N);
   end;
end;
terminede1=terminede1*2/(M^2);
for i=1:M
   terminede2=terminede2+SEGNAfft(i+1,k+1).*exp(j*2*pi*k*position(i)*teta/N).*2*pi.*k.*position(i)./N;
end;
terminede2=2*imag(conj(SEGNAfft(1,k+1)).*terminede2)/M;
de1=2/N*sum(terminede1+terminede2);

% Compute the second derivative as the summation of two terms, defined in [1]
for i=1:M
   for u=i+1:M
         terminede12=terminede12-real(SEGNAfft(i+1,k+1).*conj(SEGNAfft(u+1,k+1)).*exp(j*2*pi*k*(position(i)-position(u))*teta/N).*((2*pi*k*(position(i)-position(u))/N).^2));
   end;
end;
terminede12=terminede12*2/(M^2);

for i=1:M
   terminede22=terminede22+SEGNAfft(i+1,k+1).*exp(j*2*pi*k*position(i)*teta/N).*((2*pi*k*position(i)/N).^2);
end;
terminede22=2*real(conj(SEGNAfft(1,k+1)).*terminede22)/M;
de2=2/N*sum(terminede12+terminede22);


%Called by Segdef.m
%Translate the signal of a specified delay
%The shift is made in the frequency domain
%
%Input:
%       seg:   signal to shift: line vector
%       teta:  delay
%       fszamp: sampling frequency
%
%Output: 
%       segt: translated signal: line vector
%
%Copyright (1999) LISiN Politecnico di Torino (lisin@polito.it)
%All right reserved
%Author(s): Dario Farina

function segt=freshift(seg,teta,fsamp)

SEG=fft(seg);

f=fftshift([-0.5:1/(length(seg)):0.5-1/(length(seg))]);

SEGt=SEG.*exp(1i*2*pi*teta*f);
segt=ifft(SEGt);

function c = xcorr(a, b, option)
%XCORR	Cross-correlation function estimates.
%	XCORR(A,B), where A and B are length M vectors, returns the
%	length 2*M-1 cross-correlation sequence in a column vector.
%	XCORR(A), when A is a vector, is the auto-correlation sequence.
%	XCORR(A), when A is an M-by-N matrix, is a large matrix with
%	2*M-1 rows whose N^2 columns contain the cross-correlation
%	sequences for all combinations of the columns of A.
%	The zeroth lag of the output correlation is in the middle of the 
%	sequence, at element or row M.
%	By default, XCORR computes a raw correlation with no normalization.
%	XCORR(A,'biased') or XCORR(A,B,'biased') returns the "biased"
%	estimate of the cross-correlation function.  The biased estimate
%	scales the raw cross-correlation by 1/M.
%	XCORR(...,'unbiased') returns the "unbiased" estimate of the
%	cross-correlation function.  The unbiased estimate scales the raw
%	correlation by 1/(M-abs(k)), where k is the index into the result.
%	XCORR(...,'coeff') normalizes the sequence so that the
%	correlations at zero lag are identically 1.0.
%	See also XCOV, CORRCOEF, CONV and XCORR2.

%	Author(s): L. Shure, 1-9-88
%		   L. Shure, 4-13-92, revised
%	Copyright (c) 1984-94 by The MathWorks, Inc.
%	$Revision: 1.7 $  $Date: 1994/01/25 18:00:07 $

%	References:
%	  [1] J.S. Bendat and A.G. Piersol, "Random Data:
%	      Analysis and Measurement Procedures", John Wiley
%	      and Sons, 1971, p.332.
%	  [2] A.V. Oppenheim and R.W. Schafer, Digital Signal 
%	      Processing, Prentice-Hall, 1975, pg 539.

onearray = 1;
if nargin == 1
	option = 'none';
	if min(size(a)) == 1	% a is a vector
		a = [a(:) a(:)];
	else
		onearray = 0;
	end
elseif nargin == 2
	if isstr(b)
		option = b; clear b
		na = max(size(a));
		if min(size(a)) == 1	% a is a vector
			a = [a(:) a(:)];
		else	% a is a matrix
			onearray = 0;
			[m,n] = size(a);
		end
    else	% b is truly a second arg
		if min(size(a)) ~= 1 & min(size(b)) ~= 1
			error('You may only specify 2 vector arrays.')
		end
		option = 'none';
		onearray = 2;
	end
else
	if max(size(a)) ~= max(size(b)) & ~strcmp(option,'none')
		error('OPTION must be ''none'' for different length vectors A and B')
	end
	onearray = 2;
end
% check validity of option
nopt = nan;
if strcmp(option, 'none')
	nopt = 0;
elseif strcmp(option, 'coeff')
	nopt = 1;
elseif strcmp(option, 'biased')
	nopt = 2;
elseif strcmp(option, 'unbiased')
	nopt = 3;
end
if isnan(nopt)
	error('Unknown OPTION')
end
if onearray == 2
	[ar,ac] = size(a);
	na = max([ar ac]);
	nb = max(size(b));
	if na > nb
		b(na) = 0;
	elseif na < nb
		a(nb) = 0;
	end
	a = [a(:) b(:)];
end
[nr, nc] = size(a);
nsq  = nc^2;
mr = 2 * nr - 1;
c = zeros(mr,nsq);
ci = zeros(1,nsq);
cj = ci;
nfft = 2^nextpow2(2*nr);
for i = 1:nc
	atmpi = a(:,i);
	if ~any(any(atmpi))
		real1 = 1;
	else
		real1 = 0;
	end
	atmpi = fft([atmpi(:); zeros(nfft-nr,1)]);
	for j = 1:i
		col = (i-1)*nc+j;
		colaux = (j-1)*nc+i;
		tmp = fft([a(:,j); zeros(nfft-nr,1)]); % pad with zeros for fft
		tmp = fftshift(ifft(atmpi.*conj(tmp)));
		c(:,colaux) = tmp((1:mr)+nfft/2-nr+1);
		ci(col) = i;
		cj(col) = j;
		ci(colaux) = j;
		cj(colaux) = i;
		if ~any(any(imag(a(:,j)))) & real1
			c(:,colaux) = real(c(:,colaux));
		end
		if i~= j
			c(:,col) = conj(c(mr:-1:1,colaux));
		end
	end
end
if nopt == 1	% return normalized by sqrt of each autocorrelation at 0 lag
% do column arithmetic to get correct autocorrelations
	cdiv = ones(mr,1)*sqrt(c(nr,1+(ci-1)*(nc+1)).*c(nr,1+(cj-1)*(nc+1)));
	c = c ./ cdiv;
elseif nopt == 2	% biased result, i.e. divide by nr for each element
	c = c / nr;
elseif nopt == 3	% unbiased result, i.e. divide by nr-abs(lag)
	c = c ./ ([1:nr (nr-1):-1:1]' * ones(1,nsq));
end
if onearray == 1
	c = c(:,1);	% just want the autocorrelation
	[am, an] = size(a);
	if am == 1
		c = c.';
	end
elseif onearray == 2	% produce only cross-correlation
	c = c(:,2);
	if ar == 1
		c = c.';
	end
end
if ~any(any(imag(a)))
	c = real(c);
end

% CV SET OF FUNCTIONS --------
