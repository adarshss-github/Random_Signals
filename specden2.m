function [Y12,f] = specden2(y1,y2,N,fs)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
% By Adarsh S., Ph.D. Candidate IIT Kanpur and Chikku, The love of my life
%[Y12,f] =  specden2(y1,y2,N,fs)
%
%Function description:
%-----------------------
%Estimates the double-sided cross spectral densities of two signals y1 and y2 by
%averaging scaled FFTs
%
%Input arguments:
%----------------
%y1 and y2 are the two signals
%N is the number of points in each segment (Should be a power of 2)
%fs is the sampling rate of the measurements
%
%Ouput arguments:
%----------------
%Y12 is the cross spectral density (unit^2/Hz)
%f is the frequency axis (Hz)
%
%Note:
%------
%y1 and y2 should be row vectors
%If y1 and y2 are same signal, it estimates the power spectral density
% Windows were not used. So, care should be taken while choosing N (suffiently large) to
% mitigate leakage
%
%Example:
%-------
%Ex: [Y12,f] =  specden2(y1,y2,1024,200) ;
%
%
%
%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; October, 2019 || *********
%                                      -----------------------------------------------
%
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================

dt = 1/fs ;
df = 1/(N*dt) ;
[my1,ny1] = size(y1) ;
[my2,ny2] = size(y2) ;
r = 2^(nextpow2(N)) ;

%============================================================================================================================


if r~=N
    
    error('N should be a power of 2');
    
end

if my1~=1 || my2~=1
    
    error('Signals not row vectors') ;
    
end

if ny1~=ny2
    
    error('Signals not of same size') ;
    
end

%============================================================================================================================

nseg = floor(ny1/N) ;
Y12dum = zeros(1,N) ;

for i = 1:1:nseg
    
    Y1 = sqrt(dt/N) * ( fft( y1( N*(i-1)+1:N*i ) ) ) ;
    Y2 =  sqrt(dt/N) * ( fft( y2( N*(i-1)+1:N*i ) ) ) ;
    Y12dum = Y12dum + (Y1).*conj(Y2) ;
    
end

Y12 = Y12dum(1:N/2+1)/nseg ;
f = (0:1:N/2)*df ;

%============================================================================================================================

end

