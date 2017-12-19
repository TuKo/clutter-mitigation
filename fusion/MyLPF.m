function b=MyLPF(N,F6dB)
%   Linear-phase filter design for FIR with order N. 
%   b=MyLPF(N,F6dB) returns a length N+1 linear phase (real, symmetric
%   coefficients) FIR filter which has the best approximation to the
%   rectangular frequency response and cutoff at F6dB.
%

N = N+1;             % filter length
L=(N-1)/2;
Nodd = rem(N,2);

    % basis vectors are cos(2*pi*m*f) (see m below)
    if ~Nodd
        m=(0:L)+.5;   % type II
    else
        m=(0:L);      % type I
    end
    k=m';
    if Nodd
        k=k(2:length(k));
        b0 = F6dB/2;       %  first entry must be handled separately (where k(1)=0)
    end;
    b = F6dB/2*sinc(k*F6dB);
    if Nodd
        b=[b0; b];
    end;

    a=4*b;
    if Nodd
        a(1) = a(1)/2;
    end
    if Nodd
        h=[a(L+1:-1:2)/2; a(1); a(2:L+1)/2].';
    else
        h=.5*[flipud(a); a].';
    end;

%windowing
Wind = hamming(N); %N here is the filter length already
% to use Kaiser window, beta must be supplied
% att = 60; % dB of attenuation desired in sidelobe
% beta = 0.1102*(att-8.7);
% wind = kaiser(L,beta);
h = h.*Wind(:)'; 

%scaling to unity gain at DC
b = h / sum(h);

end


% function w = sym_window(n,window)
% %SYM_WINDOW   Symmetric generalized cosine window.
% %   SYM_WINDOW Returns an exactly symmetric N point generalized cosine 
% %   window by evaluating the first half and then flipping the same samples
% %   over the other half.
% 
% if ~rem(n,2)
%     % Even length window
%     half = n/2;
%     w = calc_window(half,n,window);
%     w = [w; w(end:-1:1)];
% else
%     % Odd length window
%     half = (n+1)/2;
%     w = calc_window(half,n,window);
%     w = [w; w(end-1:-1:1)];
% end

% function w = calc_window(m,n,window)
% %CALC_WINDOW   Calculate the generalized cosine window samples.
% %   CALC_WINDOW Calculates and returns the first M points of an N point
% %   generalized cosine window determined by the 'window' string.
% 
% % For the hamming and blackman windows we force rounding in order to achieve
% % better numerical properties.  For example, the end points of the hamming 
% % window should be exactly 0.08.
% 
% switch window
% case 'hann'
%     % Hann window
%     %    w = 0.5 * (1 - cos(2*pi*(0:m-1)'/(n-1))); 
%     a0 = 0.5;
%     a1 = 0.5;
%     a2 = 0;
%     a3 = 0;
%     a4 = 0;
% case 'hamming'
%     % Hamming window
%     %    w = (54 - 46*cos(2*pi*(0:m-1)'/(n-1)))/100;
%     a0 = 0.54;
%     a1 = 0.46;
%     a2 = 0;
%     a3 = 0;
%     a4 = 0;
% case 'blackman'
%     % Blackman window
%     %    w = (42 - 50*cos(2*pi*(0:m-1)/(n-1)) + 8*cos(4*pi*(0:m-1)/(n-1)))'/100;
%     a0 = 0.42;
%     a1 = 0.5;
%     a2 = 0.08;
%     a3 = 0;
%     a4 = 0;
% case 'flattopwin'
%     % Flattop window
%     % Coefficients as defined in the reference [1] (see flattopwin.m)
%     a0 = 0.21557895;
%     a1 = 0.41663158;
%     a2 = 0.277263158;
%     a3 = 0.083578947;
%     a4 = 0.006947368;
% end
% 
% x = (0:m-1)'/(n-1);
% w = a0 - a1*cos(2*pi*x) + a2*cos(4*pi*x) - a3*cos(6*pi*x) + a4*cos(8*pi*x);
% end
