% 3b. Function to get power in swarm script
function attn_pow_biowolf_v1(folder_name, ev_type)
% originally: FUNCTION [phase,pow] = multiphasevec3(f,S,Fs,width, silent)
% Returns the phase and power as a function of time for a range of
% frequencies (f).
%
% INPUT ARGS:
%   folder_name = file with signal to process for given patient (trials X samples)
%   output = output directory (saved power)
% HARD CODED:
%   Fs = 1000;      % Sampling frequency
%   width = 6;     % Width of Morlet wavelet (>= 5 suggested).
%
% OUTPUT ARGS:
%   power- Power data [trials,freqs,time]

f = logspace(log10(2),log10(32),200);
load(fullfile(folder_name, "sess_struct_pc.mat"));
Fs = 1000;
width = 6;

    switch ev_type
        case 'pres'
            elec = sess_struct.raw_pres; % shape: elec, events, time
        case 'recog'
            elec = sess_struct.raw_recog;
        case 'fr'
            elec = sess_struct.raw_fr;
    end
    
    nF = length(f);
    
    dt = 1/Fs;
    st = 1./(2*pi*(f/width));
    power = NaN(size(elec,1),size(elec,2),length(f));
    
for el  = 1:size(elec,1)
    for ev = 1:size(elec,2)
    
        S = squeeze(elec(el,ev,:))';
    
        if isvector(S)
            S = S(:)';
        end
        nS = size(S);
    
        % Pre-allocate output for efficiency
        curWaves = cell(1, nF);
    
        % Loop through each frequency index
        for i = 1:nF
            % Define time range based on st(i)
            t = -3.5 * st(i):dt:3.5 * st(i);
            sf = f(i)/width;
            st_t = 1/(2*pi*sf);
            A = 1/sqrt(st_t*sqrt(pi));
            
            curWaves{i} = A*exp(-t.^2/(2*st_t^2)).*exp(1i*2*pi*f(i).*t);
        end
    
        nCurWaves = cellfun( @(w) length(w), curWaves );
        
        Lys = nS(2) + nCurWaves - 1;    % length of convolution of S and curWaves{i}
        % Ly2s = pow2(nextpow2(Lys));     % next power of two (for fft)
        % Changing this to make Ly2s a vector of length nF
        
        for i=1:length(Lys)
          Ly2s(i)=pow2(nextpow2(Lys(i))); % next power of two (for fft)
        end
        ind1 = ceil(nCurWaves/2);       % start index of signal after convolution
        
        pow = zeros(nS(1), nF, nS(2));
        
        for i = 1:nF
            Ly2 = Ly2s(i);
            
            %%% Perform the convolution of curWaves{i} with every row of S
            % take the fft of S and curWaves{i}, multiply them, and take the ifft
            
            % Sfft = fft(S,Ly2,2);
            % curWaveFFT = fft(curWaves{i},Ly2);
            % Y = bsxfun(@times, Sfft, curWaveFFT);
            % y1 = ifft(Y,Ly2,2);
            
            % (EH - it's much quicker to do it in one line)
            y1 = ifft( bsxfun( @times, fft(double(S),Ly2,2), fft(curWaves{i},Ly2) ) ,Ly2,2);
        
            y1 = y1( :, ind1(i):(ind1(i)+nS(2)-1) );
            
            % find phase and power (do it inside this loop to save memory)
            pow(:,i,:) = abs(y1).^2;
        end
        
        power(el,ev,:) = nanmean(log10(pow),3);
    end
end

out_dir = fullfile(folder_name, 'power');

if ~exist(out_dir, 'dir') % make folder if doesn't exist already
    mkdir(out_dir);
end

save(fullfile(out_dir, ev_type), 'power', '-v7.3');

clear;
pause(5);

