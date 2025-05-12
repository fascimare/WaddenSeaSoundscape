%% Pulse rate from 5s clips

% Specify folder with wav files and sample rate
wavpath = 'G:\My Drive\Wadden fish sounds\Data analysis\Spectrograms\Ti\';

%Output path
outpath = 'G:\My Drive\Wadden fish sounds\Data analysis\Pulse rate results\';

%Spectrogram settings
%dur = 3.5; %spectrogram length (s)
filt_low = 50; %low end of band pass filter (Hz)
filt_high = 11000; %high end of band pass filter (Hz)
peak_thresh = 0.2; %threshold above which peaks are noted.
diff_top = 10; %blanking time between found peaks.
%nwindow = 500; %window size
nfft = 3500; %FFT size
window=hanning(nfft);
noverlap = 0.9*nfft;
% overlap = 0.9; %window overlap (0-1)
%freq_low = 0; %Lower frequency limit on spectrogram (kHz);
%freq_high = 2; %Upper frequency limit on spectrogram (kHz);
%climval = [50 80]; %dB values for colorbar
calfile= readtable('G:\Shared drives\WaddenSea fish sounds\Deployment_sensitivity.xlsx'); 
%gain = "SensitivityHighGain";
%calgain = calfile(:,1:2); %Adjust column numbers based on which gain you used.

%%  Load sound clips and measure pulse rate
filelist = dir(fullfile(wavpath,'*.wav'));

%Create empty variable for pulse rate
mean_pulserate = zeros(length(filelist),1);
mean_pulsetime = zeros(length(filelist),1);
allpulsetimes = cell(length(filelist),1);
soundfiles = cell(length(filelist),1);
calltype = cell(length(filelist),1);
peakFrequency = cell(length(filelist),1);
minFrequency = cell(length(filelist),1);
maxFrequency = cell(length(filelist),1);
SLres = cell(length(filelist),1);
RLres = cell(length(filelist),1);
down_limit=cell(length(filelist),1);
upper_limit=cell(length(filelist),1);
%%
for n = 76:length(filelist)
    tic
    %Find appropriate wav file
    soundfile = filelist(n).name;
    fileinfo = audioinfo(fullfile(wavpath,soundfile));
    fs = fileinfo.SampleRate;
    fishsound = audioread(fullfile(wavpath,soundfile));
    
    %Calibration
    soundfilesecs = strsplit(soundfile,'_');
    if contains(soundfilesecs{2},'OEST')
        dep = [soundfilesecs{2},'_',soundfilesecs{3}];
    else
        dep = [soundfilesecs{2},'_',soundfilesecs{3},'_',soundfilesecs{4}];
    end
    
    cal = calfile.Sensitivity(strcmp(dep,calfile.Deployment));
    fishsoundcal = fishsound*power(10,cal/20);
    t_a = 0:1/fs:(length(fishsoundcal)-1)/fs;

    %Filter data to make figure look nicer
    [b,a] = ellip(4,0.1,40,[filt_low,filt_high]*2/fs);
    fishsoundcal_filt = filter(b,a,fishsoundcal);

    %Normalize wave form
    idmin = find(fishsoundcal_filt<0);
    idmax = find(fishsoundcal_filt>0);

    fishnorm = fishsoundcal_filt;
    fishnorm(idmin) = fishnorm(idmin)/min(fishsoundcal_filt)*-1;
    fishnorm(idmax) = fishnorm(idmax)/max(fishsoundcal_filt);
    
    %Calculate slope to find where the slope changes from positive to
    %negative.
    fishder = [0;fishnorm(2:end)-fishnorm(1:end-1)]; %Derivative gives you the slope
    fishder2=fishder(2:end).*fishder(1:(end-1)); %if the slope turns from positive to negative or vice versa, the product of the two succeeding values in the column will be negative.
    tops = find(fishder2<0); %Find all tops and valleys.
    
    peaks = find(fishnorm(tops)>peak_thresh); % Select only tops that are above a certain threshold.
    noclippeak = peaks(peaks>1000); % Remove the first 1000 tops, they are from the initial clipping of the sound file (because sound cannot suddenly start).
   
    %figure;histogram(tops(noclippeak(2:end))-tops(noclippeak(1:end-1)),'BinWidth',1)
    %histogram to find optimal blanking time

    gaps = tops(noclippeak(2:end))-tops(noclippeak(1:end-1)); %Calculate gaps (in samples) between found tops.
    blktime = find(gaps>30); % Remove tops that are too close to each other (within blanking time).
    
    %Create a new index value with just the tops that we want to keep
    corrtops = tops(noclippeak(blktime));
    
%     %Plot waveform with final tops in one window
%     nexttile
%     plot(t_a,fishsoundcal_filt)
%     hold on
%         plot(t_a(corrtops),fishsoundcal_filt(corrtops),'x')
%         title('Final tops graph')
%     hold off
%     nexttile
%     histogram(t_a(corrtops(2:end))-t_a(corrtops(1:end-1)),'Binwidth',0.001)
%     title('Slope histogram')
% 
%     %Plot waveform with final tops in one window individually
%     figure;plot(t_a,fishsoundcal_filt)
%     hold on
%     plot(t_a(corrtops),fishsoundcal_filt(corrtops),'x')
%     title(filelist(n).name)
%     hold off
%     figure;histogram(t_a(corrtops(2:end))-t_a(corrtops(1:end-1)),'Binwidth',0.001)

%Measure peakF and dB
 %nfft = fs; % fft of fs/2 equals bin size of 1 Hz.
% fs = 2000;
soundname=soundfilesecs{1};

% Bandwidth determination
if strcmp(soundname,'water croacking')
   downlim=700;
   uplim=2000;
elseif strcmp(soundname,'coned fish grunt')
   downlim=150;
   uplim=2000;
elseif strcmp(soundname,'Gr')
   downlim=100;
   uplim=2000;
elseif strcmp(soundname,'solid fish grunt')||strcmp(soundname,'striped fish grunt')
   downlim=80;
   uplim=2000;
elseif strcmp(soundname,'Chu')||strcmp(soundname,'Sn')||strcmp(soundname,'Ti')
    downlim=50;
    uplim=11000;
else
    downlim=11;
    uplim=2000;
end


%w= window length. Default is nfft.
%nov= number of samples to overlap. Default is half of window length.
%nfft = 2000; %Will have to be lower than fs/2 because of cutting down selection to just call (window not long enough for larger nfft).
%noverlap = 0.9*nfft;
%window=hanning(nfft);
call = min(corrtops):max(corrtops); % getting timeseries for just the call to improve accuracy of spectra.
if length(call)<nfft
    call = round((min(corrtops)-(nfft/2-length(call)/2))):round(max(corrtops)+(nfft/2-length(call)/2));
end
[SL,f] = pwelch(fishsoundcal(call),window,noverlap,nfft,fs);
SL = 10*log10(SL); % counts^2/Hz
bindownlim = round(downlim*nfft/fs)+1;
binuplim = round(uplim*nfft/fs)+1;
[RLmax,idx] = max(SL(bindownlim:binuplim)); %Find peak intensity
peakF = f(idx + bindownlim-1); %Find peak frequency
peakidx = idx+bindownlim-1;
%[RLmin,idx] = min(SL(downlim:uplim)); %Select Hz bandwidth
%minF = idx + downlim-1;
RL10=find(SL(bindownlim:peakidx)>(RLmax-20));
minF = f(RL10(1)+bindownlim-1);
RL90=find(SL(peakidx:binuplim)>(RLmax-20));
maxF = f(RL90(1)+peakidx-1);
% [medianRL,idx] = median(SL(downlim:uplim)); %Select Hz bandwidth
% medianF = idx + downlim-1;
spectra{n,1} = [f,SL];

    %Save calculated pulse time and pulse rate.
    pulsetime = t_a(corrtops(2:end))-t_a(corrtops(1:end-1));
    pulserate = 1./pulsetime;
    allpulsetimes{n} = pulsetime;
    mean_pulsetime(n) = mean(pulsetime);
    mean_pulserate(n) = mean(pulserate);
    soundfiles{n} = soundfile;
    calltype{n} = soundfilesecs{1};
    peakFrequency{n} = peakF;
    minFrequency{n} = minF;
    maxFrequency{n} = maxF;
    %medianFrequency{n} = medianF;
    RLres{n} = RLmax;
    down_limit{n}=downlim;
    upper_limit{n}=uplim;
    %SLres{n} = SL;

    toc
end
%     medianFrequency=medianFrequency.';
%     calltype= strrep(calltype,'water croacking','w_c');
%     calltype= strrep(calltype,'coned fish grunt','c_f_g');
%     calltype= strrep(calltype,'solid fish grunt','sol_f_g');
%     calltype= strrep(calltype,'striped fish grunt','str_f_g');
    
save SL.mat;
pulses = table(soundfiles,calltype,mean_pulsetime,mean_pulserate,RLres,minFrequency,maxFrequency,down_limit,upper_limit);
writetable(pulses,[outpath,'Ti_pulses.csv'],'Delimiter',',','QuoteStrings','none')

