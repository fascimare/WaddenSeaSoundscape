%% Pulse rate from 5s clips

% Specify folder with wav files and sample rate
wavpath = 'G:\Shared drives\WaddenSea fish sounds\Sound library\striped fish grunt_2s';

%Output path
outpath = 'G:\My Drive\Studenten\Vassilis Galanos\pulse analysis';


%Spectrogram settings
%dur = 3.5; %spectrogram length (s)
filt_low = 50; %low end of band pass filter (Hz)
filt_high = 2000; %high end of band pass filter (Hz)
peak_thresh = 0.2; %threshold above which peaks are noted.
diff_top = 10; %blanking time between found peaks.
%nwindow = 500; %window size
%nfft = 500; %FFT size
%overlap = 0.9; %window overlap (0-1)
%freq_low = 0; %Lower frequency limit on spectrogram (kHz);
%freq_high = 2; %Upper frequency limit on spectrogram (kHz);
%climval = [50 80]; %dB values for colorbar
calfile= readtable('G:\My Drive\Sound library\Sound types Wadden Sea\Deployment_sensitivity.xlsx'); 
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

for n = 1:length(filelist)
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
    
    %Plot waveform with final tops
    figure;plot(t_a,fishsoundcal_filt)
    hold on
    plot(t_a(corrtops),fishsoundcal_filt(corrtops),'x')
    hold off
    
    %Plot final time gaps
    figure;histogram(t_a(corrtops(2:end))-t_a(corrtops(1:end-1)),'Binwidth',0.001)
    
    %Save calculated pulse time and pulse rate.
    pulsetime = t_a(corrtops(2:end))-t_a(corrtops(1:end-1));
    pulserate = 1./pulsetime;
    allpulsetimes{n} = pulsetime;
    mean_pulsetime(n) = mean(pulsetime);
    mean_pulserate(n) = mean(pulserate);
    soundfiles{n} = soundfile;
    calltype{n} = soundfilesecs{1};

    toc
end

pulses = table(soundfiles,calltype,mean_pulsetime,mean_pulserate);
writetable(pulses,'C:\Users\bill1\Documents\Work\1. University\5. Erasmus\Groningen\Work\Data analysis\Pulse analysis\Matlab results\table_results_striped fish grunt.csv')
