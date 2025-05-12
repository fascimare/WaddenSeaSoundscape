%% Script to calculate pulse rate from fish grunts

% ----------Start with input settings:------------

% Specify folder with wav files and sample rate
wavpath = 'D:\WaddenSea fish sounds\Lauwersoog recordings\February\';

%Output path
outpath = 'C:\Users\P310512\Documents\Groningen\Sound library\Fish sound clustering\TPWS library\';

% %Decide which logged call type you want to calculate the pulse rate for
calltype = 'drum vibrations';
%mkdir(outpath,calltype)

%Spectrogram settings
timebef = 0; %time before call (s)
dur = 3.5; %spectrogram length (s)
filt_low = 75; %low end of band pass filter (Hz)
filt_high = 500; %high end of band pass filter (Hz)
nwindow = 500; %window size
nfft = 500; %FFT size
overlap = 0.9; %window overlap (0-1)
freq_low = 0; %Lower frequency limit on spectrogram (kHz);
freq_high = 2; %Upper frequency limit on spectrogram (kHz);
climval = [50 80]; %dB values for colorbar
calfile= readtable('G:\My Drive\Sound library\Sound types Wadden Sea\Hydrophone calibration.xlsx'); 
%gain = "SensitivityHighGain";
calgain = calfile(:,1:2); %Adjust column numbers based on which gain you used.

%% Read in Triton log file
analysismethod = 'Triton';
[infile,inpath]=uigetfile('*.xlsx','Select a file with manual picks');
if isequal(infile,0)
    disp('Cancelled button pushed');
    return
end
logs = readtable([inpath,infile]);
%% Read in Raven detection files
analysismethod = 'Raven';
[infile,inpath]=uigetfile('*.txt','Select a Raven file with manual picks');
if isequal(infile,0)
    disp('Cancelled button pushed');
    return
end
Rlogs = readtable([inpath,infile]);
Rlogs.StartTime = datetime(Rlogs.BeginDateTime,"InputFormat","uuuu/MM/dd HH:mm:ss.SSS");
Rlogs.dur = Rlogs.EndTime_s_-Rlogs.BeginTime_s_;
numend = datenum(Rlogs.StartTime) + Rlogs.dur/3600/24;
Rlogs.EndTime = datetime(numend,"ConvertFrom","datenum");
Rlogs.Comments = Rlogs.SoundCode; %sound_code for other file, consistency needed.
depinfo = split(infile,'.');
Rlogs.Location = repmat(string(depinfo{2}),length(Rlogs.StartTime),1);
Rlogs.Deployment = repmat(string(depinfo{1}),length(Rlogs.StartTime),1);
Rlogssub = Rlogs(strcmp(Rlogs.View,'Spectrogram 1'),:);
logs = Rlogssub;

%%  Find appropriate call type and make figure
Index = find(contains(logs.Comments,calltype));
calltypelogs = logs(Index,:);

%Convert start and end times so they can be used by Matlab
starttime = datetime(calltypelogs.StartTime);
startnum = datenum(starttime);
endtime = datetime(calltypelogs.EndTime);
endnum = datenum(endtime);

%Create empty variable for pulse rate
mean_pulserate = zeros(length(startnum),1);

for n = 1:length(startnum)
    tic
    %Find appropriate wav file and calculate start time in samples
    if strcmp(analysismethod, 'Raven')
        fullwavpath = wavpath;
    else
    fullwavpath = [wavpath,calltypelogs.Location{n},'_',calltypelogs.Deployment{n}];
    end
    %Make a list of the wav files in the folder
    wavfiles = dir(fullfile(fullwavpath,'*.wav'));
    filelistwav = struct2cell(wavfiles);
    wavinfo = audioinfo(fullfile(fullwavpath,wavfiles(1).name));
    RegDate = '(?<yr>\d\d)(?<mon>\d\d)(?<day>\d\d)(?<hr>\d\d)(?<min>\d\d)(?<s>\d\d)';
    fileDates = dateregexp(filelistwav(1,:),RegDate);
    wavidx = find(fileDates<= startnum(n) & (fileDates+ datenum([0 0 0 0 0 wavinfo.Duration])) > startnum(n));

     soundfile = wavfiles(wavidx).name;
    wavpathfile = [fullwavpath,'\',wavfiles(wavidx).name];
    fileinfo = audioinfo(wavpathfile);
    fs = fileinfo.SampleRate;

    if strcmp(analysismethod, 'Triton')
        soundDate = dateregexp(soundfile,RegDate);
        startsample = (startnum(n)-soundDate)*3600*24*fs;
        startS = round(startsample - timebef*fs);
        %endS = round(startS+dur*fs);
        endsample = (endnum(n)-soundDate)*3600*24*fs;
        endS = round(endsample);
    elseif strcmp(analysismethod,'Raven')
        startS = round(calltypelogs.BegFileSamp_samples_(n)-timebef*fs);
        endS = round(startS+dur*fs);
    end

    if endS>fileinfo.TotalSamples
        timeleft = endS-fileinfo.TotalSamples;
        endS = fileinfo.TotalSamples;
    end
    fishsound = audioread(wavpathfile,[startS,endS]);
    if exist("timeleft","var")
        %idx = find(strcmp({wavfiles.name},soundfile)==1);
        nextfile = wavfiles(wavidx+1).name;
        nextfilepath = [fullwavpath,'\',nextfile];
        remainsound = audioread(nextfilepath,[1,timeleft]);
        fishsound = vertcat(fishsound,remainsound);
        clear timeleft
    end

    %Calibration
    soundfilesecs = strsplit(soundfile,'.');
    hyd = str2num(soundfilesecs{1});
    cal = table2array(calgain(calgain.Hydrophone == hyd,2));
    fishsoundcal = fishsound*power(10,cal/20);
    t_a = 0:1/fs:(length(fishsoundcal)-1)/fs;

    %Filter data to make figure look nicer
    [b,a] = ellip(4,0.1,40,[filt_low,filt_high]*2/fs);
    fishsoundcal_filt = filter(b,a,fishsoundcal);

    %Settings for spectrogram
    window=hanning(nfft);
    noverlap = round(nfft*overlap);

    %Make figure
    %h= tiledlayout(2,1);
    %nexttile
%     figure;
%     spectrogram(fishsoundcal_filt,window,noverlap,nfft,fs,'yaxis')
%     clim(climval);
%     ylim([freq_low freq_high])
%     ylabel('Frequency (kHz)')
%     fontsize(gca,16, "points")
%     title(calltype,'FontSize',24)
% 
%     %nexttile
%     figure;
%     plot(t_a,fishsoundcal_filt)
    peaks = find(fishsoundcal_filt>100000);
    timedif = peaks(2:end)-peaks(1:end-1);
    gaps = find(timedif>1);
    newpeak = [];
    for g=2:length(gaps)
        sub = peaks(gaps(g-1):gaps(g)-1);
        maxpeak = find(fishsoundcal_filt(sub)==max(fishsoundcal_filt(sub)));
        newpeak(g) = sub(maxpeak);
    end
%     hold on
%     plot(t_a(newpeak(2:end)),fishsoundcal_filt(newpeak(2:end)),'*')
    pulsetime = t_a(newpeak(3:end))-t_a(newpeak(2:end-1));
    pulserate = 1./pulsetime;
    mean_pulserate(n) = mean(pulserate);

end

%% Pulse rate from 5s clips

% Specify folder with wav files and sample rate
wavpath = 'C:\Users\P310512\Documents\Groningen\Sound library\Fish sound clustering\test for pulse rate\';

%Output path
outpath = 'C:\Users\P310512\Documents\Groningen\Sound library\Fish sound clustering\TPWS library\';


%Spectrogram settings
%dur = 3.5; %spectrogram length (s)
filt_low = 100; %low end of band pass filter (Hz)
filt_high = 2000; %high end of band pass filter (Hz)
peak_thresh = 0.2; %threshold above which peaks are noted.
diff_top = 5; %blanking time in samples between found peaks.
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
    
    peaks = find(fishnorm>peak_thresh);
    noclippeak = peaks(peaks>1000);
    fishder = [0;fishnorm(2:end)-fishnorm(1:end-1)];
    tops = find(fishder(noclippeak)<0.01 & fishder(noclippeak)>0);
    topdiff = tops(2:end) - tops(1:end-1);
    truetop = find(topdiff>diff_top);
    corrtops = noclippeak(tops(truetop));

%     timedif = peaks(2:end)-peaks(1:end-1);
%     gaps = find(timedif>1);
%     newpeak = [];
%     for g=2:length(gaps)
%         sub = peaks(gaps(g-1):gaps(g)-1);
%         maxpeak = find(fishsoundcal_filt(sub)==max(fishsoundcal_filt(sub)));
%         newpeak(g) = sub(maxpeak);
%     end
%     hold on
%     plot(t_a(newpeak(2:end)),fishsoundcal_filt(newpeak(2:end)),'*')
    %pulsetime = t_a(newpeak(3:end))-t_a(newpeak(2:end-1));
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