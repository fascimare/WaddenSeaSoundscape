%% Script to edit fish sound categories

% ----------Start with input settings:------------

% Specify folder with wav files and sample rate
wavpath = 'F:\WaddenSea fish sounds\Lauwersoog recordings\May1\'; 

%Output path
outpath = 'F:\Sound library\';

%Decide which logged call type you want to make a library for
calltype = 'woodpecker sound';
mkdir(outpath,calltype)

%Spectrogram settings
timebef = 2; %time before call (s)
dur = 5; %spectrogram length (s)
filt_low = 50; %low end of band pass filter (Hz)
filt_high = 1000; %high end of band pass filter (Hz)
nfft = 2400; %FFT size
overlap = 0.9; %window overlap (0-1)
freq_low = 0; %Lower frequency limit on spectrogram (kHz);
freq_high = 1; %Upper frequency limit on spectrogram (kHz);
cal = 177.1; %optional: calibration value in dB (see SoundTrap website for calibration value of hydrophone)
calfile= readtable('G:\My Drive\Sound library\Sound types Wadden Sea\Hydrophone calibration.xlsx'); 
%gain = "SensitivityHighGain";
calgain = calfile(:,1:2); %Adjust column numbers based on which gain you used.
%% Read in log file
[infile,inpath]=uigetfile('*.xlsx','Select a file with manual picks');
if isequal(infile,0)
    disp('Cancelled button pushed');
    return
end
logs = readtable([inpath,infile]);
%
%fs = 24000;
%%  Find appropriate call type and make figure
Index = find(contains(logs.Comments,calltype));
calltypelogs = logs(Index,:);

%Convert start times so they can be used by Matlab
starttime = datetime(calltypelogs.StartTime);
startnum = datenum(starttime);

for n = 1:length(startnum)
    tic
    %Find appropriate wav file and calculate start time in samples
    fullwavpath = [wavpath,calltypelogs.Location{n},'_',calltypelogs.Deployment{n}];

    %Make a list of the wav files in the folder
    wavfiles = dir(fullfile(fullwavpath,'*.wav'));
    filelistwav = struct2cell(wavfiles);
    wavinfo = audioinfo(fullfile(fullwavpath,wavfiles(1).name));
    RegDate = '(?<yr>\d\d)(?<mon>\d\d)(?<day>\d\d)(?<hr>\d\d)(?<min>\d\d)(?<s>\d\d)';
    fileDates = dateregexp(filelistwav(1,:),RegDate);
    wavidx = find(fileDates<= startnum(n) & (fileDates+ datenum([0 0 0 0 0 wavinfo.Duration])) > startnum(n));

% soundpath = calltypelogs.InputFile(n);
% soundpathsections = strsplit(soundpath{:},'\');
% soundfile = soundpathsections{end};
% wavpathfile = [wavpath,'\',soundfile];

    soundfile = wavfiles(wavidx).name;
    wavpathfile = [fullwavpath,'\',wavfiles(wavidx).name];
    fileinfo = audioinfo(wavpathfile);
    fs = fileinfo.SampleRate;

    soundDate = dateregexp(soundfile,RegDate);
    startsample = (startnum(n)-soundDate)*3600*24*fs;
    startS = round(startsample - timebef*fs);
    endS = round(startS+dur*fs);
    if endS>fileinfo.TotalSamples
        timeleft = endS-fileinfo.TotalSamples;
        endS = TotalSamples;
    end
    fishsound = audioread(wavpathfile,[startS,endS]);
    if exist("timeleft")
        %idx = find(strcmp({wavfiles.name},soundfile)==1);
        nextfile = wavfiles(wavidx+1).name;
        nextfilepath = [wavpath,'\',nextfile];
        remainsound = audioread(nextfilepath,[1,timeleft]);
        fishsound = vertcat(fishsound,remainsound);
        clear timeleft
    end

    %Calibration
    soundfilesecs = strsplit(soundfile,'.');
    hyd = str2num(soundfilesecs{1});
    cal = table2array(calgain(calgain.Hydrophone == hyd,2));
    fishsoundcal = fishsound*power(10,cal/20);
    t_a = 0:1/fs:dur;

    %Filter data to make figure look nicer
    [b,a] = ellip(4,0.1,40,[filt_low,filt_high]*2/fs);
    fishsoundcal_filt = filter(b,a,fishsoundcal);

    %Settings for spectrogram
    window=hanning(nfft);
    noverlap = round(nfft*overlap);

    %Make figure
    h=figure;
    spectrogram(fishsoundcal_filt,window,noverlap,nfft,fs,'yaxis')
    clim([60 90]);
    ylim([freq_low freq_high])
    ylabel('Frequency (kHz)')
    fontsize(gca,16, "points")
    title(calltype,'FontSize',24)

    %save spectrogram and soundclip
    timestart = dbSerialDateToISO8601(startnum(n));
    timestr = strrep(timestart, ':', '');
    filename = [outpath,calltype,'\',calltype,'_',calltypelogs.Location{n},'_',calltypelogs.Deployment{n},'_',timestr];
    figname = [filename,'.png'];
    wavname = [filename,'.wav'];
    saveas(gcf,figname)
    audiowrite(wavname,fishsound,fs)
    close(h)
    toc
end
close all

%n = 41;
%callname = logs.Comments(n);
%disp(callname{:});