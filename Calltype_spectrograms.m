%% Script to edit fish sound categories

% ----------Start with input settings:------------

% Specify folder with wav files and sample rate
wavpath = 'F:\WaddenSea fish sounds\Lauwersoog recordings\May2\LAUW_OFFREEF_E_6049_May2'; 

%Decide which logged call type you want to plot
calltype = 'striped fish grunt';

%Decide the name and path of the audio cilp that you save
wavname = '';
%Pick one of the calls in the subset
n = 2;

%Spectrogram settings
timebef = 2.5; %time before call (s)
dur = 5; %spectrogram length (s)
filt_low = 50; %low end of band pass filter (Hz)
filt_high = 1000; %high end of band pass filter (Hz)
nfft = 2400; %FFT size
overlap = 0.9; %window overlap (0-1)
freq_low = 0; %Lower frequency limit on spectrogram (kHz);
freq_high = 1; %Upper frequency limit on spectrogram (kHz);
cal = 177.1; %optional: calibration value in dB (see SoundTrap website for calibration value of hydrophone)
%% Read in log file
[infile,inpath]=uigetfile('*.xls','Select a xls file with manual picks');
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
%Find appropriate wav file and calculate start time in samples
soundpath = calltypelogs.InputFile(n);
soundpathsections = strsplit(soundpath{:},'\');
soundfile = soundpathsections{end};
wavpathfile = [wavpath,'\',soundfile];
fileinfo = audioinfo(wavpathfile);
fs = fileinfo.SampleRate;
RegDate = '(?<yr>\d\d)(?<mon>\d\d)(?<day>\d\d)(?<hr>\d\d)(?<min>\d\d)(?<s>\d\d)';
fileDates = dateregexp(soundfile,RegDate);
startsample = (startnum(n)-fileDates)*3600*24*fs;
startS = round(startsample - timebef*fs);
endS = round(startS+dur*fs);
fishsound = audioread(wavpathfile,[startS,endS]);

%Optional: calibration
if ~isempty(cal)
    fishsoundcal = fishsound*power(10,cal/20);
    t_a = 0:1/fs:dur;
    figure;plot(t_a,fishsoundcal)
else
    fishsoundcal = fishsound;
end

%Filter data to make figure look nicer
[b,a] = ellip(4,0.1,40,[filt_low,filt_high]*2/fs);
fishsoundcal_filt = filter(b,a,fishsoundcal);

%Settings for spectrogram
window=hanning(nfft);
noverlap = round(nfft*overlap);

%Make figure
figure;
spectrogram(fishsoundcal_filt,window,noverlap,nfft,fs,'yaxis')
clim([60 100]);
ylim([freq_low freq_high])
ylabel('Frequency (kHz)')
fontsize(gca,16, "points")
title(calltype,'FontSize',24)

%To save audio file
[b,a] = ellip(4,0.1,40,[filt_low,filt_high]*2/fs);
fishsound_filt = filter(b,a,fishsound);
audiowrite(wavname,fishsound_filt,fs)

%n = 41;
%callname = logs.Comments(n);
%disp(callname{:});