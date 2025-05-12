%% Script to create nice spectrogram figures

% Specify folder with wav files and sample rate
wavpath = 'F:\WaddenSea fish sounds\Lauwersoog recordings\'; 

%Output path
outpath = 'C:\Users\P310512\Documents\Groningen\Sound library\Fish sounds 2sec clips\';

%Call type (specify which call type name should be given as title of the
%spectrogram
calltype = 'Burst';

%Spectrogram settings
timebef = 0.5; %time before call (s)
dur = 2; %spectrogram length (s)
filt_low = 50; %low end of band pass filter (Hz)
filt_high = 2000; %high end of band pass filter (Hz)
nfft = 5000; %FFT size
overlap = 0.9; %window overlap (0-1)
freq_low = 0; %Lower frequency limit on spectrogram (kHz);
freq_high = 2; %Upper frequency limit on spectrogram (kHz);
climval = [50 90]; %dB values for colorbar

calfile= readtable('G:\My Drive\Sound library\Sound types Wadden Sea\Hydrophone calibration.xlsx'); 

calgain = calfile(:,1:2); %Adjust column numbers based on which gain you used.

%% Read in Triton log file
analysismethod = 'Triton';
[infile,inpath]=uigetfile('*.xlsx','Select a file with manual picks');
if isequal(infile,0)
    disp('Cancelled button pushed');
    return
end
logs = readtable([inpath,infile]);

%Decide which annotation you want to make a spectrogram of
rownr = 5; % put in row number of specific log

%%  Find appropriate call type and make figure

%Convert start times so they can be used by Matlab
starttime = datetime(logs(rownr).StartTime);
startnum = datenum(starttime);

%Make a list of the wav files in the folder
wavfiles = dir(fullfile(wavpath,'*.wav'));
filelistwav = struct2cell(wavfiles);
wavinfo = audioinfo(fullfile(wavpath,wavfiles(1).name));
RegDate = '(?<yr>\d\d)(?<mon>\d\d)(?<day>\d\d)(?<hr>\d\d)(?<min>\d\d)(?<s>\d\d)';
fileDates = dateregexp(filelistwav(1,:),RegDate);
wavidx = find(fileDates<= startnum(n) & (fileDates+ datenum([0 0 0 0 0 wavinfo.Duration])) > startnum(n));

soundfile = wavfiles(wavidx).name;
wavpathfile = [wavpath,'\',wavfiles(wavidx).name];
    fileinfo = audioinfo(wavpathfile);
    fs = fileinfo.SampleRate;

    if strcmp(analysismethod, 'Triton')
        soundDate = dateregexp(soundfile,RegDate);
        startsample = (startnum(n)-soundDate)*3600*24*fs;
        startS = round(startsample - timebef*fs);
        endS = round(startS+dur*fs);
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
    clim(climval);
    ylim([freq_low freq_high])
    ylabel('Frequency (kHz)')
    fontsize(gca,16, "points")
    title(calltype,'FontSize',24)

    %save spectrogram and soundclip
    timestart = dbSerialDateToISO8601(startnum(n));
    timestr = strrep(timestart, ':', '');
    if strcmp(analysismethod, 'Raven')
        filename = [outpath,'/',calltype,'/',calltype,'_',calltypelogs.Location{n},'_',calltypelogs.Deployment{n},'_',timestr];
    else
    filename = [outpath,calltype,'\',calltype,'_',hyd,'_',timestr]; %Mac computers might need forward slash here.
    end
    figname = [filename,'.png'];
    wavname = [filename,'.wav'];
    saveas(gcf,figname)
    audiowrite(wavname,fishsound,fs)
    close(h)
    toc
end
