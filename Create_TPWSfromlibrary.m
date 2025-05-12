%% Script to edit fish sound categories

% ----------Start with input settings:------------

% Specify folder with wav files and sample rate
wavpath = 'F:\WaddenSea fish sounds\Lauwersoog recordings\May2\'; 

%Output path
outpath = 'C:\Users\P310512\Documents\Groningen\Sound library\Fish sound clustering\TPWS library\';

% %Decide which logged call type you want to make a library for
%calltype = 'very low frequency drumming';
%mkdir(outpath,calltype)

%Spectrogram settings
timebef = 0.3; %time before call (s)
dur = 3.5; %spectrogram length (s)
filt_low = 75; %low end of band pass filter (Hz)
filt_high = 2000; %high end of band pass filter (Hz)
nwindow = 3500; %window size
nfft = 3500; %FFT size
overlap = 0.9; %window overlap (0-1)
freq_low = 0; %Lower frequency limit on spectrogram (kHz);
freq_high = 2; %Upper frequency limit on spectrogram (kHz);
climval = [50 80]; %dB values for colorbar
%cal = 177.1; %optional: calibration value in dB (see SoundTrap website for calibration value of hydrophone)
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
%
%fs = 24000;
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
%Convert start times so they can be used by Matlab
 Index = find(contains(logs.Comments,calltype));
 calltypelogs = logs(Index,:);
 
starttime = datetime(calltypelogs.StartTime);
startnum = datenum(starttime);
endtime = datetime(calltypelogs.EndTime);
endnum = datenum(endtime);
calltypelogs.dur = (endnum-startnum)*24*3600;


%% Trying something new
MPP = [];
MSN = [];
MSP = [];
MTT = [];
%%
%Read in all sound clips from one sound type

wavfiles = dir(fullfile("G:\Shared drives\Wadden Sea Sound Library\Sound library\solid fish grunt",'*.wav'));
wavinfo = audioinfo(fullfile(wavfiles(1).folder,wavfiles(1).name));
fs = wavinfo.SampleRate;
cal = 177;

%Create empty output variables
mpp = ones(length(wavfiles),1); %Peak-to-peak sound level, not used
msn = zeros(length(wavfiles),24001);
msp = zeros(length(wavfiles),282);
mtt = ones(length(wavfiles),1);

for n = 1:length(wavfiles)
    tic
    fishsound = audioread(fullfile(wavfiles(n).folder,wavfiles(n).name),[2*fs,3*fs]);
    fishsoundcal = fishsound*power(10,cal/20);
     [b,a] = ellip(4,0.1,40,[filt_low,filt_high]*2/fs);
    fishsound_filt = filter(b,a,fishsoundcal);
    [fishpsd,F] = pwelch(fishsound_filt,nwindow,overlap,nfft,fs);

    %Now we start filling the output
    msn(n,:) = fishsound_filt';
    msp(n,:) = fishpsd(12:293)';
end

MPP = vertcat(MPP,mpp);
MSN = vertcat(MSN,msn);
MSP = vertcat(MSP,msp);
MTT = vertcat(MTT,mtt);
%% save output
f = round(F(12:293)');
save("C:\Users\P310512\Documents\Groningen\Sound library\Fish sound clustering\LAUW_fishgrunts_TPWS1.mat","MTT","MSP","MSN","MPP","f");

%%  
%--------- New version of TPWS creation, includes centering of wave form-------------

%First, create output variables
MPP = [];
MSN = [];
MSP = [];
MTT = [];
allcallID = {};
%%
calltypes = unique(logs.Comments);

for type = 4:length(calltypes)
    calltype = calltypes{type};
Index = find(contains(logs.Comments,calltype));
calltypelogs = logs(Index,:);

%Convert start times so they can be used by Matlab
starttime = datetime(calltypelogs.StartTime);
startnum = datenum(starttime);
endtime = datetime(calltypelogs.EndTime);
endnum = datenum(endtime);
calltypelogs.duration = (endnum-startnum)*24*3600;

% Create empty output variables
mpp = ones(length(calltypelogs.StartTime),1); %Peak-to-peak sound level, not used
msn = zeros(length(calltypelogs.StartTime),dur*24000+1); %waveform
msp = zeros(length(calltypelogs.StartTime),(nfft/2+1)); %spectrum
mtt = zeros(length(calltypelogs.StartTime),1);
callID = cell(length(calltypelogs.StartTime),1);

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

% soundpath = calltypelogs.InputFile(n);
% soundpathsections = strsplit(soundpath{:},'\');
% soundfile = soundpathsections{end};
% wavpathfile = [wavpath,'\',soundfile];

    soundfile = wavfiles(wavidx).name;
    wavpathfile = [fullwavpath,'\',wavfiles(wavidx).name];
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
    t_a = 0:1/fs:dur;

    %Filter data to make figure look nicer
    [b,a] = ellip(4,0.1,40,[filt_low,filt_high]*2/fs);
    fishsoundcal_filt = filter(b,a,fishsoundcal);
   
%     %Find exact duration of call
%     esound = cumsum(fishsoundcal_filt.^2);
%     emsound = esound/max(esound);
%     ricosound = emsound(2:end)-emsound(1:end-1);
%     ricosoundsm = movmean(ricosound,2400);
%     a1 = find(ricosoundsm(1:20000)<1e-5,1);
%     a2 = find(ricosoundsm(20000:end)<1e-5,1);
%     
%     %Pad with zeros so duration stays the same
%     if isempty(a1)
%         a1 = 1;
%     end
%     fishsoundcal_filtsub = fishsoundcal_filt(a1:(20000+a2));
%     edgedur = (dur*fs - length(fishsoundcal_filtsub))/2;
%     fishsoundcal_filtpad = zeros(length(fishsoundcal),1);
%     fishsoundcal_filtpad(edgedur:(edgedur-1+length(fishsoundcal_filtsub))) = fishsoundcal_filtsub;

    %Save to TPWS file. Make sure to include original category in a
    %separate variable.
    msn(n,:) = fishsoundcal_filt;
    [fishpsd,F] = pwelch(fishsoundcal_filt,nwindow,overlap,nfft,fs);
    fishpsd_dB = 10*log10(fishpsd);
    msp(n,:) = fishpsd_dB';
    mtt(n) = startnum(n);
    callID{n} = calltype;

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
    filename = [outpath,calltype,'_',calltypelogs.Location{n},'_',calltypelogs.Deployment{n},'_',timestr];
    end
    figname = [filename,'.png'];
    wavname = [filename,'.wav'];
    saveas(gcf,figname)
    audiowrite(wavname,fishsound,fs)
    close(h)
    toc
end

allcallID = vertcat(allcallID,callID);
MPP = vertcat(MPP,mpp);
MSN = vertcat(MSN,msn);
MSP = vertcat(MSP,msp);
MTT = vertcat(MTT,mtt);
f = round(F)';
end
%% Save TPWS file per deployment
TPWSfile = [outpath,'LAUW_',calltypelogs.Deployment{1},'_alltypes_TPWS1.mat'];
save(TPWSfile,"MTT","MSP","MSN","MPP","f","allcallID");

%% Correcting PSD to dB values instead of Pa

MSP = 10*log10(MSP);
TPWSfile = 'F:\TPWS library fish sounds\Padded\TPWS\LAUW_07_fishgrunts_dB_TPWS1.mat';
save(TPWSfile,"MTT","MSP","MSN","MPP","f","allcallID");

%% only saving relevant frequencies

f = f(12:293);
MSP = MSP(:,12:293);
TPWSfile = 'C:\Users\P310512\Documents\Groningen\Sound library\Fish sound clustering\TPWS library\LAUW_07_alltypes_dB_75-2002Hz_TPWS1.mat';
save(TPWSfile,"MTT","MSP","MSN","MPP","f","allcallID");
