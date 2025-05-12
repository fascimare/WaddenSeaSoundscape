%% Script to edit fish sound categories

% ----------Start with input settings:------------

% Specify folder with wav files and sample rate
wavpath = 'E:\WaddenSea fish sounds\Lauwersoog recordings\LAUW_06\LAUW1_ON_06\'; 

%Output path
outpath = 'C:\Users\P310512\Documents\Groningen\Sound library\MatchingSounds\';

%Decide which logged call type you want to make a library for
calltype = 'Drum';
mkdir(outpath,calltype)

%Spectrogram settings
timebef = 0.5; %time before call (s)
dur = 5; %spectrogram length (s)
filt_low = 50; %low end of band pass filter (Hz)
filt_high = 2000; %high end of band pass filter (Hz)
nfft = 3500; %FFT size
overlap = 0.9; %window overlap (0-1)
freq_low = 0; %Lower frequency limit on spectrogram (kHz);
freq_high = 2; %Upper frequency limit on spectrogram (kHz);
climval = [50 90]; %dB values for colorbar
%cal = 177.1; %optional: calibration value in dB (see SoundTrap website for calibration value of hydrophone)
calfile= readtable('G:\My Drive\Sound library\Sound types Wadden Sea\Hydrophone calibration.xlsx'); 
%gain = "SensitivityHighGain";
calgain = calfile(:,1:2); %Adjust column numbers based on which gain you used.
%% Read in Triton log file
analysismethod = 'Triton';
[infile,inpath]=uigetfile('*.xls','Select a file with manual picks');
if isequal(infile,0)
    disp('Cancelled button pushed');
    return
end
logs = readtable([inpath,infile]);
%
%fs = 24000;
logsall = logs;
siteidx = find(contains(logsall.Location,"LAUW1_ON"));
logs = logsall(siteidx,:);
calltypelist = unique(logs.Comments);
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
Rlogs.Comments = Rlogs.soundcode; %sound_code for other file, consistency needed.
depinfo = split(infile,'.');
Rlogs.Location = repmat(string(depinfo{2}),length(Rlogs.StartTime),1);
Rlogs.Deployment = repmat(string(depinfo{1}),length(Rlogs.StartTime),1);
Rlogssub = Rlogs(strcmp(Rlogs.View,'Spectrogram 1'),:);
logs = Rlogssub;
calltypelist = unique(logs.soundcode);
%%  Find appropriate call type and make figure
%Loop through call types, making sure to skip "noise period" types.
for ct = 6:7
    calltype = calltypelist{ct};
    Index = find(contains(logs.Comments,calltype));
    calltypelogs = logs(Index,:);

    %Convert start times so they can be used by Matlab
    starttime = datetime(calltypelogs.StartTime);
    startnum = datenum(starttime);

    for n = 1:length(startnum)
        tic
        %Find appropriate wav file and calculate start time in samples
        %if strcmp(analysismethod, 'Raven')
            fullwavpath = wavpath;
        %else
        %fullwavpath = [wavpath,calltypelogs.Location{n},'_',calltypelogs.Deployment{n}];
        %end
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
        if contains(calltype,"?")
            calltype = calltype(1:end-1);
        elseif contains(calltype,"/")
            calltype = char(strrep(calltype,"/","-"));
        end
        %if strcmp(analysismethod, 'Raven')
            %filename = [outpath,'/',calltype,'/',calltype,'_',calltypelogs.Location{n},'_',calltypelogs.Deployment{n},'_',timestr];
            %filename = [outpath,calltypelogs.Location{n},'_',calltypelogs.Deployment{n},'_',timestr,'_',calltype];
        filename = [outpath,'LAUW1_ON_06_',timestr,'_',calltype];
            %else
        %filename = [outpath,calltype,'\',calltype,'_',calltypelogs.Location{n},'_',calltypelogs.Deployment{n},'_',timestr];
        %end
        figname = [filename,'.png'];
        wavname = [filename,'.wav'];
        saveas(gcf,figname)
        audiowrite(wavname,fishsound,fs)
        close(h)
        toc
    end
end

%n = 41;
%callname = logs.Comments(n);
%disp(callname{:});