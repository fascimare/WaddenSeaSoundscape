%% Script to edit fish sound categories

% ----------Start with input settings:------------
clear all
% Specify folder with wav files and sample rate
wavpath = 'F:\WaddenSea fish sounds\Lauwersoog recordings\June\'; 

%Output path
%outpath = 'F:\Sound library\';


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
%%  Change file path

new_inputfile = {};
%Convert start times so they can be used by Matlab
starttime = datetime(logs.StartTime);
startnum = datenum(starttime);

for n = 1:length(startnum)
    tic
    %Find appropriate wav file and calculate start time in samples
    if strcmp(analysismethod, 'Raven')
        fullwavpath = wavpath;
    else
    fullwavpath = [wavpath,logs.Location{n},'_',logs.Deployment{n}];
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
    new_inputfile{n} = wavpathfile;
    

    toc
end

logs.new_inputfile = new_inputfile';
%%
filename = split(infile,'.xlsx');
csvname = [inpath,filename{1},'_inputmod.xlsx'];
writetable(logs,csvname)
%n = 41;
%callname = logs.Comments(n);
%disp(callname{:});