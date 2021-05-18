% Codi McKee
% Texas A&M University
% First Created: 05-Feb-2019
% Last Modified: 30-May-2019

% Run OpenSees model for selected frame from database

%**************************************************************************
% Initialize
%
close all; fclose all; clc; diary off; format short g; 

runTotaltemp = tic;

% Intitial Workspace
%   Set Parallel Model Parameters in A01_initializeWorkspace.m file

A01_initializeWorkspace_2; % Calls external .m file

% Check loopRange variable
% [outTag ] = check_loopRange; % Make sure that you didn't out something special for the loopRange previously

% switch outTag
%     case 'return'
%         return
% end


% Model Toggles

tog.runOpenSees = 1;
tog.runEigen    = 1;
tog.runGravity  = 1;
tog.runPushover = 0;
tog.runDynamic  = 1;
tog.Scaled      = 1;
tog.Display     = 0;
tog.runParallel = 1;
tog.stpNum      = 4;

tog.buildRunFiles  = 1;
tog.buildInput     = 1;
tog.buildRecorders = 1;
tog.buildPushover  = 1;
tog.buildDynamic   = 1;

tog.importPushoverfiles = 0;
tog.importDynamicfiles  = 0;

tog.plotResults     = 0;
tog.NotifyTxt       = 0;
tog.NotifyEmail     = 1;

tog.saveVal         = 0;
tog.saveVal_end     = 0;

tog.plotPushover    = 0;
tog.plotDynamic     = 0;


%**************************************************************************
% Run Model
%
countVal = 1;

for scaleFac = scaleRange
    
    if scaleFac == 1
        scaledStr = '_DE';
        
    elseif scaleFac == 2
        scaledStr = '_MCE';
        
    end
    
    for frameNum = frameNumRange
        for t = tRange
            for groundMotion_num = numRange
                frameNumval(countVal) = frameNum;
                time(countVal) = t;
                GM(countVal) = groundMotion_num;
                scaleVal{countVal} = scaledStr;
                NPTS(countVal) = GroundMotions(groundMotion_num).NPTS * tog.stpNum ;
                NPtsTot(countVal) = GroundMotions(groundMotion_num).TotalTime * tog.stpNum ;
                countVal = countVal + 1;
            end
        end
    end
    
end

totPts = sum(NPtsTot);
avgPts = mean(NPtsTot);

if tog.buildRunFiles; fprintf('\n\tCreating Run Files...\n');end

for countVal = 1:length(time)
    
    frameNum = frameNumval(countVal);
    frameID = Frames(frameNum).ID;
    t = time(countVal);
    groundMotion_num = GM(countVal);
    scaledStr = scaleVal{countVal};
    
    %**************************************************************************
    % Create beamName.tcl File
    %
    
    
    filename = 'frameName.tcl';
    framefilename = [pwd '\ModelFiles\runFile_' num2str(countVal) '_D.tcl'];
    
    frameName{countVal} = [ num2str(t) 'yr_Corr' num2str(corrosion_level) '_' num2str(frameID) ];
    
    if tog.buildRunFiles
%         createFramefile(framefilename, tog, corrosion_level, t, frameID, countVal,GM(countVal),scaledStr,GroundMotions)
        buildFrameinputPARALLEL
        
%         buildDynamicPARALLEL; %save GroundMotions.mat GroundMotions
        buildRecordersPARALLEL2
        addTime_val(countVal) = GroundMotions(GM(countVal)).AddTime;
        %         countVal
    end
          
end

if tog.buildRunFiles; fprintf('\n\tRun Files COMPLETE! \n'); end

% loopRange = [1];%24 27 30 32 34 36 37 40 41];
loopRange = 1:length(time);
% return

oldDir = cd; cd([pwd '\ModelFiles']); newDir = cd;


if tog.runParallel
    parpool('local',6);
end

if ~isempty(gcp('nocreate'))
    currentPool = gcp;
    numWorks = currentPool.NumWorkers;
else
    numWorks = 1;
end

cd(oldDir);


if exist('tracker.csv', 'file') == 2
    delete('tracker.csv');
end

if exist('timeTracker.csv', 'file') == 2
    delete('timeTracker.csv');
end

timeTic = NaN(length(loopRange),1);
timeToc = NaN(length(loopRange),1);
itrTot  = length(loopRange);
dispNum = 1;

for jj = 1:length(loopRange)
    
    if tog.runParallel == 0 && ~isempty(gcp('nocreate'))
        error('Not Supposed to be running Parallel!!')
    end
    
    timeTic = tic;
    
    scaledStr = scaleVal{loopRange(jj)};
    
    fprintf('\n\tStarting OpenSees Analysis...\n');
    
    %     runTag = ['!openSees_Apr_19.exe RunFile_' num2str(loopRange(jj)) '_D.tcl'];
    runTag = ['!OpenSees.exe RunFile_' num2str(loopRange(jj)) '_D.tcl'];
    
    cd(newDir);
    runTimetic(jj) = tic;
    
    doEval(runTag)
    pause(2)
    fprintf('\n\tOpenSees Analysis COMPLETE !!\n');
    runTimetoc(jj) = toc(runTimetic(jj));
    
    cd(oldDir);
    dirOut = sprintf('')
    [status,cmdout] = system(['tar -xzf .\0yr_Corr2_1001\Dynamic\AT_1_DE.tar.gz .\0yr_Corr2_1001\Dynamic\AT_1_DE']);
       
    % Progress Update
    
    timeToc = toc(timeTic);
        
    try
                              
        filename = [pwd '\timeTracker.csv'];
        fid = fopen(filename,'a');
        fprintf(fid,'%s\n%s\n',num2str(round(timeToc,0)),num2str(NPtsTot(jj)));
        fclose(fid);
        
        filename = [pwd '\timeTracker.csv'];
        fid = fopen(filename,'r');
        timeTrackerout = textscan(fid, '%s', 'Delimiter', '\n');
        fclose(fid);
        
        timeTrackerout = str2double([timeTrackerout{:}]);
        trackerTime   = timeTrackerout(1:2:end);
        trackerNPTS   = timeTrackerout(2:2:end);
        
        trackerWeight = trackerNPTS./avgPts;
        
        weightTime = trackerTime.*trackerWeight;
        
        avgTime = mean(weightTime);
        
        itrNum = length(trackerTime);
           
        remainTime = (itrTot - itrNum) * avgTime/numWorks;
        
        numRem = rem(itrNum,dispNum);
        
        progressTime(remainTime,numRem,itrNum,itrTot);
        catchVal = 0;
    catch        
        catchVal = 1;
        fprintf('\n\n\t** Unable to read/write timeTracker.csv File!!!! **\n\n')
    end
    
    % Update Tracker File
    
    try
        numPts = GroundMotions(GM(jj)).NPTS;
        filename = [pwd '\tracker.csv'];
        fid = fopen(filename,'a');
        if catchVal == 0
            fprintf(fid,'%s %s %s\n',num2str(loopRange(jj)),num2str(round(runTimetoc(jj),0)),num2str(round(remainTime,0)));
        else
            fprintf(fid,'%s %s\n',num2str(loopRange(jj)),num2str(round(runTimetoc(jj),0)));
        end
        fclose(fid);
    catch
        fprintf('\n\n\t** Unable to read/write Tracker.csv File!!!! **\n\n')
    end
    
end

for kk = 1:length(time)
    
    scaledStr = scaleVal{kk};
    DT = GroundMotions(GM(kk)).DT / tog.stpNum;
    NPTS = GroundMotions(GM(kk)).NPTS * tog.stpNum;
    addTime = addTime_val(kk) * tog.stpNum;
    
    timeFile = [pwd '\Output\' num2str(t) 'yr_Corr' num2str(corrosion_level) '_' num2str(frameID) '\Dynamic\AT_' num2str(groundMotion_num) scaledStr '\TimeStep.txt'];
    
    if exist(timeFile, 'file') == 2
        
        fid = fopen(timeFile,'r');
        timeFileout = textscan(fid, '%s', 'Delimiter', ',');
        fclose(fid);
        timeFileout = [timeFileout{:}];
        currentTime(kk) = str2double(cell2mat(timeFileout(1)));
        totalTime(kk) = str2double(cell2mat(timeFileout(2)));
        
        Frames(frameNum).Output.(['yr_' num2str(time(kk))]).Dynamic.(['AT_' num2str(GM(kk)) scaleVal{kk}]).FinalTime = currentTime(kk);
        Frames(frameNum).Output.(['yr_' num2str(time(kk))]).Dynamic.(['AT_' num2str(GM(kk)) scaleVal{kk}]).TotalTime = totalTime(kk);
    end
    
end


cd(oldDir);
saveTimetemp = tic;
fprintf('\n\n\tSaving Frames.mat...\n');
% save('Frames.mat','Frames','-v7.3')
timeSave = toc(saveTimetemp);
fprintf('\n\tSaving COMPLETE!!\n');

if timeSave > 60
    timeMin = timeSave/60;
    timeSec = mod(timeSave,60);
    fprintf('\n\t Save Time: %0.0f Min. %0.1f Sec.\n\n',timeMin,timeSec);
else
    fprintf('\n\t Save Time: %0.1f Sec.\n\n',timeSave);
end


%**************************************************************************
% Internal Functions
%

function [outTag] = check_loopRange

% answer = questdlg('Check that loopRange is correct!','Check loopRange','Correct','Incorrect','Correct');
answer = questdlg_timer(10,'Check that loopRange is correct!','Check loopRange','Correct','Incorrect','Correct');

switch answer
    case 'Correct'
        fprintf('\n\tloopRange CORRECT... continue\n\n');
        outTag = [''];
        
    case 'Incorrect'
        fprintf('\n\tChange loopRange... Canceling Analysis\n\n');
        outTag = ['return'];
end

end

function [] = createFramefile(framefilename, tog, corrosion_level, t, frameID, countVal, groundMotion_num,scaledStr, GroundMotions)
% Create frameFile

fid = fopen(framefilename,'w');

mainDir = ['{' pwd '}'];

fprintf(fid,'# Codi McKee');
fprintf(fid,'\n# Texas A&M University');
fprintf(fid,'\n# %s',char(datetime('today')));
fprintf(fid,'\n\n# This file was written via Matlab %s',version);
fprintf(fid,'\n# File: %s\\%s.m\n\n',pwd,mfilename);


fprintf(fid,'\n set N %s',num2str(countVal));
fprintf(fid,'\n set frameName "%s"',[ num2str(t) 'yr_Corr' num2str(corrosion_level) '_' num2str(frameID) ]);
fprintf(fid,'\n set mainDir %s ;',mainDir);
fprintf(fid,'\n cd $mainDir');

fprintf(fid,'\n set frameFile [join [list "Input_" $frameName ".tcl"] ""];');
fprintf(fid,'\n set pushoverFile [join [list "Pushover_" $frameName ".tcl"] ""];');
% fprintf(fid,['\n set dynamicFile [join [list "Dynamic_" $frameName "_AT_' num2str(groundMotion_num) scaledStr '.tcl" ] ""];']);
fprintf(fid,'\n set recorderFile [join [list "Recorders_" $frameName ".tcl"] ""];');
fprintf(fid,['\n set togPush ' num2str(tog.runPushover) ';']);
fprintf(fid,['\n set togGravity ' num2str(tog.runGravity) ';']);
fprintf(fid,['\n set togDynamic ' num2str(tog.runDynamic) ';']);
fprintf(fid,['\n set togEigen ' num2str(tog.runEigen) ';']);
fprintf(fid,['\n set togDisplay ' num2str(tog.Display) ';']);
fprintf(fid,['\n set togDisplaymode ' num2str(tog.Display) ';']);

fprintf(fid,'\n\n# Dynamic Loads: Frame ID %s',num2str(frameID));

%**************************************************************************
% Ground Motion File
%
DT = GroundMotions(groundMotion_num).DT ;
NPTS = GroundMotions(groundMotion_num).NPTS;

fprintf(fid,['\n set GroundMotionFile AT_' num2str(groundMotion_num) scaledStr '.tcl;']);
fprintf(fid,['\n set GroundMotion AT_' num2str(groundMotion_num) scaledStr]);

fprintf(fid,['\n set DT ' num2str(DT) ';']);
fprintf(fid,['\n set NPTS ' num2str(NPTS) ';']);
fprintf(fid,['\n set stpNum ' num2str(tog.stpNum) ';']);
addTime = ceil(NPTS/10);
GroundMotions(groundMotion_num).AddTime = addTime;
GroundMotions(groundMotion_num).TotalTime = NPTS + addTime;

fprintf(fid,'\n set addTime %0.0f;', addTime);

fprintf(fid,['\n cd "OpenSees Files" ;']);
fprintf(fid,'\n source ModelDB.tcl');
fclose(fid);

end

function progressTime(timeValue,numRem,itrNum,itrTot)
if numRem == 0 || itrNum == 1
    
    if timeValue > 3600
        timeHr  = timeValue/3600;
        timeMin = mod(timeValue,3600)/60;
        timeSec = mod(timeMin,1)*60;
        
        if timeMin > 1
            fprintf('\n\t\t%1.0f\tof\t%1.0f:\t%1.2f%%\tApprox. %1.0f Hr. %1.0f Min. %1.0f Sec. Remaining...\n',itrNum,itrTot,itrNum*100/itrTot,timeHr,timeMin,timeSec);
        else
            fprintf('\n\t\t%1.0f\tof\t%1.0f:\t%1.2f%%\tApprox. %1.0f Hr. %1.0f Min. %1.0f Sec. Remaining...\n',itrNum,itrTot,itrNum*100/itrTot,timeHr,0,timeSec);
        end
    elseif timeValue > 60
        timeMin = timeValue/60;
        timeSec = mod(timeValue,60);
        fprintf('\n\t\t%1.0f\tof\t%1.0f:\t%1.2f%%\tApprox. %1.0f Min. %1.0f Sec. Remaining...\n',itrNum,itrTot,itrNum*100/itrTot,timeMin,timeSec);
    else
        fprintf('\n\t\t%1.0f\tof\t%1.0f:\t%1.2f%%\tApprox. %1.0f Sec. Remaining...\n',itrNum,itrTot,itrNum*100/itrTot,timeValue);
    end
end

end

