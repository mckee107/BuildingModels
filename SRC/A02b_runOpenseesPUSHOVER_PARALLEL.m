% Codi McKee
% Texas A&M University
% First Created: 05-Feb-2019
% Last Modified: 30-May-2019

% Run OpenSees model for selected frame from database

%**************************************************************************
% Initialize
%


close all; fclose all; clc; diary off; format short g

runTotaltemp = tic;

% Intitial Workspace
%   Set Parallel Model Parameters in A01_initializeWorkspace.m file

A01_initializeWorkspace_2; % Calls external .m file

% Model Toggles

tog.runOpenSees = 1;
tog.runEigen    = 1;
tog.runGravity  = 1;
tog.runPushover = 1;
tog.runDynamic  = 0;
tog.Scaled      = 1;
tog.Display     = 1;

tog.buildInput     = 0;
tog.buildRecorders = 0;
tog.buildPushover  = 0;
tog.buildDynamic   = 0;

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

for frameNum = frameNumRange
    
    for t = tRange
        frameNumval(countVal) = frameNum;
        time(countVal) = t;
        countVal = countVal + 1;
    end
    
end

fprintf('\n Creating Run Files...\n');

for countVal = 1:length(time)
     
    frameNum = frameNumval(countVal);
    frameID = Frames(frameNum).ID;    
    t = time(countVal);
    
    %**************************************************************************
    % Create beamName.tcl File
    %
    
    filename = 'frameName.tcl';
    framefilename = [pwd '\ModelFiles\runFile_' num2str(countVal) '_P.tcl'];
    
    frameName{countVal} = [ num2str(t) 'yr_Corr' num2str(corrosion_level) '_' num2str(frameID) ];
    
%     createFramefile(framefilename, tog, corrosion_level, t, frameID, countVal)
    
    buildFrameinputPARALLEL
    
%     buildRecordersPARALLEL2
    
    buildPushoverPARALLEL
    
    
end

fprintf('\n Run Files COMPLETE! \n');

% parpool('local',4);

cd(oldDir);

loopRange = [1];

if exist('tracker_pushover.csv', 'file') == 2
    delete('tracker_pushover.csv');
end

% parfor jj = 1:length(time)
for jj = 1:length(loopRange)
    
    fprintf('\n Starting OpenSees Analysis...\n');
    
    runTag = ['!OpenSees.exe RunFile_' num2str(loopRange(jj)) '_P.tcl'];
    
    cd(newDir);
    doEval(runTag)
    fprintf('\n OpenSees Analysis COMPLETE !!\n');
    cd(oldDir);
    
    filename = [pwd '\tracker_pushover.csv'];
    fid = fopen(filename,'a');
    fprintf(fid,'%s\n',num2str(loopRange(jj)));
    fclose(fid);
 
    cd(newDir);
end


cd(oldDir);
saveTimetemp = tic;
fprintf('\n\n Saving Frames.mat...\n');
save('FramesPushover.mat','Frames','-v7.3')
timeSave = toc(saveTimetemp);
fprintf('\n Saving COMPLETE!!\n');

if timeSave > 60
    timeMin = timeSave/60;
    timeSec = mod(timeSave,60);
    fprintf('\n\t Save Time: %0.0f Min. %0.1f Sec.\n\n',timeMin,timeSec);
else
    fprintf('\n\t Save Time: %0.1f Sec.\n\n',timeSave);
end
cd(oldDir);

%**************************************************************************
% Internal Functions
%

function [] = createFramefile(framefilename, tog, corrosion_level, t, frameID, countVal)
% Create frameFile
if tog.runPushover
    pushStr = '_P';
elseif tog.runDynamic
    pushStr = '';
end

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
fprintf(fid,['\n set recorderFile [join [list "Recorders_" $frameName "' pushStr '.tcl"] ""];']);
fprintf(fid,['\n set togPush ' num2str(tog.runPushover) ';']);
fprintf(fid,['\n set togGravity ' num2str(tog.runGravity) ';']);
fprintf(fid,['\n set togDynamic ' num2str(tog.runDynamic) ';']);
fprintf(fid,['\n set togEigen ' num2str(tog.runEigen) ';']);
fprintf(fid,['\n set togDisplay ' num2str(tog.Display) ';']);
fprintf(fid,'\n source model.tcl');
fclose(fid);

end

