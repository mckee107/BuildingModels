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


% Model Toggles

tog.runOpenSees = 1;
tog.runEigen    = 0;
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

% SET WALL TIME

tog.walltime = 5; % Hours
wtH = floor(tog.walltime); wtM = floor((tog.walltime - wtH)*60);

% countVal = 1;

% DYNAMIC

% Sort order of ground motions
for GM = numRange
    NPTS(GM) = GroundMotions(GM).NPTS;
    
end

[NPTS, I] = sort(NPTS);
% numRange = numRange(I);

if tog.runDynamic
    countVal = 1;
    for scaleFac = scaleRange
        
        if scaleFac == 1
            scaledStr = '_DE';
            
        elseif scaleFac == 2
            scaledStr = '_MCE';
            
        end
        
        for frameNum = frameNumRange
            
            for groundMotion_num = numRange
                frameNumval(countVal) = frameNum;
                %                     time(countVal) = t;
                GM(countVal) = groundMotion_num;
                scaleVal{countVal} = scaledStr;
                NPTS(countVal) = GroundMotions(groundMotion_num).NPTS * tog.stpNum ;
                NPtsTot(countVal) = GroundMotions(groundMotion_num).TotalTime * tog.stpNum ;
                countVal = countVal + 1;
            end
            
        end
        
    end
    
    
    if tog.buildRunFiles; fprintf('\n\tCreating Dynamic Run Files...\n');end
    for t = tRange
        for countVal = 1:length(GM)
            
            frameNum = frameNumval(countVal);
            frameID = Frames(frameNum).ID;
            %         t = time(countVal);
            groundMotion_num = GM(countVal);
            scaledStr = scaleVal{countVal};
            
            %**************************************************************************
            % Create beamName.tcl File
            %
            framefilename = [pwd '\RunFiles\RF_' num2str(t), 'yr_Corr' num2str(corrosion_level) '_' num2str(frameID) '_AT_' num2str(GM(countVal)) scaledStr nameSuffix '.tcl'];
            if tog.buildRunFiles
%                 createFramefile(framefilename, tog, corrosion_level,t, frameID, countVal,GM(countVal),scaledStr,GroundMotions, nameSuffix)
            end
        end
    end
    for frameNum = frameNumRange
        frameID = Frames(frameNum).ID;
        frameName = [ 'Corr' num2str(corrosion_level) '_' num2str(frameID) ];
        
        if tog.buildRunFiles
            %                 createFramefile(framefilename, tog, corrosion_level, t, frameID, countVal,GM(countVal),scaledStr,GroundMotions)
            buildFrameinputPARALLEL
            buildParameterUpdate
            %         buildDynamicPARALLEL; %save GroundMotions.mat GroundMotions
            buildRecordersPARALLEL2
            addTime_val(countVal) = GroundMotions(GM(countVal)).AddTime;
            %         countVal
        end
    end
end




% PUSHOVER

if tog.runPushover
    countVal = 1;
    
    for frameNum = frameNumRange
        
        for t = tRange
            frameNumval(countVal) = frameNum;
            time(countVal) = t;
            countVal = countVal + 1;
        end
        
    end
    
    fprintf('\n Creating Pushover Run Files...\n');
    
    for countVal = 1:length(time)
        
        frameNum = frameNumval(countVal);
        frameID = Frames(frameNum).ID;
        t = time(countVal);
        
        %**************************************************************************
        % Create beamName.tcl File
        %
        
        filename = 'frameName.tcl';
        framefilename = [pwd '\RunFiles\Pushover\RF_' num2str(t) 'yr_Corr' num2str(corrosion_level) '_' num2str(frameID) '_Push' nameSuffix '.tcl'];
        
        frameName{countVal} = [ num2str(t) 'yr_Corr' num2str(corrosion_level) '_' num2str(frameID) ];
        
        createFramefilePUSH(framefilename, tog, corrosion_level, t, frameID, countVal, nameSuffix)
        
        buildFrameinputPARALLEL
        
        buildRecordersPARALLEL2
        
        buildPushoverPARALLEL
        
        
    end
    
    fprintf('\n Pushover Run Files COMPLETE! \n');
    
end

% Command File

fprintf('\n\tCreating Job File...\n');
for t = tRange
    for frameNum = frameNumRange
        nNodes = length([Frames(frameNum).Node(:).X])/3.0;
        %        nNodes  = 43;
        nNodes = nNodes - rem(nNodes,2);
        nP = nNodes/2.0;
        
        frameID = Frames(frameNum).ID;
        for GM = numRange
            
            nPts = GroundMotions(GM).NPTS * 1.15 * 4;
            tog.walltime = ceil(nPts/2000); % Hours
            wtH = floor(tog.walltime); wtM = floor((tog.walltime - wtH)*60);
            for scaleFac = scaleRange
                
                if scaleFac == 1
                    scaledStr = '_DE';
                    scaledStr2 = 'DE';
                elseif scaleFac == 2
                    scaledStr = '_MCE';
                    scaledStr2 = 'MCE';
                end
                
                fmt     = sprintf(['_GM' repmat('_%u',1,numel(GM)) ],GM);
                
                % Create Job File
                jobfile     = sprintf([ '%uyr_%u%s%s%s_2'],t,frameID,scaledStr,nameSuffix,fmt);
                jobfilename = sprintf([jobfile '.lsf']);
                fid = fopen(fullfile([pwd '\Job Files\Ada\'],jobfilename),'w');
                
                fprintf(fid,'##NECESSARY JOB SPECIFICATIONS');
                fprintf(fid,'\n#BSUB -J %0.0fyr_%0.0f%s%s%s   #Set the job name to "ExampleJob1"',t,frameID,scaledStr,nameSuffix,fmt);
                fprintf(fid,'\n#BSUB -L /bin/bash             #Uses the bash login shell to initialize the jobs execution environment.');
                fprintf(fid,'\n#BSUB -W %02d:%02d             #Set the wall clock limit to 2hr',wtH,wtM);
                fprintf(fid,'\n#BSUB -n %0.0f                 #Request nodes',nNodes);
                fprintf(fid,'\n#BSUB -R "span[ptile=%0.0f]"     #Request cores per node.', nNodes/ceil(nNodes/20));
                fprintf(fid,'\n#BSUB -R "rusage[mem=150]"   #Request RAM per process (CPU) for the job');
                fprintf(fid,'\n#BSUB -M 150                 #Set the per process enforceable memory limit.');
                fprintf(fid,'\n#BSUB -o %0.0fyr_%0.0f%s%s%s_out.%%J    #Send stdout and stderr to "Example1Out.[jobID]" ',t,frameID,scaledStr,nameSuffix,fmt);
                fprintf(fid,'\n#BSUB -app resizable         #Allow dynamic release of resources');
                fprintf(fid,'\n##OPTIONAL JOB SPECIFICATIONS');
                fprintf(fid,'\n#BSUB -u cmckee@tamu.edu       #Send all emails to email_address');
                fprintf(fid,'\n#BSUB -B -N');
                
                fprintf(fid,'\n#First Executable Line:');
                % fprintf(fid,'\ncd "$SCRATCH/HPRC Frames/RunFiles"');
                fprintf(fid,'\nml myEB');
                fprintf(fid,'\nml OpenSees/3.2.0-foss-2020a-parallel');
                fprintf(fid,'\ncd "$SCRATCH/Models"');
%                 mpiexec -np 12 OpenSeesSP RunFile.tcl 1001 Corr2 "Model_Age.tcl" -dir -lib "./OpenSees Files" -runAge 5 -runGrav -runEigen -runDyn 34 -DE 4 0.15 -gmDir "./Data/Ground Motions"
                fprintf(fid,'\nmpiexec -np %u OpenSeesSP RunFile.tcl %u Corr%u "Model_Age2.tcl" -dir -lib "./OpenSees Files" -output "./Output_Presentation2" -runAge %u %0.1f -runGrav -runEigen -runDyn %u -%s 4 0.15 -gmDir "./Data/Ground Motions"   ',nNodes,frameID,corrosion_level,t,52.0, GM,scaledStr2);
                fprintf(fid,'\ncd $SCRATCH/Models/Output/Corr%u_%u/%u/Dynamic',corrosion_level,frameID,t);
                fprintf(fid,'\ntar -czf "AT_%0.0f%s.tar.gz" "./AT_%0.0f%s" ',GM,scaledStr,GM,scaledStr);
                fprintf(fid,'\nrm -rf "./AT_%0.0f%s" ',GM,scaledStr);
                fclose(fid);
            end
        end
    end
end

% Job File

for frameNum = frameNumRange
    frameID = Frames(frameNum).ID;
    
        for scaleFac = scaleRange
            
            if scaleFac == 1
                scaledStr = '_DE';
            elseif scaleFac == 2
                scaledStr = '_MCE';
            end
            
%             fmt     = sprintf(['_GM' repmat('_%u',1,numel(GM)) ],GM);
            % START FILE
            
            startfile     = sprintf([ 'START_' num2str(frameID) scaledStr nameSuffix]);
            startfilename = sprintf([startfile '.tcl']);
            fid = fopen(fullfile([pwd '\Job Files\Ada\'],startfilename),'w');
            fprintf(fid,'\n#!/usr/bin/bash           #Uses the bash login shell to initialize the jobs execution environment.');
            fprintf(fid,'\n\n ');
            for t = tRange
                for GM = numRange
                    fmt     = sprintf(['_GM' repmat('_%u',1,numel(GM)) ],GM);
                    
                    fprintf(fid,'\n\ndos2unix "./%uyr_%u%s%s%s.lsf" && bsub < "./%uyr_%u%s%s%s.lsf" && ',t,frameID,scaledStr,nameSuffix,fmt,t,frameID,scaledStr,nameSuffix,fmt);
                end
                fprintf(fid,'\n\necho "End of Age Group"');
                
            end
            fprintf(fid,'\n\necho "End of File"');
        fclose(fid);
        
        
%         % tamulauncher File
%         
%         tlfile     = sprintf([ 'TAMULAUNCHER_' num2str(frameID) scaledStr nameSuffix]);
%         tlfilename = sprintf([tlfile '.tcl']);
%         fid = fopen(fullfile([pwd '\Job Files\Ada\'],tlfilename),'w');
%         fprintf(fid,'\n#!/usr/bin/bash           #Uses the bash login shell to initialize the jobs execution environment.');
%         fprintf(fid,'\n\n ');
%         for t = tRange
%             
%             for job = 1:length(jobRange)-1
%                 
%                 GMrange = numRange(jobRange(job):jobRange(job+1)-1);
%                 fmt     = sprintf(['_GM' repmat('_%u',1,numel(GMrange)) ],GMrange);
%                 fprintf(fid,'\n\ntamulauncher --status %uyr_%u%s%s%s_commands.in && ',t,frameID,scaledStr,nameSuffix,fmt);
%                 
%                 
%                 
%                 
%             end
%             fprintf(fid,'\n\necho "End of Age Group"');
%         end
%         fprintf(fid,'\n\necho "End of File"');
%         fclose(fid);
        
        
    end
end


function [] = createFramefile(framefilename, tog, corrosion_level, t, frameID, countVal, GM,scaledStr, GroundMotions, nameSuffix)
% Create frameFile

fid = fopen(framefilename,'w');
%puts  [format $cmdOut [clock format [clock seconds] -format "%d-%b-%Y %T:"] "\t" $frameName $GroundMotion "Dynamic Analysis FAILED!"];

fprintf(fid,'#!/usr/bin/tclsh');
fprintf(fid,'\n#after [expr int(%0.0f*4)];',GM);
fprintf(fid,'\n# Codi McKee');
fprintf(fid,'\n# Texas A&M University');
fprintf(fid,'\n# %s',char(datetime('today')));
fprintf(fid,'\n\n# This file was written via Matlab %s',version);
fprintf(fid,'\n# File: %s\\%s.m\n\n',pwd,mfilename);
fprintf(fid,'\n set mainDir ../ ;');
fprintf(fid,'\n cd $mainDir ;');
fprintf(fid,'\n set N %s ;',num2str(countVal));
fprintf(fid,'\n set frameName "%s" ;',[ 'Corr' num2str(corrosion_level) '_' num2str(frameID) ]);

fprintf(fid,'\n set frameFile [join [list "Input_" $frameName %s ".tcl"] ""] ;', nameSuffix);
fprintf(fid,'\n set pushoverFile [join [list "Pushover_" $frameName ".tcl"] ""] ;');
% fprintf(fid,['\n set dynamicFile [join [list "Dynamic_" $frameName "_AT_' num2str(groundMotion_num) scaledStr '.tcl" ] ""];']);
fprintf(fid,'\n set recorderFile [join [list "Recorders_" $frameName %s ".tcl"] ""] ;', nameSuffix);
fprintf(fid,'\n set parameterFile [join [list "Parameters_" $frameName %s ".tcl"] ""] ;', nameSuffix);
fprintf(fid,['\n set togPush ' num2str(0) ' ;']);
fprintf(fid,['\n set togGravity ' num2str(tog.runGravity) ' ;']);
fprintf(fid,['\n set togDynamic ' num2str(1) ' ;']);
fprintf(fid,['\n set togAge ' num2str(1) ' ;']);
fprintf(fid,['\n set Age ' num2str(t) ' ;']);
if GM == 1
    fprintf(fid,['\n set togEigen ' num2str(1) ' ;']);
else
    fprintf(fid,['\n set togEigen ' num2str(1) ' ;']);
end
fprintf(fid,['\n set togDisplay ' num2str(tog.Display) ' ;']);
fprintf(fid,['\n set togDisplaymode ' num2str(tog.Display) ' ;']);

fprintf(fid,'\n\n# Dynamic Loads: Frame ID %s',num2str(frameID));

%**************************************************************************
% Ground Motion File
%
DT = GroundMotions(GM).DT ;
NPTS = GroundMotions(GM).NPTS;

fprintf(fid,['\n set GroundMotionFile AT_' num2str(GM) scaledStr '.tcl ;']);
fprintf(fid,['\n set GroundMotion AT_' num2str(GM) scaledStr ' ;']);
fprintf(fid,['\n set GM ' num2str(GM) ';']);
fprintf(fid,['\n set DT ' num2str(DT) ';']);
fprintf(fid,['\n set NPTS ' num2str(NPTS) ';']);
fprintf(fid,['\n set stpNum ' num2str(tog.stpNum) ';']);
fprintf(fid,['\n set nameSuffix "%s" ;'],nameSuffix);
addTime = ceil(NPTS*0.15);
GroundMotions(GM).AddTime = addTime;
GroundMotions(GM).TotalTime = NPTS + addTime;

fprintf(fid,['\n set walltime %0.0f ;  # Seconds'],tog.walltime*60*60);
fprintf(fid,'\n set addTime %0.0f;', addTime);
fprintf(fid,'\n\n# Set up Command Window Output');
fprintf(fid,'\n set cmdOut "\\n%%s %%s%%syr%%s_%%s: %%s" ;');

fprintf(fid,'\n\n\n cd "OpenSees Files" ;');
fprintf(fid,'\n puts  [format $cmdOut [clock format [clock seconds] -format "%%d-%%b-%%Y %%T:"] "" $frameName $Age  $GroundMotion "Starting Analysis..."];');

fprintf(fid,'\n source Model_Age.tcl ;');
fprintf(fid,'\n\n puts  [format $cmdOut [clock format [clock seconds] -format "%%d-%%b-%%Y %%T:"] "" $frameName $Age  $GroundMotion "Analysis Complete!"];');

fprintf(fid,'\n puts  [format $cmdOut [clock format [clock seconds] -format "%%d-%%b-%%Y %%T:"] "" $frameName $Age  $GroundMotion "Exiting...\\n\\n"];');
fprintf(fid,'\n quit ;');
fclose(fid);

end

function [] = createFramefilePUSH(framefilename, tog, corrosion_level, t, frameID, countVal, nameSuffix)
% Create frameFile
if tog.runPushover
    pushStr = '_P';
elseif tog.runDynamic
    pushStr = '';
end

fid = fopen(framefilename,'w');

% mainDir = ['{' pwd '}'];
fprintf(fid,'#!/usr/bin/tclsh');
fprintf(fid,'# Codi McKee');
fprintf(fid,'\n# Texas A&M University');
fprintf(fid,'\n# %s',char(datetime('today')));
fprintf(fid,'\n\n# This file was written via Matlab %s',version);
fprintf(fid,'\n# File: %s\\%s.m\n\n',pwd,mfilename);
fprintf(fid,'\n set mainDir "/scratch/user/mckee107/HPRC Frames/";');
fprintf(fid,'\n cd $mainDir ;');

fprintf(fid,'\n set N %s',num2str(countVal));
fprintf(fid,'\n set frameName "%s"',[ num2str(t) 'yr_Corr' num2str(corrosion_level) '_' num2str(frameID) ]);

fprintf(fid,'\n set frameFile [join [list "Input_" $frameName %s ".tcl"] ""];', nameSuffix);
fprintf(fid,'\n set pushoverFile [join [list "Pushover_" $frameName ".tcl"] ""];');
fprintf(fid,['\n set recorderFile [join [list "Recorders_" $frameName "' pushStr nameSuffix '.tcl"] ""];']);
fprintf(fid,['\n set togPush ' num2str(1) ';']);
fprintf(fid,['\n set togGravity ' num2str(tog.runGravity) ';']);
fprintf(fid,['\n set togDynamic ' num2str(0) ';']);
fprintf(fid,['\n set togEigen ' num2str(tog.runEigen) ';']);
fprintf(fid,['\n set togDisplay ' num2str(tog.Display) ';']);
fprintf(fid,['\n set GroundMotion "PUSHOVER" ;']);
fprintf(fid,['\n set GM 1 ;']);
fprintf(fid,'\n set nameSuffix "%s" ;',nameSuffix);
fprintf(fid,['\n set walltime %0.0f ;  # Seconds'],tog.walltime*60*60);

fprintf(fid,'\n\n# Set up Command Window Output');
fprintf(fid,'\n set cmdOut "\\n%%s %%s%%s_%%s: %%s" ;');

fprintf(fid,'\n\n\n cd "OpenSees Files" ;');
fprintf(fid,'\n puts  [format $cmdOut [clock format [clock seconds] -format "%%d-%%b-%%Y %%T:"] "" $frameName $GroundMotion "Starting Analysis..."];');

fprintf(fid,'\n source ModelDB_Push.tcl ;');
fprintf(fid,'\n\n puts  [format $cmdOut [clock format [clock seconds] -format "%%d-%%b-%%Y %%T:"] "" $frameName $GroundMotion "Analysis Complete!"];');

fprintf(fid,'\n puts  [format $cmdOut [clock format [clock seconds] -format "%%d-%%b-%%Y %%T:"] "" $frameName $GroundMotion "Exiting...\\n\\n"];');
fprintf(fid,'\n quit ;');
fclose(fid);

end
