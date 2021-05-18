% Codi McKee
% Texas A&M University
% First Created: 05-Feb-2019
% Last Modified: 30-May-2019

% Run OpenSees model for selected frame from database

%**************************************************************************
% Initialize
%
close all; fclose all; clc; diary off; format short g; 

A01_initializeWorkspace_2; % Calls external .m file


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
                countVal = countVal + 1;
            end
        end
    end
    
end

fprintf('\n\tCreating Command File...\n');

for t = tRange
    for frameNum = frameNumRange
        frameID = Frames(frameNum).ID;
        for scaleFac = scaleRange
            
            if scaleFac == 1
                scaledStr = '_DE';
            elseif scaleFac == 2
                scaledStr = '_MCE';
            end
            
            nameSuffix = '_Steel02';
            
            % Create Command file
            cmdfilename = [num2str(t) 'yr_' num2str(frameID) scaledStr nameSuffix '_commands.in'];
            fid = fopen(fullfile([pwd '\Job Files\Terra\'],cmdfilename),'w');
            for GM = numRange
                
                part0 = sprintf('cd "$SCRATCH/HPRC Frames/RunFiles"');
                part1 = sprintf('OpenSees < RF_%0.0fyr_Corr2_%0.0f_AT_%0.0f%s%s.tcl',t,frameID,GM,scaledStr, nameSuffix);
                part2 = sprintf('cd $SCRATCH');
                part3 = sprintf('tar -czf "./HPRC Frames/Output%s/%0.0fyr_Corr2_%0.0f/Dynamic/AT_%0.0f%s.tar.gz"',nameSuffix,t,frameID,GM,scaledStr);
                part4 = sprintf('"./HPRC Frames/Output%s/%0.0fyr_Corr2_%0.0f/Dynamic/AT_%0.0f%s" ',nameSuffix,t,frameID,GM,scaledStr);
                part5 = sprintf('rm -rf "./HPRC Frames/Output%s/%0.0fyr_Corr2_%0.0f/Dynamic/AT_%0.0f%s" ', nameSuffix,t,frameID,GM,scaledStr);
                fprintf(fid,'(%s&& %s&& %s&& %s %s&& %s)\n', part0, part1, part2, part3, part4, part5);
                
            end
            fclose(fid);
            
            % Create Job File
            jobfile     = sprintf([ num2str(t) 'yr_' num2str(frameID) scaledStr nameSuffix]);
            jobfilename = sprintf([jobfile '.slurm']);
            fid = fopen(fullfile([pwd '\Job Files\Terra\'],jobfilename),'w');
            fprintf(fid,'#!/bin/bash');
            fprintf(fid,'\n##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION');
            fprintf(fid,'\n#SBATCH --export=NONE                #Do not propagate environment');
            fprintf(fid,'\n#SBATCH --get-user-env=L             #Replicate login environment');
            
            fprintf(fid,'\n##NECESSARY JOB SPECIFICATIONS');
            fprintf(fid,'\n#SBATCH --job-name=%0.0fyr_%0.0f%s%s           #Set the job name to "ExampleJob1"',t,frameID,scaledStr,nameSuffix);
            fprintf(fid,'\n#SBATCH --time=10:00:00              #Set the wall clock limit to 1hr and 30min');
            fprintf(fid,'\n#SBATCH --ntasks=44                   #Request 1 task');
            fprintf(fid,'\n#SBATCH --ntasks-per-node=2          #Request 2 tasks/cores per node');
            fprintf(fid,'\n#SBATCH --mem=350M                    #Request 2560MB (2.5GB) per node');
            fprintf(fid,'\n#SBATCH --output=%0.0fyr_%0.0f%s%s_out.%%j      #Send stdout/err to "Example1Out.[jobID]"',t,frameID,scaledStr,nameSuffix);
            
            fprintf(fid,'\n##OPTIONAL JOB SPECIFICATIONS');
            fprintf(fid,'\n#SBATCH --mail-type=ALL              #Send email on all job events');
            fprintf(fid,'\n#SBATCH --mail-user=cmckee@tamu.edu    #Send all emails to email_address');
            
            fprintf(fid,'\n#First Executable Line');

            fprintf(fid,'\nml foss/2019b');
            fprintf(fid,'\nml Tcl/8.6.9-GCCcore-8.3.0');
            fprintf(fid,'\nml GCCcore/8.3.0');
            fprintf(fid,'\nml EasyBuild-terra-SCRATCH');
            fprintf(fid,'\nml OpenSees/3.2.0-CDM-foss-2019b');
            fprintf(fid,'\ntamulauncher --commands-pernode 1 $SCRATCH/JobFiles/%s',cmdfilename);
                       
            fclose(fid);
        end
    end
end

% fid = fopen(cmdfilename,'w');
%             % fprintf(fid,'\n# Codi McKee');
%             % fprintf(fid,'\n# Texas A&M University');
%             % fprintf(fid,'\n# %s',char(datetime('today')));
%             % fprintf(fid,'\n# This file was written via Matlab %s',version);
%             % fprintf(fid,'\n# File: %s\\%s.m\n\n',pwd,mfilename);
%             % fprintf(fid,'\n# Executable Commands\n');
%             
%             
%             for countVal = 1:length(time)
%                 frameNum = frameNumval(countVal);
%                 frameID = Frames(frameNum).ID;
%                 t = time(countVal);
%                 groundMotion_num = GM(countVal);
%                 scaledStr = scaleVal{countVal};
%                 part0 = sprintf('cd "$SCRATCH/HPRC Frames/RunFiles/"');
%                 part1 = sprintf('OpenSees < RF_%0.0fyr_Corr2_%0.0f_AT_%0.0f%s.tcl',t,frameID,groundMotion_num,scaledStr);
%                 part2 = sprintf('cd $SCRATCH');
%                 part3 = sprintf('tar -czf "./HPRC Frames/Output/%0.0fyr_Corr2_%0.0f/Dynamic/AT_%0.0f%s.tar.gz"',t,frameID,groundMotion_num,scaledStr);
%                 part4 = sprintf('"./HPRC Frames/Output/%0.0fyr_Corr2_%0.0f/Dynamic/AT_%0.0f%s" ',t,frameID,groundMotion_num,scaledStr);
%                 part5 = sprintf('rm -rf "./HPRC Frames/Output/%0.0fyr_Corr2_%0.0f/Dynamic/AT_%0.0f%s" ',t,frameID,groundMotion_num,scaledStr);
%                 fprintf(fid,'(%s; %s; %s; %s %s; %s)\n', part0, part1, part2, part3, part4, part5);
%                 
%                 
%             end
%             fclose(fid);

% 
% function [] = createFramefile(framefilename, tog, corrosion_level, t, frameID, countVal, GM,scaledStr, GroundMotions)
% % Create frameFile
% 
% fid = fopen(framefilename,'w');
% %puts  [format $cmdOut [clock format [clock seconds] -format "%d-%b-%Y %T:"] "\t" $frameName $GroundMotion "Dynamic Analysis FAILED!"];
% 
% fprintf(fid,'#!/usr/bin/tclsh');
% fprintf(fid,'\n# Codi McKee');
% fprintf(fid,'\n# Texas A&M University');
% fprintf(fid,'\n# %s',char(datetime('today')));
% fprintf(fid,'\n\n# This file was written via Matlab %s',version);
% fprintf(fid,'\n# File: %s\\%s.m\n\n',pwd,mfilename);
% fprintf(fid,'\n set mainDir ../ ;');
% fprintf(fid,'\n cd $mainDir ;');
% fprintf(fid,'\n set N %s ;',num2str(countVal));
% fprintf(fid,'\n set frameName "%s" ;',[ num2str(t) 'yr_Corr' num2str(corrosion_level) '_' num2str(frameID) ]);
% 
% fprintf(fid,'\n set frameFile [join [list "Input_" $frameName ".tcl"] ""] ;');
% fprintf(fid,'\n set pushoverFile [join [list "Pushover_" $frameName ".tcl"] ""] ;');
% % fprintf(fid,['\n set dynamicFile [join [list "Dynamic_" $frameName "_AT_' num2str(groundMotion_num) scaledStr '.tcl" ] ""];']);
% fprintf(fid,'\n set recorderFile [join [list "Recorders_" $frameName ".tcl"] ""] ;');
% fprintf(fid,['\n set togPush ' num2str(tog.runPushover) ' ;']);
% fprintf(fid,['\n set togGravity ' num2str(tog.runGravity) ' ;']);
% fprintf(fid,['\n set togDynamic ' num2str(tog.runDynamic) ' ;']);
% if GM == 1
%     fprintf(fid,['\n set togEigen ' num2str(1) ' ;']);
% else
%     fprintf(fid,['\n set togEigen ' num2str(1) ' ;']);
% end
% fprintf(fid,['\n set togDisplay ' num2str(tog.Display) ' ;']);
% fprintf(fid,['\n set togDisplaymode ' num2str(tog.Display) ' ;']);
% 
% fprintf(fid,'\n\n# Dynamic Loads: Frame ID %s',num2str(frameID));
% 
% %**************************************************************************
% % Ground Motion File
% %
% DT = GroundMotions(GM).DT ;
% NPTS = GroundMotions(GM).NPTS;
% 
% fprintf(fid,['\n set GroundMotionFile AT_' num2str(GM) scaledStr '.tcl ;']);
% fprintf(fid,['\n set GroundMotion AT_' num2str(GM) scaledStr ' ;']);
% fprintf(fid,['\n set GM ' num2str(GM) ';']);
% fprintf(fid,['\n set DT ' num2str(DT) ';']);
% fprintf(fid,['\n set NPTS ' num2str(NPTS) ';']);
% fprintf(fid,['\n set stpNum ' num2str(tog.stpNum) ';']);
% addTime = ceil(NPTS/15);
% GroundMotions(GM).AddTime = addTime;
% GroundMotions(GM).TotalTime = NPTS + addTime;
% 
% fprintf(fid,'\n set addTime %0.0f;', addTime);
% fprintf(fid,'\n\n# Set up Command Window Output');
% fprintf(fid,'\n set cmdOut {"\\n%%s %%s%%s_%%s: %%s"} ;');
% 
% fprintf(fid,'\n\n\n cd "OpenSees Files" ;');
% fprintf(fid,'\n puts  [format $cmdOut [clock format [clock seconds] -format "%%d-%%b-%%Y %%T:"] "" $frameName $GroundMotion "Starting Analysis..."];');
% 
% fprintf(fid,'\n source ModelDB.tcl ;');
% fprintf(fid,'\n\n puts  [format $cmdOut [clock format [clock seconds] -format "%%d-%%b-%%Y %%T:"] "" $frameName $GroundMotion "Analysis Complete!"];');
% 
% fprintf(fid,'\n puts  [format $cmdOut [clock format [clock seconds] -format "%%d-%%b-%%Y %%T:"] "" $frameName $GroundMotion "Exiting..."];');
% fprintf(fid,'\n quit ;');
% fclose(fid);
% 
% end
