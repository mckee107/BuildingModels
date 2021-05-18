% Codi McKee
% Texas A&M University
% First Created: 22-Feb-2019
% Last Modified: 24-Jan-2020
% TO DO:
%

% Build list of recorders for Model.tcl


%**************************************************************************
% Write Recorders File
%

if tog.runPushover
    pushStr = '_P';
elseif tog.runDynamic
    pushStr = '';
end

filename = ['Input\Recorders_Corr' num2str(corrosion_level) '_' num2str(frameID) pushStr nameSuffix '.tcl'];

fid = fopen(filename,'w');

% Header Information
fprintf(fid,'#!/usr/bin/tclsh');
fprintf(fid,'# Codi McKee');
fprintf(fid,'\n# Texas A&M University');
fprintf(fid,'\n# %s',char(datetime('today')));
fprintf(fid,'\n\n# This file was written via Matlab %s',version);
fprintf(fid,'\n# File: %s\\%s.m',pwd,mfilename);

fprintf(fid,'\n\n# Recorders Frame ID %s',num2str(frameID));
fprintf(fid,'\ncd $currentDir');


%**************************************************************************
% Make Recorders
%

fprintf(fid,'\n\n# ------------------------------------------------------------------------------');
fprintf(fid,'\n# Recorders\n');
fprintf(fid,'\n#\tNodes\n');

% Node Recorders

for colLine = 1:Frames(frameNum).Bays+1
    
    nodeSupp{1,colLine} = sprintf('10%.0f0',colLine);
    
    for flr = 2:Frames(frameNum).Stories+1
        
        nodeFrame{flr-1,colLine} = sprintf('%.0f0%.0f0',flr,colLine);
        
        if colLine == 1
            
            nodei{flr-1,colLine} = sprintf('%.0f010',flr-1);
            nodej{flr-1,colLine} = sprintf('%.0f0%.0f0',flr,colLine);
            
        end
        
    end
end

if tog.runDynamic
    fprintf(fid,'\n\t set dynStr "/$GroundMotion";');
elseif tog.runPushover
    fprintf(fid,'\n\t set dynStr "";');
end


fprintf(fid,'\n\trecorder Node -file $outputDir/$frameName/$Age/$analysisType$dynStr/RBase.out -time -node %s -dof 1 2 3 reaction',strjoin(nodeSupp));
fprintf(fid,'\n\trecorder Node -file $outputDir/$frameName/$Age/$analysisType$dynStr/DBase.out -time -node %s -dof 1 2 3 disp',strjoin(nodeSupp));
fprintf(fid,'\n\trecorder Node -file $outputDir/$frameName/$Age/$analysisType$dynStr/DFree.out -time -node %s -dof 1 2 3 disp',strjoin(nodeFrame));
fprintf(fid,'\n\trecorder Node -file $outputDir/$frameName/$Age/$analysisType$dynStr/VBase.out -time -node %s -dof 1 reaction',strjoin(nodeSupp));

% Drift Recorders

fprintf(fid,'\n\n#\tDrift Recorders\n');

% fprintf(fid,'\n\trecorder Drift -file $outputDir/$frameName/$Age/$analysisType$dynStr/DrNode.out -time -iNode %s  -jNode %s  -dof 1 -perpDirn 2;',strjoin(nodei),strjoin(nodej));
% fprintf(fid,'\n\trecorder Drift -file $outputDir/$frameName/$Age/$analysisType$dynStr/DrRoof.out -time -iNode %s  -jNode %s  -dof 1 -perpDirn 2;',nodei{1},nodej{end});


% Stress Strain Recorders

fprintf(fid,'\n\n#\tStress Strain and Force Recorders\n');
fprintf(fid,'\n#\t\tColumns\n');

for story = 1:Frames(frameNum).Stories
    
    for colLine = 1:Frames(frameNum).Bays+1
        
        if colLine == 1 || colLine == Frames(frameNum).Bays+1
            locType = 'Ext';
            h = Frames(frameNum).Columns(story).Ext.Depth;
            b = Frames(frameNum).Columns(story).Ext.Width;
            cover = Frames(frameNum).Columns(story).Ext.Cover;
            lc = Frames(frameNum).Columns(story).Ext.Width;
            R1 = Frames(frameNum).Columns(story).Ext.Width;
            R2 = Frames(frameNum).Columns(story).Int.Width;
            
        else
            locType = 'Int';
            h = Frames(frameNum).Columns(story).Int.Depth;
            b = Frames(frameNum).Columns(story).Int.Width;
            cover = Frames(frameNum).Columns(story).Int.Cover;
            lc = Frames(frameNum).Columns(story).Int.Width;
            R1 = Frames(frameNum).Columns(story).Int.Width;
            R2 = Frames(frameNum).Columns(story).Int.Width;
            
        end
        
        elemID  = sprintf('11%u%u',story,colLine);
        SmatID  = sprintf('11%u%u',story,colLine);
        UCmatID = sprintf('21%u%u',story,colLine);
        CEmatID = sprintf('31%u%u',story,colLine);
        CMmatID = sprintf('41%u%u',story,colLine);
        
        x_core = (b/2) - cover;
        y_core = (h/2) - cover;
        
        sy = y_core*2 / (n_h - 1);
        %         numIntgPts_c = 1
        if story == 1
            fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColSteel_%s_B_Li_TH.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,1,y_core,x_core,SmatID);
            fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColSteel_%s_T_Li_TH.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,numIntgPts_c,y_core,x_core,SmatID);
            fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColSteel_%s_B_Ri_TH.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,1,-y_core,x_core,SmatID);
            fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColSteel_%s_T_Ri_TH.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,numIntgPts_c,-y_core,x_core,SmatID);
            
            if tog.Steel02 ~= 1
                fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/Delta_%s_B_Li_TH.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s Delta',elemID,elemID,1,y_core,x_core,SmatID);
                fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/Delta_%s_T_Li_TH.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s Delta',elemID,elemID,numIntgPts_c,y_core,x_core,SmatID);
                fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/Delta_%s_B_Ri_TH.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s Delta',elemID,elemID,1,-y_core,x_core,SmatID);
                fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/Delta_%s_T_Ri_TH.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s Delta',elemID,elemID,numIntgPts_c,-y_core,x_core,SmatID);
            end
        end
        %         fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColSteel_%s_Ri.out -time -ele %s section 1 fiber %0.3f %0.3f %s stressStrain',elemID,elemID,x_core,-y_core,SmatID);
        fprintf(fid,'\n\trecorder EnvelopeElement  -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColSteel_%s_B_Li.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,1,y_core,x_core,SmatID);
        fprintf(fid,'\n\trecorder EnvelopeElement  -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColSteel_%s_B_Ri.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,1,-y_core,x_core,SmatID);
        fprintf(fid,'\n\trecorder EnvelopeElement  -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColSteel_%s_T_Li.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,numIntgPts_c,y_core,x_core,SmatID);
        fprintf(fid,'\n\trecorder EnvelopeElement  -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColSteel_%s_T_Ri.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,numIntgPts_c,-y_core,x_core,SmatID);
        
        %         fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColSteel_%s_Lj.out -time -ele %s section %s fiber %0.3f %0.3f %s stressStrain',elemID,elemID,num2str(numIntgPts_c),y_core,x_core,SmatID);
        %         fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColSteel_%s_Rj.out -time -ele %s section %s fiber %0.3f %0.3f %s stressStrain',elemID,elemID,num2str(numIntgPts_c),-y_core,x_core,SmatID);
        %         fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColSteel_%s_Lm.out -time -ele %s section %s fiber %0.3f %0.3f %s stressStrain',elemID,elemID,num2str((numIntgPts_c+1)/2),y_core,x_core,SmatID);
        %         fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColSteel_%s_Rm.out -time -ele %s section %s fiber %0.3f %0.3f %s stressStrain',elemID,elemID,num2str((numIntgPts_c+1)/2),-y_core,x_core,SmatID);
        
        fprintf(fid,'\n\trecorder EnvelopeElement -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColCover_%s_B_Li.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,1,0,b/2,UCmatID);
        fprintf(fid,'\n\trecorder EnvelopeElement -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColCover_%s_B_Ri.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,1,0,-b/2,UCmatID);
        fprintf(fid,'\n\trecorder EnvelopeElement -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColCover_%s_T_Li.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,numIntgPts_c,0,b/2,UCmatID);
        fprintf(fid,'\n\trecorder EnvelopeElement -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColCover_%s_T_Ri.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,numIntgPts_c,0,-b/2,UCmatID);
        %         fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColCover_%s_Lj.out -time -ele %s section %s fiber %0.3f %0.3f %s stressStrain',elemID,elemID,num2str(numIntgPts_c),0,b/2,UCmatID);
        %         fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColCover_%s_Rj.out -time -ele %s section %s fiber %0.3f %0.3f %s stressStrain',elemID,elemID,num2str(numIntgPts_c),0,-b/2,UCmatID);
        %         fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColCover_%s_Lm.out -time -ele %s section %s fiber %0.3f %0.3f %s stressStrain',elemID,elemID,num2str((numIntgPts_c+1)/2),0,b/2,UCmatID);
        %         fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColCover_%s_Rm.out -time -ele %s section %s fiber %0.3f %0.3f %s stressStrain',elemID,elemID,num2str((numIntgPts_c+1)/2),0,-b/2,UCmatID);
        
        fprintf(fid,'\n\trecorder EnvelopeElement -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColCore_%s_B_Li.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,1,0,x_core,CEmatID);
        fprintf(fid,'\n\trecorder EnvelopeElement -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColCore_%s_B_Ri.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,1,0,-x_core,CEmatID);
        fprintf(fid,'\n\trecorder EnvelopeElement -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColCore_%s_T_Li.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,numIntgPts_c,0,x_core,CEmatID);
        fprintf(fid,'\n\trecorder EnvelopeElement -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColCore_%s_T_Ri.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,numIntgPts_c,0,-x_core,CEmatID);
        %         fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColCore_%s_Lj.out -time -ele %s section %s fiber %0.3f %0.3f %s stressStrain',elemID,elemID,num2str(numIntgPts_c),0,x_core,CEmatID);
        %         fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColCore_%s_Rj.out -time -ele %s section %s fiber %0.3f %0.3f %s stressStrain',elemID,elemID,num2str(numIntgPts_c),0,-x_core,CEmatID);
        %         fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColCore_%s_Lm.out -time -ele %s section %s fiber %0.3f %0.3f %s stressStrain',elemID,elemID,num2str((numIntgPts_c+1)/2),0,x_core,CEmatID);
        %         fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/ColCore_%s_Rm.out -time -ele %s section %s fiber %0.3f %0.3f %s stressStrain',elemID,elemID,num2str((numIntgPts_c+1)/2),0,-x_core,CEmatID);
        
        fprintf(fid,'\n\trecorder EnvelopeElement -file $outputDir/$frameName/$Age/$analysisType$dynStr/IShear_%s.out -time -ele %s globalForce',elemID,elemID);
        
    end
end


% Beams

fprintf(fid,'\n\n#\tBeams\n');

for flr = 2:Frames(frameNum).Stories+1
    
    for bay = 1:Frames(frameNum).Bays
        
        if bay == 1 || bay == Frames(frameNum).Bays
            locType = 'Ext';
            h = Frames(frameNum).Beams(flr-1).Ext.Depth;
            b = Frames(frameNum).Beams(flr-1).Ext.Width;
            cover = Frames(frameNum).Beams(flr-1).Ext.Cover;
            ln = Frames(frameNum).Beams(flr-1).Ext.Length*12 - (Frames(frameNum).Columns(flr-1).Ext.Width + Frames(frameNum).Columns(flr-1).Int.Width )/2;
            b_eff = min(8*h,ln/2); % Only for interior frame
            n_slab = ceil(((b_eff-b)/2)/6);
            d_slab = tf - 2.5;
        else
            locType = 'Int';
            h = Frames(frameNum).Beams(flr-1).Int.Depth;
            b = Frames(frameNum).Beams(flr-1).Int.Width;
            cover = Frames(frameNum).Beams(flr-1).Int.Cover;
            ln = Frames(frameNum).Beams(flr-1).Int.Length*12 - (Frames(frameNum).Columns(flr-1).Int.Width + Frames(frameNum).Columns(flr-1).Int.Width )/2;
            b_eff = min(8*h,ln/2);
            n_slab = ceil(((b_eff-b)/2)/6);
            d_slab = tf - 2.5;
        end
        
        elemID  = sprintf('12%u%u',flr,bay);
        sectID  = sprintf('12%u%u',flr,bay);
        sectIDmid  = sprintf('22%u%u',flr,bay);
        SmatID  = sprintf('12%u%u',flr,bay);
        SmatID_c = sprintf('52%u%u',flr,bay);
        UCmatID = sprintf('22%u%u',flr,bay);
        CEmatID = sprintf('32%u%u',flr,bay);
        CMmatID = sprintf('42%u%u',flr,bay);
        
        x_core = (b/2) - cover;
        y_core = (h/2) - cover;
        
        
        fprintf(fid,'\n\trecorder EnvelopeElement -file $outputDir/$frameName/$Age/$analysisType$dynStr/BeamSteel_%s_L_Ti.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,1,y_core,x_core,SmatID_c);
        fprintf(fid,'\n\trecorder EnvelopeElement -file $outputDir/$frameName/$Age/$analysisType$dynStr/BeamSteel_%s_L_Bi.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,1,-y_core,x_core,SmatID);
        fprintf(fid,'\n\trecorder EnvelopeElement -file $outputDir/$frameName/$Age/$analysisType$dynStr/BeamSteel_%s_R_Ti.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,numIntgPts_b,y_core,x_core,SmatID_c);
        fprintf(fid,'\n\trecorder EnvelopeElement -file $outputDir/$frameName/$Age/$analysisType$dynStr/BeamSteel_%s_R_Bi.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,numIntgPts_b,-y_core,x_core,SmatID);
        
        %         fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/BeamSteel_%s_Tj.out -time -ele %s section %s fiber %0.3f %0.3f %s stressStrain',elemID,elemID,num2str(numIntgPts_b),x_core,y_core,SmatID_c);
        %         fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/BeamSteel_%s_Bj.out -time -ele %s section %s fiber %0.3f %0.3f %s stressStrain',elemID,elemID,num2str(numIntgPts_b),x_core,-y_core,SmatID);
        %         fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/BeamSteel_%s_Tm.out -time -ele %s section %s fiber %0.3f %0.3f %s stressStrain',elemID,elemID,num2str((numIntgPts_b+1)/2),x_core,y_core,SmatID_c);
        %         fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/BeamSteel_%s_Bm.out -time -ele %s section %s fiber %0.3f %0.3f %s stressStrain',elemID,elemID,num2str((numIntgPts_b+1)/2),x_core,-y_core,SmatID);
        
        fprintf(fid,'\n\trecorder EnvelopeElement -file $outputDir/$frameName/$Age/$analysisType$dynStr/BeamCover_%s_L_Ti.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,1,0,b/2,UCmatID);
        fprintf(fid,'\n\trecorder EnvelopeElement -file $outputDir/$frameName/$Age/$analysisType$dynStr/BeamCover_%s_L_Bi.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,1,0,-b/2,UCmatID);
        fprintf(fid,'\n\trecorder EnvelopeElement -file $outputDir/$frameName/$Age/$analysisType$dynStr/BeamCover_%s_R_Ti.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,numIntgPts_b,0,b/2,UCmatID);
        fprintf(fid,'\n\trecorder EnvelopeElement -file $outputDir/$frameName/$Age/$analysisType$dynStr/BeamCover_%s_R_Bi.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,numIntgPts_b,0,-b/2,UCmatID);
        %         fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/BeamCover_%s_Tj.out -time -ele %s section %s fiber %0.3f %0.3f %s stressStrain',elemID,elemID,num2str(numIntgPts_b),b/2,0,UCmatID);
        %         fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/BeamCover_%s_Bj.out -time -ele %s section %s fiber %0.3f %0.3f %s stressStrain',elemID,elemID,num2str(numIntgPts_b),-b/2,0,UCmatID);
        %         fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/BeamCover_%s_Tm.out -time -ele %s section %s fiber %0.3f %0.3f %s stressStrain',elemID,elemID,num2str((numIntgPts_b+1)/2),b/2,0,UCmatID);
        %         fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/BeamCover_%s_Bm.out -time -ele %s section %s fiber %0.3f %0.3f %s stressStrain',elemID,elemID,num2str((numIntgPts_b+1)/2),-b/2,0,UCmatID);
        
        fprintf(fid,'\n\trecorder EnvelopeElement -file $outputDir/$frameName/$Age/$analysisType$dynStr/BeamCore_%s_L_Ti.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,1,0,x_core,CEmatID);
        fprintf(fid,'\n\trecorder EnvelopeElement -file $outputDir/$frameName/$Age/$analysisType$dynStr/BeamCore_%s_L_Bi.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,1,0,-x_core,CEmatID);
        fprintf(fid,'\n\trecorder EnvelopeElement -file $outputDir/$frameName/$Age/$analysisType$dynStr/BeamCore_%s_R_Ti.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,numIntgPts_b,0,x_core,CEmatID);
        fprintf(fid,'\n\trecorder EnvelopeElement -file $outputDir/$frameName/$Age/$analysisType$dynStr/BeamCore_%s_R_Bi.out -time -ele %s section %0.0f fiber %0.3f %0.3f %s stressStrain',elemID,elemID,numIntgPts_b,0,-x_core,CEmatID);
        %         fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/BeamCore_%s_Tj.out -time -ele %s section %s fiber %0.3f %0.3f %s stressStrain',elemID,elemID,num2str(numIntgPts_b),x_core,0,CEmatID);
        %         fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/BeamCore_%s_Bj.out -time -ele %s section %s fiber %0.3f %0.3f %s stressStrain',elemID,elemID,num2str(numIntgPts_b),-x_core,0,CEmatID);
        %         fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/BeamCore_%s_Tm.out -time -ele %s section %s fiber %0.3f %0.3f %s stressStrain',elemID,elemID,num2str((numIntgPts_b+1)/2),x_core,0,CEmatID);
        %         fprintf(fid,'\n\trecorder Element -file $outputDir/$frameName/$Age/$analysisType$dynStr/BeamCore_%s_Bm.out -time -ele %s section %s fiber %0.3f %0.3f %s stressStrain',elemID,elemID,num2str((numIntgPts_b+1)/2),-x_core,0,CEmatID);
        %
        
    end
end

fprintf(fid,'\n\n\n\n\n\n\n\n\n\n');
fclose(fid);








