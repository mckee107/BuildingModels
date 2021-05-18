% Codi McKee
% Texas A&M University
% First Created: 05-Feb-2019
% Last Modified: 30-May-2019
% TO DO:
%

% Build frame geometry Opensees input file for Haselton frames

%**************************************************************************
% Initialize
%
% clearvars; close all; fclose all; clc;
% addpath([pwd '\Input'])
% addpath([pwd '\Material Models']);
% addpath([pwd '\Other Functions']);
%
% load Frames.mat

% plotVal = 0;

% frameID = 1001;
%
% frameNum = find([Frames.ID] == frameID);

% for frameNum = 1:29
%     frameID = Frames(frameNum).ID;


%**************************************************************************
% Input Parameters
%

% Reinforcing steel
% uniaxialMaterial Steel02 tag? fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?> siginit?> <-deter eps1? eps2? rf_min?>
% uniaxialMaterial NMBuckling tag? As? lb? fy?	 E?	   esh?	    delta_0? k? b? gamma? N1? N2? b_dam? c_dam? edam?

ageVal = 1;
Es = 29000; % [ksi] Elastic Modulus
Esi0 = Es;

if tog.Buckling
    delta_0_m   = 0.0001;
else
    delta_0_m = 0;
end

k_sh        = 5;
b_sh        = 1;
gamma       = 34;
N1          = 7.8;
N2          = 23;
b_dam       = 11.3;
c_dam       = 1.78;
edam        = 0.137;
b_dam_stff  = 0;
c_dam_stff  = 1.78;
edam_stff   = 0.0;

numPts = 17;
maxItr = 100;


% Modeling Parameters
numIntgPts_b = 25;
numIntgPts_c = 21;
integrType = 'Simpson';
elemType  = 'GI';
% elemType  = 'FBC';
maxTol = '1E-6';
minTol = '1E-8';
numIntGI = 100;

% Slab Parameters
tf = 8;         % [in.]
barSize_slab = 8;
As_slab = 0.79;

% Rigid Link Properties
E_link = 100*4500 ; % [ksi]
G_link = E_link/2.3 ; % [ksi]
% J_link = 9999 ; % [ksi]
% Iz_link = 9999 ; % [in4]
% Iy_link = 9999 ; % [in4]

% Loading Parameters

wc = 150;   % [pcf]  Density of RC
g  = 32.2;  % [ft/s2] Acceleration due to gravity
LL = 50;   % [psf] Live Load
DL = 175;  % [psf] Dead Load


%**************************************************************************
% Write Input File
%

filename = ['Input\Input_Corr' num2str(corrosion_level) '_' num2str(frameID) nameSuffix '.tcl'];

fid = fopen(filename,'w');

% Header Information
fprintf(fid,'#!/usr/bin/tclsh');
fprintf(fid,'# Codi McKee');
fprintf(fid,'\n# Texas A&M University');
fprintf(fid,'\n# %s',char(datetime('today')));
fprintf(fid,'\n\n# This file was written via Matlab %s',version);
fprintf(fid,'\n# File: %s\\%s.m',pwd,mfilename);

fprintf(fid,'\n\n# Input file for Frame ID %s',num2str(frameID));
fprintf(fid,'\n#\t\t Number of Stories: %u',Frames(frameNum).Stories);
fprintf(fid,'\n#\t\t Number of Bays: %u',Frames(frameNum).Bays);

% Frame Configuration

fprintf(fid,'\n\n# ------------------------------------------------------------------------------');
fprintf(fid,'\n# Frame Configuration');
fprintf(fid,'\n# \n');
fprintf(fid,'\n set NStories %u',Frames(frameNum).Stories);
fprintf(fid,'\n set NBays %u',Frames(frameNum).Bays);
LBuilding = 0;
for stor = 1:Frames(frameNum).Stories
    LBuilding = LBuilding + Frames(frameNum).Columns(stor).Ext.Length*12;
end
fprintf(fid,'\n set LBuilding [expr %u*$in]',LBuilding);

% Column Geometry

fprintf(fid,'\n\n# ------------------------------------------------------------------------------');
fprintf(fid,'\n# Column Geometry');
fprintf(fid,'\n#');

for stor = 1:Frames(frameNum).Stories
    
    fprintf(fid,'\n\n#\tStory %u',stor);
    
    fprintf(fid,'\n set colLength_%u [expr %u *$in]',stor,Frames(frameNum).Columns(stor).Ext.Length*12);
    fprintf(fid,'\n');
    
    for loc = 1:2%Frames(frameNum).Bays+1
        if loc == 1 %|| loc == Frames(frameNum).Bays+1
            locType = 'Ext';
        else
            locType = 'Int';
        end
        
        fprintf(fid,'\n set colDepth_%u_%s [expr %u *$in]',stor,locType,Frames(frameNum).Columns(stor).(locType).Depth);
        fprintf(fid,'\n set colWidth_%u_%s [expr %u *$in]',stor,locType,Frames(frameNum).Columns(stor).(locType).Width);
        
        fprintf(fid,'\n');
        
    end
end


% Beam Geometry

fprintf(fid,'\n\n# ------------------------------------------------------------------------------');
fprintf(fid,'\n# Beam Geometry');
fprintf(fid,'\n#');

for flr = 1:Frames(frameNum).Stories
    
    fprintf(fid,'\n\n#\tFloor %u',flr+1);
    
    fprintf(fid,'\n set beamLength_%u [expr %u *$in]',flr+1,Frames(frameNum).Beams(flr).Ext.Length*12);
    fprintf(fid,'\n');
    
    for colLine = 1:2%Frames(frameNum).Bays+1
        if colLine == 1 %|| loc == Frames(frameNum).Bays+1
            locType = 'Ext';
        else
            locType = 'Int';
        end
        
        fprintf(fid,'\n set beamDepth_%u_%s [expr %u *$in]',flr+1,locType,Frames(frameNum).Beams(flr).(locType).Depth);
        fprintf(fid,'\n set beamWidth_%u_%s [expr %u *$in]',flr+1,locType,Frames(frameNum).Beams(flr).(locType).Width);
        
        fprintf(fid,'\n');
        
    end
end

% Nodal Coordinates

fprintf(fid,'\n\n# ------------------------------------------------------------------------------');
fprintf(fid,'\n# Nodal Coordinates');
fprintf(fid,'\n# \n');
fprintf(fid,'\n#\tBeam/Column Joints');
fprintf(fid,'\n#\tnode "ABCD" X-coord Y-coord');
fprintf(fid,'\n#\t\tA: "Floor"');
fprintf(fid,'\n#\t\tB: 0');
fprintf(fid,'\n#\t\tC: "Column"');
fprintf(fid,'\n#\t\tD: 0\n');

% Beam/Column Joints

for flr = 1:Frames(frameNum).Stories+1
    
    if flr == 1
        y_coord = 0;
    else
        y_coord = y_coord + Frames(frameNum).Columns(flr-1).Ext.Length*12;
    end
    
    for colLine = 1:Frames(frameNum).Bays+1
        
        if colLine == 1
            x_coord = 0;
        elseif flr == 1
            x_coord = x_coord + Frames(frameNum).Beams(flr).Ext.Length*12;
        else
            x_coord = x_coord + Frames(frameNum).Beams(flr-1).Ext.Length*12;
        end
        fprintf(fid,'\nnode %u0%u0 [expr %u*$in] [expr %u*$in]',flr,colLine,x_coord,y_coord);
        nodeNum = str2double(sprintf('%u0%u0',flr,colLine));
        Frames(frameNum).Node(nodeNum).X = x_coord;
        Frames(frameNum).Node(nodeNum).Y = y_coord;
        
        
    end
end


% Column Rigid Link Nodes

fprintf(fid,'\n\n#\tColumn Rigid Links');
fprintf(fid,'\n#\tnode "ABCD" X-coord Y-coord');
fprintf(fid,'\n#\t\tA: "Floor"');
fprintf(fid,'\n#\t\tB: 1 = "Lower" rigid link; 2 = "Upper" rigid link');
fprintf(fid,'\n#\t\tC: "Column"');
fprintf(fid,'\n#\t\tD: 0\n');


for flr = 1:Frames(frameNum).Stories
    
    for colLine = 1:Frames(frameNum).Bays+1
        
        if flr == 1
            
            if colLine == 1 || colLine == Frames(frameNum).Bays+1
                colLinklength2 = Frames(frameNum).Beams(flr).Ext.Depth/2;
            else
                colLinklength2 = max(Frames(frameNum).Beams(flr).Ext.Depth/2,Frames(frameNum).Beams(flr).Int.Depth/2);
            end
            
        elseif flr > 1 && flr ~= Frames(frameNum).Stories+1
            
            if colLine == 1 || colLine == Frames(frameNum).Bays+1
                colLinklength1 = Frames(frameNum).Beams(flr-1).Ext.Depth/2;
                colLinklength2 = Frames(frameNum).Beams(flr).Ext.Depth/2;
            else
                colLinklength1 = max(Frames(frameNum).Beams(flr-1).Ext.Depth/2,Frames(frameNum).Beams(flr-1).Int.Depth/2);
                colLinklength2 = max(Frames(frameNum).Beams(flr).Ext.Depth/2,Frames(frameNum).Beams(flr).Int.Depth/2);
            end
            
        end
        
        for colLink = 1:2
            
            if flr == 1
                
                if colLink == 1
                    continue
                else
                    nodeNum = str2double(sprintf('%u0%u0',flr,colLine));
                    colLinknode_x = Frames(frameNum).Node(nodeNum).X;
                    colLinknode_y = Frames(frameNum).Node(nodeNum).Y + Frames(frameNum).Columns(flr).Ext.Length*12 - colLinklength2;
                    linkNodenum = str2double(sprintf('%u%u%u0',flr,colLink,colLine));
                    fprintf(fid,'\nnode %u%u%u0 [expr %u*$in] [expr %u*$in]',flr,colLink,colLine,colLinknode_x,colLinknode_y);
                    Frames(frameNum).Node(linkNodenum).X = colLinknode_x;
                    Frames(frameNum).Node(linkNodenum).Y = colLinknode_y;
                end
                
            else
                nodeNum = str2double(sprintf('%u0%u0',flr,colLine));
                colLinknode_x = Frames(frameNum).Node(nodeNum).X;
                
                if colLink == 1
                    colLinknode_y = Frames(frameNum).Node(nodeNum).Y + colLinklength1;
                else
                    colLinknode_y = Frames(frameNum).Node(nodeNum).Y + Frames(frameNum).Columns(flr).Ext.Length*12 - colLinklength2;
                    
                end
                
                linkNodenum = str2double(sprintf('%u%u%u0',flr,colLink,colLine));
                fprintf(fid,'\nnode %u%u%u0 [expr %u*$in] [expr %u*$in]',flr,colLink,colLine,colLinknode_x,colLinknode_y);
                Frames(frameNum).Node(linkNodenum).X = colLinknode_x;
                Frames(frameNum).Node(linkNodenum).Y = colLinknode_y;
                
            end
            
        end
    end
end


% Beam Rigid Link Nodes

fprintf(fid,'\n\n#\tBeam Rigid Links');
fprintf(fid,'\n#\tnode "ABCD" X-coord Y-coord');
fprintf(fid,'\n#\t\tA: "Floor"');
fprintf(fid,'\n#\t\tB: 0');
fprintf(fid,'\n#\t\tC: "Bay"');
fprintf(fid,'\n#\t\tD: 1 = "Left" rigid link; "Right" rigid link\n');


for colLine = 1:Frames(frameNum).Bays
    
    for flr = 2:Frames(frameNum).Stories+1
        
        if flr == Frames(frameNum).Stories+1
            
            if colLine == 1 || colLine == Frames(frameNum).Bays
                beamLinklength1 = Frames(frameNum).Columns(flr-1).Ext.Width/2;
                beamLinklength2 = Frames(frameNum).Columns(flr-1).Int.Width/2;
            else
                beamLinklength1 = Frames(frameNum).Columns(flr-1).Int.Width/2;
                beamLinklength2 = Frames(frameNum).Columns(flr-1).Int.Width/2;
            end
            
        else
            
            if colLine == 1 || colLine == Frames(frameNum).Bays
                beamLinklength1 = max(Frames(frameNum).Columns(flr-1).Ext.Width/2,Frames(frameNum).Columns(flr).Ext.Width/2);
                beamLinklength2 = max(Frames(frameNum).Columns(flr-1).Int.Width/2,Frames(frameNum).Columns(flr).Int.Width/2);
            else
                beamLinklength1 = max(Frames(frameNum).Columns(flr-1).Int.Width/2,Frames(frameNum).Columns(flr).Int.Width/2);
                beamLinklength2 = max(Frames(frameNum).Columns(flr-1).Int.Width/2,Frames(frameNum).Columns(flr).Int.Width/2);
            end
            
        end
        
        for beamLink = 1:2
            nodeNum = str2double(sprintf('%u0%u0',flr,colLine));
            beamLinknode_y = Frames(frameNum).Node(nodeNum).Y;
            
            if beamLink == 1
                beamLinknode_x = Frames(frameNum).Node(nodeNum).X + beamLinklength1;
            else
                if colLine == 1 || colLine == Frames(frameNum).Bays
                    beamLinknode_x = Frames(frameNum).Node(nodeNum).X + Frames(frameNum).Beams(flr-1).Ext.Length*12 - beamLinklength2;
                else
                    beamLinknode_x = Frames(frameNum).Node(nodeNum).X + Frames(frameNum).Beams(flr-1).Int.Length*12 - beamLinklength2;
                end
                
            end
            
            linkNodenum = str2double(sprintf('%u0%u%u',flr,colLine,beamLink));
            fprintf(fid,'\nnode %u0%u%u [expr %u*$in] [expr %u*$in]',flr,colLine,beamLink,beamLinknode_x,beamLinknode_y);
            Frames(frameNum).Node(linkNodenum).X = beamLinknode_x;
            Frames(frameNum).Node(linkNodenum).Y = beamLinknode_y;
            
            
        end
    end
end


% Boundary Conditions

fprintf(fid,'\n\n# ------------------------------------------------------------------------------');
fprintf(fid,'\n# Boundary Conditions');
fprintf(fid,'\n#\tAssume 1st story columns are fixed at foundation \n');

for colLine = 1:Frames(frameNum).Bays+1
    
    nodeNum = str2double(sprintf('10%u0',colLine));
    
    fprintf(fid,'\nfix 10%u0 1 1 1',colLine);
    
end


% Uniaxial Materials

fprintf(fid,'\n\n# ------------------------------------------------------------------------------');
fprintf(fid,'\n# Uniaxial Materials');
fprintf(fid,'\n\n#\t uniaxialMaterial "Model" "ABCD"');
fprintf(fid,'\n#\t\t A: 1 = Steel 2 = Unconfined Concrete 3 = Confined Concrete (near ends) 4 = Confined Concrete (near midspan)');
fprintf(fid,'\n#\t\t B: 1 = Column 2 = Beam');
fprintf(fid,'\n#\t\t C: Story/Floor (Column/Beam)');
fprintf(fid,'\n#\t\t D: Column Line/Bay (Column/Beam)');

% Steel Models
fprintf(fid,'\n\n#\tSteel Models');

%   Columns
fprintf(fid,'\n#\t\tColumns\n');
% fprintf(fid,'\n# uniaxialMaterial Steel02 tag? fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?> siginit?> <-deter eps1? eps2? rf_min?>');
fprintf(fid,'\n# uniaxialMaterial NMBuckling tag? As? lb? fy?	 E?	   esh?	    delta_0? k? b? gamma? N1? N2? b_dam? c_dam? edam? numPts? maxItr?');

for story = 1:Frames(frameNum).Stories
    
    for colLine = 1:Frames(frameNum).Bays+1
        
        % Determine Rebar Layout
        
        if colLine == 1 || colLine == Frames(frameNum).Bays+1
            locType = 'Ext';
            fyi = Frames(frameNum).Columns(story).Ext.fy ;
            x = Frames(frameNum).Columns(story).Ext.Cover ;
            fc = Frames(frameNum).Columns(story).Ext.fc ;
            h = Frames(frameNum).Columns(story).Ext.Depth;
            b = Frames(frameNum).Columns(story).Ext.Width;
            cover = Frames(frameNum).Columns(story).Ext.Cover;
            rho = Frames(frameNum).Columns(story).Ext.rho_tot;
            rho_c = NaN;
            rho_sh = Frames(frameNum).Columns(story).Ext.rho_sh;
            s = Frames(frameNum).Columns(story).Ext.s;
            
        else
            locType = 'Int';
            fyi = Frames(frameNum).Columns(story).Int.fy ;
            x = Frames(frameNum).Columns(story).Int.Cover;
            fc = Frames(frameNum).Columns(story).Int.fc ;
            h = Frames(frameNum).Columns(story).Int.Depth;
            b = Frames(frameNum).Columns(story).Int.Width;
            cover = Frames(frameNum).Columns(story).Int.Cover;
            rho = Frames(frameNum).Columns(story).Int.rho_tot;
            rho_c = NaN;
            rho_sh = Frames(frameNum).Columns(story).Int.rho_sh;
            s = Frames(frameNum).Columns(story).Int.s;
        end
        
        [n,n_h,n_b,n_c, a_barSize, a_barArea, rho_actual,a_barSize_c, a_barArea_c, rho_actual_c, n_sh, a_barSize_sh, a_barArea_sh, a_rho_sh] = rebarLayout(h, b, cover, rho, rho_c, rho_sh,s,colLine,story,[], 'column');
        
        Frames(frameNum).Columns(story).(locType).Rebar.n = n;
        Frames(frameNum).Columns(story).(locType).Rebar.n_h = n_h;
        Frames(frameNum).Columns(story).(locType).Rebar.n_b = n_b;
        Frames(frameNum).Columns(story).(locType).Rebar.n_sh = n_sh;
        Frames(frameNum).Columns(story).(locType).Rebar.barSize = a_barSize;
        Frames(frameNum).Columns(story).(locType).Rebar.barArea = a_barArea;
        Frames(frameNum).Columns(story).(locType).Rebar.rho = rho_actual;
        Frames(frameNum).Columns(story).(locType).Rebar.barSize_sh = a_barSize_sh;
        Frames(frameNum).Columns(story).(locType).Rebar.rho_sh = a_rho_sh;
        
        % Longitudinal Rebar Properties
        As  = Frames(frameNum).Columns(story).(locType).Rebar.barArea;
        D   = Rebar(find([Rebar.Area] == As)).Diameter ;
        As  = As ;
        Ai  = As ;
        Di  = D  ;
        Esi = Esi0 ;
        fy  = fyi;
        Es  = Esi;
        
        Frames(frameNum).Columns(story).(locType).Rebar.As = As;
        fu = 1.25 * fy;
        lb = Frames(frameNum).Columns(story).(locType).s;
        
        
        esh         = 0.02;
        matID = sprintf('11%u%u',story,colLine);
        % uniaxialMaterial Steel02 tag? fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?> siginit?> <-deter eps1? eps2? rf_min?>
        % uniaxialMaterial NMBuckling tag? As? lb? fy?	 E?	   esh?	    delta_0? k? b? gamma? N1? N2? b_dam? c_dam? edam?
        if tog.Steel02
            fprintf(fid,'\nuniaxialMaterial Steel02 %s %0.2f %0.0f %0.6f %u %0.3f %0.3f -deter %0.3f %0.3f %0.3f',num2str(matID),fy,Es,bs,R0,cR1,cR2,eps1,eps2,rf_min);
        else
            fprintf(fid,'\nuniaxialMaterial NMBuckling %s %0.4f %0.4f %0.2f %0.0f %0.8f %0.9f %0.1f %0.1f %0.2f %0.2f %0.2f %0.2f %0.2f %0.8f %0.2f %0.2f %0.8f %u %u -corr %0.3e %0.3e %0.3e %0.2f %0.2f -fiber; ',num2str(matID),As,lb,fy,Es,esh,delta_0_m*lb,k_sh,b_sh,gamma,N1,N2,b_dam,c_dam,edam, b_dam_stff,c_dam_stff,edam_stff,numPts, maxItr, Cs, Dc, Ccr, fc, x);
        end
        
        % Transverse Rebar Properties
        As_sh = Rebar(find([Rebar.Size] == Frames(frameNum).Columns(story).(locType).Rebar.barSize_sh)).Area ;
        D_sh  = Rebar(find([Rebar.Size] == Frames(frameNum).Columns(story).(locType).Rebar.barSize_sh)).Diameter ;
        x_sh  = x - Di/2 - D_sh/2;
        Ai_sh = As_sh;
        Di_sh = D_sh;
        fy_sh = fyi;
        Es_sh = Esi;
        
    end
    
end

%   Beams
fprintf(fid,'\n\n#\t\tBeams\n');

for flr = 2:Frames(frameNum).Stories+1
    
    for bay = 1:Frames(frameNum).Bays
        
        % Determine Rebar Layout
        
        if bay == 1 || bay == Frames(frameNum).Bays
            locType = 'Ext';
            fc = Frames(frameNum).Beams(flr-1).Ext.fc ;
            h = Frames(frameNum).Beams(flr-1).Ext.Depth;
            b = Frames(frameNum).Beams(flr-1).Ext.Width;
            cover = Frames(frameNum).Beams(flr-1).Ext.Cover;
            fyi = Frames(frameNum).Beams(flr-1).Ext.fy ;
            x = Frames(frameNum).Beams(flr-1).Ext.Cover ;
            rho = Frames(frameNum).Beams(flr-1).Ext.rho;
            rho_c = Frames(frameNum).Beams(flr-1).Ext.rho_c;
            rho_sh = Frames(frameNum).Beams(flr-1).Ext.rho_sh;
            s = Frames(frameNum).Beams(flr-1).Ext.s;
        else
            locType = 'Int';
            fc = Frames(frameNum).Beams(flr-1).Int.fc ;
            h = Frames(frameNum).Beams(flr-1).Int.Depth;
            b = Frames(frameNum).Beams(flr-1).Int.Width;
            cover = Frames(frameNum).Beams(flr-1).Int.Cover;
            fyi = Frames(frameNum).Beams(flr-1).Int.fy ;
            x = Frames(frameNum).Beams(flr-1).Int.Cover ;
            rho = Frames(frameNum).Beams(flr-1).Int.rho;
            rho_c = Frames(frameNum).Beams(flr-1).Int.rho_c;
            rho_sh = Frames(frameNum).Beams(flr-1).Int.rho_sh;
            s = Frames(frameNum).Beams(flr-1).Int.s;
        end
        
        [n,n_h,n_b,n_c, a_barSize, a_barArea, rho_actual,a_barSize_c, a_barArea_c, rho_actual_c, n_sh, a_barSize_sh, a_barArea_sh, a_rho_sh] = rebarLayout(h, b, cover, rho, rho_c, rho_sh,s,bay,flr,[], 'beam');
        
        Frames(frameNum).Beams(flr-1).(locType).Rebar.n_b = n_b;
        Frames(frameNum).Beams(flr-1).(locType).Rebar.n_c = n_c;
        Frames(frameNum).Beams(flr-1).(locType).Rebar.n_sh = n_sh;
        Frames(frameNum).Beams(flr-1).(locType).Rebar.barSize = a_barSize;
        Frames(frameNum).Beams(flr-1).(locType).Rebar.barArea = a_barArea;
        Frames(frameNum).Beams(flr-1).(locType).Rebar.rho = rho_actual;
        Frames(frameNum).Beams(flr-1).(locType).Rebar.barSize_c = a_barSize_c;
        Frames(frameNum).Beams(flr-1).(locType).Rebar.barArea_c = a_barArea_c;
        Frames(frameNum).Beams(flr-1).(locType).Rebar.rho_c = rho_actual_c;
        Frames(frameNum).Beams(flr-1).(locType).Rebar.barSize_sh = a_barSize_sh;
        Frames(frameNum).Beams(flr-1).(locType).Rebar.rho_sh = a_rho_sh;
        
        % Longitudinal Rebar Properties
        As  = Frames(frameNum).Beams(flr-1).(locType).Rebar.barArea;
        As_c = Frames(frameNum).Beams(flr-1).(locType).Rebar.barArea_c;
        Di   = Rebar(find([Rebar.Area] == As)).Diameter ;
        Di_c   = Rebar(find([Rebar.Area] == As_c)).Diameter ;
        Ai  = As;
        Ai_c  = As_c ;
        Esi = Esi0 ;
        fy = fyi;
        fy_c = fyi;
        Es = Esi;
        Es_c = Esi;
        D = Di;
        D_c = Di_c;
        
        Frames(frameNum).Beams(flr-1).(locType).Rebar.As = As;
        Frames(frameNum).Beams(flr-1).(locType).Rebar.As_c = As_c;
        lb = Frames(frameNum).Beams(flr-1).(locType).s;
        
        fu = 1.25 * fy;
        fu_c = 1.25 * fy_c;
        esh         = 0.02;
        
        matID = sprintf('12%u%u',flr,bay);
        matID_c = sprintf('52%u%u',flr,bay);
        matID_s = sprintf('62%u%u',flr,bay);
        % uniaxialMaterial Steel02 tag? fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?> siginit?> <-deter eps1? eps2? rf_min?>
        if tog.Steel02
            fprintf(fid,'\nuniaxialMaterial Steel02 %s %0.2f %0.0f %0.6f %u %0.3f %0.3f -deter %0.3f %0.3f %0.3f',num2str(matID),fy,Es,bs,R0,cR1,cR2,eps1,eps2,rf_min);
            fprintf(fid,'\nuniaxialMaterial Steel02 %s %0.2f %0.0f %0.6f %u %0.3f %0.3f -deter %0.3f %0.3f %0.3f',num2str(matID_c),fy_c,Es_c,bs_c,R0,cR1,cR2,eps1,eps2,rf_min);
        else
            fprintf(fid,'\nuniaxialMaterial NMBuckling %s %0.4f %0.4f %0.2f %0.0f %0.8f %0.9f %0.1f %0.1f %0.2f %0.2f %0.2f %0.2f %0.2f %0.8f %0.2f %0.2f %0.8f %u %u -corr %0.3e %0.3e %0.3e %0.2f %0.2f -fiber;',num2str(matID),As,lb,fy,Es,esh,delta_0_m*lb,k_sh,b_sh,gamma,N1,N2,b_dam,c_dam,edam,b_dam_stff,c_dam_stff,edam_stff,numPts, maxItr, Cs, Dc, Ccr, fc, x);
            fprintf(fid,'\nuniaxialMaterial NMBuckling %s %0.4f %0.4f %0.2f %0.0f %0.8f %0.9f %0.1f %0.1f %0.2f %0.2f %0.2f %0.2f %0.2f %0.8f %0.2f %0.2f %0.8f %u %u -corr %0.3e %0.3e %0.3e %0.2f %0.2f -fiber;',num2str(matID_c),As_c,lb,fy_c,Es_c,esh,delta_0_m*lb,k_sh,b_sh,gamma,N1,N2,b_dam,c_dam,edam,b_dam_stff,c_dam_stff,edam_stff,numPts, maxItr, Cs, Dc, Ccr, fc, x);
        end
        
        % Slab
        %         fprintf(fid,'\nuniaxialMaterial Steel02 %s %0.2f %0.0f %0.6f %u %0.3f %0.3f -deter %0.3f %0.3f %0.3f',num2str(matID_s),fy,Es,bs,R0,cR1,cR2,eps1,eps2,rf_min);
        fprintf(fid,'\nuniaxialMaterial NMBuckling %s %0.4f %0.4f %0.2f %0.0f %0.8f %0.9f %0.1f %0.1f %0.2f %0.2f %0.2f %0.2f %0.2f %0.8f %0.2f %0.2f %0.8f %u %u;',num2str(matID_s),As_slab,5,fyi,Esi,esh,0,k_sh,b_sh,gamma,N1,N2,b_dam,c_dam,edam,b_dam_stff,c_dam_stff,edam_stff,numPts, maxItr);
        
        % Transverse Rebar Properties
        As_sh = Rebar(find([Rebar.Size] == Frames(frameNum).Beams(flr-1).(locType).Rebar.barSize_sh)).Area ;
        D_sh  = Rebar(find([Rebar.Size] == Frames(frameNum).Beams(flr-1).(locType).Rebar.barSize_sh)).Diameter ;
        x_sh  = x - ((Di+Di_c)/2)/2 - D_sh/2;
        Ai_sh = As_sh;
        Di_sh = D_sh;
        
    end
    
end


% Concrete Models
fprintf(fid,'\n\n#\tConcrete Models');

%   Columns
fprintf(fid,'\n#\t\tColumns\n');
fprintf(fid,'\n# uniaxialMaterial ConfinedConcrete04 tag? fpc? h? b? s? cover? -column n_h? n_b? A_long? A_sh? fy_sh? Xp? Xn? < -corr? Cs? Dc? Ccr? < -time? t? > -loc?');
% # 									tag? fpc? h? b? s? cover? -column n_h? n_b? A_long? A_sh? fy_sh? Xp? Xn? < -corr? Cs? Dc? Ccr? < -time? t? > -loc?
% uniaxialMaterial ConfinedConcrete04 1 5.0 22.0 22.0 5.0 2.5 -column 4 4 0.44 0.20 60.0 2.0 30.0 -corr 2.95E-9 1.29E2 0.90E-9 -end;

for story = 1:Frames(frameNum).Stories
    
    for colLine = 1:Frames(frameNum).Bays+1
        
        if colLine == 1 || colLine == Frames(frameNum).Bays+1
            locType = 'Ext';
            fc = Frames(frameNum).Columns(story).Ext.fc;
            h = Frames(frameNum).Columns(story).Ext.Depth;
            b = Frames(frameNum).Columns(story).Ext.Width;
            cover = Frames(frameNum).Columns(story).Ext.Cover;
            fy = Frames(frameNum).Columns(story).Ext.fy;
            rho = Frames(frameNum).Columns(story).Ext.rho_tot;
            rho_c = NaN;
            rho_sh = Frames(frameNum).Columns(story).Ext.rho_sh;
            s = Frames(frameNum).Columns(story).Ext.s;
            
        else
            locType = 'Int';
            fc = Frames(frameNum).Columns(story).Int.fc;
            h = Frames(frameNum).Columns(story).Int.Depth;
            b = Frames(frameNum).Columns(story).Int.Width;
            cover = Frames(frameNum).Columns(story).Int.Cover;
            fy = Frames(frameNum).Columns(story).Int.fy;
            rho = Frames(frameNum).Columns(story).Int.rho_tot;
            rho_c = NaN;
            rho_sh = Frames(frameNum).Columns(story).Int.rho_sh;
            s = Frames(frameNum).Columns(story).Int.s;
        end
        
        % Confined (near end)
        matID = sprintf('31%u%u',story,colLine);
        
        [n,n_h,n_b,n_c, a_barSize, a_barArea, rho_actual,a_barSize_c, a_barArea_c, rho_actual_c, n_sh, a_barSize_sh, a_barArea_sh, a_rho_sh] = rebarLayout(h, b, cover, rho, rho_c, rho_sh,s,colLine,story,[], 'column');
        
        %         [fcc_e, eps_cc_e, eps_cu_e, ~, ~, Ec_e] = C_Concrete04(fc,h,b,s,cover,n_h,n_b,n_c,a_barSize,a_barSize_c,a_barSize_sh,fy_sh,2,30,'end','column',corrRebar);
        %         tag? fpc? h? b? s? cover? -column n_h? n_b? A_long? A_sh? fy_sh? Xp? Xn? < -corr? Cs? Dc? Ccr? < -time? t? > -loc?
        fprintf(fid,'\nuniaxialMaterial ConfinedConcrete04 %s %.2f %.2f %.2f %.2f %.2f -res %0.2f %0.2f %s %u %u %.3f %.3f %.2f %u %u -corr %.3e %.3e %.3e %s',num2str(matID),fc,h,b,s,cover,0.1,2,'-column',n_h,n_b,a_barArea,a_barArea_sh,fy_sh,2,30,Cs,Dc,Ccr,'-end');%,ft,eps_tu);
        %         fprintf(fid,'\nuniaxialMaterial Concrete04 %s -6 %.6f %.6f %.0f %.4f %.7f',num2str(matID),-eps_cc,-eps_cu,Ec);%,ft,eps_tu);
        
        % Confined (near mid)
        matID = sprintf('41%u%u',story,colLine);
        
        %         fprintf(fid,'\nuniaxialMaterial Concrete04 %s %.2f %.6f %.6f %.0f %.4f %.7f',num2str(matID),-fcc_m,-eps_cc_m,-eps_cu_m,Ec_m);%,ft,eps_tu);
        fprintf(fid,'\nuniaxialMaterial ConfinedConcrete04 %s %.2f %.2f %.2f %.2f %.2f -res %0.2f %0.2f %s %u %u %.3f %.3f %.2f %u %u -corr %.3e %.3e %.3e %s',num2str(matID),fc,h,b,s,cover,0.1, 2, '-column',n_h,n_b,a_barArea,a_barArea_sh,fy_sh,2,30,Cs,Dc,Ccr,'-mid');%,ft,eps_tu);
        
        % Unconfined
        matID = sprintf('21%u%u',story,colLine);
        
        % uniaxialMaterial Concrete04 tag? fpc? epsc0? epscu? Ec0? <ft? etu? <beta?>> <-res rf? epscf?>
        
        [eps_c0,eps_cu,~,~,Ec] = UC_Concrete04(fc,2,2.3);
        
        %         if ageVal == 1
        %             Di   = Rebar(find([Rebar.Size] == a_barSize)).Diameter;
        %             D = corrRebar.D;
        %             X = (Di - D)/2;
        %             [fc_star_h] = CoverCorrosion(fc, k, v_rs, eps_c0, n_h, X, h );
        %             [fc_star_b] = CoverCorrosion(fc, k, v_rs, eps_c0, n_b, X, b );
        %             fc_star = min(fc_star_b,fc_star_h);
        %             fprintf(fid,'\nuniaxialMaterial Concrete04 %s %.2f %.6f %.6f %.0f %.4f %.7f',num2str(matID),-fc_star,-eps_c0,-eps_cu,Ec);%,ft,eps_tu);
        %             Frames(frameNum).Columns(story).(locType).Rebar.corr.(['yr_' num2str(t)]).fc_star = -fc_star;
        %             Frames(frameNum).Columns(story).(locType).Rebar.corr.(['yr_' num2str(t)]).eps_c0 = -eps_c0;
        %             Frames(frameNum).Columns(story).(locType).Rebar.corr.(['yr_' num2str(t)]).eps_cu = -eps_cu;
        %             Frames(frameNum).Columns(story).(locType).Rebar.corr.(['yr_' num2str(t)]).Ec = Ec;
        %
        %             Frames(frameNum).Columns(story).(locType).Rebar.corr.(['yr_' num2str(t)]).fcc_e = -fcc_e;
        %             Frames(frameNum).Columns(story).(locType).Rebar.corr.(['yr_' num2str(t)]).eps_cc_e = -eps_cc_e;
        %             Frames(frameNum).Columns(story).(locType).Rebar.corr.(['yr_' num2str(t)]).eps_cu_e = -eps_cu_e;
        %             Frames(frameNum).Columns(story).(locType).Rebar.corr.(['yr_' num2str(t)]).Ec_e = Ec_e;
        %
        %             Frames(frameNum).Columns(story).(locType).Rebar.corr.(['yr_' num2str(t)]).fcc_m = -fcc_m;
        %             Frames(frameNum).Columns(story).(locType).Rebar.corr.(['yr_' num2str(t)]).eps_cc_m = -eps_cc_m;
        %             Frames(frameNum).Columns(story).(locType).Rebar.corr.(['yr_' num2str(t)]).eps_cu_m = -eps_cu_m;
        %             Frames(frameNum).Columns(story).(locType).Rebar.corr.(['yr_' num2str(t)]).Ec_m = Ec_m;
        %         else
        fprintf(fid,'\nuniaxialMaterial Concrete04 %s %.2f %.6f %.6f %.0f %.4f %.7f',num2str(matID),-fc,-eps_c0,-eps_cu,Ec);%,ft,eps_tu);
        %         end
    end
    
end

%   Beams
fprintf(fid,'\n\n#\t\tBeams\n');

for flr = 2:Frames(frameNum).Stories+1
    
    for bay = 1:Frames(frameNum).Bays
        
        if bay == 1 || bay == Frames(frameNum).Bays
            locType = 'Ext';
            fc = Frames(frameNum).Beams(flr-1).Ext.fc;
            h = Frames(frameNum).Beams(flr-1).Ext.Depth;
            b = Frames(frameNum).Beams(flr-1).Ext.Width;
            cover = Frames(frameNum).Beams(flr-1).Ext.Cover;
            fy = Frames(frameNum).Beams(flr-1).Ext.fy;
            rho = Frames(frameNum).Beams(flr-1).Ext.rho;
            rho_c = Frames(frameNum).Beams(flr-1).Ext.rho_c;
            rho_sh = Frames(frameNum).Beams(flr-1).Ext.rho_sh;
            s = Frames(frameNum).Beams(flr-1).Ext.s;
        else
            locType = 'Int';
            fc = Frames(frameNum).Beams(flr-1).Int.fc;
            h = Frames(frameNum).Beams(flr-1).Int.Depth;
            b = Frames(frameNum).Beams(flr-1).Int.Width;
            cover = Frames(frameNum).Beams(flr-1).Int.Cover;
            fy = Frames(frameNum).Beams(flr-1).Int.fy;
            rho = Frames(frameNum).Beams(flr-1).Int.rho;
            rho_c = Frames(frameNum).Beams(flr-1).Int.rho_c;
            rho_sh = Frames(frameNum).Beams(flr-1).Int.rho_sh;
            s = Frames(frameNum).Beams(flr-1).Int.s;
        end
        
        fy_sh = fy;
        
        
        % Confined (near end)
        matID = sprintf('32%u%u',flr,bay);
        
        [n,n_h,n_b,n_c, a_barSize, a_barArea, rho_actual,a_barSize_c, a_barArea_c, rho_actual_c, n_sh, a_barSize_sh,a_barArea_sh, a_rho_sh] = rebarLayout(h, b, cover, rho, rho_c, rho_sh,s,bay,flr,[], 'beam');
        
%         [fcc_e, eps_cc_e, eps_cu_e, ~, ~, Ec_e] = C_Concrete04(fc,h,b,s,cover,n_h,n_b,n_c,a_barSize,a_barSize_c,a_barSize_sh,fy_sh,2,30,'end','beam',corrRebar);
        
%         fprintf(fid,'\nuniaxialMaterial Concrete04 %s %.2f %.6f %.6f %.0f %.4f %.7f',num2str(matID),-fcc_e,-eps_cc_e,-eps_cu_e,Ec_e);%,ft,eps_tu);
        fprintf(fid,'\nuniaxialMaterial ConfinedConcrete04 %s %.2f %.2f %.2f %.2f %.2f -res %0.2f %0.2f %s %u %u %u %.3f %.3f %.3f %.2f %u %u -corr %.3e %.3e %.3e %s',num2str(matID),fc,h,b,s,cover,0.1,2,'-beam',0,n_b,n_c,a_barArea,a_barArea_c,a_barArea_sh,fy_sh,2,30,Cs,Dc,Ccr,'-end');%,ft,eps_tu);
        
        % Confined (near mid)
        matID = sprintf('42%u%u',flr,bay);
        
%         [fcc_m, eps_cc_m, eps_cu_m, ~, ~, Ec_m] = C_Concrete04(fc,h,b,s,cover,n_h,n_b,n_c,a_barSize,a_barSize_c,a_barSize_sh,fy_sh,2,30,'mid','beam',corrRebar);
        
%         fprintf(fid,'\nuniaxialMaterial Concrete04 %s %.2f %.6f %.6f %.0f %.4f %.7f',num2str(matID),-fcc_m,-eps_cc_m,-eps_cu_m,Ec_m);%,ft,eps_tu);
        fprintf(fid,'\nuniaxialMaterial ConfinedConcrete04 %s %.2f %.2f %.2f %.2f %.2f -res %0.2f %0.2f %s %u %u %u %.3f %.3f %.3f %.2f %u %u -corr %.3e %.3e %.3e %s',num2str(matID),fc,h,b,s,cover,0.1,2,'-beam',0,n_b,n_c,a_barArea,a_barArea_c,a_barArea_sh,fy_sh,2,30,Cs,Dc,Ccr,'-mid');%,ft,eps_tu);
        
        % Unconfined
        matID = sprintf('22%u%u',flr,bay);
        [eps_c0,eps_cu,~,~,Ec] = UC_Concrete04(fc,2,2.3);
        
%         if ageVal == 1
%             Di   = Rebar(find([Rebar.Size] == a_barSize)).Diameter;
%             D = corrRebar.D;
%             X = (Di - D)/2;
%             
%             Di_c   = Rebar(find([Rebar.Size] == a_barSize_c)).Diameter;
%             D_c = corrRebar.D_c;
%             X_c = (Di_c - D_c)/2 ;
%             
%             [fc_star] = CoverCorrosion(fc, k, v_rs, eps_c0, n_b, X, b );
%             [fc_star_c] = CoverCorrosion(fc, k, v_rs, eps_c0, n_c, X_c, b );
%             fc_star = min(fc_star,fc_star_c);
%             
%             fprintf(fid,'\nuniaxialMaterial Concrete04 %s %.2f %.6f %.6f %.0f %.4f %.7f',num2str(matID),-fc_star,-eps_c0,-eps_cu,Ec);%,ft,eps_tu);
%             
%             Frames(frameNum).Beams(flr-1).(locType).Rebar.corr.(['yr_' num2str(t)]).fc_star = -fc_star;
%             Frames(frameNum).Beams(flr-1).(locType).Rebar.corr.(['yr_' num2str(t)]).eps_c0 = -eps_c0;
%             Frames(frameNum).Beams(flr-1).(locType).Rebar.corr.(['yr_' num2str(t)]).eps_cu = -eps_cu;
%             Frames(frameNum).Beams(flr-1).(locType).Rebar.corr.(['yr_' num2str(t)]).Ec = Ec;
%             
%             Frames(frameNum).Beams(flr-1).(locType).Rebar.corr.(['yr_' num2str(t)]).fcc_e = -fcc_e;
%             Frames(frameNum).Beams(flr-1).(locType).Rebar.corr.(['yr_' num2str(t)]).eps_cc_e = -eps_cc_e;
%             Frames(frameNum).Beams(flr-1).(locType).Rebar.corr.(['yr_' num2str(t)]).eps_cu_e = -eps_cu_e;
%             Frames(frameNum).Beams(flr-1).(locType).Rebar.corr.(['yr_' num2str(t)]).Ec_e = Ec_e;
%             
%             Frames(frameNum).Beams(flr-1).(locType).Rebar.corr.(['yr_' num2str(t)]).fcc_m = -fcc_m;
%             Frames(frameNum).Beams(flr-1).(locType).Rebar.corr.(['yr_' num2str(t)]).eps_cc_m = -eps_cc_m;
%             Frames(frameNum).Beams(flr-1).(locType).Rebar.corr.(['yr_' num2str(t)]).eps_cu_m = -eps_cu_m;
%             Frames(frameNum).Beams(flr-1).(locType).Rebar.corr.(['yr_' num2str(t)]).Ec_m = Ec_m;
%         else
            fprintf(fid,'\nuniaxialMaterial Concrete04 %s %.2f %.6f %.6f %.0f %.4f %.7f', num2str(matID),-fc,-eps_c0,-eps_cu,Ec);%,ft,eps_tu);
%         end
        
    end
    
end

%   Rigid Link Material
fprintf(fid,'\n\n#\t\tRigid Link Material\n');
fprintf(fid,'\nuniaxialMaterial Elastic %.0f %0.0f',1,999999);

% Sections

fprintf(fid,'\n\n# ------------------------------------------------------------------------------');
fprintf(fid,'\n# Sections');
fprintf(fid,'\n#\tsection "ABCD"');
fprintf(fid,'\n#\t\tA: 1 = End Section; 2 = Interior Section');
fprintf(fid,'\n#\t\tB: 1 = Column; 2 = Beam');
fprintf(fid,'\n#\t\tC: Story/Floor');
fprintf(fid,'\n#\t\tD: Column/Bay\n');
fprintf(fid,'\n#\tColumns\n');

for story = 1:Frames(frameNum).Stories
    
    for colLine = 1:Frames(frameNum).Bays+1
        
        if colLine == 1 || colLine == Frames(frameNum).Bays+1
            locType = 'Ext';
            h = Frames(frameNum).Columns(story).Ext.Depth;
            b = Frames(frameNum).Columns(story).Ext.Width;
            cover = Frames(frameNum).Columns(story).Ext.Cover;
            
        else
            locType = 'Int';
            h = Frames(frameNum).Columns(story).Int.Depth;
            b = Frames(frameNum).Columns(story).Int.Width;
            cover = Frames(frameNum).Columns(story).Int.Cover;
            
        end
        
        sectID  = sprintf('11%u%u',story,colLine);
        sectIDmid  = sprintf('21%u%u',story,colLine);
        SmatID  = sprintf('11%u%u',story,colLine);
        UCmatID = sprintf('21%u%u',story,colLine);
        CEmatID = sprintf('31%u%u',story,colLine);
        CMmatID = sprintf('41%u%u',story,colLine);
        
        n_h =  Frames(frameNum).Columns(story).(locType).Rebar.n_h;
        n_b =  Frames(frameNum).Columns(story).(locType).Rebar.n_b;
        As  =  Frames(frameNum).Columns(story).(locType).Rebar.As;
        
        x_core = (b/2) - cover;
        y_core = (h/2) - cover;
        
        %         numCov_bt_y = floor(cover/1);
        numCov_bt_y = 10;
        numCov_bt_x = 10;
        %         numCov_bt_x = floor(b/1);
        
        %         numCov_side_y = floor((h - 2*cover)/1);
        numCov_side_y = 30;
        numCov_side_x = 10;
        %         numCov_side_x = floor(cover/0.5);
        
        %         numCore_y = floor((h - 2*cover)/1);
        numCore_y = 30;
        numCore_x = 30;
        %         numCore_x = floor((b - 2*cover)/1);
        %         numCore_x = floor((b - 2*cover)/1);
        
        sy = y_core*2 / (n_h - 1);
        
        % End Section
        fprintf(fid,'\n\n#\tStory %u Column Line %u',story,colLine);
        fprintf(fid,'\nsection\tFiber %s\t{ ;', num2str(sectID));
        
        % Core Concrete
        fprintf(fid,'\n#\t\tCore');
        fprintf(fid,'\n\tpatch rect %s %.0f %.0f %.2f %.2f %.2f %.2f ;',num2str(CEmatID),numCore_x,numCore_y,-y_core,-x_core,y_core,x_core);
        
        % Cover Concrete
        fprintf(fid,'\n#\t\tCover');
        fprintf(fid,'\n\tpatch rect %s %.0f %.0f %.2f %.2f %.2f %.2f ;\t# Left',num2str(UCmatID),numCov_side_x,numCov_side_y,-y_core,-b/2,y_core,-x_core);
        fprintf(fid,'\n\tpatch rect %s %.0f %.0f %.2f %.2f %.2f %.2f ;\t# Right',num2str(UCmatID),numCov_side_x,numCov_side_y,-y_core,x_core,y_core,b/2);
        fprintf(fid,'\n\tpatch rect %s %.0f %.0f %.2f %.2f %.2f %.2f ;\t# Top',num2str(UCmatID),numCov_bt_x,numCov_bt_y,y_core,-b/2,h/2,b/2);
        fprintf(fid,'\n\tpatch rect %s %.0f %.0f %.2f %.2f %.2f %.2f ;\t# Bot',num2str(UCmatID),numCov_bt_x,numCov_bt_y,-h/2,-b/2,-y_core, b/2);
        
        % Reinforcing Steel
        fprintf(fid,'\n#\t\tSteel');
        
        y = -y_core;
        
        for layST = 1:n_h
            
            if layST == 1
                fprintf(fid,'\n\tlayer straight %s %.0f %.3f %.2f %.2f %.2f %.2f ;', num2str(SmatID),n_b,As,y,-x_core,y,x_core);
            elseif layST == n_h
                y = y + sy;
                fprintf(fid,'\n\tlayer straight %s %.0f %.3f %.2f %.2f %.2f %.2f ;', num2str(SmatID),n_b,As,y,-x_core,y,x_core);
            else
                y = y + sy;
                fprintf(fid,'\n\tlayer straight %s %.0f %.3f %.2f %.2f %.2f %.2f ;', num2str(SmatID),2,As,y,-x_core,y,x_core);
            end
        end
        
        fprintf(fid,'\n };');
        
        % Mid Section
        
        fprintf(fid,'\nsection\tFiber %s\t{ ;', num2str(sectIDmid));
        
        % Core Concrete
        fprintf(fid,'\n#\t\tCore');
        fprintf(fid,'\n\tpatch rect %s %.0f %.0f %.2f %.2f %.2f %.2f ;',num2str(CMmatID),numCore_x,numCore_y,-y_core,-x_core,y_core,x_core);
        
        % Cover Concrete
        fprintf(fid,'\n#\t\tCover');
        fprintf(fid,'\n\tpatch rect %s %.0f %.0f %.2f %.2f %.2f %.2f ;\t# Left',num2str(UCmatID),numCov_side_x,numCov_side_y,-y_core,-b/2,y_core,-x_core);
        fprintf(fid,'\n\tpatch rect %s %.0f %.0f %.2f %.2f %.2f %.2f ;\t# Right',num2str(UCmatID),numCov_side_x,numCov_side_y,-y_core,x_core,y_core,b/2);
        fprintf(fid,'\n\tpatch rect %s %.0f %.0f %.2f %.2f %.2f %.2f ;\t# Top',num2str(UCmatID),numCov_bt_x,numCov_bt_y,y_core,-b/2,h/2,b/2);
        fprintf(fid,'\n\tpatch rect %s %.0f %.0f %.2f %.2f %.2f %.2f ;\t# Bot',num2str(UCmatID),numCov_bt_x,numCov_bt_y,-h/2,-b/2,-y_core, b/2);
        
        % Reinforcing Steel
        fprintf(fid,'\n#\t\tSteel');
        
        y = -y_core;
        
        for layST = 1:n_h
            
            if layST == 1
                fprintf(fid,'\n\tlayer straight %s %.0f %.3f %.2f %.2f %.2f %.2f ;', num2str(SmatID),n_b,As,y,-x_core,y,x_core);
            elseif layST == n_h
                y = y + sy;
                fprintf(fid,'\n\tlayer straight %s %.0f %.3f %.2f %.2f %.2f %.2f ;', num2str(SmatID),n_b,As,y,-x_core,y,x_core);
            else
                y = y + sy;
                fprintf(fid,'\n\tlayer straight %s %.0f %.3f %.2f %.2f %.2f %.2f ;', num2str(SmatID),2,As,y,-x_core,y,x_core);
            end
        end
        
        fprintf(fid,'\n };');
        
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
        
        sectID  = sprintf('12%u%u',flr,bay);
        sectIDmid  = sprintf('22%u%u',flr,bay);
        SmatID  = sprintf('12%u%u',flr,bay);
        SmatID_c = sprintf('52%u%u',flr,bay);
        SmatID_s = sprintf('62%u%u',flr,bay);
        UCmatID = sprintf('22%u%u',flr,bay);
        CEmatID = sprintf('32%u%u',flr,bay);
        CMmatID = sprintf('42%u%u',flr,bay);
        
        n_c =  Frames(frameNum).Beams(flr-1).(locType).Rebar.n_c;
        n_b =  Frames(frameNum).Beams(flr-1).(locType).Rebar.n_b;
        As  =  Frames(frameNum).Beams(flr-1).(locType).Rebar.As;
        As_c  =  Frames(frameNum).Beams(flr-1).(locType).Rebar.barArea_c;
        
        
        x_core = (b/2) - cover;
        y_core = (h/2) - cover;
        
        %         numCov_bt_y = floor(cover/1);
        numCov_bt_y = 10;
        numCov_bt_x = 10;
        %         numCov_bt_x = floor(b/1);
        
        %         numCov_side_y = floor((h - 2*cover)/1);
        numCov_side_y = 30;
        numCov_side_x = 10;
        %         numCov_side_x = floor(cover/0.5);
        
        %         numCore_y = floor((h - 2*cover)/1);
        numCore_y = 30;
        numCore_x = 30;
        %         numCore_x = floor((b - 2*cover)/1);
        %         numCore_x = floor((b - 2*cover)/1);
        
        %         numSlab_y = floor(tf/1);
        numSlab_y = 20;
        numSlab_x = 10;
        
        % End Section
        
        fprintf(fid,'\n\n#\tFloor %u Bay %u',flr,bay);
        fprintf(fid,'\nsection\tFiber %s\t{ ;', num2str(sectID));
        
        % Core Concrete
        fprintf(fid,'\n#\t\tCore');
        fprintf(fid,'\n\tpatch rect %s %.0f %.0f %.2f %.2f %.2f %.2f ;',num2str(CEmatID),numCore_x,numCore_y,-y_core,-x_core,y_core,x_core);
        
        % Cover Concrete
        fprintf(fid,'\n#\t\tCover');
        fprintf(fid,'\n\tpatch rect %s %.0f %.0f %.2f %.2f %.2f %.2f ;\t# Left',num2str(UCmatID),numCov_side_x,numCov_side_y,-y_core,-b/2,y_core,-x_core);
        fprintf(fid,'\n\tpatch rect %s %.0f %.0f %.2f %.2f %.2f %.2f ;\t# Right',num2str(UCmatID),numCov_side_x,numCov_side_y,-y_core,x_core,y_core,b/2);
        fprintf(fid,'\n\tpatch rect %s %.0f %.0f %.2f %.2f %.2f %.2f ;\t# Top',num2str(UCmatID),numCov_bt_x,numCov_bt_y,y_core,-b/2,h/2,b/2);
        fprintf(fid,'\n\tpatch rect %s %.0f %.0f %.2f %.2f %.2f %.2f ;\t# Bot',num2str(UCmatID),numCov_bt_x,numCov_bt_y,-h/2,-b/2,-y_core,b/2);
        
        % Reinforcing Steel
        fprintf(fid,'\n#\t\tSteel');
        
        fprintf(fid,'\n\tlayer straight %s %.0f %.3f %.2f %.2f %.2f %.2f ;', num2str(SmatID),n_b,As,-y_core,-x_core,-y_core,x_core);
        fprintf(fid,'\n\tlayer straight %s %.0f %.3f %.2f %.2f %.2f %.2f ;', num2str(SmatID_c),n_c,As_c,y_core,-x_core,y_core,x_core);
        
        % Slab Sections
        fprintf(fid,'\n#\t\tSlab');
        fprintf(fid,'\n\tpatch rect %s %.0f %.0f %.2f %.2f %.2f %.2f ;\t\t# Left',num2str(UCmatID),numSlab_x,numSlab_y, (h/2 - tf),-(b_eff-b)/2,h/2,-b/2);
        fprintf(fid,'\n\tpatch rect %s %.0f %.0f %.2f %.2f %.2f %.2f ;\t\t\t# Right',num2str(UCmatID),numSlab_x,numSlab_y, (h/2 - tf),b/2,h/2,(b_eff-b)/2);
        fprintf(fid,'\n\tlayer straight %s %.0f %.3f %.2f %.2f %.2f %.2f ;', num2str(SmatID_s),n_slab,As_slab,(h/2)-d_slab,-(b/2 + 3),(h/2)-d_slab,-((b_eff-b)/2 - 3));
        fprintf(fid,'\n\tlayer straight %s %.0f %.3f %.2f %.2f %.2f %.2f ;', num2str(SmatID_s),n_slab,As_slab,(h/2)-d_slab,(b/2 + 3),(h/2)-d_slab,((b_eff-b)/2 - 3));
        
        fprintf(fid,'\n };');
        
        % Mid Section
        
        fprintf(fid,'\nsection\tFiber %s\t{ ;', num2str(sectIDmid));
        
        % Core Concrete
        fprintf(fid,'\n#\t\tCore');
        fprintf(fid,'\n\tpatch rect %s %.0f %.0f %.2f %.2f %.2f %.2f ;',num2str(CMmatID),numCore_x,numCore_y,-y_core,-x_core,y_core,x_core);
        
        % Cover Concrete
        fprintf(fid,'\n#\t\tCover');
        fprintf(fid,'\n\tpatch rect %s %.0f %.0f %.2f %.2f %.2f %.2f ;\t# Left',num2str(UCmatID),numCov_side_x,numCov_side_y,-y_core,-b/2,y_core,-x_core);
        fprintf(fid,'\n\tpatch rect %s %.0f %.0f %.2f %.2f %.2f %.2f ;\t# Right',num2str(UCmatID),numCov_side_x,numCov_side_y,-y_core,x_core,y_core,b/2);
        fprintf(fid,'\n\tpatch rect %s %.0f %.0f %.2f %.2f %.2f %.2f ;\t# Top',num2str(UCmatID),numCov_bt_x,numCov_bt_y,y_core,-b/2,h/2,b/2);
        fprintf(fid,'\n\tpatch rect %s %.0f %.0f %.2f %.2f %.2f %.2f ;\t# Bot',num2str(UCmatID),numCov_bt_x,numCov_bt_y,-h/2,-b/2,-y_core,b/2);
        
        % Reinforcing Steel
        fprintf(fid,'\n#\t\tSteel');
        
        fprintf(fid,'\n\tlayer straight %s %.0f %.3f %.2f %.2f %.2f %.2f ;', num2str(SmatID),n_b,As,-y_core,-x_core,-y_core,x_core);
        fprintf(fid,'\n\tlayer straight %s %.0f %.3f %.2f %.2f %.2f %.2f ;', num2str(SmatID_c),n_c,As_c,y_core,-x_core,y_core,x_core);
        
        % Slab Sections
        fprintf(fid,'\n#\t\tSlab');
        fprintf(fid,'\n\tpatch rect %s %.0f %.0f %.2f %.2f %.2f %.2f ;\t\t# Left',num2str(UCmatID),numSlab_x,numSlab_y, (h/2 - tf),-(b_eff-b)/2,h/2,-b/2);
        fprintf(fid,'\n\tpatch rect %s %.0f %.0f %.2f %.2f %.2f %.2f ;\t\t\t# Right',num2str(UCmatID),numSlab_x,numSlab_y, (h/2 - tf),b/2,h/2,(b_eff-b)/2);
        fprintf(fid,'\n\tlayer straight %s %.0f %.3f %.2f %.2f %.2f %.2f ;', num2str(SmatID_s),n_slab,As_slab,(h/2)-d_slab,-(b/2 + 3),(h/2)-d_slab,-((b_eff-b)/2 - 3));
        fprintf(fid,'\n\tlayer straight %s %.0f %.3f %.2f %.2f %.2f %.2f ;', num2str(SmatID_s),n_slab,As_slab,(h/2)-d_slab,(b/2 + 3),(h/2)-d_slab,((b_eff-b)/2 - 3));
        
        fprintf(fid,'\n };');
    end
end

% Elements

fprintf(fid,'\n\n# ------------------------------------------------------------------------------');
fprintf(fid,'\n# Elements\n');

% Columns
fprintf(fid,'\n#\tColumns');
fprintf(fid,'\n#\telement "ABCD"');
fprintf(fid,'\n#\t\tA: 1');
fprintf(fid,'\n#\t\tB: 1');
fprintf(fid,'\n#\t\tC: Story');
fprintf(fid,'\n#\t\tD: Column\n');

switch elemType
    case 'GI'
        fprintf(fid,'\n#element gradientInelasticBeamColumn eleTag? iNode? jNode? numIntgrPts? endSecTag1? intSecTag? endSecTag2? lengthR1? lengthR2? lc? transfTag? <-constH> <-integration integrType?> <-iter maxIter? minTol? maxTol?> <-corControl auto/maxEpsInc? maxPhiInc?>');
    case 'FBC'
        fprintf(fid,'\n#element forceBeamColumn $eleTag $iNode $jNode $numIntgrPts $secTag $transfTag <-mass $massDens> <-iter $maxIters $tol> <-integration $intType>');
        
end
for story = 1:Frames(frameNum).Stories
    
    for colLine = 1:Frames(frameNum).Bays+1
        
        if colLine == 1 || colLine == Frames(frameNum).Bays + 1
            lc = Frames(frameNum).Columns(story).Ext.Width;
            R1 = Frames(frameNum).Columns(story).Ext.Width/(Frames(frameNum).Columns(story).Ext.Width*12);
            R2 = Frames(frameNum).Columns(story).Ext.Width/(Frames(frameNum).Columns(story).Ext.Width*12);
            %             lc = Frames(frameNum).Columns(story).Ext.Width;
            %             R1 = Frames(frameNum).Columns(story).Ext.Width;
            %             R2 = Frames(frameNum).Columns(story).Int.Width;
        else
            %             lc = Frames(frameNum).Columns(story).Int.Width;
            %             R1 = Frames(frameNum).Columns(story).Int.Width;
            %             R2 = Frames(frameNum).Columns(story).Int.Width;
            lc = Frames(frameNum).Columns(story).Int.Width;
            R1 = Frames(frameNum).Columns(story).Int.Width/(Frames(frameNum).Columns(story).Int.Width*12);
            R2 = Frames(frameNum).Columns(story).Int.Width/(Frames(frameNum).Columns(story).Int.Width*12);
        end
        
        elemID  = sprintf('11%u%u',story,colLine);
        sectID  = sprintf('11%u%u',story,colLine);
        sectIDmid  = sprintf('21%u%u',story,colLine);
        
        if story == 1
            iNode = sprintf('%u0%u0',story,colLine);
        else
            iNode = sprintf('%u1%u0',story,colLine);
        end
        jNode = sprintf('%u2%u0',story,colLine);
        
        switch elemType
            case 'GI'
                % element gradientInelasticBeamColumn eleTag? iNode? jNode? numIntgrPts? firstSecTag? intSecTag? lastSecTag? secOrder? lc? transfTag? <-integration integrType?> <-iter maxIter? minTol? maxTol?> <-corControl maxEpsInc? maxPhiInc?> <-doPF>
                %                 fprintf(fid,'\nelement gradientInelasticBeamColumn %s %s %s %.0f %s %s %s %.0f %.0f %.0f -integration %s -iter %0.0f %s %s',elemID,iNode,jNode,numIntgPts_c,sectID,sectID,sectID,2,lc,1,integrType,numIntGI,minTol, maxTol);
                %                 fprintf(fid,'\nelement gradientInelasticBeamColumn %s %s %s %.0f %s %s %s %.0f %.0f %.0f %.0f %.0f -integration %s -iter %0.0f %s %s',elemID,iNode,jNode,numIntgPts_c,sectID,sectIDmid,sectID,R1,R2,2,lc,1,integrType,numIntGI,minTol, maxTol);
                fprintf(fid,'\nelement gradientInelasticBeamColumn %s %s %s %.0f %s %s %s %.4f %.4f %.0f %.0f -integration %s -iter %0.0f %s %s',elemID,iNode,jNode,numIntgPts_c,sectID,sectIDmid,sectID,R1,R2,lc,1,integrType,numIntGI,minTol, maxTol);
                
            case 'FBC'
                % element forceBeamColumn $eleTag $iNode $jNode $numIntgrPts $secTag $transfTag <-mass $massDens> <-iter $maxIters $tol> <-integration $intType>');
                fprintf(fid,'\nelement forceBeamColumn %s %s %s %.0f %s %.0f ',elemID,iNode,jNode,numIntgPts_c,sectID,1);
                
        end
        
    end
end

% Beams
fprintf(fid,'\n\n#\tBeams');
fprintf(fid,'\n#\telement "ABCD"');
fprintf(fid,'\n#\t\tA: 1');
fprintf(fid,'\n#\t\tB: 2');
fprintf(fid,'\n#\t\tC: Floor');
fprintf(fid,'\n#\t\tD: Bay\n');

for flr = 2:Frames(frameNum).Stories+1
    
    for bay = 1:Frames(frameNum).Bays
        
        if bay == 1 || bay == Frames(frameNum).Bays
            %             lc = Frames(frameNum).Beams(flr-1).Ext.Depth;
            %             R1 = Frames(frameNum).Beams(flr-1).Ext.Width;
            %             R2 = Frames(frameNum).Beams(flr-1).Int.Width;
            lc = Frames(frameNum).Beams(flr-1).Ext.Depth;
            R1 = Frames(frameNum).Beams(flr-1).Ext.Width/(Frames(frameNum).Beams(flr-1).Ext.Length*12);
            R2 = Frames(frameNum).Beams(flr-1).Int.Width/(Frames(frameNum).Beams(flr-1).Ext.Length*12);
        else
            %             lc = Frames(frameNum).Beams(flr-1).Int.Depth;
            %             R1 = Frames(frameNum).Beams(flr-1).Int.Width;
            %             R2 = Frames(frameNum).Beams(flr-1).Int.Width;
            lc = Frames(frameNum).Beams(flr-1).Int.Depth;
            R1 = Frames(frameNum).Beams(flr-1).Int.Width/(Frames(frameNum).Beams(flr-1).Int.Length*12);
            R2 = Frames(frameNum).Beams(flr-1).Int.Width/(Frames(frameNum).Beams(flr-1).Int.Length*12);
        end
        
        elemID  = sprintf('12%u%u',flr,bay);
        sectID  = sprintf('12%u%u',flr,bay);
        sectIDmid  = sprintf('22%u%u',flr,bay);
        
        iNode = sprintf('%u0%u1',flr,bay);
        jNode = sprintf('%u0%u2',flr,bay);
        
        switch elemType
            case 'GI'
                % element gradientInelasticBeamColumn eleTag? iNode? jNode? numIntgrPts? firstSecTag? intSecTag? lastSecTag? secOrder? lc? transfTag? <-integration integrType?> <-iter maxIter? minTol? maxTol?> <-corControl maxEpsInc? maxPhiInc?> <-doPF>
                %                 fprintf(fid,'\nelement gradientInelasticBeamColumn %s %s %s %.0f %s %s %s %.0f %.0f %.0f -integration %s -iter %0.0f %s %s',elemID,iNode,jNode,numIntgPts_b,sectID,sectID,sectID,2,lc,1,integrType,numIntGI,minTol, maxTol);
                %                 fprintf(fid,'\nelement gradientInelasticBeamColumn %s %s %s %.0f %s %s %s %.0f %.0f %.0f %.0f %.0f -integration %s -iter %0.0f %s %s',elemID,iNode,jNode,numIntgPts_c,sectID,sectIDmid,sectID,R1,R2,2,lc,1,integrType,numIntGI,minTol, maxTol);
                fprintf(fid,'\nelement gradientInelasticBeamColumn %s %s %s %.0f %s %s %s %.4f %.4f %.0f %.0f -integration %s -iter %0.0f %s %s',elemID,iNode,jNode,numIntgPts_b,sectID,sectIDmid,sectID,R1,R2,lc,1,integrType,numIntGI,minTol, maxTol);
            case 'FBC'
                % element forceBeamColumn $eleTag $iNode $jNode $numIntgrPts $secTag $transfTag <-mass $massDens> <-iter $maxIters $tol> <-integration $intType>');
                fprintf(fid,'\nelement forceBeamColumn %s %s %s %.0f %s %.0f ',elemID,iNode,jNode,numIntgPts_b,sectID,1);
                
        end
    end
end

% Rigid Links

fprintf(fid,'\n\n#\tRigid Links\n');
%
% % Column Links
fprintf(fid,'\n#\t\tColumn Links');
% fprintf(fid,'\n#\telement "ABCD"');
% fprintf(fid,'\n#\t\tA: Story');
% fprintf(fid,'\n#\t\tB: 1 = "Lower" rigid link; 2 = "Upper" rigid link');
% fprintf(fid,'\n#\t\tC: Column');
% fprintf(fid,'\n#\t\tD: 0\n');

for story = 1:Frames(frameNum).Stories
    
    for colLine = 1:Frames(frameNum).Bays+1
        
        %         elemID1  = sprintf('%u1%u0',story,colLine);
        %         elemID2  = sprintf('%u2%u0',story,colLine);
        %
        %         if colLine == 1 || colLine == Frames(frameNum).Bays+1
        %             b = Frames(frameNum).Columns(story).Ext.Width;
        %             h = Frames(frameNum).Columns(story).Ext.Depth;
        %         else
        %             b = Frames(frameNum).Columns(story).Ext.Width;
        %             h = Frames(frameNum).Columns(story).Ext.Depth;
        %         end
        %
        %         A_link = b*h;
        %         Iz_link = 1/12 * b*h^3;
        %         J_link = 2*Iz_link;
        
        if story == 1
            iNode1 = NaN;
            iNode2 = sprintf('%u0%u0',story+1,colLine);
            jNode1 = NaN;
            jNode2 = sprintf('%u2%u0',story,colLine);
            
        else
            iNode1 = sprintf('%u0%u0',story,colLine);
            iNode2 = sprintf('%u0%u0',story+1,colLine);
            
            jNode1 = sprintf('%u1%u0',story,colLine);
            jNode2 = sprintf('%u2%u0',story,colLine);
            
        end
        
        if story >1
            fprintf(fid,'\nrigidLink beam %s %s ;',iNode1,jNode1 );
            %             fprintf(fid,'\nelement elasticBeamColumn %s %s %s %.0f %.0f %.0f %.0f %.0f %u ',elemID1,iNode1,jNode1,A_link,E_link,G_link,J_link,Iz_link,1);
            %             fprintf(fid,'\nelement elasticBeamColumn %s %s %s %.0f %.0f %.0f %u ',elemID1,iNode1,jNode1,A_link,E_link,Iz_link,1);
        end
        fprintf(fid,'\nrigidLink beam %s %s ;',iNode2,jNode2 );
        %         fprintf(fid,'\nelement elasticBeamColumn %s %s %s %.0f %.0f %.0f %.0f %.0f %u ',elemID2,iNode2,jNode2,A_link,E_link,G_link,J_link,Iz_link,1);
        %         fprintf(fid,'\nelement elasticBeamColumn %s %s %s %.0f %.0f %.0f %u ',elemID2,iNode2,jNode2,A_link,E_link,Iz_link,1);
        
    end
end

% Beam Links
fprintf(fid,'\n\n#\t\tBeam Links');
% fprintf(fid,'\n#\telement "ABCD"');
% fprintf(fid,'\n#\t\tA: Floor');
% fprintf(fid,'\n#\t\tB: 0');
% fprintf(fid,'\n#\t\tC: Bay');
% fprintf(fid,'\n#\t\tD: 1 = "Left" rigid link; 2 = "Right" rigid link\n');

for flr = 2:Frames(frameNum).Stories+1
    
    for bay = 1:Frames(frameNum).Bays
        
        %         elemID1  = sprintf('%u0%u1',flr,bay);
        %         elemID2  = sprintf('%u0%u2',flr,bay);
        %
        %         if bay == 1 || bay == Frames(frameNum).Bays
        %             b = Frames(frameNum).Beams(flr-1).Ext.Width;
        %             h = Frames(frameNum).Beams(flr-1).Ext.Depth;
        %         else
        %             b = Frames(frameNum).Beams(flr-1).Ext.Width;
        %             h = Frames(frameNum).Beams(flr-1).Ext.Depth;
        %         end
        %
        %         A_link = b*h;
        %         Iz_link = 1/12 * b*h^3;
        %         J_link = 2*Iz_link;
        
        iNode1 = sprintf('%u0%u0',flr,bay);
        iNode2 = sprintf('%u0%u0',flr,bay+1);
        
        jNode1 = sprintf('%u0%u1',flr,bay);
        jNode2 = sprintf('%u0%u2',flr,bay);
        
        fprintf(fid,'\nrigidLink beam %s %s ;',iNode1,jNode1 );
        fprintf(fid,'\nrigidLink beam %s %s ;',iNode2,jNode2 );
        %         fprintf(fid,'\nelement elasticBeamColumn %s %s %s %.0f %.0f %.0f %.0f %.0f %u ',elemID1,iNode1,jNode1,A_link,E_link,G_link,J_link,Iz_link,1);
        %         fprintf(fid,'\nelement elasticBeamColumn %s %s %s %.0f %.0f %.0f %u ',elemID1,iNode1,jNode1,A_link,E_link,Iz_link,1);
        
        %         fprintf(fid,'\nelement elasticBeamColumn %s %s %s %.0f %.0f %.0f %.0f %.0f %u ',elemID2,iNode2,jNode2,A_link,E_link,G_link,J_link,Iz_link,1);
        %         fprintf(fid,'\nelement elasticBeamColumn %s %s %s %.0f %.0f %.0f %u ',elemID2,iNode2,jNode2,A_link,E_link,Iz_link,1);
        
    end
end

% Material Parameters parameter 10001 element 1111 1111 t;

fprintf(fid,'\n\n# ------------------------------------------------------------------------------');
fprintf(fid,'\n# Parameters\n');

% Columns
fprintf(fid,'\n#\tColumns');
fprintf(fid,'\n#\telement "ABCD"');
fprintf(fid,'\n#\t\tA: 1');
fprintf(fid,'\n#\t\tB: 1');
fprintf(fid,'\n#\t\tC: Story');
fprintf(fid,'\n#\t\tD: Column\n');


for story = 1:Frames(frameNum).Stories
    
    for colLine = 1:Frames(frameNum).Bays+1
        
        elemID  = sprintf('11%u%u',story,colLine);
        SmatID  = sprintf('11%u%u',story,colLine);
        UCmatID = sprintf('21%u%u',story,colLine);
        CEmatID = sprintf('31%u%u',story,colLine);
        CMmatID = sprintf('41%u%u',story,colLine);
        
        fprintf(fid,'\nparameter 10%s element %s %s t;',SmatID, elemID, SmatID );
        fprintf(fid,'\nparameter 10%s element %s %s t;',CEmatID, elemID, CEmatID );
        fprintf(fid,'\nparameter 10%s element %s %s t;',CMmatID, elemID, CMmatID );
    end
end

% Beams
fprintf(fid,'\n\n#\tBeams');
fprintf(fid,'\n#\telement "ABCD"');
fprintf(fid,'\n#\t\tA: 1');
fprintf(fid,'\n#\t\tB: 2');
fprintf(fid,'\n#\t\tC: Floor');
fprintf(fid,'\n#\t\tD: Bay\n');

for flr = 2:Frames(frameNum).Stories+1
    
    for bay = 1:Frames(frameNum).Bays 
        elemID  = sprintf('12%u%u',flr,bay);
        SmatID  = sprintf('12%u%u',flr,bay);
        SmatID_c = sprintf('52%u%u',flr,bay);
        SmatID_s = sprintf('62%u%u',flr,bay);
        UCmatID = sprintf('22%u%u',flr,bay);
        CEmatID = sprintf('32%u%u',flr,bay);
        CMmatID = sprintf('42%u%u',flr,bay);
        
        fprintf(fid,'\nparameter 10%s element %s %s t;',SmatID, elemID, SmatID );
        fprintf(fid,'\nparameter 10%s element %s %s t;',SmatID_c, elemID, SmatID_c );
%         fprintf(fid,'\nparameter %s element %s %s t;',SmatID_s, elemID, SmatID_s );
        fprintf(fid,'\nparameter 10%s element %s %s t;',CEmatID, elemID, CEmatID );
        fprintf(fid,'\nparameter 10%s element %s %s t;',CMmatID, elemID, CMmatID );
    end
end





% Nodal Mass

fprintf(fid,'\n\n# ------------------------------------------------------------------------------');
fprintf(fid,'\n# Nodal Mass\n');


for flr = 1:Frames(frameNum).Stories+1
    
    for colLine = 1:Frames(frameNum).Bays+1
        % Columns
        if flr == 1
            if colLine == 1
                h_col_t = Frames(frameNum).Columns(flr).Ext.Depth;
                b_col_t = Frames(frameNum).Columns(flr).Ext.Width;
                l_col_t = Frames(frameNum).Columns(flr).Ext.Length;
                
                h_col_b = 0;
                b_col_b = 0;
                l_col_b = 0;
            else
                h_col_t = Frames(frameNum).Columns(flr).Int.Depth;
                b_col_t = Frames(frameNum).Columns(flr).Int.Width;
                l_col_t = Frames(frameNum).Columns(flr).Int.Length;
                
                h_col_b = 0;
                b_col_b = 0;
                l_col_b = 0;
            end
            
        elseif  flr == Frames(frameNum).Stories+1
            if colLine == 1
                h_col_b = Frames(frameNum).Columns(flr-1).Ext.Depth;
                b_col_b = Frames(frameNum).Columns(flr-1).Ext.Width;
                l_col_b = Frames(frameNum).Columns(flr-1).Ext.Length;
                
                h_col_t = 0;
                b_col_t = 0;
                l_col_t = 0;
            else
                h_col_b = Frames(frameNum).Columns(flr-1).Int.Depth;
                b_col_b = Frames(frameNum).Columns(flr-1).Int.Width;
                l_col_b = Frames(frameNum).Columns(flr-1).Int.Length;
                
                h_col_t = 0;
                b_col_t = 0;
                l_col_t = 0;
            end
            
        else
            
            if colLine == 1 || colLine == Frames(frameNum).Bays+1
                h_col_t = Frames(frameNum).Columns(flr).Ext.Depth;
                b_col_t = Frames(frameNum).Columns(flr).Ext.Width;
                l_col_t = Frames(frameNum).Columns(flr).Ext.Length;
                
                h_col_b = Frames(frameNum).Columns(flr-1).Ext.Depth;
                b_col_b = Frames(frameNum).Columns(flr-1).Ext.Width;
                l_col_b = Frames(frameNum).Columns(flr-1).Ext.Length;
            else
                h_col_t = Frames(frameNum).Columns(flr).Int.Depth;
                b_col_t = Frames(frameNum).Columns(flr).Int.Width;
                l_col_t = Frames(frameNum).Columns(flr).Int.Length;
                
                h_col_b = Frames(frameNum).Columns(flr-1).Int.Depth;
                b_col_b = Frames(frameNum).Columns(flr-1).Int.Width;
                l_col_b = Frames(frameNum).Columns(flr-1).Int.Length;
            end
            
        end
        
        m_col = ((h_col_t/12)*(b_col_t/12)*(l_col_t/2) * wc + (h_col_b/12)*(b_col_b/12)*(l_col_b/2) * wc)/g;
        
        
        % Beams & Slab
        
        if flr == 1
            DL_beam = 0;
            DL_slab = 0;
            m_beam  = 0;
            m_slab  = 0;
        else
            
            if colLine == 1
                h_beam_l = 0;
                b_beam_l = 0;
                l_beam_l = 0;
                
                h_beam_r = Frames(frameNum).Beams(flr-1).Ext.Depth;
                b_beam_r = Frames(frameNum).Beams(flr-1).Ext.Width;
                l_beam_r = Frames(frameNum).Beams(flr-1).Ext.Length;
                
            elseif colLine == 2
                h_beam_l = Frames(frameNum).Beams(flr-1).Ext.Depth;
                b_beam_l = Frames(frameNum).Beams(flr-1).Ext.Width;
                l_beam_l = Frames(frameNum).Beams(flr-1).Ext.Length;
                
                h_beam_r = Frames(frameNum).Beams(flr-1).Int.Depth;
                b_beam_r = Frames(frameNum).Beams(flr-1).Int.Width;
                l_beam_r = Frames(frameNum).Beams(flr-1).Int.Length;
                
            elseif colLine == Frames(frameNum).Bays
                h_beam_l = Frames(frameNum).Beams(flr-1).Int.Depth;
                b_beam_l = Frames(frameNum).Beams(flr-1).Int.Width;
                l_beam_l = Frames(frameNum).Beams(flr-1).Int.Length;
                
                h_beam_r = Frames(frameNum).Beams(flr-1).Ext.Depth;
                b_beam_r = Frames(frameNum).Beams(flr-1).Ext.Width;
                l_beam_r = Frames(frameNum).Beams(flr-1).Ext.Length;
                
            elseif colLine == Frames(frameNum).Bays +1
                h_beam_l = Frames(frameNum).Beams(flr-1).Ext.Depth;
                b_beam_l = Frames(frameNum).Beams(flr-1).Ext.Width;
                l_beam_l = Frames(frameNum).Beams(flr-1).Ext.Length;
                
                h_beam_r = 0;
                b_beam_r = 0;
                l_beam_r = 0;
            else
                h_beam_l = Frames(frameNum).Beams(flr-1).Int.Depth;
                b_beam_l = Frames(frameNum).Beams(flr-1).Int.Width;
                l_beam_l = Frames(frameNum).Beams(flr-1).Int.Length;
                
                h_beam_r = Frames(frameNum).Beams(flr-1).Int.Depth;
                b_beam_r = Frames(frameNum).Beams(flr-1).Int.Width;
                l_beam_r = Frames(frameNum).Beams(flr-1).Int.Length;
            end
            
            DL_beam = ((h_beam_l/12)*(b_beam_l/12)*(l_beam_l/2) * wc + (h_beam_r/12)*(b_beam_r/12)*(l_beam_r/2) * wc);
            DL_slab = ((tf/12)*(l_beam_l/2)*(l_beam_l/2) * wc + (tf/12)*(l_beam_r/2)*(l_beam_r/2) * wc) + ((l_beam_l/2)*(l_beam_l/2) * DL + (l_beam_r/2)*(l_beam_r/2) * DL);
            LL_slab = ((l_beam_l/2)*(l_beam_r/2) * LL);
            m_beam  = DL_beam/g;
            m_slab  = (1.0*DL_slab + 0.2*LL_slab)/g;
        end
        
        m_node = (m_col + m_beam + m_slab);
        
        nodeNum = str2double(sprintf('%u0%u0',flr,colLine));
        fprintf(fid,'\nmass %u0%u0 [expr %.1f*$lbf*$sec*$sec/$ft] [expr %.1f*$lbf*$sec*$sec/$ft] 0.0000001',flr,colLine,m_node,m_node);
        Frames(frameNum).Mass(flr,colLine) = m_node;
        
        wt_node(flr,colLine) = m_node*g;
        
    end
    % Lateral Load Dist. Cals
    Wi(flr) = sum(wt_node(flr,:))/1000;
    nodeNum = sprintf('%u0%u0',flr,colLine);
    Hi(flr) = Frames(frameNum).Node(str2num(nodeNum)).Y;
    WiHi(flr) = Wi(flr) * Hi(flr);
end
sumWiHi = sum(WiHi(:));
W_tot = sum(Wi(:));
Fj = (WiHi./sumWiHi) * W_tot;
% Fi = Fj ./ (Frames(frameNum).Bays+1);
Fi = Wi./sum(Wi);
Frames(frameNum).Height = Hi(end);

% Nodal Loads

fprintf(fid,'\n\n# ------------------------------------------------------------------------------');
fprintf(fid,'\n# Gravity Loads\n');
fprintf(fid,'\npattern Plain $patternTag Linear {');

for flr = 1:Frames(frameNum).Stories+1
    
    for colLine = 1:Frames(frameNum).Bays+1
        % Columns
        if flr == 1
            if colLine == 1 || colLine == Frames(frameNum).Bays+1
                h_col = Frames(frameNum).Columns(flr).Ext.Depth;
                b_col = Frames(frameNum).Columns(flr).Ext.Width;
                l_col = Frames(frameNum).Columns(flr).Ext.Length;
                
            else
                h_col = Frames(frameNum).Columns(flr).Int.Depth;
                b_col = Frames(frameNum).Columns(flr).Int.Width;
                l_col = Frames(frameNum).Columns(flr).Int.Length;
            end
            
        elseif  flr == Frames(frameNum).Stories+1
            h_col = 0;
            b_col = 0;
            l_col = 0;
            
        else
            
            if colLine == 1 || colLine == Frames(frameNum).Bays+1
                h_col = Frames(frameNum).Columns(flr).Ext.Depth;
                b_col = Frames(frameNum).Columns(flr).Ext.Width;
                l_col = Frames(frameNum).Columns(flr).Ext.Length;
                
            else
                h_col = Frames(frameNum).Columns(flr).Int.Depth;
                b_col = Frames(frameNum).Columns(flr).Int.Width;
                l_col = Frames(frameNum).Columns(flr).Int.Length;
                
            end
            
        end
        
        w_col = ((h_col/12)*(b_col/12)*(l_col) * wc )/1000;
        
        
        % Beams & Slab
        
        if flr == 1
            w_beam = 0;
            w_slab = 0;
        else
            
            if colLine == 1
                h_beam_l = 0;
                b_beam_l = 0;
                l_beam_l = 0;
                
                h_beam_r = Frames(frameNum).Beams(flr-1).Ext.Depth;
                b_beam_r = Frames(frameNum).Beams(flr-1).Ext.Width;
                l_beam_r = Frames(frameNum).Beams(flr-1).Ext.Length;
                
            elseif colLine == 2
                h_beam_l = Frames(frameNum).Beams(flr-1).Ext.Depth;
                b_beam_l = Frames(frameNum).Beams(flr-1).Ext.Width;
                l_beam_l = Frames(frameNum).Beams(flr-1).Ext.Length;
                
                h_beam_r = Frames(frameNum).Beams(flr-1).Int.Depth;
                b_beam_r = Frames(frameNum).Beams(flr-1).Int.Width;
                l_beam_r = Frames(frameNum).Beams(flr-1).Int.Length;
                
            elseif colLine == Frames(frameNum).Bays
                h_beam_l = Frames(frameNum).Beams(flr-1).Int.Depth;
                b_beam_l = Frames(frameNum).Beams(flr-1).Int.Width;
                l_beam_l = Frames(frameNum).Beams(flr-1).Int.Length;
                
                h_beam_r = Frames(frameNum).Beams(flr-1).Ext.Depth;
                b_beam_r = Frames(frameNum).Beams(flr-1).Ext.Width;
                l_beam_r = Frames(frameNum).Beams(flr-1).Ext.Length;
                
            elseif colLine == Frames(frameNum).Bays +1
                h_beam_l = Frames(frameNum).Beams(flr-1).Ext.Depth;
                b_beam_l = Frames(frameNum).Beams(flr-1).Ext.Width;
                l_beam_l = Frames(frameNum).Beams(flr-1).Ext.Length;
                
                h_beam_r = 0;
                b_beam_r = 0;
                l_beam_r = 0;
            else
                h_beam_l = Frames(frameNum).Beams(flr-1).Int.Depth;
                b_beam_l = Frames(frameNum).Beams(flr-1).Int.Width;
                l_beam_l = Frames(frameNum).Beams(flr-1).Int.Length;
                
                h_beam_r = Frames(frameNum).Beams(flr-1).Int.Depth;
                b_beam_r = Frames(frameNum).Beams(flr-1).Int.Width;
                l_beam_r = Frames(frameNum).Beams(flr-1).Int.Length;
            end
            
            DL_beam = ((h_beam_l/12)*(b_beam_l/12)*(l_beam_l/2) * wc + (h_beam_r/12)*(b_beam_r/12)*(l_beam_r/2) * wc)/1000;
            DL_slab = ((tf/12)*(l_beam_l/2)*(l_beam_l/2) * wc + (tf/12)*(l_beam_r/2)*(l_beam_r/2) * wc)/1000 + ((l_beam_l/2)*(l_beam_l/2) * DL + (l_beam_r/2)*(l_beam_r/2) * DL)/1000;
            LL_slab = ((l_beam_l/2)*(l_beam_r/2) * LL)/1000;
            w_beam  = DL_beam;
            w_slab  = 1.0*DL_slab + 0.2*LL_slab;
        end
        
        w_node = (w_col + w_beam + w_slab);
        
        nodeNum = str2double(sprintf('%u0%u0',flr,colLine));
        fprintf(fid,'\nload %u0%u0 0 -[expr %.1f*$kip]  0',flr,colLine,w_node);
        Frames(frameNum).Load(flr,colLine) = w_node;
        
    end
end
fprintf(fid,'\n}');

% % Gravity Loads
%
% fprintf(fid,'\n\n# ------------------------------------------------------------------------------');
% fprintf(fid,'\n# Gravity Loads\n');
%
% fprintf(fid,'\npattern Plain 101 Linear {');
%
% % DEAD LOAD ONLY
% % Columns
% for story = 1:Frames(frameNum).Stories
%
%     for colLine = 1:Frames(frameNum).Bays+1
%
%         if colLine == 1 || colLine == Frames(frameNum).Bays+1
%             h_col = Frames(frameNum).Columns(story).Ext.Depth;
%             b_col = Frames(frameNum).Columns(story).Ext.Width;
%             l_col = Frames(frameNum).Columns(story).Ext.Length;
%
%         else
%             h_col = Frames(frameNum).Columns(story).Int.Depth;
%             b_col = Frames(frameNum).Columns(story).Int.Width;
%             l_col = Frames(frameNum).Columns(story).Int.Length;
%
%         end
%
%         w_col = ((h_col/12)*(b_col/12)*wc);
%
%         fprintf(fid,'\n\teleLoad -ele 11%u%u -type -beamUniform 0 -[expr %.0f*$lbf/$ft]',story,colLine,w_col);
%
%     end
% end
%
%
% % Beams
% for flr = 2:Frames(frameNum).Stories+1
%
%     for bay = 1:Frames(frameNum).Bays
%
%         if flr == 1
%             m_beam = 0;
%             m_slab = 0;
%         else
%
%             if bay == 1 || bay == Frames(frameNum).Bays
%
%                 h_beam = Frames(frameNum).Beams(flr-1).Ext.Depth;
%                 b_beam = Frames(frameNum).Beams(flr-1).Ext.Width;
%                 l_beam_l = Frames(frameNum).Beams(flr-1).Ext.Length;
%                 l_beam_r = Frames(frameNum).Beams(flr-1).Ext.Length;
%
%             elseif bay == 2 || bay == Frames(frameNum).Bays-1
%                 h_beam = Frames(frameNum).Beams(flr-1).Int.Depth;
%                 b_beam = Frames(frameNum).Beams(flr-1).Int.Width;
%                 l_beam_r = Frames(frameNum).Beams(flr-1).Ext.Length;
%                 l_beam_l = Frames(frameNum).Beams(flr-1).Int.Length;
%             else
%                 h_beam = Frames(frameNum).Beams(flr-1).Int.Depth;
%                 b_beam = Frames(frameNum).Beams(flr-1).Int.Width;
%                 l_beam_r = Frames(frameNum).Beams(flr-1).Int.Length;
%                 l_beam_l = Frames(frameNum).Beams(flr-1).Int.Length;
%             end
%
%             w_beam = (h_beam/12)*(b_beam/12) * wc;
%             w_slab = (tf/12)*((l_beam_r+l_beam_l)/2) * wc ;
%
%             w_tot = (w_beam + w_slab);
%             fprintf(fid,'\n\teleLoad -ele 12%u%u -type -beamUniform -[expr %.0f*$lbf/$ft] 0 ',flr,bay,w_tot);
%
%         end
%     end
% end

% fprintf(fid,'\n}');


fprintf(fid,'\n\n\n\n\n\n\n\n\n\n');
fclose(fid);