% Codi McKee
% Texas A&M University
% First Created: 23-Feb-2019
% Last Modified: 23-Feb-2019
% TO DO:
%   

% Build Pushover Analysis File


%**************************************************************************
% Write Pushover File
%

filename = ['Input\Pushover_'  num2str(t) 'yr_Corr' num2str(corrosion_level) '_' num2str(frameID) '.tcl'];


fid = fopen(filename,'w');

% Header Information

fprintf(fid,'# Codi McKee');
fprintf(fid,'\n# Texas A&M University');
fprintf(fid,'\n# %s',char(datetime('today')));
fprintf(fid,'\n\n# This file was written via Matlab %s',version);
fprintf(fid,'\n# File: %s\\%s.m',pwd,mfilename);

fprintf(fid,'\n\n# Pushover Loads: Frame ID %s',num2str(frameID));

%**************************************************************************
% Pushover Loads
%
fprintf(fid,'\n\nset IDctrlNode %u010 ;',Frames(frameNum).Stories+1);
fprintf(fid,'\nset IDctrlDOF 1;');


%**************************************************************************
% Pushover Loads
%

fprintf(fid,'\n\n# ------------------------------------------------------------------------------');
fprintf(fid,'\n# Pushover Loads');
fprintf(fid,'\n#\n');
fprintf(fid,'\npattern Plain $patternTag Linear {');

for flr = 2:Frames(frameNum).Stories+1
    
    for colLine = 1:Frames(frameNum).Bays+1
        
        nodeNum = sprintf('%u0%u0',flr,colLine);
        
        fprintf(fid,'\n\tload %u0%u0 [expr %.5f*$kip] 0  0',flr,colLine,Fi(flr));
        
    end
end

fprintf(fid,'\n}');








fprintf(fid,'\n\n\n\n\n\n\n\n\n\n');
fclose(fid);