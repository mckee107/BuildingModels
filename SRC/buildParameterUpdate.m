% Codi McKee
% Texas A&M University
% First Created: 05-Feb-2019
% Last Modified: 30-May-2019
% TO DO:
%

%**************************************************************************
% Write File
%

filename = ['Input\Parameters_Corr' num2str(corrosion_level) '_' num2str(frameID) nameSuffix '.tcl'];

fid = fopen(filename,'w');

% Header Information
fprintf(fid,'#!/usr/bin/tclsh');
fprintf(fid,'# Codi McKee');
fprintf(fid,'\n# Texas A&M University');
fprintf(fid,'\n# %s',char(datetime('today')));
fprintf(fid,'\n\n# This file was written via Matlab %s',version);
fprintf(fid,'\n# File: %s\\%s.m',pwd,mfilename);

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
        
        fprintf(fid,'\nupdateParameter 10%s $TVal;',SmatID );
        fprintf(fid,'\nupdateParameter 10%s $TVal;',CEmatID );
        fprintf(fid,'\nupdateParameter 10%s $TVal;',CMmatID );
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
        
        fprintf(fid,'\nupdateParameter 10%s $TVal;',SmatID);
        fprintf(fid,'\nupdateParameter 10%s $TVal;',SmatID_c);
%         fprintf(fid,'\nupdateParameter 10%s $Age;',SmatID_s);
        fprintf(fid,'\nupdateParameter 10%s $TVal;',CEmatID);
        fprintf(fid,'\nupdateParameter 10%s $TVal;',CMmatID);
    end
end


fprintf(fid,'\n\n\n\n\n\n\n\n\n\n');
fclose(fid);