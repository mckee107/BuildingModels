% Codi McKee
% Texas A&M University
% First Created: 30-Feb-2019
% Last Modified:
% TO DO:
%

% Create Ground Motion Time Series Files for Dynamic Analysis


for ii = 1:2
    
    if ii == 1
        scaledStr = '_DE';
        scaleFac = 1;
    elseif ii == 2
        scaledStr = '_MCE';
        scaleFac = 1.5;
    end
    
    for GroundMotion_num = [1:44]
        
        filename = ['Data\Ground Motions\AT_' num2str(GroundMotion_num ) scaledStr '.tcl'];
        % return
        
        fid = fopen(filename,'w');
        
        fprintf(fid,'%0.6E\n', GroundMotions(GroundMotion_num).TimeSeries.*scaleFac);
        
        fclose(fid);
        
    end
end

% save GroundMotions.mat GroundMotions 