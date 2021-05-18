function [n,n_h,n_b,n_c, a_barSize, a_barArea, rho_actual,a_barSize_c, a_barArea_c, rho_actual_c, n_sh, a_barSize_sh, a_barArea_sh, a_rho_sh] = rebarLayout(h, b, cover, rho, rho_c, rho_sh,s,col,beam ,plotVal,memType)
%Determines a likely rebar size based on member geometery and given
%reinforcement ratios
%
%[n,n_h,n_b,n_c, a_barSize, a_barArea, rho_actual,a_barSize_c, a_barArea_c, rho_actual_c, n_sh, a_barSize_sh, a_rho_sh] = rebarLayout(h, b, cover, rho, rho_c, rho_sh,s,col,beam ,plotVal,memType)

% Find "Most Probable" rebar layout from given dimensions and reinforcement
% ratios


% clearvars; close all; clc

load Rebar.mat Rebar

% memType = 'Column'; %'Beam'

% h       = 22;   % [in] Section depth
% b       = 22;   % [in] Section width
% cover   = 2.5;  % [in] Cover, to center of longitudinal rebar

% rho     = 0.0110;
% rho_c   = NaN;
% rho_sh  = 0.007;

switch memType
    
    case 'column'
        
        % Longitudinal
        
        As      = rho*b *h;
        barArea = nan(length(Rebar),1);
        nBar = nan(length(Rebar),1);
        
        for jj = 1:length(Rebar)
            barArea(jj) = Rebar(jj).Area;
            nBar(jj) = As/barArea(jj);
            nBar(jj) = nBar(jj) - mod(nBar(jj),2);
            
            if nBar(jj) < 4 || barArea(jj) < Rebar(5).Area || barArea(jj) > Rebar(10).Area || nBar(jj) > 20
                nBar(jj) = NaN;
            end
            
            rho_a(jj) = (nBar(jj)*barArea(jj))/(b*h);
            
            diffNum(jj) = abs(rho_a(jj)-rho);
            
        end
        
        n = nBar(diffNum == min(diffNum));
        a_barSize = Rebar(diffNum == min(diffNum)).Size;
        a_barArea = Rebar(diffNum == min(diffNum)).Area;
        
        if length(n) >1
            n = min(n);
            a_barSize = max(Rebar(diffNum == min(diffNum)).Size);
            a_barArea = max(Rebar(diffNum == min(diffNum)).Area);
        end
        
        rho_actual  = (n.*a_barArea)/(b*h);
        
        if mod(n,4) > 0
            
            if b>=h
                n_b = (n - mod(n,4))/4 + 2;
                n_h = (n - mod(n,4))/4 + 1;
            elseif h>b
                n_b = (n - mod(n,4))/4 + 1;
                n_h = (n - mod(n,4))/4 + 2;
            end
            
            
        else
            
            n_h = n/4 +1;
            n_b = n/4 +1;
        end
        
        
        % Shear
        
        As_sh       = rho_sh*b*s;
        n_sh        = (n + 4)/4 ;
        
        for jj = 1:length(Rebar)
            
            barAreash(jj) = Rebar(jj).Area;
            
            if barAreash(jj) < Rebar(2).Area || barAreash(jj) > Rebar(4).Area
                barAreash(jj) = NaN;
            end
            
            rho_ash(jj) = (n_sh*barAreash(jj))/(b*s);
            
            diffNum(jj) = abs(rho_ash(jj)-rho_sh);
        end
        a_barSize_sh = Rebar(diffNum == min(diffNum)).Size;
        a_barArea_sh = Rebar(diffNum == min(diffNum)).Area;
        a_rho_sh = n_sh * a_barArea_sh/ (b*s);
        
        
        n_c = NaN; a_barSize_c = NaN; a_barArea_c = NaN; rho_actual_c = NaN;
        
        
        % return
        
    case 'beam'
        % Longitudinal Tension
        
        As      = rho * b * (h-cover);
        barArea = nan(length(Rebar),1);
        nBar = nan(length(Rebar),1);
        
        for jj = 1:length(Rebar)
            barArea(jj) = Rebar(jj).Area;
            nBar(jj) = As/barArea(jj);
            nBar(jj) = nBar(jj) - mod(nBar(jj),2);
            
            if nBar(jj) < 2 || barArea(jj) < Rebar(5).Area || barArea(jj) > Rebar(10).Area || nBar(jj) > 10
                nBar(jj) = NaN;
            end
            
            rho_a(jj) = (nBar(jj)*barArea(jj))/(b*(h-cover));
            
            diffNum(jj) = abs(rho_a(jj)-rho);
            
        end
        
        n_b = nBar(diffNum == min(diffNum));
        a_barSize = Rebar(diffNum == min(diffNum)).Size;
        a_barArea = Rebar(diffNum == min(diffNum)).Area;
        
        if length(n_b) >1
            n_b = min(n_b);
            a_barSize = max(Rebar(diffNum == min(diffNum)).Size);
            a_barArea = max(Rebar(diffNum == min(diffNum)).Area);
        end
        
        rho_actual  = (n_b.*a_barArea)/(b*(h-cover));
        
        % Longitudinal Compression
        
        As_c      = rho_c * b * (h-cover);
        barArea = nan(length(Rebar),1);
        nBar = nan(length(Rebar),1);
        clearvars diffNum
        
        for jj = 1:length(Rebar)
            barArea(jj) = Rebar(jj).Area;
            nBar(jj) = As_c/barArea(jj);
            nBar(jj) = nBar(jj) - mod(nBar(jj),2);
            
            if nBar(jj) < 2 || barArea(jj) < Rebar(2).Area || barArea(jj) > Rebar(10).Area || nBar(jj) > 10
                nBar(jj) = NaN;
            end
            
            rho_a(jj) = (nBar(jj)*barArea(jj))/(b*(h-cover));
            
            diffNum(jj) = abs(rho_a(jj)-rho_c);
            
        end
        
        n_c = nBar(diffNum == min(diffNum));
        a_barSize_c = Rebar(diffNum == min(diffNum)).Size;
        a_barArea_c = Rebar(diffNum == min(diffNum)).Area;
        
        if length(n_c) >1
            n_c = min(n_c);
            a_barSize_c = max(Rebar(diffNum == min(diffNum)).Size);
            a_barArea_c = max(Rebar(diffNum == min(diffNum)).Area);
        end
        
        rho_actual_c  = (n_c.*a_barArea_c)/(b*(h-cover));
        
        
        % Shear
        
        As_sh       = rho_sh*b*s;
        n_sh        = min(n_b,n_c);
        
        for jj = 1:length(Rebar)
            
            barAreash(jj) = Rebar(jj).Area;
            
            if barAreash(jj) < Rebar(2).Area || barAreash(jj) > Rebar(4).Area
                barAreash(jj) = NaN;
            end
            
            rho_ash(jj) = (n_sh*barAreash(jj))/(b*s);
            
            diffNum(jj) = abs(rho_ash(jj)-rho_sh);
        end
        a_barSize_sh = Rebar(diffNum == min(diffNum)).Size;
        a_barArea_sh = Rebar(diffNum == min(diffNum)).Area;
        a_rho_sh = n_sh * a_barArea_sh/ (b*s);
        
        
        n_h = NaN; n = NaN;
end
% 
% if plotVal == 1
%     
%     % Draw Section
%     
%     f = figure(); hold all; axis('equal');
%     set(f,'Color','w');
%     xlim([-(max(h,b)/2 +5) (max(h,b)/2 +5)])
%     ylim([-(max(h,b)/2 +5) (max(h,b)/2 +5)])
%     
%     title({[memType], ['Column Line/Bay: ' num2str(col), '     Story/Floor: ', num2str(beam)]})
%     
%     x = -b/2;
%     y = -h/2;
%     
%     % Plot member outline
%     rectangle('Position',[x y b h],'FaceColor',[0.75 0.75 0.75])
%     
%     % Plot Rebar
%     switch memType
%         case 'column'
%             
%             diam = a_barSize/8;
%             sx = (b - 2*cover)/(n_b-1);
%             sy = (h - 2*cover)/(n_h-1);
%             
%             
%             
%             for jj = 1:n_h
%                 
%                 if jj == 1
%                     y(jj) = -(h/2) + cover;
%                 elseif jj >1
%                     y(jj) = y(jj-1) + sy;
%                 end
%                 
%                 for ii = 1:n_b
%                     
%                     if ii == 1
%                         x(jj,ii) = -(b/2) + cover;
%                     elseif ii == n_b
%                         x(jj,ii) = (b/2) - cover;
%                     else
%                         if jj == 1 || jj == n_h
%                             x(jj,ii) = x(jj,ii-1) + sx;
%                         else
%                             x(jj,ii) = NaN;
%                         end
%                     end
%                     
%                     rebar(x(jj,ii),y(jj),diam/2);
%                     
%                     
%                 end
%             end
%             
%                 % Plot Dimensions
%             
%                 for jj = 1:n_h-1
%                     % Veritcal arrow
%                     arrow([-(b/2)-2 ((y(jj+1)+y(jj)))/2],[-(b/2)-2 y(jj)],'Length',8)
%                     arrow([-(b/2)-2 ((y(jj+1)+y(jj)))/2],[-(b/2)-2 y(jj+1)],'Length',8)
%                 end
%             
%                 for ii = 1:n_b-1
%             
%                     % Horizontal arrow
%                     arrow([((x(1,ii+1)+x(1,ii)))/2 -(h/2)-2],[x(1,ii) -(h/2)-2],'Length',8)
%                     arrow([((x(1,ii+1)+x(1,ii)))/2 -(h/2)-2],[x(1,ii+1) -(h/2)-2],'Length',8)
%             
%             
%                 end
%         case 'beam'
%             diam = a_barSize/8;
%             diam_c = a_barSize_c/8;
%             sb = (b - 2*cover)/(n_b-1);
%             sc = (b - 2*cover)/(n_c-1);
%             sy = (h - 2*cover);
%             
%             
%             for ii = 1:n_b
%                 
%                 if ii == 1
%                     x(ii) = -(b/2) + cover;
%                 elseif ii == n_b
%                     x(ii) = (b/2) - cover;
%                 else
%                     x(ii) = x(ii-1) + sb;
%                    
%                 end
%                 
%                 y = -h/2 + cover;
%                 
%                 rebar(x(ii),y,diam/2);
%                 
%             end
%             
%             clearvars x
%             
%             for ii = 1:n_c
%                 
%                 if ii == 1
%                     x(ii) = -(b/2) + cover;
%                 elseif ii == n_c
%                     x(ii) = (b/2) - cover;
%                 else
%                     x(ii) = x(ii-1) + sc;
%                    
%                 end
%                 y = h/2 - cover;
%                 
%                 rebar(x(ii),y,diam_c/2);
%                 
%             end
%             
%     end
%     
%     
%     
%     
% end
end

