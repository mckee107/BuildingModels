function [units] = setUnits(BaseUnits)
%setUnits Function to set unit system for use in calculations
%   BaseUnits: 'US' or 'SI'

switch BaseUnits
    
    case 'US'
        % length
        units.in = 1;
        units.ft = 12 * units.in;
        units.yd = 3 * units.ft;
        units.mi = 5280*units.ft;
        
        units.mm = 0.03937 * units.in;
        units.cm = 10 * units.mm;
        units.m  = 100 * units.cm;
        units.km  = 1000 * units.m;
        
        % Force
        units.kip = 1;
        units.lbf = units.kip/1000;
        
        units.N   = 0.224809 * units.lbf;
        units.kN  = 1000 * units.N;
        
        
        % Stress
        units.ksi = 1;
        units.psi = units.ksi/1000;
        
        units.Pa  = 0.000145038 * units.psi;
        units.kPa = 1000 * units.Pa;
        units.mPa = 1000000 * units.Pa;
        units.gPa = 1000000000 * units.Pa;
        
        
        % Mass
        units.slug = 1;
        
        units.g    = 1/14593.903 * units.slug;
        units.kg   = 1000 * units.g;
        
        % Area
        units.in2 = units.in * units.in;
        units.ft2 = units.ft * units.ft;
        
        units.mm2 = units.mm * units.mm;
        units.cm2 = units.cm * units.cm;
        units.m2  = units.m * units.m;
        
        % Volume
        units.in3 = units.in * units.in * units.in;
        units.ft3 = units.ft * units.ft * units.ft;
        
        units.mm3 = units.mm * units.mm * units.mm;
        units.cm3 = units.cm * units.cm * units.cm;
        units.m3  = units.m * units.m * units.m;
        
        % Time
        units.yr  = 1;
        units.min = 365*24*60;
        units.s   = 60 * units.min;
        
    case 'SI'
        % length
        %         units.m  = 100 * units.cm;
        %         units.cm = units.m / 100;
        units.mm = 1;
        units.cm = 10 * units.mm;
        units.m  = 100 * units.cm;
        units.km = 1000 * units.m;
        
        units.in = 2.54 * units.cm;
        units.ft = 12 * units.in;
        units.yd = 3 * units.ft;
        units.mi = 5280 * units.ft;
        
        
        
        % Force
        units.N   = 1;
        units.kN  = 1000 * units.N;
        
        units.lbf = 4.44822 * units.N;
        units.kip = 1000 * units.lbf;
        
        % Stress
        units.mPa = 1;
        units.Pa  = units.mPa / 1000000;
        units.kPa = units.mPa / 1000;
        units.gPa = 1000 * units.mPa;
        
        units.ksi = 6.895 * units.mPa;
        units.psi = units.ksi/1000;
        
        
        % Mass
        units.kg   = 1;
        units.g    = units.kg / 1000;
        
        units.slug = 14593.9 * units.g;
        
        % Area
        units.in2 = units.in * units.in;
        units.ft2 = units.ft * units.ft;
        
        units.mm2 = units.mm * units.mm;
        units.cm2 = units.cm * units.cm;
        units.m2  = units.m * units.m;
        
        % Volume
        units.in3 = units.in * units.in * units.in;
        units.ft3 = units.ft * units.ft * units.ft;
        
        units.mm3 = units.mm * units.mm * units.mm;
        units.cm3 = units.cm * units.cm * units.cm;
        units.m3  = units.m * units.m * units.m;
        
        % Time
        units.yr  = 1;
        units.min = 365*24*60;
        units.s   = 60 * units.min;
        
    case ''
        
end
















end

