function H = fmcwParametersRMB2(dataFile)

%Extract from the hex codes the actual paramaters used by RMB2
%The contents of config.ini are copied into a data header.
%Note this script assumes that the format of the hex codes have quotes
%e.g. Reg02="0D1F41C8"

%Checks for a sampling frequency of 40 or 80 KHz.  Apart from Lai Bun's
%variant (WDC + Greenland) it could be hard coded to 40 KHz.

%However, there is no check made by the system that the N_ADC_SAMPLES
%matches the requested chirp length

%NOT COMPLETE - needs a means of checking for profile mode, where multiple sweeps
%per period are transmitted- see last line

H = struct;

fsysclk = 1e9;
H.fs = 4e4;
% fprintf(1,'ASSUMPTIONS:DDS clock = %d\n',fsysclk);

if nargin == 0
     [fileName, pathName] = uigetfile('*.dat','Choose data file');
     dataFile = [pathName fileName];
end
fid = fopen(dataFile,'rt');
A = fread(fid,2000,'*char');
A = A';
fclose(fid);
loc1 = strfind(A,'Reg0');
loc2 = strfind(A,'="');

for k = 1:length(loc1)
   switch(A(loc1(k):loc2(k)))
       case 'Reg01=' %Control Function Register 2 (CFR2)—Address 0x01 Four bytes
        %Bit 19 (Digital ramp enable)= 1 = Enables digital ramp generator functionality.
        %Bit 18 (Digital ramp no-dwell high) 1 = enables no-dwell high functionality.
        %Bit 17 (Digital ramp no-dwell low) 1 = enables no-dwell low functionality.
        %With no-dwell high, a positive transition of the DRCTL pin initiates a positive slope ramp, which
        %continues uninterrupted (regardless of any activity on the DRCTL pin) until the upper limit is reached.
        %Setting both no-dwell bits invokes a continuous ramping mode of operation;
        loc3 = strfind(A(loc2(k)+2:end),'"');
        val = A((loc2(k)+2:loc2(k)+loc3(1)));
        val = dec2bin(hex2dec(val)); val = fliplr(val);
        H.noDwellHigh = str2num(val(18+1));
        H.noDwellLow = str2num(val(17+1));
 
%        case 'Reg08' %Phase offset word Register (POW)—Address 0x08. 2 Bytes dTheta = 360*POW/2^16.
%         val = char(reg{1,2}(k)); 
%         H.phaseOffsetDeg = hex2dec(val(1:4))*360/2^16;

       case 'Reg0B=' %Digital Ramp Limit Register—Address 0x0B
        %63:32 Digital ramp upper limit 32-bit digital ramp upper limit value.
        %31:0 Digital ramp lower limit 32-bit digital ramp lower limit value.
         loc3 = strfind(A(loc2(k)+2:end),'"');
        val = A((loc2(k)+2:loc2(k)+loc3(1)));
        H.startFreq = hex2dec(val(9:end))*fsysclk/2^32;
        H.stopFreq = hex2dec(val(1:8))*fsysclk/2^32;
        
       case 'Reg0C='  %Digital Ramp Step Size Register—Address 0x0C
        %63:32 Digital ramp decrement step size 32-bit digital ramp decrement step size value.
        %31:0 Digital ramp increment step size 32-bit digital ramp increment step size value.
        loc3 = strfind(A(loc2(k)+2:end),'"');
        val = A((loc2(k)+2:loc2(k)+loc3(1)));
        H.rampUpStep = hex2dec(val(9:end))*fsysclk/2^32;
        H.rampDownStep = hex2dec(val(1:8))*fsysclk/2^32;
        
       case 'Reg0D='  %Digital Ramp Rate Register—Address 0x0D
        %31:16 Digital ramp negative slope rate 16-bit digital ramp negative slope value that defines the time interval between decrement values.
        %15:0 Digital ramp positive slope rate 16-bit digital ramp positive slope value that defines the time interval between increment values.
        loc3 = strfind(A(loc2(k)+2:end),'"');
        val = A((loc2(k)+2:loc2(k)+loc3(1)));
        H.tstepUp = hex2dec(val(5:end))*4/fsysclk;
        H.tstepDown = hex2dec(val(1:4))*4/fsysclk;       
   end
end

loc = strfind(A,'SamplingFreqMode=');
searchCR = strfind(A(loc(1):end),[char(10)]);
H.fs = sscanf(A(loc(1)+length(['SamplingFreqMode=']):searchCR(1)+loc(1)),'%d\n');
if H.fs == 1
    H.fs = 8e4;
else
    H.fs = 4e4;
end
% if(H.fs > 70e3)
%     H.fs = 80e3;
% else
%     H.fs = 40e3;
% end

loc = strfind(A,'N_ADC_SAMPLES=');
searchCR = strfind(A(loc(1):end),[char(10)]);
H.Nsamples = sscanf(A(loc(1)+length(['N_ADC_SAMPLES=']):searchCR(1)+loc(1)),'%d\n');

H.nstepsDDS = round(abs((H.stopFreq - H.startFreq)/H.rampUpStep));%abs as ramp could be down
H.chirpLength = H.nstepsDDS * H.tstepUp;
H.nchirpSamples = round(H.chirpLength * H.fs);

% If number of ADC samples collected is less than required to collect
% entire chirp, set chirp length to length of series actually collected
if H.nchirpSamples > H.Nsamples
    H.chirpLength = H.Nsamples / H.fs;
end

H.K = 2*pi*(H.rampUpStep/H.tstepUp); % chirp gradient (rad/s/s)
if(H.stopFreq > 400e6)
    H.rampDir = 'down';
else
    H.rampDir = 'up';
end

if(H.noDwellHigh && H.noDwellLow)
    H.rampDir = 'upDown';
    H.nchirpsPerPeriod = NaN;% H.nchirpSamples/(H.chirpLength);
end
