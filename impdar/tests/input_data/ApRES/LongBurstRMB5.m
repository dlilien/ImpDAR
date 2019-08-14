function vdat = LongBurstRMB5(Filename, Burst, SamplingFrequency, MaxNChirp)

% vdat = LongBurstRMB5(Filename, Burst, SamplesPerChirp, MaxNChirp)
%
% Load data from RMB2 burst with 1 attenuator setting. Burst average by stacking to limit memory use
% Used for averaging very long bursts.

MaxHeaderLen = 1500;
BytesPerWord = 2;
burstpointer = 0;
vdat.Code = 0;
NChirp = 0;

% Find length of chirp to allocate space for average
fid = fopen(Filename,'r');
A = fread(fid,MaxHeaderLen,'*char');
A = A';
SearchString = 'N_ADC_SAMPLES=';
searchind = strfind(A,SearchString);
searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
vdat.Nsamples = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'%d');
WperChirpCycle = vdat.Nsamples;
fseek(fid,0,'bof');

% Allocate space for average v
vdat.v = zeros(1,vdat.Nsamples);

if fid >= 0
    fseek(fid,0,'eof');
    filelength = ftell(fid);
    BurstCount = 1;
    while BurstCount <= Burst && burstpointer <= filelength - MaxHeaderLen
        fseek(fid,burstpointer,'bof');
        A = fread(fid,MaxHeaderLen,'*char');
        A = A';
        try
            SearchString = 'NSubBursts=';
            searchind = strfind(A,SearchString);
            searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
            vdat.SubBurstsInBurst = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'%d');
            
            SearchString = 'Average=';
            searchind = strfind(A, SearchString);
            if isempty(searchind)
                vdat.Average = 0; %cls 9/jan/14 -average not included in mooring deploy
            else
                searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
                vdat.Average = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'%d');
            end
            
            SearchString = 'nAttenuators=';
            searchind = strfind(A, SearchString);
            searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
            vdat.NAttenuators = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'%d',1);
            
            SearchString = 'Attenuator1=';
            searchind = strfind(A, SearchString);
            searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
            vdat.Attenuator_1 = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'%f',vdat.NAttenuators);
            
            SearchString = 'AFGain=';
            searchind = strfind(A, SearchString);
            searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
            vdat.Attenuator_2 = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'%f',vdat.NAttenuators);
            
            SearchString = 'TxAnt=';
            searchind = strfind(A, SearchString);
            searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
            vdat.TxAnt = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'%d',8);
            
            SearchString = 'RxAnt=';
            searchind = strfind(A, SearchString);
            searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
            vdat.RxAnt = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'%d',8);
            
            ind = find(vdat.TxAnt~=1);
            vdat.TxAnt(ind) = [];
            ind = find(vdat.RxAnt~=1);
            vdat.RxAnt(ind) = [];
            
            if vdat.Average
                vdat.ChirpsInBurst = 1;
            else
                vdat.ChirpsInBurst = vdat.SubBurstsInBurst * length(vdat.TxAnt) * ...
                    length(vdat.RxAnt) * vdat.NAttenuators;
            end
            
            SearchString = '*** End Header ***';
            searchind = strfind(A, SearchString);
            
            burstpointer = burstpointer + searchind(1) + length(SearchString);
        catch
            vdat.Code = -2;
            vdat.Burst = BurstCount;
            return
        end
        WordsPerBurst = vdat.ChirpsInBurst * WperChirpCycle;
        if BurstCount < Burst && burstpointer <= filelength - MaxHeaderLen
            burstpointer = burstpointer + vdat.ChirpsInBurst * WperChirpCycle * BytesPerWord;
        end
        BurstCount = BurstCount + 1;
    end
    
    % Extract remaining information from header
    SearchString = 'Time stamp=';
    searchind = strfind(A, SearchString);
    if isempty(searchind)
        vdat.Code = -4;
        return
    end
    try
        searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
        td = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),...
            '%d-%d-%d %d:%d:%d');
        vdat.TimeStamp = datenum(td(1),td(2),td(3),td(4),td(5),td(6));
    catch err
        vdat.Code = 1;
    end
    
    SearchString = 'Temp1=';
    searchind = strfind(A, SearchString);
    try
        searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
        vdat.Temperature_1 = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'%f');
    catch err
        vdat.Code = 1;
    end
    
    SearchString = 'Temp2=';
    searchind = strfind(A, SearchString);
    try
        searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
        vdat.Temperature_2 = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'%f');
    catch err
        vdat.Code = 1;
    end
    
    SearchString = 'BatteryVoltage=';
    searchind = strfind(A, SearchString);
    try
        searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
        vdat.BatteryVoltage = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'%f');
    catch err
        vdat.Code = 1;
    end
    
    fseek(fid,burstpointer-1,'bof');
    if BurstCount == Burst + 1
        %while vdat.Code == 0 && NChirp < MaxNChirp
        while vdat.Code == 0 && NChirp < vdat.ChirpsInBurst && NChirp < MaxNChirp
            if vdat.Average == 2
                [tmp count] = fread(fid,vdat.Nsamples,'*uint32','ieee-le');
            else
                [tmp count] = fread(fid,vdat.Nsamples,'*uint16','ieee-le');
            end
            if count < vdat.Nsamples
                vdat.Code = 2;
            else
                vdat.v = vdat.v + double(tmp');
                NChirp = NChirp + 1;
                fseek(fid,burstpointer-1+vdat.Nsamples * NChirp * BytesPerWord,'bof');
            end
        end
        %vdat.v = vdat.v/NChirp;
        vdat.v = vdat.v * 2.5 / 2^16;
        vdat.Startind = 1;
        vdat.Endind = vdat.Startind + WperChirpCycle - 1;
        vdat.Burst = Burst;
        vdat.NChirp = NChirp; % save number of chirps so we can average multiple bursts
    else
        % Too few bursts in file
        vdat.Burst = BurstCount - 1;
        vdat.Code = -4;
    end
    fclose(fid);
else
    % Unknown file
    vdat.Code = -1;
end
