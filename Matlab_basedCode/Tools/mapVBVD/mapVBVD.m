function [image_obj noise_obj phasecor_obj refscan_obj refscanPC_obj...
    RTfeedback_obj phasestab_obj] = mapVBVD(filename,varargin)

%  Reads Siemens raw .dat file from VB/VD MRI raw data.
%
%  Requires twix_map_obj.m
%
%  Author: Philipp Ehses (philipp.ehses@tuebingen.mpg.de)
%
%
%  Philipp Ehses 11.02.07, original version
%  [..]
%  Philipp Ehses 22.03.11, port to VD11
%  Felix Breuer  31.03.11, added support for noise & ref scans, speed fixes
%  Philipp Ehses 19.08.11, complete reorganization of code, added
%                          siemens_data_obj class to improve readability
%  Wolf Blecher  15.05.12, readout of slice position parameters for VB Data sets
%  Wolf Blecher  11.06.12, added distinction between PATREF and PATREF PHASCOR
%  Philipp Ehses 02.08.12, again massive code reorganization. The new class
%                          twix_map_obj.m now stores the memory position of
%                          each dataset (coils are included - size: NCol*NCha)
%                          The actual data is not read until it is demanded
%                          by a "data_obj()" call. This makes it possible
%                          to selectively read in only parts of the data
%                          (e.g. to preserve memory).
%  Philipp Ehses 27.08.12, speed optimizations (avoiding of .-subsref calls)
%  07.09.12 Thanks to Stephen Yutzy for implementing support for raw data 
%           correction (currently only supported for VB software version).
% 
% Input:
%
% filename or simply meas. id, e.g. mapVBVD(122) (if file is in same path)
% optional arguments (see below)
%
% Output:
%
% image_obj:       object for image scan
% noise_obj:       object for noise scan
% phasecor_obj:    object for phase correction scan
% refscan_obj:     object for reference scan
% refscanPC_obj:   object for phase correction scan for reference data
% RTfeedback_obj:  object for realtime feedback data
% phasestab_obj:   object for phase stabilization scan
%
%
% The raw data can be obtained by calling e.g. image_obj() or for squeezed 
% data image_obj{''} (the '' are needed due to a limitation of matlab's 
% overloading capabilities). 
% Slicing is supported as well, e.g. image_obj(:,:,1,:) will return only 
% the first line of the full data set (all later dimensions are squeezed 
% into one). Thus, slicing of the "memory-mapped" data objects works 
% exactly the same as regular matlab array slicing - with one exception: 
% The keyword 'end' is not supported.
% Overloading of the '()' and '{}' operators works by overloading matlab's
% built-in 'subsref' function. Matlab calls subsref whenever the operators
% '()', '{}', or '.' are called. In the latter case, the overloaded subsref
% just calls the built-in subsref function since we don't want to change
% the behaviour for '.'-calls. However, this has one minor consequence:
% There's no way (afaik) to know whether the original call was terminated
% with a semicolon. Thus, a call to e.g. image_obj.NLin will produce no
% output with or without semicolon termination. 'a = image_obj.NLin" will 
% however produce the expected result.
%
%
% Order of raw data:
%  1) Columns
%  2) Channels/Coils
%  3) Lines
%  4) Partitions
%  5) Slices
%  6) Averages
%  7) (Cardiac-) Phases
%  8) Contrasts/Echoes
%  9) Measurements
% 10) Sets
% 11) Segments
% 12) Ida
% 13) Idb
% 14) Idc
% 15) Idd
% 16) Ide
%
%
% Optional parameters/flags:
% 
% removeOS:          removes oversampling (factor 2) in read direction
% doAverage:         performs average (resulting avg-dim has thus size 1)
% ignoreSeg:         ignores segment mdh index (works basically the same as 
%                    the average flag)
% noWeightedAverage: disables proper weighting of averages (i.e. data is 
%                    summed up instead of averaged)
% doRawDataCorrect:  enables raw data correction if used in the acquisition
%                    (only works for VB atm)
%
% These flags can also be set/unset later, e.g "image_obj.flagRemoveOS = 1"
%
%
% Examples:
%   [image_obj noise_obj phasecor_obj] = mapVBVD(measID);
%   
%   % return all image-data:
%   image_data = image_obj();
%   % return all image-data with all singular dimensions removed/squeezed:
%   image_data = image_obj{''}; % '' necessary due to a matlab limitation
%   % return only data for line numbers 1 and 5; all dims higher than 4 are
%   % grouped into dim 5):
%   image_data = image_obj(:,:,[1 5],:,:);
%   % return only data for coil channels 2 to 6; all dims higher than 4 are
%   % grouped into dim 5); but work with the squeezed data order 
%   % => use '{}' instead of '()':                                      
%   image_data = image_obj{:,2:6,:,:,:};
%
%   So basically it works like regular matlab array slicing (but 'end' is
%   not supported).
%

if ischar(filename) 
    % assume that complete path is given
    if  ~strcmpi(filename(end-3:end),'.dat');
        filename=[filename '.dat'];   %% adds filetype ending to file
    end
else 
    
    % filename not a string, so assume that it is the MeasID
    filename=ls(['meas_MID' num2str(filename) '_*.dat']);
    if isunix
        filename=filename(1:end-1); % crops line break
    end
end


% add absolute path, when no path is given
[pathstr, name, ext] = fileparts(filename);

if isempty(pathstr)
    pathstr  = pwd;
    filename = fullfile(pathstr, [name ext]);
end

%%%%% Parse varargin %%%%%

    % Definition of default parameters
    arg.removeOS          = false;
    arg.doAverage         = false;
    arg.ignoreSeg         = false;
    arg.noWeightedAverage = false;
    arg.bReadImaScan      = true;
    arg.bReadNoiseScan    = false;
    arg.bReadPCScan       = false;
    arg.bReadRefScan      = false;
    arg.bReadRefPCScan    = false;
    arg.bReadRTfeedback   = false;
    arg.bReadPhaseStab    = false;
    arg.doRawDataCorrect  = false; %SRY
    
    if (nargout>1)
        arg.bReadNoiseScan  = true;
    end
    if (nargout>2)
        arg.bReadPCScan     = true;
    end
    if (nargout>3)
        arg.bReadRefScan    = true;
    end
    if (nargout>4)
        arg.bReadRefPCScan  = true;
    end
    if (nargout>5)
        arg.bReadRTfeedback = true;
    end
    if (nargout>6)
        arg.bReadPhaseStab  = true;
    end
    
    k=1;

    while k <= numel(varargin)
        
        if ~ischar(varargin{k})
            error('string expected');
        end
        
        switch lower(varargin{k})
            case {'removeos','rmos'}
                if numel(varargin) > k && ~ischar(varargin{k+1})
                    arg.removeOS = logical(varargin{k+1});
                    k = k+2;
                else
                    arg.removeOS = true;
                    k = k+1;
                end
            case {'doaverage','doave','ave','average'}
                if numel(varargin) > k && ~ischar(varargin{k+1})
                    arg.doAverage = logical(varargin{k+1});
                    k = k+2;
                else
                    arg.doAverage = true;
                    k = k+1;
                end
            case {'ignseg','ignsegments','ignoreseg','ignoresegments'}
                if numel(varargin) > k && ~ischar(varargin{k+1})
                    arg.ignoreSeg = logical(varargin{k+1});
                    k = k+2;
                else
                    arg.ignoreSeg = true;
                    k = k+1;
                end
            case {'noweightedaverage','noweightedave','noweighting'}
                if numel(varargin) > k && ~ischar(varargin{k+1})
                    arg.noWeightedAverage = logical(varargin{k+1});
                    k = k+2;
                else
                    arg.noWeightedAverage = true;
                    k = k+1;
                end
            case {'rawdatacorrect','dorawdatacorrect'}
                %SRY: handle raw data correct arguments
                if numel(varargin) > k && ~ischar(varargin{k+1})
                    arg.doRawDataCorrect = logical(varargin{k+1});
                    k = k+2;
                else
                    arg.doRawDataCorrect = true;
                    k = k+1;
                end
            otherwise
                error('Argument not recognized.');
        end
    end
    clear varargin
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    tic;
    fid = fopen(filename,'r','l','US-ASCII');
    fseek(fid,0,'eof');
    fileSize = ftell(fid);
    
    % start of actual measurment data (sans header)
    fseek(fid,0,'bof');
    
    firstInt  = fread(fid,1,'uint32');
    secondInt = fread(fid,1,'uint32');
    
    % lazy software version check (VB or VD?)
    if and(firstInt < 10000, secondInt <= 64)
        version = 'vd';
        disp('Software version: VD (!?)');
        
        % number of different scans in file stored in 2nd int (wip, only 
        % one supported for now)
        NScans = secondInt;
        measID = fread(fid,1,'uint32');
        fileID = fread(fid,1,'uint32');
        % measOffset: points to beginning of header, usually at 10240 bytes
        measOffset = fread(fid,1,'uint64');
        measLength = fread(fid,1,'uint64');
        fseek(fid,measOffset,'bof');
        hdrLength  = fread(fid,1,'uint32');
        datStart   = measOffset + hdrLength;
    else
        % in VB versions, the first 4 bytes indicate the beginning of the
        % raw data part of the file
        version  = 'vb';
        disp('Software version: VB (!?)');
        datStart = firstInt;
        NScans   = 1; % VB does not support multiple scans in one file
    end
    
    %SRY read data correction factors
    % do this for all VB datasets, so that the factors are available later
    % in the image_obj if the user chooses to set the correction flag
    if (strcmp(version, 'vb')) % not implemented/tested for vd, yet
      frewind(fid);
      while ( (ftell(fid) < datStart) && ~exist('rawfactors', 'var'))
         line = fgetl(fid);
         %find the section of the protocol
         %note: the factors are also available in <ParamArray."CoilSelects">
         %along with element name and FFT scale, but the parsing is
         %significantly more difficult
         if (~isempty(strfind(line, '<ParamArray."axRawDataCorrectionFactor">')))
            while (ftell(fid) < datStart)
               line = fgetl(fid);
               %find the line with correction factors
               %the factors are on the first line that begins with this
               %pattern
               if (~isempty(strfind(line, '{ {  { ')))
                  line = strrep(line, '}  { } ', '0.0');
                  line = strrep(line, '{', '');
                  line = strrep(line, '}', '');
                  rawfactors = textscan(line, '%f');
                  rawfactors = rawfactors{1}; %textscan returns a 1x1 cell array
                  % this does not work in this location because the MDHs
                  % have not been parsed yet
                  %                    if (length(rawfactors) ~= 2*max(image_obj.NCha))
                  %                       error('Number of raw factors (%f) does not equal channel count (%d)', length(rawfactors)/2, image_obj.NCha);
                  %                    end;
                  if (mod(length(rawfactors),2) ~= 0)
                     error('Error reading rawfactors');
                  end;
                  %note the transpose, this makes the vector
                  %multiplication during the read easier
                  arg.rawDataCorrectionFactors = rawfactors(1:2:end).' + 1i*rawfactors(2:2:end).';
                  break;
               end
            end
         end
      end
      disp('Read raw data correction factors');
    end
   
   
    % data will be read in two steps (two while loops):
    %   1) reading all MDHs to find maximum line no., partition no.,... for
    %      ima, ref,... scan
    %   2) reading the data
    
    % declare data objects:
    image_obj      = twix_map_obj(arg,'image',filename,version);
    noise_obj      = twix_map_obj(arg,'noise',filename,version);
    phasecor_obj   = twix_map_obj(arg,'phasecor',filename,version);
    refscan_obj    = twix_map_obj(arg,'refscan',filename,version);
    refscanPC_obj  = twix_map_obj(arg,'refscan_phasecor',filename,version);
    RTfeedback_obj = twix_map_obj(arg,'rtfeedback',filename,version);
    phasestab_obj  = twix_map_obj(arg,'phasestab',filename,version);
    
    tic;
    mask.MDH_ACQEND = 0;
    cPos            = datStart;
    percentFinished = 0;
    frewind(fid);
    while ftell(fid)+128 < fileSize % fail-safe; in case we miss MDH_ACQEND

        switch version
            case 'vb'
                [mdh mask nBytes] = evalMDHvb(fid,cPos);
            case 'vd'
                [mdh mask nBytes] = evalMDHvd(fid,cPos);
            otherwise
                disp('error: only vb/vd software versions supported');                    
        end
        
        if mask.MDH_ACQEND
            break;
        end
        
        if (mask.MDH_IMASCAN && arg.bReadImaScan)
            image_obj.readMDH(mdh,cPos);
        end
        
        if (mask.MDH_NOISEADJSCAN && arg.bReadNoiseScan)
            noise_obj.readMDH(mdh,cPos);
        end

        if (and(mask.MDH_PHASCOR,~mask.MDH_PATREFSCAN) && arg.bReadPCScan)
            phasecor_obj.readMDH(mdh,cPos);
        end
        
        if (and(~mask.MDH_PHASCOR,(mask.MDH_PATREFSCAN || mask.MDH_PATREFANDIMASCAN)) && arg.bReadRefScan)
            refscan_obj.readMDH(mdh,cPos);
        end

        if (and(mask.MDH_PATREFSCAN,mask.MDH_PHASCOR) && arg.bReadRefPCScan)
            refscanPC_obj.readMDH(mdh,cPos);
        end

        if ((mask.MDH_RTFEEDBACK || mask.MDH_HPFEEDBACK) && arg.bReadRTfeedback)
            RTfeedback_obj.readMDH(mdh,cPos);
        end

        if ((mask.MDH_PHASESTABSCAN || mask.MDH_REFPHASESTABSCAN) && arg.bReadPhaseStab)
            phasestab_obj.readMDH(mdh,cPos);
        end
        
        % jump to mdh of next scan
        cPos = cPos + nBytes;
        
        if (cPos/fileSize*100 > percentFinished + 1)
            percentFinished = floor(cPos/fileSize*100);
            elapsed_time  = toc;
            time_left     = (fileSize/cPos-1) * elapsed_time;

            if ~exist('progress_str','var')
                prevLength = 0;
            else
                prevLength = numel(progress_str);
            end

            progress_str = sprintf('%3.0f %% parsed in %4.0f s; estimated time left: %4.0f s \n',...
            percentFinished,elapsed_time, time_left);

            fprintf([repmat('\b',1,prevLength) '%s'],progress_str);
        end
        
    end % while         
    
    if arg.bReadImaScan
        image_obj.clean();
    end
   
    if arg.bReadNoiseScan
        noise_obj.clean();
    end
    
    if arg.bReadPCScan
        phasecor_obj.clean();
    end
    
    if arg.bReadRefScan
        refscan_obj.clean();
    end
    
    if arg.bReadRefPCScan
        refscanPC_obj.clean();
    end
    
    if arg.bReadRTfeedback
        RTfeedback_obj.clean();
    end
    
    if arg.bReadPhaseStab
        phasestab_obj.clean();
    end
   
    elapsed_time = toc;
    progress_str = sprintf('100 %% parsed in %4.0f s; estimated time left:    0 s \n', elapsed_time);
    fprintf([repmat('\b',1,prevLength) '%s'],progress_str);
end


function [mdh mask nBytes cur_read] = evalMDHvb(fid,cPos)

    % no difference between 'scan' and 'channel' header in VB
    szMDH = 128; % [bytes]

    % inlining of readMDH
    fseek(fid,cPos+20,'bof');
    mdh.aulEvalInfoMask            = fread(fid,  [1 2], 'uint32');
    dummy                          = fread(fid,      2, 'uint16');
    mdh.ushSamplesInScan           = dummy(1);
    mdh.ushUsedChannels            = dummy(2);
    mdh.sLC                        = fread(fid, [1 14], 'ushort');  
    dummy                          = fread(fid,     18, 'uint16');
%     mdh.sCutOff                    = dummy(1:2);
    mdh.ushKSpaceCentreColumn      = dummy(3);
%     mdh.ushCoilSelect              = dummy(4);
    mdh.ushKSpaceCentreLineNo      = dummy(9);
    mdh.ushKSpaceCentrePartitionNo = dummy(10);
    mdh.aushIceProgramPara         = dummy(11:14);
    mdh.aushFreePara               = dummy(15:18);
    mdh.SlicePos                   = fread(fid,      7, 'float');
    
    % inlining of evalInfoMask
    mask.MDH_ACQEND             = min(bitand(mdh.aulEvalInfoMask(1), 2^0),1);
    mask.MDH_RTFEEDBACK         = min(bitand(mdh.aulEvalInfoMask(1), 2^1),1);
    mask.MDH_HPFEEDBACK         = min(bitand(mdh.aulEvalInfoMask(1), 2^2),1);
    mask.MDH_REFPHASESTABSCAN   = min(bitand(mdh.aulEvalInfoMask(1), 2^14),1);
    mask.MDH_PHASESTABSCAN      = min(bitand(mdh.aulEvalInfoMask(1), 2^15),1);
    mask.MDH_PHASCOR            = min(bitand(mdh.aulEvalInfoMask(1), 2^21),1);
    mask.MDH_PATREFSCAN         = min(bitand(mdh.aulEvalInfoMask(1), 2^22),1);
    mask.MDH_PATREFANDIMASCAN   = min(bitand(mdh.aulEvalInfoMask(1), 2^23),1);
    mask.MDH_REFLECT            = min(bitand(mdh.aulEvalInfoMask(1), 2^24),1);
    mask.MDH_RAWDATACORRECTION  = min(bitand(mdh.aulEvalInfoMask(1), 2^10),1);
    mask.MDH_SIGNREV            = min(bitand(mdh.aulEvalInfoMask(1), 2^17),1);
    mask.MDH_NOISEADJSCAN       = min(bitand(mdh.aulEvalInfoMask(1), 2^25),1);
    mask.MDH_IMASCAN            = 1;
    
    if (mask.MDH_ACQEND || mask.MDH_RTFEEDBACK || mask.MDH_HPFEEDBACK...
            || mask.MDH_REFPHASESTABSCAN || mask.MDH_PHASESTABSCAN...
            || mask.MDH_PHASCOR          || mask.MDH_NOISEADJSCAN)
        mask.MDH_IMASCAN = 0; 
    end
    
    if ( mask.MDH_PATREFSCAN && ~mask.MDH_PATREFANDIMASCAN )
        mask.MDH_IMASCAN = 0;
    end
    
    % size of current data set (2*4 because of complex + float)
    nBytes = mdh.ushUsedChannels * (szMDH + 2*4*mdh.ushSamplesInScan);

    % nothing to skip
    cur_read.skip = 0;
    
    % size for fread
    cur_read.sz   = [2 nBytes/(2*4)];
    
    % reshape size
    cur_read.shape = [mdh.ushSamplesInScan + szMDH/8, mdh.ushUsedChannels];
    
    % we need to cut MDHs from fread data
    cur_read.cut   = szMDH/8+1:mdh.ushSamplesInScan+szMDH/8;

end
    

function [mdh mask nBytes cur_read] = evalMDHvd(fid,cPos)

    % we need to differentiate between 'scan header' and 'channel header'
    % since these are used in VD versions:
    szScanHeader    = 192; % [bytes]
    szChannelHeader = 32;  % [bytes]

    % inlining of readScanHeader
    fseek(fid,cPos+40,'bof');
    mdh.aulEvalInfoMask            = fread(fid,  [1 2], 'uint32');
    dummy                          = fread(fid,      2, 'uint16');
    mdh.ushSamplesInScan           = dummy(1);
    mdh.ushUsedChannels            = dummy(2);
    mdh.sLC                        = fread(fid, [1 14], 'ushort');
    dummy                          = fread(fid,     10, 'uint16');
%     mdh.sCutOff                    = dummy(1:2);
    mdh.ushKSpaceCentreColumn      = dummy(3);
%     mdh.ushCoilSelect              = dummy(4);
    mdh.ushKSpaceCentreLineNo      = dummy(9);
    mdh.ushKSpaceCentrePartitionNo = dummy(10);
    mdh.SlicePos                   = fread(fid,      7, 'float');
    dummy                          = fread(fid,     28, 'uint16');
    mdh.aushIceProgramPara         = dummy(1:24);
    mdh.aushFreePara               = dummy(25:28); % actually aushReservedPara;
                                                   % there's no freePara in VD
  
    % inlining of evalInfoMask
    mask.MDH_ACQEND             = min(bitand(mdh.aulEvalInfoMask(1), 2^0),1);
    mask.MDH_RTFEEDBACK         = min(bitand(mdh.aulEvalInfoMask(1), 2^1),1);
    mask.MDH_HPFEEDBACK         = min(bitand(mdh.aulEvalInfoMask(1), 2^2),1);
    mask.MDH_REFPHASESTABSCAN   = min(bitand(mdh.aulEvalInfoMask(1), 2^14),1);
    mask.MDH_PHASESTABSCAN      = min(bitand(mdh.aulEvalInfoMask(1), 2^15),1);
    mask.MDH_PHASCOR            = min(bitand(mdh.aulEvalInfoMask(1), 2^21),1);
    mask.MDH_PATREFSCAN         = min(bitand(mdh.aulEvalInfoMask(1), 2^22),1);
    mask.MDH_PATREFANDIMASCAN   = min(bitand(mdh.aulEvalInfoMask(1), 2^23),1);
    mask.MDH_REFLECT            = min(bitand(mdh.aulEvalInfoMask(1), 2^24),1);
    mask.MDH_RAWDATACORRECTION  = min(bitand(mdh.aulEvalInfoMask(1), 2^10),1);
    mask.MDH_SIGNREV            = min(bitand(mdh.aulEvalInfoMask(1), 2^17),1);
    mask.MDH_NOISEADJSCAN       = min(bitand(mdh.aulEvalInfoMask(1), 2^25),1);
    mask.MDH_IMASCAN            = 1;
    
    if (mask.MDH_ACQEND || mask.MDH_RTFEEDBACK || mask.MDH_HPFEEDBACK...
        || mask.MDH_REFPHASESTABSCAN || mask.MDH_PHASESTABSCAN...
        || mask.MDH_PHASCOR          || mask.MDH_NOISEADJSCAN)
        mask.MDH_IMASCAN = 0; 
    end
    
    if ( mask.MDH_PATREFSCAN && ~mask.MDH_PATREFANDIMASCAN )
        mask.MDH_IMASCAN = 0;
    end
    
    nBytes = szScanHeader + mdh.ushUsedChannels * (szChannelHeader + 2*4*mdh.ushSamplesInScan);
    
    % skip line header
	cur_read.skip = szScanHeader;
    
    % size for fread
    cur_read.sz   = [2 (nBytes-szScanHeader)/8];
    
    % reshape size
    cur_read.shape = [mdh.ushSamplesInScan + szChannelHeader/8, mdh.ushUsedChannels];
    
    % we need to cut MDHs from fread data
    cur_read.cut   = szChannelHeader/8+1:mdh.ushSamplesInScan+szChannelHeader/8;
    
end
