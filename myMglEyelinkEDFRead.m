% myMglEyelinkEDFRead.m
% slightly modified version of mglEyelinkEDFRead.m
%
%      usage: mglEyelinkEDFRead(filename,<verbose>)
%         by: justin gardner
%       date: 04/04/10
%    purpose: Function to read EyeLink eye-tracker files into matlab
%
function retval = myMglEyelinkEDFRead(filename,verbose)

% default return argument
retval = [];

% check arguments
if ~any(nargin == [1 2])
  help mglEyelinkEDFRead
  return
end

% default arguments
if nargin < 2,verbose = 1;end

% check for compiled file
if exist('mglPrivateEyelinkEDFRead')~=3
  disp(sprintf('(mglEyelinkEDFRead) You must compile the eyelink files: mglMake(''Eyelink'')'));
  return
end

[p,n,e] = fileparts(filename);
if isempty(e)
    filename = fullfile(p, [n '.edf']);
end
if isfile(filename)
  % mglPrivateEleyinkReadEDF returns a structre
  retval = mglPrivateEyelinkEDFRead(filename,verbose);
  if isempty(retval),return,end
else
  disp(sprintf('(mglEyelinkEDFRead) Could not open file %s',filename));
end

%% let's parse some additional info
% this could be parsed to provide information about the calibration quality
retval.cal = char(strtrim({retval.messages(strmatch('!CAL',{retval.messages.message})).message}));
% mode
retval.mode = strtrim({retval.messages(strmatch('!MODE',{retval.messages.message})).message});
[t,m] = strtok(retval.mode); % should be !MODE
[t,m] = strtok(m); % should be RECORD
if ~strcmp(t,'RECORD')
    warning('mglEyelinkEDFRead:UnknownMode', 'Unknown mode encountered in edf file.');
end
[retval.trackmode,m] = strtok(m); % will be CR or P? (pupil only)
[t,m] = strtok(m); % sample rate
% this is the true sample rate.
retval.samplerate = str2double(t);
[t,m] = strtok(m); % filer mode
retval.filter = str2double(t);
[t,m] = strtok(m); % number of eyes
retval.numeye = str2double(t);
[t,m] = strtok(m); % which eye


if strcmp(t,'LR')
    retval.whicheye = 'Both';
    retval.gaze = retval.gazeLeft;
    retval.gazeBoth{1} = retval.gazeLeft;
    retval.gazeBoth{2} = retval.gazeRight;
    retval = rmfield(retval,'gazeLeft');
    retval = rmfield(retval,'gazeRight');
    disp(sprintf('(mglEyelinkEDFRead) !!! Both eyes were recorded. Setting gaze variable to left eye !!!'));
elseif retval.numeye == 2
    retval.whicheye = 'Both';
    disp(sprintf('(mglEyelinkEDFRead) !!! Both eyes were recorded. Setting gaze variable to left eye !!!'));
%     retval.gaze = retval.gazeLeft;
    retval.gaze = retval.gazeRight;
elseif strcmp(t,'R')
    retval.whicheye = 'Right';
    retval.gaze = retval.gazeRight;
    % remove left and right gaze
    retval = rmfield(retval,'gazeLeft');
    retval = rmfield(retval,'gazeRight');
else
    retval.whicheye = 'Left';
    retval.gaze = retval.gazeLeft;
    % remove left and right gaze
    retval = rmfield(retval,'gazeLeft');
    retval = rmfield(retval,'gazeRight');
end

return


