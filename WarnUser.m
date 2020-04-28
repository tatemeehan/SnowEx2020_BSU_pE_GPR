
function WarnUser(errorMessage,filename,workingDirectory)
	% Alert user via the command window and a popup message.
% 	fprint(1, '%s\n', errorMessage); % To command window.
% 	uiwait(warndlg(errorMessage));
	
	% Open the Error Log file for appending.
	fullFileName = fullfile(workingDirectory,filename);
	fid = fopen(fullFileName, 'at');
	fprintf(fid, '%s\n', errorMessage); % To file
	fclose(fid);
	return; % from WarnUser()
end