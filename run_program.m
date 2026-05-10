function run_program(program_name)
% Run a MATLAB script from the programs folder.
% Example:
%   run_program('MTE_1d')

if nargin < 1 || isempty(program_name)
    error('Please provide a program name, e.g. run_program(''MTE_1d'').');
end

root_dir = fileparts(mfilename('fullpath'));
programs_dir = fullfile(root_dir, 'programs');
addpath(programs_dir);

if isempty(regexp(program_name, '\\.m$', 'once'))
    program_name = [program_name, '.m'];
end

program_path = fullfile(programs_dir, program_name);
if exist(program_path, 'file') ~= 2
    error('Program not found: %s', program_path);
end

run(program_path);
end
