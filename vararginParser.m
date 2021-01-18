function vararg = vararginParser(vararg, vararg_in)
%VARARGINPARSER Input parser.

% How to use:
% % Set default 
% vararg = {'firstparameter', 1, 'secondparameter', magic(3)};
% % Overwrite by input
% vararg = vararginParser(vararg, varargin);
% % Generate variables
% for pair = reshape(vararg, 2, []) % pair is {propName;propValue}
%    eval([pair{1} '= pair{2};']);
% end

% Copyright Chong You @ Johns Hopkins University, 2016
% chong.you1987@gmail.com

% count arguments
if mod(length(vararg_in), 2) ~= 0
    error('varargin needs propertyName/propertyValue pairs')
end
% overwrite vararg.
optionNames = vararg(1:2:end);
for pair = reshape(vararg_in, 2, []) % pair is {propName;propValue}
    optName = pair{1};
    index = find( strcmpi(optName, optionNames) );
    if ~isempty(index)
        vararg{index * 2} = pair{2};
    else
        error('%s is not a recognized parameter name', optName)
    end
end

end