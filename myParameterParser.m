classdef myParameterParser < inputParser
    % inputParser expects the required arguments in be explicitly passed,
    % and not as a parameter/value pair. I think this is annoying, so this
    % is a wrapper / work around... explained better in the link.
    %
    % Code modified from: dumbmatter
    % REF: http://stackoverflow.com/questions/9221690/matlab-using-inputparser-with-varargin
    properties
        required = {};
    end

    methods
        function obj = myParameterParser(varargin)
            obj = obj@inputParser(varargin{:});
        end

        function addRequired(obj, argname, validator)
            obj.required = [obj.required, {argname}];
            obj.addOptional(argname, NaN, validator);
        end

        function parse(obj, varargin)
            params_input = varargin(1:2:end);
            % Look for all the values of the required parameters
            for i = 1:length(obj.required)
                if isempty(validatestring(obj.required{i}, params_input))
                    error( 'Required named parameter %s was not passed to function', obj.required{i} );
                end
            end
            parse@inputParser(obj, varargin{2},varargin{3:end});
            
        end

    end
end