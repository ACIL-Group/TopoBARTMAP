% -----------------------------------------------------------------------------	%
% Description:									%
% 										%
% this file is an implementation of fuzzyNorm1 function. Which calculates city	%
% block distance for the given arguments.					%
% 										%
% test?	 ignore									%
% Takes a column vector as input. 						%
%------------------------------------------------------------------------------	%

function output = fuzzyNorm1(varargin)
	
	switch nargin
		
		case 1
			output= sum(abs(varargin{1}));
		case 2
			if size(varargin{1}) == size(varargin{2})
				output = sum(abs(min(varargin{1},varargin{2})));
            else
				%check if size of one of the two vectors is of size 1 across second dimension
				if size(varargin{1},2)==1
					temp=size(varargin{2},2);
					output = sum(abs(min(varargin{1}*ones(1,temp),varargin{2})));
				elseif size(varargin{2},2)==1
					temp=size(varargin{1},2);
					output = sum(abs(min(varargin{2}*ones(1,temp), varargin{1})));
                else	
                    disp('I am here');
					error('Error:The given input to the fuzzyNorm1 function is not supported');
				end
			end
		case 0
			disp('no input given to fuzzyNorm1');
			output=0;
		otherwise
			disp('more than 2 inputs not supported');
			output=0;	
	end
end

% EOF
