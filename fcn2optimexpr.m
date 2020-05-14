function varargout = fcn2optimexpr(func, varargin)
%FCN2OPTIMEXPR Create an optimization expression from a function
% 
%   EXPR = FCN2OPTIMEXPR(FCN, X1, ..., XN) creates an optimization
%   expression from a function. This is an optimization expression that
%   wraps the function given by the handle or function name, FCN, at the
%   given arguments, X1, ..., XN. Xi can be optimization variables,
%   expressions or any other MATLAB datatype.
% 
%   EXPR = FCN2OPTIMEXPR(FCN, X1, ..., XN, PARAM, VAL, ...) creates the
%   optimization expression as above with specified values of optional
%   parameters:
%
%   'OutputSize', size of the expression, EXPR. 
%
%      When the optimization function expression is created, FCN will be
%      evaluated with dummy data to determine the size of EXPR. If you know
%      the size of EXPR, you can set OutputSize to avoid this evaluation.
%
%   'ReuseEvaluation', specifies whether to reuse common function evaluations,
%                      Boolean, default false 
%
%      If you know that FCN is expensive, set ReuseEvaluation to true. This
%      allows common evaluations of FCN to be reused where possible.
%      
%   [EXPR1, ..., EXPRM] = FCN2OPTIMEXPR(FCN, ...) returns multiple
%   expressions from the specified function.
%
% Examples: 
%    
%      % Given a function, myfun(x, y), create an optimization expression:
%      x = optimvar('x');
%      y = optimvar('y');
%      expr = fcn2optimexpr(@myfun, x, y);
% 
%      % Create an optimization expression with data as the second argument and set the 
%      % output size:
%      x = optimvar('x', 2, 2);
%      expr = fcn2optimexpr(@peaks, x, ones(2, 2), 'OutputSize', [2 2]);
% 
%      % Create an optimization expression from an expensive function, 
%      myexpensive(x, y, z)
%      x = optimvar('x', 2, 1); y = optimvar('y', 2, 1);
%      z = 2*x + 3*y;
%      expr = fcn2optimexpr(@myexpensive, x, y, z, 'ReuseEvaluation', true);
%  
%      % Create an optimization expression from the first two output arguments of eig
%      x = optimvar('x', 3, 3);
%      [e1, e2] = fcn2optimexpr(@eig, x);
%
%   See also OPTIMPROBLEM, OPTIMVAR, OPTIMEXPR

%   Copyright 2018 The MathWorks, Inc.

%%% Start validation %%%

% Check that func is a function handle or anonymous function.
if ~isa(func, 'function_handle')
    error(message('optim_problemdef:fcn2optimexpr:FirstArgNotFcnHdl'));
end

% Return at least one output
nout = max(1, nargout);

% Parse the optional params from the right end. The number of
% data inputs is what's left on the left end. A char row vector
% or a scalar string that is intended to be a data input may
% be interpreted as a parameter name in unlucky cases.
pnames = {'OutputSize' 'ReuseEvaluation'};
dflts =  {    []        false};
partialMatchPriority = [0 0];
[NumVars, OutputSize, ReuseEvaluation, supplied] ...
    = matlab.internal.datatypes.reverseParseArgs(pnames,dflts,partialMatchPriority,varargin{:});

% Validate the function inputs 
fcnInputs = varargin(1:NumVars);
% Determine the number of inputs to the function. If nFcnInputs >= 0, we
% can additonally check that the correct number of inputs has been passed
% to fcn2optimexpr
try    
    nFcnInputs = nargin(func);
catch ME
    % It is possible for valid functions to throw an error when passed to
    % nargin, e.g. MEX functions. In this case, we'll return nFcnInputs =
    % -1 and skip the number of inputs check. Otherwise we'll throw the
    % nargin error - note we throw the error, as opposed to rethrow, to
    % make the error appear to originate from fcn2optimexpr.
    if strcmp(ME.identifier, 'MATLAB:narginout:doesNotApply')
        nFcnInputs = -1;
    else
        throw(ME);
    end
end
validateFcnInputs(fcnInputs, nFcnInputs);

% Validate the optional parameters
OutputSize = validateOptionalParameters(OutputSize, ReuseEvaluation, supplied, nout);

% Warn if a parameter name has appeared to be used as an input
warnIfAmbiguousParamName(fcnInputs, pnames)

% Create an expression from the function
[varargout{1:nout}] = createFunctionExpression(func, ...
    fcnInputs, OutputSize, supplied, ReuseEvaluation);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Helper function to check that the input sizes specified are valid
function cellSz = checkSizes(cellSz, nout)

if ~iscell(cellSz)
    % The user provided a single size
    % Convert to cell for the common checks below
    cellSz = {cellSz};
end
nSizes = numel(cellSz);
if nSizes < 1 % empty cell input
    throwAsCaller(MException('optim_problemdef:fcn2optimexpr:OutputSizeEmptyCell', ...
        getString(message('optim_problemdef:fcn2optimexpr:OutputSizeEmptyCell'))));
elseif nSizes == 1
    % If there is a single size for multiple outputs, expand
    cellSz = repmat(cellSz, 1, nout);
elseif nSizes >  nout % too many size inputs in cell
    throwAsCaller(MException('optim_problemdef:fcn2optimexpr:OutputSizeCellTooManyElts', ...
        getString(message('optim_problemdef:fcn2optimexpr:OutputSizeCellTooManyElts'))));
elseif nSizes < nout % too few size inputs
    throwAsCaller(MException('optim_problemdef:fcn2optimexpr:OutputSizeCellTooFewElts', ...
        getString(message('optim_problemdef:fcn2optimexpr:OutputSizeCellTooFewElts'))));
end

% Checks for every size array
for i = 1:nSizes
    szI = cellSz{i};

    % Must be numeric
    if ~isnumeric(szI)
        throwAsCaller(MException('optim_problemdef:fcn2optimexpr:OutputSizeNotNumeric', ...
            getString(message('optim_problemdef:fcn2optimexpr:OutputSizeNotNumeric'))));            
    end
    
    % Must be a size vector, i.e have at least two elements
    if numel(szI) < 2
        throwAsCaller(MException('optim_problemdef:fcn2optimexpr:OutputSizeNotEnoughElts', ...
            getString(message('MATLAB:getReshapeDims:sizeVector'))));
    end
    
    % Size must be a row vector
    if ~isrow(szI)
        throwAsCaller(MException('optim_problemdef:fcn2optimexpr:OutputSizeNotRowVector', ...
            getString(message('MATLAB:checkDimRow:rowSize'))));
    end
    
    % All elements of size must be nonnegative
    if any(szI < 0)
        throwAsCaller(MException('optim_problemdef:fcn2optimexpr:OutputSizeMustBeNonnegative', ...
            getString(message('MATLAB:checkDimCommon:nonnegativeSize'))));
    end
    
    % All elements of size vector must be real
    if ~isreal(szI)
        throwAsCaller(MException('optim_problemdef:fcn2optimexpr:OutputSizeNotReal', ...
            getString(message('MATLAB:checkDimCommon:complexSize'))));
    end
    
    % All elements of size must be real, finite, integers
    if any(~isfinite(szI)) || any(floor(szI) ~= szI)
        throwAsCaller(MException('optim_problemdef:fcn2optimexpr:OutputSizeNotInteger', ...
            getString(message('MATLAB:checkDimRow:rowSize'))));
    end
    
end

function varargout = createFunctionExpression(func, inputs, outSize, supplied, ReuseEvaluation)
%CREATEFUNCTIONEXPRESSION Create N OptimizationExpressions for a black-box
%                         function.
%
% EOUT = CREATEFUNCTIONEXPRESSION(FUNC, INPUTS) creates N
% OptimizationExpressions representing the operation EOUT = FUNC(INPUTS{:})
%        FUNC: Function handle
%      INPUTS: Cell array of inputs to the function


% Determine the size of each output.
if ~supplied.OutputSize
    evalOut = cell(1, nargout);
    try
        % Generate an initial point to evaluate the user-function at so
        % that we can get the size information.
        evalInputs = cellfun(@generateInputPoint, inputs, 'UniformOutput', false);
        [evalOut{:}] = func(evalInputs{:});        
    catch
        throwAsCaller(MException('optim_problemdef:fcn2optimexpr:FcnError', ...
            getString(message('optim_problemdef:fcn2optimexpr:FcnError'))))
    end
    outSize = cellfun(@size, evalOut, 'UniformOutput', false);
end

% Wrap the function handle and its inputs and returns the variables vars in
% the black-box expression.
[optimFunc, vars, depth] = optim.internal.problemdef.FunctionWrapper.createFunctionWrapper(func, inputs, nargout, ReuseEvaluation);

% Create the outputs
varargout = cell(1,nargout);
for i = 1:nargout
    outputi = optim.problemdef.OptimizationExpression();
    % Create a nonlinear expression node for the ith output
    outputi = createFunction(outputi, optimFunc, vars, depth, outSize{i}, i);
    varargout{i} = outputi;
end

% Helper function to generate an initial point for every input of the user
% function.
function evalArg = generateInputPoint(arg)
if ~isa(arg, 'optim.problemdef.OptimizationExpression')
    % If the input is standard MATLAB data, use it directly
    evalArg = arg;
else
    % Otherwise, the input is of type OptimizationExpression, so we need to
    % generate a double array to pass to the user function.
    % Get the expression's variables
    vars = getVariables(arg);
    % Generate an initial point structure to evaluate the expression at.
    initPt = vars;
    varnames = fieldnames(vars);
    % Loop through the variables in the expression to generate a valid
    % point.
    for i = 1:numel(varnames)
        % Current variable
        name = varnames{i};
        thisVar = vars.(name);
        % Generate a value for the variable based on its bounds.
        % If there are no bounds, we select 1 + eps as the default value.
        % If there are only lower bounds, we select lb + max(1,abs(lb))*eps.
        % If there are only upper bounds, we select ub - max(1,abs(ub))*eps.
        % If there are both types of bounds, we select (lb + ub)/2 + ((ub - lb)/2)*eps.
        varValue = ones(size(thisVar)) + eps;
        % Get the bounds on the variable
        lb = thisVar.LowerBound;
        ub = thisVar.UpperBound;
        finiteLb = isfinite(lb);
        finiteUb = isfinite(ub);
        bothFinite = finiteLb & finiteUb;
        finiteLb = finiteLb & ~bothFinite;
        finiteUb = finiteUb & ~bothFinite;
        % Update the value for the finite bounds.
        varValue(finiteLb) = lb(finiteLb) + max(1,abs(lb(finiteLb)))*eps;
        varValue(finiteUb) = ub(finiteUb) - max(1,abs(ub(finiteUb)))*eps;
        varValue(bothFinite) = (lb(bothFinite) + ub(bothFinite))/2 + ((ub(bothFinite)-lb(bothFinite))/2)*eps;
        % Round the point if the variable is integer
        if strcmpi(thisVar.Type, 'integer')
            varValue = ceil(varValue);
        end
        % Add the variable value to the initial point structure.
        initPt.(name) = varValue;
    end
    evalArg = evaluate(arg, initPt);
end

function validateFcnInputs(inps, nFcnInputs)

FCNHASFIXEDINPUTS = nFcnInputs > -1;
if FCNHASFIXEDINPUTS && length(inps) > nFcnInputs
    throwAsCaller(MException('optim_problemdef:fcn2optimexpr:TooManyFcnInputs', ...
        getString(message('optim_problemdef:fcn2optimexpr:TooManyFcnInputs', nFcnInputs, length(inps)))));    
end


function OutputSize = validateOptionalParameters(OutputSize, ReuseEvaluation, supplied, nOut)

if supplied.OutputSize   
    try
        OutputSize = checkSizes(OutputSize, nOut);
    catch ME
        throwAsCaller(ME);
    end
end

if supplied.ReuseEvaluation
    % ReuseEvaluation must be boolean
    try
        validateattributes(ReuseEvaluation, {'logical'}, {'scalar'});
    catch 
        throwAsCaller(MException('optim_problemdef:fcn2optimexpr:ReuseEvaluationNotBoolean', ...
            getString(message('optim_problemdef:fcn2optimexpr:ReuseEvaluationNotBoolean'))));        
    end    
end

function warnIfAmbiguousParamName(vars, pnames)

% Check to see if a function input matches a parameter name. If so,
% warn about the ambiguity.
for k = 1:numel(vars)
    if matlab.internal.datatypes.isScalarText(vars{k})
        pNameMatch = strcmpi(vars{k},pnames);
        if any(pNameMatch)
            warning(message('optim_problemdef:fcn2optimexpr:ParamNameAmbiguity', pnames{find(pNameMatch,1)}));
            break;
        end
    end
end

