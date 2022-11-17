% MOVCORR Moving Pearson product-moment correlation coefficient r. This is
% basically a wrapper to MOVSUM and the low-memory overhead computation of r.
% See second formula in 
%   https://en.wikipedia.org/wiki/Pearson_correlation_coefficient#For_a_sample
%
% @Synopsis
%   r = MOVCORR(x, y, k) computes r for x and y over a centered window specified
%   by k.
%
%   MOVCORR(x, y, k, missing) specifies how to treats NaNs.
%
%   MOVSUM(...,'Endpoints',endpoints) controls how the start/end of the vectors
%   are handled.
%
%   [r, p] = MOVSUM(...) also returns the p-values for each computed r.
%
%   [r, p] = MOVSUM(...,'Tail',tail) controls what kind of hypothesis test is
%   used to compute the p-values.
%
% @Input
%   x, y (required)
%     m x 1 column vectors with the data to correlate.
%
%   k (required)
%      Positive integer scalar with the window size or 1x2 array of non-negative
%      integer scalars [nB, nF] with backward & forward samples for the window.
%      If k is a scalar and is even, the window is centered on the current &
%      previous sample, e.g. ranges from i-k/2 to i+k/2-1. If k is a vector the
%      window ranges from i-nB to i+nF.
%
%   missing (optional)
%     Flag of what to do with missing values. Either 'includenan' or 'omitnan'.
%     For 'omitnan' the r values are computed pairwise, resembling the behavior
%     of CORR(...,'Rows','pairwise') - but of course much faster.
%     DEFAULT: 'includenan'
%
%   'Endpoints',endpoints (Name/value)
%     Specifies how to handle the start/end of the vectors. Either a numeric or
%     logical scalar or 'shrink', 'fill' or 'discard'. Be aware that using
%     arbitrary numeric values does not really make sense in this setting (as
%     the correlation coefficient is bound between -1 and 1). Moreover, when
%     using numeric/logical endpoints, the according p-values are set to NaN.
%     This result in unexpected NaNs: If missing=='omitnan' and endpoints is
%     numeric, then p will contain NaNs at the endpoints, despite r beeing
%     finite.
%     DEFAULT: 'shrink'
%
%   'Tail',tail (Name/value)
%     Specifies the hypothesis with which the p-values are computed. Either
%     'both', 'left' or 'right'. See CORR for more details.
%     DEFAULT: 'both'
%
%   'MinN',nMin (Name/value)
%     Positive integer scalar with the minimum number of non-NaN values in a
%     local window required for computation of r. If a window does not have the
%     required number of values, r is not computed (e.g. set to NaN). 'MinN' has
%     to be equal or smaller than the window length specified by k. This
%     option is ignored if missing=='includenan'. Setting 'MinN' > 1 when using 
%     'Endpoints'=='shrink' may cause the values at the borders to be NaN,
%     despite the 'shrink' option.
%     DEFAULT: 1 (always compute r)
%
% @Output
%   r
%     m x 1 column vector with the local Pearson correlation coefficients for
%     each mutual sample in x and y.
%
%   p
%     m x 1 column vector with the p-values for the local Pearson correlation
%     coefficients.
%
%   n
%     The number of non-NaN values in each local window. If
%     missing=='includenan', this is a scalar equal to the window size. If
%     missing=='omitnan', this is an m x 1 column vector with the number of
%     non-NaN elements in each local window. This can be used to compute the
%     degrees-of-freedom for partial correlations.
%
% @Remarks
%   - For more detailed info on the inputs & syntax see MOVSUM.
%
% @Dependencies
%   none
%
% See also CORR, MOVSUM
%
% ML-FEX ID:65342
%
% @Changelog
%   2017-12-07 (DJM): 
%     - [RELEASE] Golden.
%     - [MOD] Computation of r with NaNs + 'omitnan' now gives same results as
%       CORR(...,'Rows','pairwise').
%   2017-12-11 (DJM): 
%     - [ADD] Name-value input 'Tail' and according p-value output.
%     - [MOD] Now uses INPUTPARSER. Changes error messages for invalid inputs.
%     - [FIX] Erroneous window size computation when using k = [nB, nF].
%     - [FIX] Complex output when using numeric endpoints.
%   2017-12-12 (DJM): 
%     - [ADD] Output n.
%   2017-12-14 (DJM): 
%     - [ADD] Name-value input 'MinN'.
%     - [MOD] Now returns NaN for zero-variances, to be consistent with CORR.
%
function [r, p, n] = movcorr(x, y, k, varargin)
%% Parse inputs
[x, y, k, missing, endpoints, tail, nMin, n] ...
    = parseInputs(x, y, k, varargin);
isNumericEndpoints = ~ischar(endpoints);
if isNumericEndpoints
	% With numeric/logical endpoints, we set the method to 'fill' and set the
	% values in the end, to avoid problems when using MOVSUM on x & x^2.
	endpointVals = endpoints;
	endpoints    = 'fill';
end
%% Assemble lambda function for MOVSUM
movingSum = @(v) movsum(v, k, missing, 'Endpoints',endpoints);
%% Compute sums
sumX = movingSum(x);
sumY = movingSum(y);
%% Compute r
% Use SQRT in denominator on x & y separately to avoid under/overflow.
r = (n .* movingSum(x.*y) - sumX.*sumY) ...
 ./ (sqrt(n .* movingSum(x.^2) - sumX.^2) .* sqrt(n .* movingSum(y.^2) - sumY.^2));
%% Handle values > +/-1
% Use same method as CORR to be consistent. This will cause INFs to be NaN...
% which might be argued upon.
isExceeding = abs(r) > 1;
r(isExceeding) = r(isExceeding) ./ abs(r(isExceeding));
%% Remove invalid windows
if strcmp(missing, 'omitnan')
	r(n < nMin) = NaN; % n is vector if missing=='omitnan'.
end
%% Compute p-values
if nargout() > 1
	p = computePvalues(r, n, tail);
end
%% Handle endpoints
if isNumericEndpoints
	ids = [1 : k(1), length(r) - k(2) + 1 : length(r)];
	r(ids) = endpointVals;
	p(ids) = NaN;
end
end % MAINFUNCTION
%% SUBFUNCTION: parseInputs
% PARSEINPUTS Parses the inputs to the main function.
%
% @Input
%   x,y,k (required)
%     Required main function inputs.
%
%   varArgIn (required)
%     Optional main function inputs (VARARGIN).
%
% @Output
%   x,y,k
%     Required main function inputs.
%
%   missing
%     Optional main function input.
%
%   endpoints,tail,nMin
%     Name-value main function inputs.
%
%   n
%     Window length (based on k). If isAnyNan==true this is a mx1 vector with
%     the local window sizes for finite-values.
%
function [x, y, k, missing, endpoints, tail, nMin, n] = parseInputs(x, y, k,...
                                                                    varArgIn)
%% >SUB: Init
% Valid enumeration vals. First is used as default.
VALID_MISS    = {'includenan', 'omitnan'};
VALID_ENDPNTS = {'shrink', 'discard', 'fill'};
VALID_TAILS   = {'both', 'left', 'right'};
% Error messages - partly resembles errors from MOVSUM.
ERR_BASE_ID = 'DJM:IOFailure';
ERR_MSG_X_Y_DIM         = '''x'' and ''y'' have to be of equal size.';
ERR_MSG_K_INVALID       = '''k'' must be a finite scalar positive integer or 2-element vector of finite nonnegative integers.';
ERR_MSG_MISS_INVALID    = '''missing'' value must be ''omitnan'' or ''includenan''.';
ERR_MSG_ENDPNTS_INVALID = '''Endpoints'' value must be ''shrink'', ''discard'', ''fill'', or a numeric or logical scalar.';
ERR_MSG_TAIL_INVALID    = '''Tail'' value must be ''both'', ''left'' or ''right''.';
ERR_MSG_MIN_N_INVALID   = '''MinN'' has to be equal or smaller than the window length given by ''k''.';
ERR_MSG_MIN_N_IGNORED   = '''MinN'' is ignored with option ''includenan''.';
ERR_MSG_MIN_N_SHRINKED  = '''MinN'' may cause unexpected NaNs at the borders with option ''shrink''.';
%% >SUB:  Create INPUTPARSER
IP = inputParser();
IP.StructExpand  = true;
IP.KeepUnmatched = false;
%% >SUB: Set up validation functions
validateNumCol          = @(v)validateattributes(v,{'numeric'},{'column'});
validateNumScaFinIntPos = @(v)validateattributes(v,{'numeric'},{'scalar','finite','integer','positive'});
validateNumVecFinIntNN  = @(v)validateattributes(v,{'numeric'},{'vector','finite','integer','nonnegative'});
validateCharOrNumRow    = @(v)validateattributes(v,{'char','numeric','logical'},{'row'});
validateCharRow         = @(v)validateattributes(v,{'char'},{'scalartext'});
%% >SUB: Define inputs
IP.addRequired( 'x'        , validateNumCol);
IP.addRequired( 'y'        , validateNumCol);
IP.addRequired( 'k'        , validateNumVecFinIntNN);
IP.addParameter('Endpoints', VALID_ENDPNTS{1}, validateCharOrNumRow);
IP.addParameter('Tail'     , VALID_TAILS{1}  , validateCharRow);
IP.addParameter('MinN'     , 1               , validateNumScaFinIntPos);
%% >SUB: Parse 'missing' if given
% This is an optional input which can be ommitted and thus cannot be handled by
% INPUTPARSER.
isDefault = isempty(varArgIn) || ~ischar(varArgIn{1}) ...
         || ismember(varArgIn{1}, {'Endpoints', 'Tail', 'MinN'});
if isDefault
	missing = VALID_MISS{1}; % Use default.
else
	missing = varArgIn{1};
	varArgIn(1) = [];
	errId = [ERR_BASE_ID ':MissingInvalid'];
	assert(ismember(missing, VALID_MISS), errId, ERR_MSG_MISS_INVALID);
end
%% >SUB: Parse & assign
IP.parse(x, y, k, varArgIn{:});
endpoints = IP.Results.Endpoints;
tail      = IP.Results.Tail;
nMin      = IP.Results.MinN;
%% >SUB: Check k
% Convert scalar k to [nB,nF] notation.
if isscalar(k)
	k = floor(k / 2) + [0, -double(mod(k, 2) == 0)];
end
errId = [ERR_BASE_ID ':KInvalid'];
assert(numel(k) == 2 && all(k >= 0), errId, ERR_MSG_K_INVALID);
winSize = sum(k) + 1; % nB + nF + current sample
%% >SUB: Check size of x & y
errId = [ERR_BASE_ID ':XYInvalid'];
assert(isequal(size(x), size(y)), errId, ERR_MSG_X_Y_DIM);
%% >SUB: Check NaNs in x & y
if strcmp(missing, 'includenan')
    n = winSize;
else
	isNanX = isnan(x);
	isNanY = isnan(y);
	isNanXorY = isNanX | isNanY;
	if any(isNanXorY)
		% Set NaNs mutually in both vectors.
		x(isNanY) = NaN;
		y(isNanX) = NaN;
		n = movsum(~isNanXorY, k);
	end
end
%% >SUB: Check endpoints
errId = [ERR_BASE_ID ':EndpointsInvalid'];
if ischar(endpoints)
	isValid = ismember(endpoints,VALID_ENDPNTS);
else % Numeric or logical
	isValid = isscalar(endpoints);
end
assert(isValid, errId, ERR_MSG_ENDPNTS_INVALID);
   
%% >SUB: Check tail
errId = [ERR_BASE_ID ':TailInvalid'];
assert(ismember(tail, VALID_TAILS), errId, ERR_MSG_TAIL_INVALID);
%% >SUB: Check MinN
errId = [ERR_BASE_ID ':MinNInvalid'];
assert(nMin <= winSize, errId, ERR_MSG_MIN_N_INVALID);
if ~ismember('MinN', IP.UsingDefaults)
    if strcmp(missing, 'includenan')
        errId = [ERR_BASE_ID ':MinNIgnored'];
        warning(errId, ERR_MSG_MIN_N_IGNORED);
    elseif nMin > 1 && strcmp(endpoints,'shrink') % & missing==omitnan
        errId = [ERR_BASE_ID ':MinNCausingNaNs'];
        warning(errId, ERR_MSG_MIN_N_SHRINKED);
    end
end
end % PARSEINPUTS
%% SUBFUNCTION: computePvalues
% COMPUTEPVALUES Computes the p-values for the r's.
%
% @Input
%   r, n, tail (required)
%     Main function variables.
%
% @Output
%   p
%     Array same size as r with p-values.
%
function p = computePvalues(r, n, tail)
%% >SUB: Compute t-statistics
% t-statistic for r from Rudolf2008 - Biostatistik, p.213 eq. 7.9.
dof = n - 2; % Degrees-of-freedom.
t   = r .* sqrt(dof ./ (1 - r.^2));
%% >SUB: Compute p-values
if strcmp(tail, 'both')
	% Two-sided.
	p = 2 * tcdf(-abs(t), dof); % == 2 * (1 - tcdf(-abs(t), dof))
else
	% One-sided: left or right.
	if strcmp(tail, 'right')
		t = -t; % 1 - tcdf(t, dof) == tcdf(-t, dof)
	end
	p = tcdf(t, dof);
end
end % COMPUTEPVALUES