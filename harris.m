%% Harris corner detection (COPYRIGHT see below)
% downloaded from http://www.csse.uwa.edu.au/~pk/Research/MatlabFns/
%%
% HARRIS - Harris corner detector
%
% Usage:                 cim = harris(im, sigma)
%                [cim, r, c] = harris(im, sigma, thresh, radius, disp)
%  [cim, r, c, rsubp, csubp] = harris(im, sigma, thresh, radius, disp)
%
% Arguments:   
%            im     - image to be processed.
%            sigma  - standard deviation of smoothing Gaussian. Typical
%                     values to use might be 1-3.
%            thresh - threshold (optional). Try a value ~1000.
%            radius - radius of region considered in non-maximal
%                     suppression (optional). Typical values to use might
%                     be 1-3.
%            disp   - optional flag (0 or 1) indicating whether you want
%                     to display corners overlayed on the original
%                     image. This can be useful for parameter tuning. This
%                     defaults to 0
%
% Returns:
%            cim    - binary image marking corners.
%            r      - row coordinates of corner points.
%            c      - column coordinates of corner points.
%            rsubp  - If five return values are requested sub-pixel
%            csubp  - localization of feature points is attempted and
%                     returned as an additional set of floating point
%                     coords. Note that you may still want to use the integer
%                     valued coords to specify centres of correlation windows
%                     for feature matching.
%
% If thresh and radius are omitted from the argument list only 'cim' is returned
% as a raw corner strength image.  You may then want to look at the values
% within 'cim' to determine the appropriate threshold value to use. Note that
% the Harris corner strength varies with the intensity gradient raised to the
% 4th power.  Small changes in input image contrast result in huge changes in
% the appropriate threshold.
%
% Note that this code computes Noble's version of the detector which does not
% require the parameter 'k'.  See comments in code if you wish to use Harris'
% original measure.
%
% See also: NONMAXSUPPTS, DERIVATIVE5

% References: 
% C.G. Harris and M.J. Stephens. "A combined corner and edge detector", 
% Proceedings Fourth Alvey Vision Conference, Manchester.
% pp 147-151, 1988.
%
% Alison Noble, "Descriptions of Image Surfaces", PhD thesis, Department
% of Engineering Science, Oxford University 1989, p45.

% Copyright (c) 2002-2010 Peter Kovesi
% Centre for Exploration Targeting
% The University of Western Australia
% http://www.csse.uwa.edu.au/~pk/research/matlabfns/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% March    2002 - Original version
% December 2002 - Updated comments
% August   2005 - Changed so that code calls nonmaxsuppts
% August   2010 - Changed to use Farid and Simoncelli's derivative filters

function [cim, r, c, rsubp, csubp] = harris(im, sigma, thresh, radius, disp)
    
    error(nargchk(2,5,nargin));
    if nargin == 4
	disp = 0;
    end
    
    if ~isa(im,'double')
	im = double(im);
    end

    subpixel = nargout == 5;

    % Compute derivatives and elements of the structure tensor.
    [Ix, Iy] = derivative5(im, 'x', 'y');
    Ix2 = gaussfilt(Ix.^2,  sigma);
    Iy2 = gaussfilt(Iy.^2,  sigma);    
    Ixy = gaussfilt(Ix.*Iy, sigma);    

    % Compute the Harris corner measure. Note that there are two measures
    % that can be calculated.  I prefer the first one below as given by
    % Nobel in her thesis (reference above).  The second one (commented out)
    % requires setting a parameter, it is commonly suggested that k=0.04 - I
    % find this a bit arbitrary and unsatisfactory. 

    cim = (Ix2.*Iy2 - Ixy.^2)./(Ix2 + Iy2 + eps); % My preferred  measure.
%    k = 0.04;
%    cim = (Ix2.*Iy2 - Ixy.^2) - k*(Ix2 + Iy2).^2; % Original Harris measure.

    if nargin > 2   % We should perform nonmaximal suppression and threshold

	if disp  % Call nonmaxsuppts to so that image is displayed
	    if subpixel
		[r,c,rsubp,csubp] = nonmaxsuppts(cim, radius, thresh, im);
	    else
		[r,c] = nonmaxsuppts(cim, radius, thresh, im);		
	    end
	else     % Just do the nonmaximal suppression
	    if subpixel
		[r,c,rsubp,csubp] = nonmaxsuppts(cim, radius, thresh);
	    else
		[r,c] = nonmaxsuppts(cim, radius, thresh);		
	    end
	end
    end
   
    
    
%%
% DERIVATIVE5 - 5-Tap 1st and 2nd discrete derivatives
%
% This function computes 1st and 2nd derivatives of an image using the 5-tap
% coefficients given by Farid and Simoncelli.  The results are significantly
% more accurate than MATLAB's GRADIENT function on edges that are at angles
% other than vertical or horizontal. This in turn improves gradient orientation
% estimation enormously.  If you are after extreme accuracy try using DERIVATIVE7.
%
% Usage:  [gx, gy, gxx, gyy, gxy] = derivative5(im, derivative specifiers)
%
% Arguments:
%                       im - Image to compute derivatives from.
%    derivative specifiers - A comma separated list of character strings
%                            that can be any of 'x', 'y', 'xx', 'yy' or 'xy'
%                            These can be in any order, the order of the
%                            computed output arguments will match the order
%                            of the derivative specifier strings.
% Returns:
%  Function returns requested derivatives which can be:
%     gx, gy   - 1st derivative in x and y
%     gxx, gyy - 2nd derivative in x and y
%     gxy      - 1st derivative in y of 1st derivative in x
%
%  Examples:
%    Just compute 1st derivatives in x and y
%    [gx, gy] = derivative5(im, 'x', 'y');  
%                                           
%    Compute 2nd derivative in x, 1st derivative in y and 2nd derivative in y
%    [gxx, gy, gyy] = derivative5(im, 'xx', 'y', 'yy')
%
% See also: DERIVATIVE7

% Reference: Hany Farid and Eero Simoncelli "Differentiation of Discrete
% Multi-Dimensional Signals" IEEE Trans. Image Processing. 13(4): 496-508 (2004)

% Copyright (c) 2010 Peter Kovesi
% Centre for Exploration Targeting
% The University of Western Australia
% http://www.csse.uwa.edu.au/~pk/research/matlabfns/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.
%
% April 2010

function varargout = derivative5(im, varargin)

    varargin = varargin(:);
    varargout = cell(size(varargin));
    
    % Check if we are just computing 1st derivatives.  If so use the
    % interpolant and derivative filters optimized for 1st derivatives, else
    % use 2nd derivative filters and interpolant coefficients.
    % Detection is done by seeing if any of the derivative specifier
    % arguments is longer than 1 char, this implies 2nd derivative needed.
    secondDeriv = false;    
    for n = 1:length(varargin)
        if length(varargin{n}) > 1
            secondDeriv = true;
            break
        end
    end
    
    if ~secondDeriv
        % 5 tap 1st derivative cofficients.  These are optimal if you are just
        % seeking the 1st deriavtives
        p = [0.037659  0.249153  0.426375  0.249153  0.037659];
        d1 =[0.109604  0.276691  0.000000 -0.276691 -0.109604];
    else         
        % 5-tap 2nd derivative coefficients. The associated 1st derivative
        % coefficients are not quite as optimal as the ones above but are
        % consistent with the 2nd derivative interpolator p and thus are
        % appropriate to use if you are after both 1st and 2nd derivatives.
        p  = [0.030320  0.249724  0.439911  0.249724  0.030320];
        d1 = [0.104550  0.292315  0.000000 -0.292315 -0.104550];
        d2 = [0.232905  0.002668 -0.471147  0.002668  0.232905];
    end

    % Compute derivatives.  Note that in the 1st call below MATLAB's conv2
    % function performs a 1D convolution down the columns using p then a 1D
    % convolution along the rows using d1. etc etc.
    gx = false;
    
    for n = 1:length(varargin)
      if strcmpi('x', varargin{n})
          varargout{n} = conv2(p, d1, im, 'same');    
          gx = true;   % Record that gx is available for gxy if needed
          gxn = n;
      elseif strcmpi('y', varargin{n})
          varargout{n} = conv2(d1, p, im, 'same');
      elseif strcmpi('xx', varargin{n})
          varargout{n} = conv2(p, d2, im, 'same');    
      elseif strcmpi('yy', varargin{n})
          varargout{n} = conv2(d2, p, im, 'same');
      elseif strcmpi('xy', varargin{n}) | strcmpi('yx', varargin{n})
          if gx
              varargout{n} = conv2(d1, p, varargout{gxn}, 'same');
          else
              gx = conv2(p, d1, im, 'same');    
              varargout{n} = conv2(d1, p, gx, 'same');
          end
      else
          error(sprintf('''%s'' is an unrecognized derivative option',varargin{n}));
      end
    end

%%
% GAUSSFILT -  Small wrapper function for convenient Gaussian filtering
%
% Usage:  smim = gaussfilt(im, sigma)
%
% Arguments:   im - Image to be smoothed.
%           sigma - Standard deviation of Gaussian filter.
%
% Returns:   smim - Smoothed image.
%
% See also:  INTEGGAUSSFILT

% Peter Kovesi
% Centre for Explortion Targeting
% The University of Western Australia
% http://www.csse.uwa.edu.au/~pk/research/matlabfns/

% March 2010

function smim = gaussfilt(im, sigma)
 
    assert(ndims(im) == 2, 'Image must be greyscale');
    
    % If needed convert im to double
    if ~strcmp(class(im),'double')
        im = double(im);  
    end
    
    sze = ceil(6*sigma);  
    if ~mod(sze,2)    % Ensure filter size is odd
        sze = sze+1;
    end
    sze = max(sze,1); % and make sure it is at least 1
    
    h = fspecial('gaussian', [sze sze], sigma);

    smim = filter2(h, im);
