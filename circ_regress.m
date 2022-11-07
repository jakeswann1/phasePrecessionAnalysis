function [m, phi0, rho, p, phiOff] = circ_regress(phase, pos, maxNCyclesPerField, minMeth)
% circ_regress - circular-linear regression according to Kempter et al. (2016)
% https://www.sciencedirect.com/science/article/pii/S016502701200101X
% 
% One general snag with using the original Kempter method is that it seems fairly common that if the phase range is not close-ish to optimal for a particular cell,
% i.e. decent precession across Gaussian-shaped field, then instead of a negative slope, the minimisation returns a positive one (despite by eye there being clear 
% precession). (I think) this is because on a linear track fields are typically asymmetric and not Gaussian-shaped and the phase-pos relationship is often weak(er) 
% in the first half where rates tend to be lower  
% So here you have the option to also use a phaseshift for the cost function minimisation which  
%
% Usage:
%
%       [beta, phi0, rho] = circ_regress(phase, pos)
%       [beta, phi0, rho] = circ_regress(phase, pos, maxNCyclesPerField)
%
% Inputs:   phase              - phase values
%           pos                - pos samples
%           maxNCyclesPerField - bounds for max/min slope for regression
%                                (this many cycles per range of pos)
%
% Outputs:  m     - slope for circ-lin regression
%           phi0  - phase offset for circ-lin regression (y intercept)
%           rho   - correlation coefficient for regression
%                  (circ-lin analogue for Person's product moment corr.) 
%
%
% LM 2020

%% To-Do:
% pVal calculation for rho

%%
if nargin < 3
    % how many cycles of theta per field are 'allowed' for precession (default = 1)
    maxNCyclesPerField = 1; % N.B. it's def. not uncommon to get better fits when using > 1
end

if nargin < 4
   minMeth = 'slope'; 
end

% remove possible nan's
nanInd = isnan(phase) | isnan(pos);
phase  = phase(~nanInd);
phase  = phase(:);
pos    = pos(~nanInd);
% also normalise pos to [0 1]
pos   = (pos(:) - min(pos)) ./ range(pos);
% ... or between -1 1?? this might acount better for anisotropy of fields??
% pos = pos(:) - mean(pos(:));
% pos = pos./max(abs(pos));

%%
% get slope - it seems fair to generally assume to allow a slope of max. 1-2 full theta cycle(s)/field
% currently allowing negative gradients only (only phase precession, not recession)
slopeBounds = [-1 0] .* (maxNCyclesPerField * 2*pi / ( max(pos) - min(pos)) ); 
% 2 options to do the regression...
switch lower(minMeth)
    case 'slope'
        % as in Kempter et al. 
        m      = fminbnd(@(in) cost(in,phase,pos,false), slopeBounds(1), slopeBounds(2)); % 
        phiOff = 0;
    case 'slope+off'
        % add a phase shift parameter to minimisation routine
        rtn    = fminsearch(@(in) cost(in,phase,pos,true), [slopeBounds(1)/2 0]); % 
        m      = rtn(1); % slope
        phiOff = rtn(2); % phase shift
end
% phase       = phase + phiOff;
% phase offset (Eq. 2)
phi0        = atan2(nansum(sin(phase - m*pos)), nansum(cos(phase - m*pos)) );


% calc. correlation coefficient (Eq. 3)
circ_pos    = mod(abs(m)*pos, 2*pi);
meanPos     = mod(circ_mean(circ_pos), 2*pi);
meanPhase   = mod(circ_mean(phase), 2*pi);
rho         = nansum( sin(phase-meanPhase).*sin(circ_pos-meanPos) ) / sqrt( nansum( sin(phase-meanPhase).^2 ).*nansum( sin(circ_pos-meanPos).^2 ) );

% p value (in Appendix, after A.17)
n           = length(phase);
L20         = nansum( sin(phase-meanPhase).^2) / n;
L02         = nansum( sin(circ_pos-meanPos).^2) / n;
L22         = nansum( sin(phase-meanPhase).^2 .* sin(circ_pos-meanPos).^2 ) / n;
z           = rho*sqrt(n*L20*L02/L22);
p           = erfc(abs(z)/sqrt(2));  

end

function R = cost(in,phase,pos,addPhaseShift)
% cost function

if addPhaseShift
    phase = mod(phase+in(2),2*pi);
end
% Equation 1 in Kempter et al. (2012) (times -1 as we can only minimise but really want to maximise)
R = -1*sqrt( ( nansum( cos(phase - in(1)*pos) )/length(phase) ).^2 + ( nansum( sin(phase - in(1)*pos) )/length(phase) ).^2 ); % minimise this



end

        % minimise this, i.e. find best phaseshift
        
%         [phi0, r]  = fminbnd(@(P) phasePosR(P,phase,pos),-pi,pi); % 
%         b = regress(mod(phase+phi0,2*pi),[ones(length(phase),1), pos(:)]);
%         m = b(2);
%         [rho, p] = deal(nan);

% function r = phasePosR(P,phase,pos)
% 
% phase = mod(phase+P,2*pi);
% cov = nansum( (phase-circ_mean(phase)) .* (pos-mean(pos)) ) / (length(phase)-1);
% 
% r = cov / ( nanstd(phase+P) * nanstd(pos) );
% 
% 
% 
% end

