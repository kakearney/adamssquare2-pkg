function [x,y] = adamssquare2(phi, lam)
%ADAMSSQUARE2 Adams projection of the world in a square II
% 
% [x,y] = adamssquare2(phi, lam)
%
% This function converts lat/lon coordinates to the Adams projection of the
% world in a square II, assuming a spherical earth, based on the following
% publication:
%
% Adams, O. S. (1929). Conformal Projection of the Sphere Within a Square.
% Washington: U.S. Coast and Geodetic Survey Special Publication 153. 
%
% This code is a direct translation of proj_guyou.c, which was part of
% release 3 of libproj4. I think. Accessed from the following url on
% 2020/3/25:   
%
% https://github.com/jeffbaumes/jeffbaumes-vtk/blob/master/Utilities/vtklibproj4/proj_guyou.c
%
% Input variables:
%
%   phi:    latitude, array
%
%   lam:    longitude, array
%
% Output variables:
%
%   x:      x-coordinate in projected coordinates
%
%   y:      y-coordinate in projected coordinates

% Copyright 2020 Kelly Kearney

    spp = tan(0.5 .* phi);
    a = cos(asin(spp)) .* sin(0.5 .* lam);
    sm = (spp + a) < 0.;
    sn = (spp - a) < 0.;
    b = acos(spp);
    a = acos(a);
    
    m = asin(sqrt(abs(1. + cos(a + b))));
    m(sm) = -m(sm);
    
    n = asin(sqrt(abs(1. - cos(a - b))));
    n(sn) = -n(sn);
    
    x = ell_int_5(m);
    y = ell_int_5(n);
    
end

% Procedure to compute elliptic integral of the first kind where k^2=0.5.
% Precision good to better than 1e-7 The approximation is performed with an
% even Chebyshev series, thus the coefficients below are the even values
% and where series evaluation must be multiplied by the argument.

function out = ell_int_5(phi)

    y = phi .* 2/pi;
    y = 2 .* y .* y - 1;
    y2 = 2 .* y;

    d1 = 0;
    d2 = 0;

    C = [...
        2.19174570831038
        0.0914203033408211
       -0.00575574836830288
       -0.0012804644680613
        5.30394739921063e-05
        3.12960480765314e-05
        2.02692115653689e-07
       -8.58691003636495e-07];

    ord = 8;
    Cp = C;

    for ii = ord:-1:2
        temp = d1;
        d1 = y2 .* d1 - d2 + Cp(ii);
        d2 = temp;
    end
    
    out = phi .* (y .* d1 - d2 + 0.5 .* Cp(1));

end

% proj_guyou.c source code, for reference

% #define GUYOU 1
% #define PEIRCE_Q 2
% #define ADAMS_HEMI 3
% #define ADAMS_WSI 4
% #define ADAMS_WSII 5
% #define TOL 1e-9
% #define RSQRT2 0.7071067811865475244008443620
% #define PROJ_LIB__
% #define PROJ_PARMS__ \
%   int mode;
% #include  <lib_proj.h>
% PROJ_HEAD(guyou, "Guyou") "\n\tMisc., Sph., NoInv.";
% PROJ_HEAD(peirce_q, "Pierce Quincuncial") "\n\tMisc., Sph., NoInv.";
% PROJ_HEAD(adams_hemi, "Adams Hemisphere in a Square") "\n\tMisc., Sph., NoInv.";
% PROJ_HEAD(adams_wsI, "Adams World in a Square I") "\n\tMisc., Sph., NoInv.";
% PROJ_HEAD(adams_wsII, "Adams World in a Square II") "\n\tMisc., Sph., NoInv.";
% #define TWO_OVER_PI 0.6366197723675813430755350534
% #define ORDER 8
% /* Procedure to compute elliptic integral of the first kind
%  * where k^2=0.5.  Precision good to better than 1e-7
%  * The approximation is performed with an even Chebyshev
%  * series, thus the coefficients below are the even values
%  * and where series evaluation  must be multiplied by the argument. */
%   double
% ell_int_5(double phi) {
%   int i = ORDER;
%   double d1 = 0., d2 = 0., y, y2, temp;
%   const double C[] = { /* even coefficients */
%     2.19174570831038,
%     0.0914203033408211,
%     -0.00575574836830288,
%     -0.0012804644680613,
%     5.30394739921063e-05,
%     3.12960480765314e-05,
%     2.02692115653689e-07,
%     -8.58691003636495e-07};
%   double const *Cp = C + ORDER - 1;
% 
%   y = phi * TWO_OVER_PI;
%   y = 2. * y * y - 1.;
%   y2 = 2. * y;
%   while (--i) {
%     temp = d1;
%     d1 = y2 * d1 - d2 + *Cp--;
%     d2 = temp;
%   }
%   return phi * (y * d1 - d2 + 0.5 * *Cp);
% }
% FORWARD(s_forward); /* spheroid */
%   double m, n, a=0., b=0.;
%   int sm=0, sn=0;
%   switch (P->mode) {
%   case GUYOU:
%     if ((fabs(lp.lam) - TOL) > HALFPI) F_ERROR;
%     if (fabs(fabs(lp.phi) - HALFPI) < TOL) {
%       xy.x = 0;
%       xy.y = lp.phi < 0 ? -1.85407 : 1.85407;
%       return xy;
%     } else {
%       double sl = sin(lp.lam);
%       double sp = sin(lp.phi);
%       double cp = cos(lp.phi);
%       a = proj_acos((cp * sl - sp) * RSQRT2);
%       b = proj_acos((cp * sl + sp) * RSQRT2);
%       sm = lp.lam < 0.;
%       sn = lp.phi < 0.;
%     }
%     break;
%   case PEIRCE_Q: {
%       double sl = sin(lp.lam);
%       double cl = cos(lp.lam);
%       double cp = cos(lp.phi);
%       a = proj_acos(cp * (sl + cl) * RSQRT2);
%       b = proj_acos(cp * (sl - cl) * RSQRT2);
%       sm = sl < 0.;
%       sn = cl > 0.;
%     }
%     break;
%   case ADAMS_HEMI: {
%       double sp = sin(lp.phi);
%       if ((fabs(lp.lam) - TOL) > HALFPI) F_ERROR;
%       a = cos(lp.phi) * sin(lp.lam);
%       sm = (sp + a) < 0.;
%       sn = (sp - a) < 0.;
%       a = proj_acos(a);
%       b = HALFPI - lp.phi;
%     }
%     break;
%   case ADAMS_WSI: {
%       double sp = tan(0.5 * lp.phi);
%       b = cos(proj_asin(sp)) * sin(0.5 * lp.lam);
%       a = proj_acos((b - sp) * RSQRT2);
%       b = proj_acos((b + sp) * RSQRT2);
%       sm = lp.lam < 0.;
%       sn = lp.phi < 0.;
%     }
%     break;
%   case ADAMS_WSII: {
%       double spp = tan(0.5 * lp.phi);
%       a = cos(proj_asin(spp)) * sin(0.5 * lp.lam);
%       sm = (spp + a) < 0.;
%       sn = (spp - a) < 0.;
%       b = proj_acos(spp);
%       a = proj_acos(a);
%     }
%     break;
%   }
%   m = proj_asin(sqrt(fabs(1. + cos(a + b))));
%   if (sm) m = -m;
%   n = proj_asin(sqrt(fabs(1. - cos(a - b))));
%   if (sn) n = -n;
%   xy.x = ell_int_5(m);
%   xy.y = ell_int_5(n);
%   if (P->mode == ADAMS_HEMI || P->mode == ADAMS_WSII) { /* rotate by 45deg. */
%     double temp = xy.x;
%     xy.x = RSQRT2 * (xy.x - xy.y);
%     xy.y = RSQRT2 * (temp + xy.y);
%   }
%   return (xy);
% }
% FREEUP; if (P) free(P); }
%   static void*
% setup(PROJ *P) {
%     P->es = 0;
%   P->fwd = s_forward;
%   return P;
% }
% ENTRY0(guyou) P->mode = GUYOU; ENDENTRY(setup(P))
% ENTRY0(peirce_q) P->mode = PEIRCE_Q; ENDENTRY(setup(P))
% ENTRY0(adams_hemi) P->mode = ADAMS_HEMI; ENDENTRY(setup(P))
% ENTRY0(adams_wsI) P->mode = ADAMS_WSI; ENDENTRY(setup(P))
% ENTRY0(adams_wsII) P->mode = ADAMS_WSII; ENDENTRY(setup(P))