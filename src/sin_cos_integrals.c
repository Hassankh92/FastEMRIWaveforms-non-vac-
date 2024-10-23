
// #include "auxiliary_sin_cos_integrals_fi_gi.c"
// #include "power_series_Cin.c"
// #include "power_series_Si.c"
////////////////////////////////////////////////////////////////////////////////
// File: sin_cos_integrals.c                                                  //
// Routine(s):                                                                //
//    Sin_Integral_Si                                                         //
//    xSin_Integral_Si                                                        //
//    Entire_Cos_Integral_Cin                                                 //
//    xEntire_Cos_Integral_Cin                                                //
//    Cos_Integral_Ci                                                         //
//    xCos_Integral_Ci                                                        //
//    Sin_Cos_Integrals_Si_Ci                                                 //
//    xSin_Cos_Integrals_Si_Ci                                                //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>        // required for fabsl(), sinl(), cosl(), and logl()
#include <float.h>       // required for LDBL_EPSILON







long double xPower_Series_Si( long double x )
{ 
   long double sum;
   long double xx = - x * x;
   int n = (int) (3.42 * xx + 7.46 * fabsl(x) + 6.95);
   int k = n + n;

       // If the argument is sufficiently small use the approximation //
                         // Si(x) = x exp(-x^2/18) //
  
   if ( fabsl(x) <= 0.0003L ) return x * expl(xx / 18.0L);

        // Otherwise evaluate the power series expansion for Si(x). //

   sum = xx / ( (k + 1) * (k + 1) * k) + 1.0L / (k - 1);
   for (k--; k >= 3; k -=2) {
      sum *= xx / ( k * (k - 1) );
      sum += 1.0L / (k - 2);
   }
   return x * sum;
}







long double xPower_Series_Cin( long double x )
{ 
   long double sum;
   long double xx = - x * x;
   int n = (int) (2.41 * xx + 7.15 * fabsl(x) + 7.00);
   int k = n + n;
  
       // If the argument is sufficiently small use the approximation //
                      // Cin(x) = x^2 exp(-x^2/24) / 4 //

   if ( fabsl(x) < 0.00025L ) return - (expl(xx / 24.0L) * xx) / 4.0L;

       // Otherwise evaluate the power series expansion for Cin(x). //

   sum = xx / ( k * k * (k - 1)) + 1.0L / (k - 2);
   for (k -= 2; k > 2; k -= 2) {
      sum *= xx / ( k * (k - 1) );
      sum += 1.0L / (k - 2);
   }
   return -xx * sum / 2.0L;
}








////////////////////////////////////////////////////////////////////////////////
// File: auxiliary_sin_cos_integrals_fi_gi.c                                  //
// Routine(s):                                                                //
//    Auxiliary_Sin_Integral_fi                                               //
//    Auxiliary_Cos_Integral_gi                                               //
//    Auxiliary_Sin_Cos_Integrals_fi_gi                                       //
//    xAuxiliary_Sin_Integral_fi                                              //
//    xAuxiliary_Cos_Integral_gi                                              //
//    xAuxiliary_Sin_Cos_Integrals_fi_gi                                      //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Description:                                                              //
//     The auxiliary sin integral, fi(x), is the integral from 0 to infinity  //
//     of the integrand                                                       //
//                            sin(t) / (t + x) dt.                            //
//     The auxiliary cosine integral, gi(x), is the integral from 0 to        //
//     infinity of the integrand                                              //
//                            cos(t) / (t + x) dt.                            //
//                                                                            //
//     In terms of the sine integral, Si, and the cosine integral, Ci, the    //
//     auxiliary sine integral fi(x) is given by                              //
//               fi(x) = sin(x) Ci(x) + cos(x) [pi/2 - Si(x)]                 //
//     and the auxiliary cosine integral gi(x) is given by                    //
//               gi(x) = sin(x) [pi/2 - Si(x)] - cos(x) Ci(x).                //
//                                                                            //
//     The domain of the auxiliary sin integral fi(x) is {x : x >= 0 } where  //
//     fi(0) = pi/2.  The domain of the auxiliary cos integral gi(x) is       //
//     {x : x > 0} where as x -> 0+, gi(x) -> inf.                            //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>           // required for fabsl(), expl(), sinl(), cosl(),
                            // logl()
#include <float.h>          // required for LDBL_EPSILON

//                         Internally Defined Routines                        //

double      Auxiliary_Sin_Integral_fi( double x );
double      Auxiliary_Cos_Integral_gi( double x );
void        Auxiliary_Sin_Cos_Integrals_fi_gi(double x, double *fi, double *gi);
long double xAuxiliary_Sin_Integral_fi( long double x );
long double xAuxiliary_Cos_Integral_gi( long double x );
void        xAuxiliary_Sin_Cos_Integrals_fi_gi(long double x, long double *fi, 
                                                               long double *gi);

static long double xCos_Integral_Ci( long double x );
static long double fi_rational_polynomial(long double x, long double a[],
                                                      long double b[], int n );
static long double Asymptotic_Series_fi( long double x );
static long double gi_rational_polynomial(long double x, long double a[],
                                                      long double b[], int n );
static long double Asymptotic_Series_gi( long double x );

//                         Externally Defined Routines                        //

extern long double xPower_Series_Si( long double x );
extern long double xPower_Series_Cin( long double x );

//                         Internally Defined Constants                       //

static const long double pi2 = 1.57079632679489661923132169L;       // pi / 2
static const long double euler_gamma = 0.577215664901532860606512090L;
static const double auxiliary_asymptotic_cutoff = 48.0;

static long double a_x_ge_1_le_4_fi[] = {
      +3.131622691136541251894e+6L,  +5.865887504115410010938e+8L,
      +1.634852375578508416146e+10L, +1.592481384106901732624e+11L,
      +7.184770514348595264787e+11L, +1.726730020205455640781e+12L,
      +2.397017133822436251930e+12L, +2.020697105077248035167e+12L,
      +1.067232555649863576986e+12L, +3.595836616885923865165e+11L,
      +7.789746108788072914678e+10L, +1.083563302486680874140e+10L,
      +9.574882063563057212637e+8L,  +5.257964657853357906628e+7L,
      +1.727886704287183044067e+6L,  +3.186889399585378551937e+4L,
      +2.926771594419498165548e+2L
   };
static long double b_x_ge_1_le_4_fi[] = {
      +4.436542812456388065099e+7L,  +3.071881739597743437918e+9L,
      +5.510695064187223810111e+10L, +4.064528338807937104680e+11L,
      +1.502383531521515631047e+12L, +3.100277228892702060035e+12L,
      +3.813360580503500561372e+12L, +2.914078460895404297552e+12L,
      +1.419825394894026616675e+12L, +4.475808471509418423489e+11L,
      +9.178793064550390753770e+10L, +1.220855825313365329368e+10L,
      +1.040609666863148007442e+9L,  +5.554659381265650120714e+7L,
      +1.786400496083247945938e+6L,  +3.243425514601407346662e+4L,
      +2.946771482142805033246e+2L
   };
static long double a_x_ge_1_le_4_gi[] = {
      +9.011634207324336137169e+5L,  +7.479818286024998460948e+8L, 
      +4.151156375407831323555e+10L, +7.527803170191763096250e+11L,
      +6.273399733237371076085e+12L, +2.814715541899249302011e+13L,
      +7.442080767131902041599e+13L, +1.227172725914716222093e+14L,
      +1.309767511841246149009e+14L, +9.271225348708999857908e+13L,
      +4.418956912530285701879e+13L, +1.429482655697021907140e+13L,
      +3.143569573598121793475e+12L, +4.678982847861465840256e+11L,
      +4.663335634051774987907e+10L, +3.055838078958224702739e+9L,
      +1.280057381534594891504e+8L,  +3.283789466836908440869e+6L,
      +4.818060733773778820102e+4L,  +3.575079810165216346615e+2L
   };
static long double b_x_ge_1_le_4_gi[] = {
      +3.473778902563924058876e+8L,  +2.845671273312673204906e+10L,
      +6.887224173494194811858e+11L, +7.375036329278632360411e+12L,
      +4.176080452260044111884e+13L, +1.381922611468670308990e+14L,
      +2.845010797960102251342e+14L, +3.797756529707299562974e+14L,
      +3.379834764627141276920e+14L, +2.042485720392467096358e+14L,
      +8.475284361332246080070e+13L, +2.427267535696371015657e+13L,
      +4.796526275835169465639e+12L, +6.502922666518397649596e+11L,
      +5.978699373743563855764e+10L, +3.657882344026889055127e+9L,
      +1.447319540468370039281e+8L,  +3.546640511226990055118e+6L,
      +5.024169863961865278657e+4L,  +3.635079182389876878272e+2L
   };
static long double a_x_ge_4_le_12_fi[] = {
      +8.629036659345232923178e+15L, +9.470743102805298529462e+16L,
      +1.568021122342358329530e+17L, +9.015832733196613551192e+16L,
      +2.373367953145819143578e+16L, +3.275410521405716571530e+15L,
      +2.556227076494300926751e+14L, +1.177702886070105437976e+13L,
      +3.270951405687038350516e+11L, +5.490274976211303931784e+9L,
      +5.472393083052247561960e+7L,  +3.092021722264748314966e+5L,
      +8.929706311321410431845e+2L
   };
static long double b_x_ge_4_le_12_fi[] = {
      +3.688863305339062824609e+16L, +1.876827085370834659310e+17L,
      +2.346441340788672968041e+17L, +1.163521165422882838284e+17L,
      +2.802875478319020095488e+16L, +3.655276206330722751898e+15L,
      +2.748086189268929173963e+14L, +1.234713082649844595139e+13L,
      +3.371469088301994839064e+11L, +5.594066018240082151795e+9L,
      +5.532510775982731500507e+7L,  +3.109681134906426986992e+5L,
      +8.949706311313908973230e+2L
   };
static long double a_x_ge_4_le_12_gi[] = {
      +9.760124389962086158256e+17L, +2.768135717060729724771e+19L,
      +7.269925460678163397319e+19L, +6.335403079477117544205e+19L,
      +2.521611356160483301958e+19L, +5.326725622049037767865e+18L,
      +6.500059887901948040470e+17L, +4.822924381737713175777e+16L,
      +2.243854020350856804468e+15L, +6.651665288514689504327e+13L,
      +1.260766706261790080221e+12L, +1.513530299892650289088e+10L,
      +1.121138426325906850959e+8L,  +4.860629732996342070790e+5L,
      +1.106668096706748177652e+3L 
   };
static long double b_x_ge_4_le_12_gi[] = {
      +2.816012818637797223215e+19L, +1.453987140070137119268e+20L,
      +2.126102863700101349915e+20L, +1.330695747874759471888e+20L,
      +4.272711621417996420464e+19L, +7.778583943891982632436e+18L,
      +8.532286887837346444555e+17L, +5.856829003446583158094e+16L,
      +2.573169752130006183553e+15L, +7.312517856321160958464e+13L,
      +1.343719194435306279692e+12L, +1.577108016388663217206e+10L,
      +1.149410763356329587152e+8L,  +4.926189818912136539829e+5L,
      +1.112668096702659763721e+3L
   };
static long double a_x_ge_12_le_48_fi[] = {
      +8.190718946165709238422e+17L, +1.209912798380869069939e+18L,
      +2.685711451753038556686e+17L, +2.031432644806673394287e+16L,
      +6.849516346373244528380e+14L, +1.167908359237227948685e+13L,
      +1.071365422608890062545e+11L, +5.395836264116777645374e+8L,
      +1.462073394608352079917e+6L,  +1.959326763594685895502e+3L
   };
static long double b_x_ge_12_le_48_fi[] = {
      +1.759376483182613052616e+18L, +1.549737809630230245083e+18L,
      +3.002314821022841548975e+17L, +2.150253471166368305136e+16L,
      +7.064600781175281798566e+14L, +1.188341971751225609460e+13L,
      +1.081876692043348699994e+11L, +5.424692186656225562683e+8L,
      +1.465972048135541454369e+6L,  +1.961326763594685895323e+3L
   };
static long double a_x_ge_12_le_48_gi[] = {
      +5.524091612614961621464e+19L, +1.284075904576105184520e+20L,
      +3.447334407523257944528e+19L, +3.121715037272484722094e+18L,
      +1.282539019600256176592e+17L, +2.740263968387649522824e+15L,
      +3.267265290103262920765e+13L, +2.245923126260050126684e+11L,
      +8.923806059854096302378e+8L,  +1.985082566703293127903e+6L,
      +2.254025115381787893881e+3L
   };
static long double b_x_ge_12_le_48_gi[] = {
      +3.247999301164088453284e+20L, +2.442688918303073183435e+20L,
      +4.767807497134760332700e+19L, +3.740845893032137972381e+18L,
      +1.425986072860589430641e+17L, +2.920317933370183472849e+15L,
      +3.395218149102856121458e+13L, +2.297881510221565965240e+11L,
      +9.041055792759368518992e+8L,  +1.998522717395583928785e+6L,
      +2.260025115381787888363e+3L
   };
static int n_x_ge_1_le_4_fi = sizeof(a_x_ge_1_le_4_fi) / sizeof(long double);
static int n_x_ge_1_le_4_gi = sizeof(a_x_ge_1_le_4_gi) / sizeof(long double);
static int n_x_ge_4_le_12_fi = sizeof(a_x_ge_4_le_12_fi) / sizeof(long double);
static int n_x_ge_4_le_12_gi = sizeof(a_x_ge_4_le_12_gi) / sizeof(long double);
static int n_x_ge_12_le_48_fi = sizeof(a_x_ge_12_le_48_fi)
                                                         / sizeof(long double);
static int n_x_ge_12_le_48_gi = sizeof(a_x_ge_12_le_48_gi)
                                                         / sizeof(long double);

////////////////////////////////////////////////////////////////////////////////
// double Auxiliary_Sin_Integral_fi( double x )                               //
//                                                                            //
//  Description:                                                              //
//     This routine returns the auxiliary sin integral fi(x) for x >= 0.      //
//                                                                            //
//  Arguments:                                                                //
//     double  x  The argument of the auxiliary sin integral fi(), x >= 0.    //
//                                                                            //
//  Return Value:                                                             //
//     The value of the auxiliary sin integral fi evaluated at x >= 0.        //
//                                                                            //
//  Example:                                                                  //
//     double y, x;                                                           //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Auxiliary_Sin_Integral_fi( x );                                    //
////////////////////////////////////////////////////////////////////////////////
double Auxiliary_Sin_Integral_fi( double x )
{
   if (x == 0.0) return (double) pi2;
   return (double) xAuxiliary_Sin_Integral_fi((long double) x);
}


////////////////////////////////////////////////////////////////////////////////
// double Auxiliary_Cos_Integral_gi( double x )                               //
//                                                                            //
//  Description:                                                              //
//     This routine returns the auxiliary cos integral, gi(x), for x > 0.     //
//     Note! As x -> 0+, gi(x) -> inf.                                        //
//                                                                            //
//  Arguments:                                                                //
//     double  x  The argument of the auxiliary cos integral gi(), x > 0.     //
//                                                                            //
//  Return Value:                                                             //
//     The value of the auxiliary cos integral gi evaluated at x > 0.         //
//     If x = 0.0, then DBL_MAX is returned.                                  //
//                                                                            //
//  Example:                                                                  //
//     double y, x;                                                           //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Auxiliary_Cos_Integral_gi( x );                                    //
////////////////////////////////////////////////////////////////////////////////
double Auxiliary_Cos_Integral_gi( double x )
{
   long double gi = xAuxiliary_Cos_Integral_gi((long double) x);
   return (gi >= DBL_MAX) ? DBL_MAX : (double) gi;
}


////////////////////////////////////////////////////////////////////////////////
// void Auxiliary_Sin_Cos_Integrals_fi_gi( double x, double *fi, double *gi ) //
//                                                                            //
//  Description:                                                              //
//     This routine returns both the auxiliary sin integral fi(x) and the     //
//     auxiliary cos integral gi(x) via the addresses in the argument list    //
//     for x > 0.  For x = 0, *fi is set to pi/2 and *gi is set to DBL_MAX.   //
//                                                                            //
//  Arguments:                                                                //
//     double  x                                                              //
//        The argument of the auxiliary integrals fi() and gi().  The argument//
//        x must be nonnegative. If x = 0, the *fi = pi/2 and *gi = DBL_MAX.  //
//     double *fi                                                             //
//        The address of the auxiliary sin integral fi() evaluated at x.      //
//     double *gi                                                             //
//        The address of the auxiliary cos integral gi() evaluated at x.      //
//                                                                            //
//  Return Value:                                                             //
//     Type void.  The results are returned via the addresses in the argument //
//     list.                                                                  //
//                                                                            //
//  Example:                                                                  //
//     double x, fi, gi;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     Auxiliary_Sin_Cos_Integrals_fi_gi( x, &fi, &gi );                      //
////////////////////////////////////////////////////////////////////////////////
void Auxiliary_Sin_Cos_Integrals_fi_gi( double x, double *fi, double *gi )
{
   long double xfi, xgi;

   xAuxiliary_Sin_Cos_Integrals_fi_gi((long double) x, &xfi, &xgi);

   *fi = (double) xfi;
   *gi = (xgi >= DBL_MAX) ? DBL_MAX : (double) xgi;
}


////////////////////////////////////////////////////////////////////////////////
// long double xAuxiliary_Sin_Integral_fi( long double x )                    //
//                                                                            //
//  Description:                                                              //
//     This routine returns the sin integral fi(x) for x >= 0.                //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of the auxiliary sin integral fi(),       //
//                     x >= 0.                                                //
//                                                                            //
//  Return Value:                                                             //
//     The value of the auxiliary sin integral fi evaluated at x >= 0.        //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//     y = xAuxiliary_Sin_Integral_fi( x );                                   //
////////////////////////////////////////////////////////////////////////////////
long double xAuxiliary_Sin_Integral_fi(long double x)
{
   long double si;
   long double ci;

   if (x == 0.0L) return pi2;
   if (x <= 1.0L) {
      si = xPower_Series_Si(x);
      ci = xCos_Integral_Ci(x);
      return sinl(x) * ci + cosl(x) * (pi2 - si);
   }
   if (x <= 4.0L) return fi_rational_polynomial(x, a_x_ge_1_le_4_fi, 
                                       b_x_ge_1_le_4_fi, n_x_ge_1_le_4_fi);
   if (x <= 12.0L) return fi_rational_polynomial(x, a_x_ge_4_le_12_fi, 
                                       b_x_ge_4_le_12_fi, n_x_ge_4_le_12_fi);
   if (x < auxiliary_asymptotic_cutoff)
                   return fi_rational_polynomial(x, a_x_ge_12_le_48_fi,
                                       b_x_ge_12_le_48_fi, n_x_ge_12_le_48_fi);
   return Asymptotic_Series_fi( x );
}


////////////////////////////////////////////////////////////////////////////////
// long double xAuxiliary_Cos_Integral_gi( long double x )                    //
//                                                                            //
//  Description:                                                              //
//     This routine returns the cos integral gi(x) for x > 0.                 //
//     Note! As x -> 0+, gi(x) -> inf.                                        //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of the auxiliary cos integral gi(),       //
//                     x >  0.                                                //
//                                                                            //
//  Return Value:                                                             //
//     The value of the auxiliary cos integral gi evaluated at x > 0.         //
//     If x = 0.0, then LDBL_MAX is returned.                                 //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//     y = xAuxiliary_Cos_Integral_gi( x );                                   //
////////////////////////////////////////////////////////////////////////////////
long double xAuxiliary_Cos_Integral_gi(long double x)
{
   long double si;
   long double ci;

   if (x == 0.0L) return LDBL_MAX;
   if (x <= 1.0L) {
      si = xPower_Series_Si(x);
      ci = xCos_Integral_Ci(x);
      return sinl(x) * (pi2 - si) - cosl(x) * ci;
   }
   if (x <= 4.0L) return gi_rational_polynomial(x, a_x_ge_1_le_4_gi, 
                                       b_x_ge_1_le_4_gi, n_x_ge_1_le_4_gi);
   if (x <= 12.0L) return gi_rational_polynomial(x, a_x_ge_4_le_12_gi, 
                                       b_x_ge_4_le_12_gi, n_x_ge_4_le_12_gi);
   if (x < auxiliary_asymptotic_cutoff)
                   return gi_rational_polynomial(x, a_x_ge_12_le_48_gi,
                                       b_x_ge_12_le_48_gi, n_x_ge_12_le_48_gi);
   return Asymptotic_Series_gi( x );
}


////////////////////////////////////////////////////////////////////////////////
// void xAuxiliary_Sin_Cos_Integral_fi_gi( long double x, long double *fi,    //
//                                                         long double *gi )  //
//                                                                            //
//  Description:                                                              //
//     This routine returns both the auxiliary sin integral fi(x) and the     //
//     auxiliary cos integral gi(x) via the addresses in the argument list    //
//     for x > 0.  For x = 0, *fi is set to pi/2 and *gi is set to DBL_MAX.   //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//        The argument of the auxiliary integrals fi() and gi().  The argument//
//        x must be nonnegative. If x = 0, the *fi = pi/2 and *gi = LDBL_MAX. //
//     long double *fi                                                        //
//        The address of the auxiliary sin integral fi() evaluated at x.      //
//     long double *gi                                                        //
//        The address of the auxiliary cos integral gi() evaluated at x.      //
//                                                                            //
//  Return Value:                                                             //
//     Type void.  The results are returned via the addresses in the argument //
//     list.                                                                  //
//                                                                            //
//  Example:                                                                  //
//     long double x, fi, gi;                                                 //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     xAuxiliary_Sin_Cos_Integrals_fi_gi( x, &fi, &gi );                     //
////////////////////////////////////////////////////////////////////////////////
void xAuxiliary_Sin_Cos_Integrals_fi_gi(long double x, long double *fi,
                                                               long double *gi)
{
   long double si;
   long double ci;
   long double sx;
   long double cx;

   if (x == 0.0L) {
      *fi = pi2;
      *gi = LDBL_MAX;
   }      
   else if (x <= 1.0L) {
      si = xPower_Series_Si(x);
      ci = xCos_Integral_Ci(x);
      sx = sinl(x);
      cx = cosl(x);
      *fi = sx * ci + cx * (pi2 - si);
      *gi = sx * (pi2 - si) - cx * ci;
   }
   else if (x <= 4.0) {
      *fi = fi_rational_polynomial(x, a_x_ge_1_le_4_fi, b_x_ge_1_le_4_fi,
                                                             n_x_ge_1_le_4_fi);
      *gi = gi_rational_polynomial(x, a_x_ge_1_le_4_gi, b_x_ge_1_le_4_gi,
                                                           n_x_ge_1_le_4_gi);
   }
   else if (x <= 12.0L) {
      *fi = fi_rational_polynomial(x, a_x_ge_4_le_12_fi, b_x_ge_4_le_12_fi,
                                                            n_x_ge_4_le_12_fi);
      *gi = gi_rational_polynomial(x, a_x_ge_4_le_12_gi, b_x_ge_4_le_12_gi,
                                                            n_x_ge_4_le_12_gi);
   }
   else if (x < auxiliary_asymptotic_cutoff) {
      *fi = fi_rational_polynomial(x, a_x_ge_12_le_48_fi, b_x_ge_12_le_48_fi,
                                                          n_x_ge_12_le_48_fi);
      *gi = gi_rational_polynomial(x, a_x_ge_12_le_48_gi, b_x_ge_12_le_48_gi,
                                                          n_x_ge_12_le_48_gi);
   }
   else { 
      *fi = Asymptotic_Series_fi( x );
      *gi = Asymptotic_Series_gi( x );
   }
   return;
}


////////////////////////////////////////////////////////////////////////////////
// static long double xCos_Integral_Ci( long double x )                       //
//                                                                            //
//  Description:                                                              //
//     The cos integral, Ci(x), is the integral with integrand                //
//                               - cos(t) / t                                 //
//     where the integral extends from x to inf.                              //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of the cos integral Ci().                 //
//                                                                            //
//  Return Value:                                                             //
//     The value of the cos integral Ci evaluated at x.                       //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = xCos_Integral_Ci( x );                                             //
////////////////////////////////////////////////////////////////////////////////
static long double xCos_Integral_Ci( long double x )
{
   if (x == 0.0L) return -LDBL_MAX;
   return logl(fabsl(x)) + euler_gamma - xPower_Series_Cin(x);
}


////////////////////////////////////////////////////////////////////////////////
// static long double fi_rational_polynomial(long double x, long double a[],  //
//                                                   long double b[], int n ) //
//                                                                            //
//  Description:                                                              //
//     This function returns a rational polynomial approximation to the       //
//     auxiliary sin integral for the argument x.  The rational polynomial    //
//     has the form p(x) / x q(x) where p(x) and q(x) are monic polynomials   //
//     of degree 2n,                                                          //
//                p(x) = x^2n + a[n-1] x^2(n-1) + ... + a[0]                  //
//          and   q(x) = x^2n + b[b-1] x^2(n-1) + ... + b[0].                 //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of the auxiliary sin integral fi().       //
//     long double a[] The arrays of coefficients of the numerator.           //
//     long double b[] The arrays of coefficients of the denominator.         //
//     int          n  The half-order of the two polynomials.                 //
//                                                                            //
//  Return Value:                                                             //
//     The value of the auxiliary sin integral fi evaluated at x.             //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//     long double a[], b[];                                                  //
//     int n;                                                                 //
//                                                                            //
//     ( code to initialize x, a, b, and n )                                  //
//                                                                            //
//     y = fi_rational_polynomial(x, a, b, n);                                //
////////////////////////////////////////////////////////////////////////////////
static long double fi_rational_polynomial(long double x, long double a[],
                                                       long double b[], int n )
{
   long double xx = x * x;
   long double numerator = xx + a[n-1];
   long double denominator = xx + b[n-1];
   int i; 

   for (i = n-2; i >= 0; i--) {
      numerator *= xx;
      numerator += a[i];
      denominator *= xx;
      denominator += b[i];
   }
   return (numerator / denominator) / x; 
}


////////////////////////////////////////////////////////////////////////////////
// static long double Asymptotic_Series_fi( long double x )                   //
//                                                                            //
//  Description:                                                              //
//     For a large argument x, the auxiliary sin integral, fi(x), can be      //
//     expressed as the asymptotic series                                     //
//      fi(x) ~ 1/x - 2/x^3 + 24/x^5 - 720/x^7 + ... + (2j)!/x(-x^2)^j + ...  //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//                The argument of the auxiliary sin integral fi().            //
//                                                                            //
//  Return Value:                                                             //
//     The value of the auxiliary sin integral fi evaluated at x.             //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Asymptotic_Series_fi( x );                                         //
////////////////////////////////////////////////////////////////////////////////

static long double Asymptotic_Series_fi( long double x )
{
   long double term = 1.0L;
   long double xx = - x * x;
   long double xn = 1.0L;
   long double factorial = 1.0L;
   long double fi = 0.0L;
   long double old_term = 0.0L;
   int j = 2;

   do {
      fi += term;
      old_term = term;
      factorial *= (long double) ( j * ( j - 1 ) );
      xn *= xx;
      term = factorial / xn;
      j += 2;
   } while (fabsl(term) < fabsl(old_term));

   return fi / x;    
}


////////////////////////////////////////////////////////////////////////////////
// static long double gi_rational_polynomial(long double x, long double a[],  //
//                                                   long double b[], int n ) //
//                                                                            //
//  Description:                                                              //
//     This function returns a rational polynomial approximation to the       //
//     auxiliary cos integral for the argument x.  The rational polynomial    //
//     has the form p(x) / x^2 q(x) where p(x) and q(x) are monic polynomials //
//     of degree 2n,                                                          //
//                p(x) = x^2n + a[n-1] x^2(n-1) + ... + a[0]                  //
//          and   q(x) = x^2n + b[n-1] x^2(n-1) + ... + b[0].                 //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of the auxiliary cos integral gi().       //
//     long double a[] The arrays of coefficients of the numerator.           //
//     long double b[] The arrays of coefficients of the denominator.         //
//     int          n  The half-order of the two polynomials.                 //
//                                                                            //
//  Return Value:                                                             //
//     The value of the auxiliary cos integral gi evaluated at x.             //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//     long double a[], b[];                                                  //
//     int n;                                                                 //
//                                                                            //
//     ( code to initialize x, a, b, and n )                                  //
//                                                                            //
//     y = gi_rational_polynomial(x, a, b, n);                                //
////////////////////////////////////////////////////////////////////////////////
static long double gi_rational_polynomial(long double x, long double a[],
                                                       long double b[], int n )
{
   long double xx = x * x;
   long double numerator = xx + a[n-1];
   long double denominator = xx + b[n-1];
   int i; 

   for (i = n-2; i >= 0; i--) {
      numerator *= xx;
      numerator += a[i];
      denominator *= xx;
      denominator += b[i];
   }
   return (numerator / denominator) / xx; 
}


////////////////////////////////////////////////////////////////////////////////
// static long double Asymptotic_Series_gi( long double x )                   //
//                                                                            //
//  Description:                                                              //
//     For a large argument x, the auxiliary cosine integral, gi(x), can be   //
//     expressed as the asymptotic series                                     //
//      gi(x) ~ 1/x^2 - 6/x^4 + 120/x^6 - 5040/x^8 + ...                      //
//                                               + (2j+1)!/x^2(-x^2)^j + ...  //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//                The argument of the auxiliary cos integral gi().            //
//                                                                            //
//  Return Value:                                                             //
//     The value of the auxiliary cos integral gi evaluated at x.             //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Asymptotic_Series_gi( x );                                         //
////////////////////////////////////////////////////////////////////////////////
static long double Asymptotic_Series_gi( long double x )
{
   long double term = 1.0L;
   long double xx = - x * x;
   long double xn = 1.0L;
   long double factorial = 1.0L;
   long double gi = 0.0L;
   long double old_term = 0.0L;
   int j = 3;

   do {
      gi += term;
      old_term = term;
      factorial *= (long double) ( j * ( j - 1 ) );
      xn *= xx;
      term = factorial / xn;
      j += 2;
   } while (fabsl(term) < fabsl(old_term));

   return -gi / xx;    
}













//                         Internally Defined Routines                        //
double      Sin_Integral_Si( double x );
double      Entire_Cos_Integral_Cin( double x );
double      Cos_Integral_Ci( double x );
void        Sin_Cos_Integrals_Si_Ci( double x, double *Si, double *Ci );
long double xSin_Integral_Si( long double x );
long double xEntire_Cos_Integral_Cin( long double x );
long double xCos_Integral_Ci( long double x );
void        xSin_Cos_Integrals_Si_Ci( long double x, long double *Si,
                                                            long double *Ci );

static long double Asymptotic_Series_Ci( long double x );

//                         Externally Defined Routines                        //
extern void xAuxiliary_Sin_Cos_Integrals_fi_gi(long double x, long double *fi, 
                                                               long double *gi);
extern long double xPower_Series_Si( long double x );
extern long double xPower_Series_Cin( long double x );

//                         Internally Defined Constants                       //
static const long double pi =  3.1415926535897932384626433832795029L; 
// static const long double pi2 = 1.5707963267948966192313L;      // pi / 2
// static const long double euler_gamma = 0.577215664901532860606512090L;
// static const double auxiliary_asymptotic_cutoff = 48.0;

////////////////////////////////////////////////////////////////////////////////
// double Sin_Integral_Si( double x )                                         //
//                                                                            //
//  Description:                                                              //
//     The sin integral, Si(x), is the integral with integrand                //
//                             sin(t) / t                                     //
//     where the integral extends from 0 to x.                                //
//                                                                            //
//  Arguments:                                                                //
//     double  x  The argument of the sin integral Si().                      //
//                                                                            //
//  Return Value:                                                             //
//     The value of the sin integral Si evaluated at x.                       //
//                                                                            //
//  Example:                                                                  //
//     double y, x;                                                           //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Sin_Integral_Si( x );                                              //
////////////////////////////////////////////////////////////////////////////////
double Sin_Integral_Si( double x )
{
   return (double) xSin_Integral_Si( (long double) x);
}


////////////////////////////////////////////////////////////////////////////////
// long double xSin_Integral_Si( long double x )                              //
//                                                                            //
//  Description:                                                              //
//     The sin integral, Si(x), is the integral with integrand                //
//                             sin(t) / t                                     //
//     where the integral extends from 0 to x.                                //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of the sin integral Si().                 //
//                                                                            //
//  Return Value:                                                             //
//     The value of the sin integral Si evaluated at x.                       //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = xSin_Integral_Si( x );                                             //
////////////////////////////////////////////////////////////////////////////////

long double xSin_Integral_Si( long double x )
{
   long double fi, gi, si;

   if ( fabsl(x) <= 1.0L ) return xPower_Series_Si(x);
   xAuxiliary_Sin_Cos_Integrals_fi_gi(fabsl(x), &fi, &gi);
   si = pi2 - cosl(x) * fi - sinl(fabs(x)) * gi;
   return (x < 0.0L) ? -si : si;
}


////////////////////////////////////////////////////////////////////////////////
// double Entire_Cos_Integral_Cin( double x )                                 //
//                                                                            //
//  Description:                                                              //
//     The entire cos integral, Cin(x), is the integral with integrand        //
//                             (1 -  cos(t)) / t                              //
//     where the integral extends from 0 to x.                                //
//                                                                            //
//  Arguments:                                                                //
//     double  x  The argument of the entire cos integral Cin().              //
//                                                                            //
//  Return Value:                                                             //
//     The value of the entire cos integral Cin() evaluated at x.             //
//                                                                            //
//  Example:                                                                  //
//     double y, x;                                                           //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Entire_Cos_Integral_Cin( x );                                      //
////////////////////////////////////////////////////////////////////////////////
double Entire_Cos_Integral_Cin( double x )
{
   return (double) xEntire_Cos_Integral_Cin( (long double) x);
}


////////////////////////////////////////////////////////////////////////////////
// long double xEntire_Cos_Integral_Cin( long double x )                      //
//                                                                            //
//  Description:                                                              //
//     The entire cos integral, Cin(x), is the integral with integrand        //
//                             (1 -  cos(t)) / t                              //
//     where the integral extends from 0 to x.                                //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of the entire cos integral Cin().         //
//                                                                            //
//  Return Value:                                                             //
//     The value of the entire cos integral Cin() evaluated at x.             //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = xEntire_Cos_Integral_Cin( x );                                     //
////////////////////////////////////////////////////////////////////////////////

long double xEntire_Cos_Integral_Cin( long double x )
{
   if ( fabsl(x) <= 1.0L ) return xPower_Series_Cin(x);
   return logl(fabsl(x)) + euler_gamma - Asymptotic_Series_Ci(x);
}


////////////////////////////////////////////////////////////////////////////////
// double Cos_Integral_Ci( double x )                                         //
//                                                                            //
//  Description:                                                              //
//     The cos integral, Ci(x), is the integral with integrand                //
//                               - cos(t) / t                                 //
//     where the integral extends from x to inf.                              //
//                                                                            //
//  Arguments:                                                                //
//     double  x  The argument of the cos integral Ci().                      //
//                                                                            //
//  Return Value:                                                             //
//     The value of the cos integral Ci() evaluated at x.                     //
//                                                                            //
//  Example:                                                                  //
//     double y, x;                                                           //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Cos_Integral_Ci( x );                                              //
////////////////////////////////////////////////////////////////////////////////
double Cos_Integral_Ci( double x )
{
   long double ci = xCos_Integral_Ci( (long double) x);
      
   return (fabsl(ci) < DBL_MAX) ? (double) ci :
                                             (ci < 0.0L) ? -DBL_MAX : DBL_MAX;
}


////////////////////////////////////////////////////////////////////////////////
// long double xCos_Integral_Ci( long double x )                              //
//                                                                            //
//  Description:                                                              //
//     The cos integral, Ci(x), is the integral with integrand                //
//                               - cos(t) / t                                 //
//     where the integral extends from x to inf.                              //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of the cos integral Ci().                 //
//                                                                            //
//  Return Value:                                                             //
//     The value of the cos integral Ci evaluated at x.                       //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = xCos_Integral_Ci( x );                                             //
////////////////////////////////////////////////////////////////////////////////

// long double xCos_Integral_Ci( long double x )
// {
//    if (x == 0.0L) return -LDBL_MAX;
//    if (fabsl(x) <= 1.0L) 
//       return logl(fabsl(x)) + euler_gamma - xPower_Series_Cin(x);
//    return Asymptotic_Series_Ci(x);
// }


////////////////////////////////////////////////////////////////////////////////
// static long double Asymptotic_Series_Ci( long double x )                   //
//                                                                            //
//  Description:                                                              //
//     For a large argument x, the cosine integral, Ci(x), can be expressed as//
//                    Ci(x) = sin(x) fi(x) - cos(x) gi(x)                     //
//     where fi(x) and gi(x) are the auxiliary sine and cosine integrals which//
//     in turn are approximated using their respective asymptotic expansions. //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//                The argument of the cosine integral Ci().                   //
//                                                                            //
//  Return Value:                                                             //
//     The value of the cosine integral Ci evaluated at x.                    //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Asymptotic_Series_Ci( x );                                         //
////////////////////////////////////////////////////////////////////////////////

static long double Asymptotic_Series_Ci( long double x )
{
   long double fi, gi;

   xAuxiliary_Sin_Cos_Integrals_fi_gi(fabsl(x), &fi, &gi);

   return sinl(fabsl(x)) * fi - cosl(x) * gi;
}


////////////////////////////////////////////////////////////////////////////////
// void Sin_Cos_Integrals_Si_Ci( double x, double *Si, double *Ci )           //
//                                                                            //
//  Description:                                                              //
//     This routine returns both the sin integral Si(x) and the cos integral  //
//     Ci(x) via the addresses in the argument list.  For x = 0, *Ci is set   //
//     - DBL_MAX.                                                             //
//                                                                            //
//  Arguments:                                                                //
//     double  x                                                              //
//       The argument of the sin integral Si() and the cos integral Ci().     //
//     double *Si                                                             //
//        The address of the sin integral Si() evaluated at x.                //
//     double *Ci                                                             //
//        The address of the cos integral Ci() evaluated at x.                //
//                                                                            //
//  Return Value:                                                             //
//     Type void.  The results are returned via the addresses in the argument //
//     list.                                                                  //
//                                                                            //
//  Example:                                                                  //
//     double x, Si, Ci;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     Sin_Cos_Integrals_Si_Ci( x, &Si, &Ci );                                //
////////////////////////////////////////////////////////////////////////////////
void Sin_Cos_Integrals_Si_Ci( double x, double *Si, double *Ci )
{
   long double xSi, xCi;

   xSin_Cos_Integrals_Si_Ci( (long double) x, &xSi, &xCi);
   *Si = (double) xSi;
   *Ci = (fabsl(xCi) < DBL_MAX) ? (double) xCi : (xCi < 0) ? -DBL_MAX:DBL_MAX;
}


////////////////////////////////////////////////////////////////////////////////
// void xSin_Cos_Integrals_Si_Ci( long double x, long double *Si,             //
//                                                          long double *Ci ) //
//                                                                            //
//  Description:                                                              //
//     This routine returns both the sin integral Si(x) and the cos integral  //
//     Ci(x) via the addresses in the argument list.  For x = 0, *Ci is set   //
//     - LDBL_MAX.                                                            //
//                                                                            //
//  Arguments:                                                                //
//     long double  x                                                         //
//       The argument of the sin integral Si() and the cos integral Ci().     //
//     long double *Si                                                        //
//        The address of the sin integral Si() evaluated at x.                //
//     long double *Ci                                                        //
//        The address of the cos integral Ci() evaluated at x.                //
//                                                                            //
//  Return Value:                                                             //
//     Type void.  The results are returned via the addresses in the argument //
//     list.                                                                  //
//                                                                            //
//  Example:                                                                  //
//     long double x, Si, Ci;                                                 //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     xSin_Cos_Integrals_Si_Ci( x, &Si, &Ci );                               //
////////////////////////////////////////////////////////////////////////////////
void xSin_Cos_Integrals_Si_Ci(long double x, long double *Si, long double *Ci)
{
   long double fi, gi, sx, cx;

   if ( x == 0.0L) {
      *Si = 0.0L;
      *Ci = -LDBL_MAX;
      return;
   }
   if ( fabsl(x) <= 1.0L ) {
      *Si = xPower_Series_Si(x);
      *Ci = logl(fabsl(x)) + euler_gamma - xPower_Series_Cin(x);
      return;
   }
   xAuxiliary_Sin_Cos_Integrals_fi_gi(fabsl(x), &fi, &gi);
   sx = sinl(fabsl(x));
   cx = cosl(x);
   *Si = pi2 - cx * fi - sx * gi;
   *Ci = sx * fi - cx * gi;
   if (x < 0.0L) *Si = - *Si;
}
