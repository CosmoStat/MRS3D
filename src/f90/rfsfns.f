*	 							        ADCK0001
*									ADCK0002
*********************************************************************** ADCK0003
*    RFSFNS package                                                   * ADCK0004
*                                                                     * ADCK0005
*    The package RFSFNS (Root-Finder of Special FuNctionS) is written * ADCK0006
*  for the numerical localization, isolation and computation of roots * ADCK0007
*  of Bessel functions.                                               * ADCK0008
*                                                                     * ADCK0009
*    It  implements the  topological  degree of a  continuous mapping * ADCK0010
*  and, especially, Kronecker-Picard integral, for the computation of * ADCK0011
*  the total number  of simple real zeros of the Bessel  functions in * ADCK0012
*  a given interval.                                                  * ADCK0013
*                                                                     * ADCK0014
*    The package gives with certainty  the total number of the  above * ADCK0015
*  zeros within a predetermined interval. Also, it isolates each  one * ADCK0016
*  of these roots  and then  it computes  them within a given desired * ADCK0017
*  accuracy.                                                          * ADCK0018
*                                                                     * ADCK0019
*                                                                     * ADCK0020
*                                 RFSFNS.  Version of 21-07-1995      * ADCK0021
*                                                                     * ADCK0022
*                                                                     * ADCK0023
* Title of paper :  RFSFNS : A portable package for the numerical     * ADCK0024
*                   determination of the number and the calculation   * ADCK0025
*                   of roots of Bessel functions                      * ADCK0026
*                                                                     * ADCK0027
* To appear in   :  Computer Physics Communications                   * ADCK0028
*                                                                     * ADCK0029
* Authors        :  M.N. Vrahatis, O. Ragos, T. Skiniotis,            * ADCK0030
*                   F.A. Zafiropoulos and T.N. Grapsa                 * ADCK0031
* Address        :  Department of Mathematics                         * ADCK0032
*                   University of Patras                              * ADCK0033
*                   GR-261.10 Patras                                  * ADCK0034
*                   Greece                                            * ADCK0035
*                   Email : vrahatis@math.upatras.gr                  * ADCK0036
*                                                                     * ADCK0037
*********************************************************************** ADCK0038
      SUBROUTINE ROOTBESSJ( XXNU,A,B,MAXRT, INR, IERR, ROOTS )          ADCK0039
      IMPLICIT REAL*8(A-H, O-Z), INTEGER*8(I-N)                         ADCK0040
      INTEGER*8 MAXRT,IERR,NRF                                          ADCK0041
      DIMENSION  ROOTIS(MAXRT,2), ROOTS(MAXRT), VAR(MAXRT)              ADCK0042
      INTEGER*8  INF                                                    ADCK0043
      REAL*8 XXNU, XNU                                                  ADCK0044
      COMMON / BLK2 / GAMMA                                             ADCK0045
*     Set the starting values :                                         ADCK0046
*         a) the considered Bessel function specified by ICASE,         ADCK0047
*         b) the order of the corresponding Bessel function XNU,        ADCK0048
*         c) the endpoints of the given interval [A, B] and             ADCK0049
*         d) the conditional variable ICON. If ICON = 4, supply         ADCK0050
*            also the desired number of roots NR.                       ADCK0051
*     WRITE(*,*) XNU,XXNU                                               ADCK0052
*     Input data for Test Runs 1-4                                      ADCK0053
*      DATA   ICASE,    XNU,   A,       B,         ICON, NR              ADCK0054
*     +     / 1,        0.5D0, 0.D0,    2000.D0,   3,    0  /            ADCK0055
*     XNU=0.5D0                                                         ADCK0056
      XNU=XXNU!+XNU                                                     ADCK0057
*     WRITE(*,*) XNU,XXNU	   		                        ADCK0058
*      A=0.D0                                                            ADCK0059
*      B=2000.D0                                                         ADCK0060
      NR=0                                                              ADCK0061
      ICASE=1                                                           ADCK0062
      ICON=3                                                            ADCK0063
*     Set the values of EPSILO and DELTA.                               ADCK0064
*                                                                       ADCK0065
      DATA   EPSILO, DELTA                                              ADCK0066
     +     / 1.D-3,  1.D-15 /                                           ADCK0067
*                                                                       ADCK0068
*     GAMMA is a parameter involved in the Kronecker-Picard             ADCK0069
*           integral. Its default value is taken equal to 1.D0          ADCK0070
*           but the user is able to change it in his own main           ADCK0071
*           program. If this value is less than the machine             ADCK0072
*           precision EPSMCH, INTSUB sets GAMMA equal to 1.D0.          ADCK0073
*                                                                       ADCK0074
      GAMMA=1.D0                                                        ADCK0075
*                                                                       ADCK0076
*      PRINT 9999, ICASE, XNU, A, B, ICON, EPSILO, DELTA                ADCK0077
*      PRINT 9999, ICASE, IXNU, A, B, ICON, EPSILO, DELTA               ADCK0078
*     Call the interface subroutine INTSUB.                             ADCK0079
      CALL INTSUB (ICASE, A, B, ICON, EPSILO, DELTA, MAXRT, NR,         ADCK0080
     +             ROOTIS, ROOTS, VAR,XNU,NRF,INF)                      ADCK0081
      IERR = INF                                                        ADCK0082
      INR = NR                                                          ADCK0083
      IF ( INF .EQ. 0 ) THEN                                            ADCK0084
         PRINT 9998                                                     ADCK0085
         GO TO 10                                                       ADCK0086
      ENDIF                                                             ADCK0087
      IF ( INF .EQ. 2 ) THEN                                            ADCK0088
         PRINT 9997                                                     ADCK0089
         GO TO 10                                                       ADCK0090
      ENDIF                                                             ADCK0091
      IF ( INF .EQ. 3 ) THEN                                            ADCK0092
         PRINT 9996                                                     ADCK0093
*        PRINT 9995, NR                                                 ADCK0094
         GO TO 10                                                       ADCK0095
      ENDIF                                                             ADCK0096
*      PRINT 9995, NR                                                   ADCK0097
      IF ( ICON .EQ. 1 ) GO TO 10                                       ADCK0098
      IF ( INF .EQ. 4 ) THEN                                            ADCK0099
         PRINT 9994                                                     ADCK0100
         GO TO 10                                                       ADCK0101
      ENDIF                                                             ADCK0102
*      PRINT 9990, NR,NRF                                               ADCK0103
*      PRINT 9993, (J, ROOTIS(J,1), ROOTIS(J,2), J = 1, NRF)            ADCK0104
      IF ( ICON .EQ. 2 ) GO TO 10                                       ADCK0105
      IF ( INF .EQ. 5 ) THEN                                            ADCK0106
*        PRINT 9992                                                     ADCK0107
         GO TO 10                                                       ADCK0108
      ENDIF                                                             ADCK0109
*      PRINT 9991, (J, ROOTS(J), VAR(J), J = 1, NRF)                    ADCK0110
*      PRINT 9989, INF                                                  ADCK0111
      RETURN                                                            ADCK0112
   10 RETURN                                                            ADCK0113
*                                                                       ADCK0114
 9999 FORMAT (/2X, ' STARTING VALUES :' /3X, 17('-'),                   ADCK0115
     +        /2X, ' ICASE  :    ', I1,                                 ADCK0116
     +        /2X, ' ORDER  : ', F20.15,                                ADCK0117
     +        /2X, ' A      : ', F20.15,                                ADCK0118
     +        /2X, ' B      : ', F20.15,                                ADCK0119
     +        /2X, ' ICON   :    ', I1,                                 ADCK0120
     +        /2X, ' EPSILO : ', F20.15,                                ADCK0121
     +        /2X, ' DELTA  : ', F20.15 )                               ADCK0122
 9998 FORMAT (/2X, '  * * * IMPROPER INPUT PARAMETERS * * *'//)         ADCK0123
 9997 FORMAT (/2X, '  * * * THE PROCEDURE FOR THE CALCULATION',         ADCK0124
     +             ' OF THE BESSEL FUNCTION FAILED * * *'//)            ADCK0125
 9996 FORMAT (/2X, '  * * * THE PROCEDURE FOR THE CALCULATION',         ADCK0126
     +             ' OF THE TOPOLOGICAL DEGREE FAILED * * *'//)         ADCK0127
 9995 FORMAT (/2X, ' THE COMPUTED TOTAL NUMBER OF ROOTS',/,             ADCK0128
     +         2X, ' WITHIN THE INTERVAL (A,B) IS    : ', I5)           ADCK0129
 9994 FORMAT (/2X, ' * * * THE PROCEDURE FOR THE ROOTS',                ADCK0130
     +             ' ISOLATION FAILED * * *'//)                         ADCK0131
 9993 FORMAT (/2X, ' INTERVALS OF THE ISOLATED ROOTS :'                 ADCK0132
     +        /3X, 33('-') /(2X, I4,')',                                ADCK0133
     +         2X, ' (', F20.15, ',', F20.15,' )'))                     ADCK0134
 9992 FORMAT (/2X, '* * * THE ROOTFINDING PORTION',                     ADCK0135
     +             ' FAILED * * *'//)                                   ADCK0136
 9991 FORMAT (/2X, ' FINAL APPROXIMATE ROOTS :',                        ADCK0137
     +         7X, ' VERIFICATION :' /3X, 25('-'), 8X, 14('-')          ADCK0138
     +               /(2X, I4,')', F20.15, 5X, F20.15))                 ADCK0139
 9990 FORMAT (/2X, ' NUMBER OF ROOTS REQUESTED       : ',I5,/           ADCK0140
     +        /2X, ' NUMBER OF ROOTS ISOLATED        : ',I5)            ADCK0141
 9989 FORMAT (/2X, ' EXIT PARAMETER :   INF = ',I2)                     ADCK0142
*                                                                       ADCK0143
*     Last statement of the main program.                               ADCK0144
*                                                                       ADCK0145
      END                                                               ADCK0146
*_______________________________________________________________        ADCK0147
                                                                        ADCK0148
      SUBROUTINE INTSUB (ICASE, A, B, ICON, EPSILO, DELTA, MAXRT,       ADCK0149
     +                   NR, ROOTIS, ROOTS, VAR, XNU, NRF,INF)          ADCK0150
*********************************************************************** ADCK0151
*                                                                     * ADCK0152
*     SUBROUTINE INTSUB                                               * ADCK0153
*                                                                     * ADCK0154
*                                                                     * ADCK0155
*       This is an interface subroutine between the  main program and * ADCK0156
*       the subroutines  TOPDEG,  ISOLAT  and  BISECT.                * ADCK0157
*                                                                     * ADCK0158
*                                                                     * ADCK0159
*     The subroutine statement is :                                   * ADCK0160
*                                                                     * ADCK0161
*      SUBROUTINE INTSUB (ICASE, A, B, ICON, EPSILO, DELTA, MAXRT,    * ADCK0162
*     +                   NR, ROOTIS, ROOTS, VAR)                     * ADCK0163
*                                                                     * ADCK0164
*     where :                                                         * ADCK0165
*                                                                     * ADCK0166
*       ICASE is an integer input variable which specifies the        * ADCK0167
*         considered function as follows :                            * ADCK0168
*                                                                     * ADCK0169
*         ICASE = 1  Bessel function of first kind, J.                * ADCK0170
*                                                                     * ADCK0171
*         ICASE = 2  The derivative of Bessel function of first kind. * ADCK0172
*                                                                     * ADCK0173
*         ICASE = 3  Bessel function of second kind, Y.               * ADCK0174
*                                                                     * ADCK0175
*         ICASE = 4  The derivative of Bessel function of second kind.* ADCK0176
*                                                                     * ADCK0177
*       A is an input variable which defines the left endpoint of the * ADCK0178
*         given interval.                                             * ADCK0179
*                                                                     * ADCK0180
*       B is an input variable which defines the right endpoint of    * ADCK0181
*         the given interval.                                         * ADCK0182
*                                                                     * ADCK0183
*       ICON is an integer input variable with values :               * ADCK0184
*                                                                     * ADCK0185
*         ICON = 1   Calculate  only the total number of roots in the * ADCK0186
*                    given interval (A,B).                            * ADCK0187
*                                                                     * ADCK0188
*         ICON = 2   Calculate the total number of roots in the given * ADCK0189
*                    interval (A,B) and isolate each one of them.     * ADCK0190
*                                                                     * ADCK0191
*         ICON = 3   Calculate the total number of roots in the given * ADCK0192
*                    interval  (A,B), isolate and compute each one of * ADCK0193
*                    them.                                            * ADCK0194
*                                                                     * ADCK0195
*         ICON = 4   Calculate NR roots in the given interval  (A,B). * ADCK0196
*                                                                     * ADCK0197
*       EPSILO is a positive  input  variable  which  determines  the * ADCK0198
*         accuracy of the computation of Kronecker-Picard integral.   * ADCK0199
*         If EPSILO is less than the machine precision EPSMCH, EPSILO * ADCK0200
*         becomes equal to 1.D5*DSQRT(EPSMCH). The value of EPSMCH is * ADCK0201
*         computed within INTSUB.                                     * ADCK0202
*                                                                     * ADCK0203
*       DELTA is a  positive  input  variable  which  determines  the * ADCK0204
*         accuracy of the computation of the roots of the function.   * ADCK0205
*         If DELTA is less than the machine precision  EPSMCH,  DELTA * ADCK0206
*         becomes equal to  EPSMCH.  The value of  EPSMCH is computed * ADCK0207
*         within INTSUB.                                              * ADCK0208
*                                                                     * ADCK0209
*       MAXRT is a positive integer input  variable  which determines * ADCK0210
*         the maximum number of roots required in the given interval. * ADCK0211
*                                                                     * ADCK0212
*       NR is a positive integer output variable which determines the * ADCK0213
*         computed total number of roots in the given interval.       * ADCK0214
*         NR must be less than MAXRT.                                 * ADCK0215
*                                                                     * ADCK0216
*       ROOTIS is  an  output  MAXRT times 2  array.  Its two columns * ADCK0217
*         provide lower and upper root bounds.                        * ADCK0218
*                                                                     * ADCK0219
*       ROOTS is  an output array  of length  NR.  It determines  the * ADCK0220
*         final approximate roots.                                    * ADCK0221
*                                                                     * ADCK0222
*       VAR is  an output array  of  length NR  which  specifies  the * ADCK0223
*         function values at the final approximate roots.             * ADCK0224
*                                                                     * ADCK0225
*                                                                     * ADCK0226
*     The COMMON-transferred variables used here are :                * ADCK0227
*                                                                     * ADCK0228
*       INF is an integer output variable set as follows :            * ADCK0229
*                                                                     * ADCK0230
*         INF = 0    The initial value of  ICASE,XNU,ICON  or  NR is  * ADCK0231
*                    out of range.                                    * ADCK0232
*                                                                     * ADCK0233
*         INF = 1    The required answer is accurate.                 * ADCK0234
*                                                                     * ADCK0235
*         INF = 2    The evaluation of the function FNC or G failed.  * ADCK0236
*                                                                     * ADCK0237
*         INF = 3    The procedure for the calculation of  Kronecker- * ADCK0238
*                    Picard integral failed because                   * ADCK0239
*                       a) the number of iterations exceeded MAXI, or * ADCK0240
*                       b) the extrapolation procedure failed.        * ADCK0241
*                    ISOLAT and BISECT will not be evoked.            * ADCK0242
*                                                                     * ADCK0243
*         INF = 4    The  procedure  for the  isolation  of the roots * ADCK0244
*                    failed because                                   * ADCK0245
*                       a) the maximum number of iterations MAXDIV or * ADCK0246
*                       b) the maximum dimension of the internal      * ADCK0247
*                          working arrays WA1,WA2 of ISOLAT was       * ADCK0248
*                          exceeded.                                  * ADCK0249
*                    BISECT will not be evoked.                       * ADCK0250
*                                                                     * ADCK0251
*         INF = 5    The  rootfinding  procedure  failed  because the * ADCK0252
*                    maximum number of steps of BISECT was exceeded.  * ADCK0253
*                                                                     * ADCK0254
*       NRF is an integer  output  variable indicating  the number of * ADCK0255
*         roots which have been isolated.                             * ADCK0256
*                                                                     * ADCK0257
*       XNU is an  input  variable  which  defines  the order  of the * ADCK0258
*         Bessel function.                                            * ADCK0259
*                                                                     * ADCK0260
*       GAMMA  is  a  parameter  involved  in  the  Kronecker-Picard  * ADCK0261
*         integral. Its default value is taken equal to 1.D0 but the  * ADCK0262
*         user is able to change it in his own main program. If this  * ADCK0263
*         value is less than the machine precision EPSMCH, INTSUB     * ADCK0264
*         sets GAMMA equal to 1.D0.                                   * ADCK0265
*                                                                     * ADCK0266
*       EPSMCH is the machine precision.                              * ADCK0267
*                                                                     * ADCK0268
*                                                                     * ADCK0269
*     Subprograms required :                                          * ADCK0270
*                                                                     * ADCK0271
*       FORTRAN-Supplied ... DMAX1, DSQRT                             * ADCK0272
*                                                                     * ADCK0273
*       RFSFNS-Supplied .... TOPDEG, ISOLAT, BISECT, FNC              * ADCK0274
*                                                                     * ADCK0275
*                                                                     * ADCK0276
*                                 RFSFNS.  Version of 21-07-1995      * ADCK0277
*                                                                     * ADCK0278
*                                 Department of Mathematics           * ADCK0279
*                                 University of Patras                * ADCK0280
*                                 GR-261.10 Patras                    * ADCK0281
*                                 Greece                              * ADCK0282
*                                 Email : vrahatis@math.upatras.gr    * ADCK0283
*                                                                     * ADCK0284
*********************************************************************** ADCK0285
      IMPLICIT REAL*8(A-H, O-Z), INTEGER*8(I-N)                         ADCK0286
      PARAMETER (ZERO = 0.D0, HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0,     ADCK0287
     +           ONTH = 0.001D0, TEN5 = 1.D5)                           ADCK0288
      INTEGER*8  INF,NRF                                                ADCK0289
      COMMON / BLK2 / GAMMA                                             ADCK0290
      COMMON / BLK3 / EPSMCH                                            ADCK0291
      DIMENSION ROOTIS(MAXRT,2), ROOTS(MAXRT), VAR(MAXRT)               ADCK0292
      INTRINSIC DMAX1, DSQRT                                            ADCK0293
      REAL*8 XNU,TMPEPS                                                 ADCK0294
                                                                        ADCK0295
*     Compute the machine precision                                     ADCK0296
      TMPEPS = ONE                                                      ADCK0297
   10 TMPEPS = HALF*TMPEPS                                              ADCK0298
         TOL = ONE+TMPEPS                                               ADCK0299
         IF (TOL .GT. ONE) GOTO 10                                      ADCK0300
      EPSMCH = TWO*TMPEPS                                               ADCK0301
                                                                        ADCK0302
*     Check the input parameters for errors.                            ADCK0303
*                                                                       ADCK0304
      INF = 0                                                           ADCK0305
*                                                                       ADCK0306
      IF (ICASE .LT. 1 .OR. ICASE .GT. 4) RETURN                        ADCK0307
      IF (XNU .LT. ZERO) RETURN                                         ADCK0308
      IF (ICON .LT. 1 .OR. ICON .GT. 4) RETURN                          ADCK0309
      IF (ICON .EQ. 4 .AND. NR .GT. MAXRT) RETURN                       ADCK0310
*                                                                       ADCK0311
      ALIMIT = DMAX1(XNU,HALF)                                          ADCK0312
      IF (A .LT. ALIMIT) A = ALIMIT                                     ADCK0313
      IF (B .LE. A) RETURN                                              ADCK0314
*                                                                       ADCK0315
*     Check if any of the endpoints of the given interval is a root.    ADCK0316
*                                                                       ADCK0317
      IF (DABS(FNC(ICASE,A,XNU,INF)) .LT. EPSMCH) A = A-TWO*DELTA       ADCK0318
      IF (DABS(FNC(ICASE,B,XNU,INF)) .LT. EPSMCH) B = B+TWO*DELTA       ADCK0319
*                                                                       ADCK0320
*     Determine if the input values of GAMMA, EPSILO and DELTA,         ADCK0321
*     are less than the machine precision.                              ADCK0322
*                                                                       ADCK0323
      IF (GAMMA .LT. ONTH) GAMMA = ONE                                  ADCK0324
      IF (EPSILO .LT. DSQRT(EPSMCH)) EPSILO = TEN5*DSQRT(EPSMCH)        ADCK0325
      IF (DELTA .LT. EPSMCH) DELTA = EPSMCH                             ADCK0326
*                                                                       ADCK0327
      IF (ICON .EQ. 4) THEN                                             ADCK0328
         ICON = 3                                                       ADCK0329
         GOTO 20                                                        ADCK0330
      ENDIF                                                             ADCK0331
*                                                                       ADCK0332
*     Call TOPDEG.                                                      ADCK0333
*                                                                       ADCK0334
      CALL TOPDEG (ICASE, A, B, EPSILO, NR,XNU,INF)                     ADCK0335
      IF (NR .EQ. 0 .OR. ICON .EQ. 1) RETURN                            ADCK0336
      IF (NR .GT. MAXRT) THEN                                           ADCK0337
         INF = 0                                                        ADCK0338
         RETURN                                                         ADCK0339
      ENDIF                                                             ADCK0340
*                                                                       ADCK0341
*     Call ISOLAT.                                                      ADCK0342
*                                                                       ADCK0343
   20 CALL ISOLAT (ICASE, A, B, NR, MAXRT, ROOTIS,XNU,NRF,INF)          ADCK0344
      CALL ISOLAT (ICASE, A, B, NR, MAXRT, ROOTIS,XNU,NRF,INF)          ADCK0345
      IF (NRF .GT. NR) NRF = NR                                         ADCK0346
      IF (ICON .EQ. 2) RETURN                                           ADCK0347
      IF (INF .NE. 1) RETURN                                            ADCK0348
*     Call BISECT.                                                      ADCK0349
*                                                                       ADCK0350
      DO I=1, NRF                                                       ADCK0351
         A = ROOTIS(I,1)                                                ADCK0352
         B = ROOTIS(I,2)                                                ADCK0353
         CALL BISECT (ICASE, A, B, DELTA, ROOTS(I), VAR(I),XNU,INF)     ADCK0354
      ENDDO                                                             ADCK0355
*                                                                       ADCK0356
      RETURN                                                            ADCK0357
*                                                                       ADCK0358
*     Last statement of the interface subroutine INTSUB.                ADCK0359
*                                                                       ADCK0360
      END                                                               ADCK0361
*---------------------------------------------------------------------* ADCK0362
                                                                        ADCK0363
      SUBROUTINE TOPDEG (ICASE, A, B, EPSILO, NR,XNU,INF)               ADCK0364
*********************************************************************** ADCK0365
*                                                                     * ADCK0366
*     SUBROUTINE TOPDEG                                               * ADCK0367
*                                                                     * ADCK0368
*                                                                     * ADCK0369
*       The purpose of this subroutine is to compute the total number * ADCK0370
*       of roots  of  a  Bessel  function  in a  given  interval.  It * ADCK0371
*       implements the topological degree of a continuous mapping and * ADCK0372
*       especially  Kronecker-Picard  integral.  The  corresponding   * ADCK0373
*       integration is performed by Romberg's method.                 * ADCK0374
*                                                                     * ADCK0375
*                                                                     * ADCK0376
*     The subroutine statement is :                                   * ADCK0377
*                                                                     * ADCK0378
*     SUBROUTINE TOPDEG (ICASE, A, B, EPSILO, NR)                     * ADCK0379
*                                                                     * ADCK0380
*     where :                                                         * ADCK0381
*                                                                     * ADCK0382
*       ICASE is  an  integer  input  variable  which  specifies  the * ADCK0383
*         considered function as follows :                            * ADCK0384
*                                                                     * ADCK0385
*         ICASE = 1  Bessel function of first kind, J.                * ADCK0386
*                                                                     * ADCK0387
*         ICASE = 2  The derivative of Bessel function of first kind. * ADCK0388
*                                                                     * ADCK0389
*         ICASE = 3  Bessel function of second kind, Y.               * ADCK0390
*                                                                     * ADCK0391
*         ICASE = 4  The derivative of Bessel function of second kind.* ADCK0392
*                                                                     * ADCK0393
*       A is an input variable which defines the left endpoint of the * ADCK0394
*         given interval.                                             * ADCK0395
*                                                                     * ADCK0396
*       B is an input variable which defines the right endpoint of    * ADCK0397
*         the given interval.                                         * ADCK0398
*                                                                     * ADCK0399
*       EPSILO is a  positive  input  variable  which  determines the * ADCK0400
*         accuracy of the computation of Kronecker-Picard integral.   * ADCK0401
*         If EPSILO is less than the machine precision EPSMCH, EPSILO * ADCK0402
*         becomes equal to 1.D5*DSQRT(EPSMCH). The value of EPSMCH is * ADCK0403
*         computed within INTSUB.                                     * ADCK0404
*                                                                     * ADCK0405
*       NR is a positive integer output variable which determines the * ADCK0406
*         computed total number of roots in the given interval.       * ADCK0407
*         NR must be less than MAXRT.                                 * ADCK0408
*                                                                     * ADCK0409
*                                                                     * ADCK0410
*     The COMMON-transferred variables used here are :                * ADCK0411
*                                                                     * ADCK0412
*       INF is an integer output variable set as follows :            * ADCK0413
*                                                                     * ADCK0414
*         INF = 1    The required answer is accurate.                 * ADCK0415
*                                                                     * ADCK0416
*         INF = 3    The procedure for the calculation of  Kronecker- * ADCK0417
*                    Picard integral failed because                   * ADCK0418
*                       a) the number of iterations exceeded MAXI, or * ADCK0419
*                       b) the extrapolation procedure failed.        * ADCK0420
*                    ISOLAT and BISECT will not be evoked.            * ADCK0421
*                                                                     * ADCK0422
*       XNU is  an  input  variable  which  defines the order  of the * ADCK0423
*         Bessel function.                                            * ADCK0424
*                                                                     * ADCK0425
*       GAMMA  is  a  parameter  involved  in  the  Kronecker-Picard  * ADCK0426
*         integral. Its default value is taken equal to 1.D0 but the  * ADCK0427
*         user is able to change it in his own main program. If this  * ADCK0428
*         value is less than the machine precision EPSMCH, INTSUB     * ADCK0429
*         sets GAMMA equal to 1.D0.                                   * ADCK0430
*                                                                     * ADCK0431
*                                                                     * ADCK0432
*     Subprograms required :                                          * ADCK0433
*                                                                     * ADCK0434
*       FORTRAN-Supplied ... DATAN, IDNINT, DBLE, DABS                * ADCK0435
*                                                                     * ADCK0436
*       RFSFNS-Supplied .... FNC, G                                   * ADCK0437
*                                                                     * ADCK0438
*                                                                     * ADCK0439
*                                 RFSFNS.  Version of 21-07-1995      * ADCK0440
*                                                                     * ADCK0441
*                                 Department of Mathematics           * ADCK0442
*                                 University of Patras                * ADCK0443
*                                 GR-261.10 Patras                    * ADCK0444
*                                 Greece                              * ADCK0445
*                                 Email : vrahatis@math.upatras.gr    * ADCK0446
*                                                                     * ADCK0447
*********************************************************************** ADCK0448
      IMPLICIT REAL*8(A-H, O-Z), INTEGER*8(I-N)                         ADCK0449
      PARAMETER (ZERO = 0.D0, HALF = 0.5D0, ONE = 1.D0,                 ADCK0450
     +           QUART = 0.25D0, FOUR = 4.D0,                           ADCK0451
     +           MAXI = 1000, K = 3)                                    ADCK0452
      INTEGER*8  INF                                                    ADCK0453
      COMMON / BLK2 / GAMMA                                             ADCK0454
      DIMENSION S(K), H(K), C(K), D(K)                                  ADCK0455
      INTRINSIC DATAN, IDNINT, DBLE, DABS                               ADCK0456
      REAL*8 XNU                                                        ADCK0457
      PI = FOUR*DATAN(ONE)                                              ADCK0458
*                                                                       ADCK0459
      ASAV = A                                                          ADCK0460
      BSAV = B                                                          ADCK0461
      KM1 = K-1                                                         ADCK0462
*                                                                       ADCK0463
      WV = B-A                                                          ADCK0464
*                                                                       ADCK0465
      IF (WV .LE. ONE) THEN                                             ADCK0466
         ITERA = 1                                                      ADCK0467
         HSTEP = WV                                                     ADCK0468
      ELSE                                                              ADCK0469
         ITERA = IDNINT(WV)                                             ADCK0470
         HSTEP = WV/DBLE(ITERA)                                         ADCK0471
      ENDIF                                                             ADCK0472
*                                                                       ADCK0473
*     Evaluate Kronecker-Picard integral.                               ADCK0474
*                                                                       ADCK0475
      SUMINT = ZERO                                                     ADCK0476
      DO ITA = 1, ITERA                                                 ADCK0477
         H(1) = ONE                                                     ADCK0478
         L = 1                                                          ADCK0479
         B = A+HSTEP                                                    ADCK0480
*                                                                       ADCK0481
         DO I = 1, MAXI-1                                               ADCK0482
*                                                                       ADCK0483
            IF (I .EQ. 1) THEN                                          ADCK0484
               FA = G(ICASE,A,XNU,INF)                                          ADCK0485
               IF (INF .EQ. 2) RETURN                                   ADCK0486
               FB = G(ICASE,B,XNU,INF)                                          ADCK0487
               IF (INF .EQ. 2) RETURN                                   ADCK0488
               S(1) = HALF*(B-A)*(FA+FB)                                ADCK0489
            ELSE                                                        ADCK0490
               IST = 2**(I-2)                                           ADCK0491
               ST = DBLE(IST)                                           ADCK0492
               STEP = (B-A)/ST                                          ADCK0493
               X = A+HALF*STEP                                          ADCK0494
               SUM = ZERO                                               ADCK0495
               DO J = 1, IST                                            ADCK0496
                  SUM = SUM+G(ICASE,X,XNU,INF)                                  ADCK0497
                  IF (INF .EQ. 2) RETURN                                ADCK0498
                  X = X+STEP                                            ADCK0499
               ENDDO                                                    ADCK0500
               S(L) = HALF*(S(L-1)+SUM*STEP)                            ADCK0501
            ENDIF                                                       ADCK0502
*                                                                       ADCK0503
            IF (I .GE. K) THEN                                          ADCK0504
               SINT = ZERO                                              ADCK0505
               NS = 1                                                   ADCK0506
               DIF = DABS(SINT-H(1))                                    ADCK0507
               DO N = 1, K                                              ADCK0508
                  DIFT = DABS(SINT-H(N))                                ADCK0509
                  IF (DIFT .LT. DIF) THEN                               ADCK0510
                     NS = N                                             ADCK0511
                     DIF = DIFT                                         ADCK0512
                  ENDIF                                                 ADCK0513
                  C(N) = S(N)                                           ADCK0514
                  D(N) = C(N)                                           ADCK0515
               ENDDO                                                    ADCK0516
               YINT = S(NS)                                             ADCK0517
               NS = NS-1                                                ADCK0518
               DO M = 1, K-1                                            ADCK0519
                  DO N = 1, K-M                                         ADCK0520
                     HO = H(N)-SINT                                     ADCK0521
                     HP = H(N+M)-SINT                                   ADCK0522
                     W = C(N+1)-D(N)                                    ADCK0523
                     DEN = HO-HP                                        ADCK0524
                     IF (DEN .EQ. ZERO) THEN                            ADCK0525
*                                                                       ADCK0526
*     The evaluation of the integral has failed.                        ADCK0527
*                                                                       ADCK0528
                        INF = 3                                         ADCK0529
                        RETURN                                          ADCK0530
                     ENDIF                                              ADCK0531
                     DEN = W/DEN                                        ADCK0532
                     D(N) = HP*DEN                                      ADCK0533
                     C(N) = HO*DEN                                      ADCK0534
                  ENDDO                                                 ADCK0535
                  IF (2*NS .LT. K-M) THEN                               ADCK0536
                     ERREST = C(NS+1)                                   ADCK0537
                  ELSE                                                  ADCK0538
                     ERREST = D(NS)                                     ADCK0539
                     NS = NS-1                                          ADCK0540
                  ENDIF                                                 ADCK0541
                  YINT = YINT+ERREST                                    ADCK0542
               ENDDO                                                    ADCK0543
               IF (DABS(ERREST) .LT. EPSILO*DABS(YINT)) THEN            ADCK0544
*                                                                       ADCK0545
*     The evaluation of the integral has been successful.               ADCK0546
*                                                                       ADCK0547
                  INF = 1                                               ADCK0548
                  SUMINT = SUMINT+YINT                                  ADCK0549
                  GOTO 10                                               ADCK0550
               ENDIF                                                    ADCK0551
               DO N = 1, KM1                                            ADCK0552
                  H(N) = H(N+1)                                         ADCK0553
                  S(N) = S(N+1)                                         ADCK0554
               ENDDO                                                    ADCK0555
               H(K) = QUART*H(KM1)                                      ADCK0556
            ELSE                                                        ADCK0557
               L = L+1                                                  ADCK0558
               H(L) = QUART*H(L-1)                                      ADCK0559
            ENDIF                                                       ADCK0560
*                                                                       ADCK0561
         ENDDO                                                          ADCK0562
*                                                                       ADCK0563
*     The evaluation of the integral has failed.                        ADCK0564
*                                                                       ADCK0565
         INF = 3                                                        ADCK0566
         RETURN                                                         ADCK0567
*                                                                       ADCK0568
   10    A = B                                                          ADCK0569
      ENDDO                                                             ADCK0570
*                                                                       ADCK0571
*     Computation of the topological degree.                            ADCK0572
*                                                                       ADCK0573
      A = ASAV                                                          ADCK0574
      B = BSAV                                                          ADCK0575
      YINT = SUMINT                                                     ADCK0576
                                                                        ADCK0577
      IF (ICASE .LT. 3) THEN                                            ADCK0578
         FA = FNC(ICASE,A,XNU,INF)                                          ADCK0579
         FB = FNC(ICASE,B,XNU,INF)                                          ADCK0580
         FAP = FNC(ICASE+2,A,XNU,INF)                                       ADCK0581
         FBP = FNC(ICASE+2,B,XNU,INF)                                       ADCK0582
         ARCTAN = DATAN(GAMMA*(FA*FBP-FB*FAP)                           ADCK0583
     +                  /(FA*FB+GAMMA*GAMMA*FAP*FBP))                   ADCK0584
      ELSE                                                              ADCK0585
         FA = FNC(ICASE-2,A,XNU,INF)                                        ADCK0586
         FB = FNC(ICASE-2,B,XNU,INF)                                        ADCK0587
         FAP = FNC(ICASE,A,XNU,INF)                                         ADCK0588
         FBP = FNC(ICASE,B,XNU,INF)                                         ADCK0589
         FAPP = -FAP/A-FA*(ONE-XNU*XNU/A/A)                             ADCK0590
         FBPP = -FBP/B-FB*(ONE-XNU*XNU/B/B)                             ADCK0591
         ARCTAN = DATAN(GAMMA*(FAP*FBPP-FBP*FAPP)                       ADCK0592
     +                  /(FAP*FBP+GAMMA*GAMMA*FAPP*FBPP))               ADCK0593
      ENDIF                                                             ADCK0594
*                                                                       ADCK0595
      NR = IDNINT((GAMMA*YINT+ARCTAN)/PI)                               ADCK0596
*                                                                       ADCK0597
      RETURN                                                            ADCK0598
*                                                                       ADCK0599
*     Last statement of the subroutine TOPDEG.                          ADCK0600
*                                                                       ADCK0601
      END                                                               ADCK0602
*---------------------------------------------------------------------* ADCK0603
                                                                        ADCK0604
      SUBROUTINE ISOLAT (ICASE, A, B, NR, MAXRT, ROOTIS,XNU,NRF,INF)    ADCK0605
*********************************************************************** ADCK0606
*                                                                     * ADCK0607
*     SUBROUTINE ISOLAT                                               * ADCK0608
*                                                                     * ADCK0609
*                                                                     * ADCK0610
*       The purpose of this subroutine is to isolate the roots in  NR * ADCK0611
*       subintervals of [A,B]. The left endpoints  of these intervals * ADCK0612
*       are stored  in the first column  of matrix  ROOTIS, while the * ADCK0613
*       right endpoints are stored in the second column.              * ADCK0614
*                                                                     * ADCK0615
*                                                                     * ADCK0616
*     The subroutine statement is :                                   * ADCK0617
*                                                                     * ADCK0618
*     SUBROUTINE ISOLAT (ICASE, A, B, NR, MAXRT, ROOTIS)              * ADCK0619
*                                                                     * ADCK0620
*     where :                                                         * ADCK0621
*                                                                     * ADCK0622
*       ICASE is  an  integer  input  variable  which  specifies  the * ADCK0623
*         considered function as follows :                            * ADCK0624
*                                                                     * ADCK0625
*         ICASE = 1  Bessel function of first kind, J.                * ADCK0626
*                                                                     * ADCK0627
*         ICASE = 2  The derivative of Bessel function of first kind. * ADCK0628
*                                                                     * ADCK0629
*         ICASE = 3  Bessel function of second kind, Y.               * ADCK0630
*                                                                     * ADCK0631
*         ICASE = 4  The derivative of Bessel function of second kind.* ADCK0632
*                                                                     * ADCK0633
*       A is an input variable which defines the left endpoint of the * ADCK0634
*         given interval.                                             * ADCK0635
*                                                                     * ADCK0636
*       B is  an input variable  which defines  the right endpoint of * ADCK0637
*         the given interval.                                         * ADCK0638
*                                                                     * ADCK0639
*       NR is a positive integer input variable  which determines the * ADCK0640
*         number of roots to be isolated in the given interval.       * ADCK0641
*         NR must be less than MAXRT.                                 * ADCK0642
*                                                                     * ADCK0643
*       MAXRT is a positive integer input  variable  which determines * ADCK0644
*         the maximum number of roots required in the given interval. * ADCK0645
*                                                                     * ADCK0646
*      ROOTIS is  an output  MAXRT times 2  array.  Its  two  columns * ADCK0647
*         provide lower and upper root bounds.                        * ADCK0648
*                                                                     * ADCK0649
*                                                                     * ADCK0650
*     The COMMON-transferred variables used here are :                * ADCK0651
*                                                                     * ADCK0652
*       INF is an integer output variable set as follows :            * ADCK0653
*                                                                     * ADCK0654
*         INF = 1    The required answer is accurate.                 * ADCK0655
*                                                                     * ADCK0656
*         INF = 4    The procedure for the isolation of the roots     * ADCK0657
*                    failed because                                   * ADCK0658
*                       a) the maximum number of iterations MAXDIV or * ADCK0659
*                       b) the maximum dimension of the internal      * ADCK0660
*                          working arrays WA1,WA2 of ISOLAT was       * ADCK0661
*                          exceeded.                                  * ADCK0662
*                    BISECT will not be evoked.                       * ADCK0663
*                                                                     * ADCK0664
*       NRF is an integer  output  variable indicating  the number of * ADCK0665
*         roots which have been isolated.                             * ADCK0666
*                                                                     * ADCK0667
*                                                                     * ADCK0668
*     Subprograms required :                                          * ADCK0669
*                                                                     * ADCK0670
*       FORTRAN-Supplied ... DBLE                                     * ADCK0671
*                                                                     * ADCK0672
*       RFSFNS-Supplied .... FNC                                      * ADCK0673
*                                                                     * ADCK0674
*                                                                     * ADCK0675
*                                 RFSFNS.  Version of 21-07-1995      * ADCK0676
*                                                                     * ADCK0677
*                                 Department of Mathematics           * ADCK0678
*                                 University of Patras                * ADCK0679
*                                 GR-261.10 Patras                    * ADCK0680
*                                 Greece                              * ADCK0681
*                                 Email : vrahatis@math.upatras.gr    * ADCK0682
*                                                                     * ADCK0683
*********************************************************************** ADCK0684
      IMPLICIT REAL*8(A-H, O-Z), INTEGER*8(I-N)                         ADCK0685
      PARAMETER (ZERO = 0.D0, HALF = 0.5D0, TWO = 2.D0,                 ADCK0686
     +           MAXIND = 2000, MAXDIV = 5000)                          ADCK0687
      INTEGER*8       INF                                               ADCK0688
      DIMENSION ROOTIS(MAXRT,2), WA1(MAXDIV), WA2(MAXDIV)               ADCK0689
      INTRINSIC DBLE                                                    ADCK0690
      REAL*8 XNU                                                        ADCK0691
      INTEGER*8 NRF                                                     ADCK0692
      NRP1 = NR+1                                                       ADCK0693
      H = (B-A)/DBLE(NR)                                                ADCK0694
      DO I = 1, NRP1                                                    ADCK0695
         WA1(I) = A+DBLE(I-1)*H                                         ADCK0696
         WA2(I) = FNC(ICASE,WA1(I),XNU,INF)                             ADCK0697
      ENDDO                                                             ADCK0698
*     Check whether all the roots have been isolated.                   ADCK0699
*                                                                       ADCK0700
      DO INDEX= 1, MAXIND                                               ADCK0701
                                                                        ADCK0702
         NRF = 0                                                        ADCK0703
         DO I = 1, NRP1-1                                               ADCK0704
            IP1 = I+1                                                   ADCK0705
            IF (WA2(I)*WA2(IP1) .LT. ZERO) THEN                         ADCK0706
               NRF = NRF+1                                              ADCK0707
               ROOTIS(NRF,1) = WA1(I)                                   ADCK0708
               ROOTIS(NRF,2) = WA1(IP1)                                 ADCK0709
            ENDIF                                                       ADCK0710
         ENDDO                                                          ADCK0711
         IF (NRF .GE. NR) THEN                                          ADCK0712
*                                                                       ADCK0713
*     Subroutine ISOLAT succeded.                                       ADCK0714
*                                                                       ADCK0715
            INF = 1                                                     ADCK0716
            RETURN                                                      ADCK0717
         ENDIF                                                          ADCK0718
*                                                                       ADCK0719
*     Rescale the arrays WA1 and WA2.                                   ADCK0720
*     Test also if this procedure violates the dimension MAXDIV.        ADCK0721
*     If this is true return the number of roots and the corresponding  ADCK0722
*     intervals already found.                                          ADCK0723
*                                                                       ADCK0724
         NNRP1 = NRP1+TWO**(INDEX-1)*NR                                 ADCK0725
         IF (NNRP1 .GT. MAXDIV) THEN                                    ADCK0726
            INF = 4                                                     ADCK0727
            RETURN                                                      ADCK0728
         ENDIF                                                          ADCK0729
         J = NRP1                                                       ADCK0730
         DO I = NNRP1, 1, -2                                            ADCK0731
            WA1(I) = WA1(J)                                             ADCK0732
            WA2(I) = WA2(J)                                             ADCK0733
            J = J-1                                                     ADCK0734
         ENDDO                                                          ADCK0735
         H = HALF*H                                                     ADCK0736
         DO I = 2, NNRP1-1, 2                                           ADCK0737
            WA1(I) = WA1(I-1)+H                                         ADCK0738
            WA2(I) = FNC(ICASE,WA1(I),XNU,INF)                              ADCK0739
         ENDDO                                                          ADCK0740
         NRP1 = NNRP1                                                   ADCK0741
                                                                        ADCK0742
      ENDDO                                                             ADCK0743
*                                                                       ADCK0744
*     Subroutine ISOLAT exceeded maximum of steps.                      ADCK0745
      PRINT*,"max step exceeded"                                           ADCK0746
      INF = 4                                                           ADCK0747
                                                                        ADCK0748
      RETURN                                                            ADCK0749
*                                                                       ADCK0750
*     Last statement of the subroutine ISOLAT.                          ADCK0751
*                                                                       ADCK0752
      END                                                               ADCK0753
*---------------------------------------------------------------------* ADCK0754
                                                                        ADCK0755
      SUBROUTINE BISECT (ICASE, A, B, DELTA, ROOT, FROOT,XNU,INF)       ADCK0756
*********************************************************************** ADCK0757
*                                                                     * ADCK0758
*     SUBROUTINE BISECT                                               * ADCK0759
*                                                                     * ADCK0760
*                                                                     * ADCK0761
*       The purpose of the subroutine  BISECT  is to compute the root * ADCK0762
*       which is isolated in the interval (A,B).                      * ADCK0763
*                                                                     * ADCK0764
*                                                                     * ADCK0765
*    The subroutine statement is :                                    * ADCK0766
*                                                                     * ADCK0767
*     SUBROUTINE BISECT (ICASE, A, B, DELTA, ROOT, FROOT)             * ADCK0768
*                                                                     * ADCK0769
*     where :                                                         * ADCK0770
*                                                                     * ADCK0771
*       ICASE  is  an integer  input  variable  which  specifies  the * ADCK0772
*         considered function as follows :                            * ADCK0773
*                                                                     * ADCK0774
*         ICASE = 1  Bessel function of first kind, J.                * ADCK0775
*                                                                     * ADCK0776
*         ICASE = 2  The derivative of Bessel function of first kind. * ADCK0777
*                                                                     * ADCK0778
*         ICASE = 3  Bessel function of second kind, Y.               * ADCK0779
*                                                                     * ADCK0780
*         ICASE = 4  The derivative of Bessel function of second kind.* ADCK0781
*                                                                     * ADCK0782
*       A is an input variable which defines the left endpoint of the * ADCK0783
*         given interval.                                             * ADCK0784
*                                                                     * ADCK0785
*       B is an input variable  which  defines  the right endpoint of * ADCK0786
*         the given interval.                                         * ADCK0787
*                                                                     * ADCK0788
*       DELTA is a  positive  input  variable  which  determines  the * ADCK0789
*         accuracy of the computation of the roots of the function.   * ADCK0790
*         If DELTA is less than the machine precision  EPSMCH,  DELTA * ADCK0791
*         becomes equal to  EPSMCH.  The value of  EPSMCH is computed * ADCK0792
*         within INTSUB.                                              * ADCK0793
*                                                                     * ADCK0794
*       ROOT  is  an  output  variable  which  determines  the  final * ADCK0795
*          approximate root.                                          * ADCK0796
*                                                                     * ADCK0797
*       FROOT is  an output  variable  which  specifies  the function * ADCK0798
*          value at ROOT.                                             * ADCK0799
*                                                                     * ADCK0800
*     The COMMON-transferred variables used here are :                * ADCK0801
*                                                                     * ADCK0802
*       INF is an integer output variable set as follows :            * ADCK0803
*                                                                     * ADCK0804
*         INF = 1    The required answer is accurate.                 * ADCK0805
*                                                                     * ADCK0806
*         INF = 5    The  rootfinding  procedure  failed  because the * ADCK0807
*                    maximum number of steps of BISECT was exceeded.  * ADCK0808
*                                                                     * ADCK0809
*                                                                     * ADCK0810
*     Subprograms required :                                          * ADCK0811
*                                                                     * ADCK0812
*       FORTRAN-Supplied ... IDNINT, DLOG10, DSIGN, DBLE, DABS        * ADCK0813
*                                                                     * ADCK0814
*       RFSFNS-Supplied .... FNC                                      * ADCK0815
*                                                                     * ADCK0816
*                                                                     * ADCK0817
*                                 RFSFNS.  Version of 21-07-1995      * ADCK0818
*                                                                     * ADCK0819
*                                 Department of Mathematics           * ADCK0820
*                                 University of Patras                * ADCK0821
*                                 GR-261.10 Patras                    * ADCK0822
*                                 Greece                              * ADCK0823
*                                 Email : vrahatis@math.upatras.gr    * ADCK0824
*                                                                     * ADCK0825
*********************************************************************** ADCK0826
      IMPLICIT REAL*8(A-H, O-Z), INTEGER*8(I-N)                         ADCK0827
      PARAMETER (ONE = 1.D0, TWO = 2.D0)                                ADCK0828
      INTEGER*8 INF                                                     ADCK0829
      INTRINSIC IDNINT, DLOG10, DSIGN, DBLE, DABS                       ADCK0830
      REAL*8 XNU                                                        ADCK0831
      STEP = B-A                                                        ADCK0832
*                                                                       ADCK0833
*     Compute the number of iterations MAXI which are required          ADCK0834
*     to obtain a root to the predetermined accuracy, DELTA.            ADCK0835
*                                                                       ADCK0836
      MAXI = IDNINT(DLOG10(STEP/DELTA)/DLOG10(TWO))+1                   ADCK0837
*                                                                       ADCK0838
      SA = DSIGN(ONE,FNC(ICASE,A,XNU,INF))*STEP                                 ADCK0839
      XOLD = A                                                          ADCK0840
*                                                                       ADCK0841
*     Begin the iterations.                                             ADCK0842
*                                                                       ADCK0843
      DO I = 1, MAXI                                                    ADCK0844
*                                                                       ADCK0845
*        Compute the new approximate root.                              ADCK0846
*                                                                       ADCK0847
         ROOT=XOLD+SA*DSIGN(ONE,FNC(ICASE,XOLD,XNU,INF))/(TWO**DBLE(I)) ADCK0848
*                                                                       ADCK0849
*        Specify a stopping test using the difference                   ADCK0850
*        between two successive iterations.                             ADCK0851
*                                                                       ADCK0852
         IF (DABS(ROOT-XOLD) .LE. DELTA) THEN                           ADCK0853
            INF = 1                                                     ADCK0854
            FROOT = FNC(ICASE,ROOT,XNU,INF)                                     ADCK0855
            RETURN                                                      ADCK0856
         ENDIF                                                          ADCK0857
         XOLD = ROOT                                                    ADCK0858
      ENDDO                                                             ADCK0859
      INF = 5                                                           ADCK0860
      RETURN                                                            ADCK0861
*                                                                       ADCK0862
*     Last statement of the subroutine BISECT.                          ADCK0863
*                                                                       ADCK0864
      END                                                               ADCK0865
*---------------------------------------------------------------------* ADCK0866
                                                                        ADCK0867
      REAL*8 FUNCTION FNC (ICASE, X, XNU, INF)                          ADCK0868
*********************************************************************** ADCK0869
*                                                                     * ADCK0870
*     FUNCTION FNC                                                    * ADCK0871
*                                                                     * ADCK0872
*                                                                     * ADCK0873
*       The purpose of this subprogram is to compute the value of the * ADCK0874
*     specified  Bessel function.  Here the Bessel functions of first * ADCK0875
*     and second kind as well as their first derivatives are examined.* ADCK0876
*                                                                     * ADCK0877
*                                                                     * ADCK0878
*     The function statement is :                                     * ADCK0879
*                                                                     * ADCK0880
*     REAL*8 FUNCTION FNC (ICASE, X)                                  * ADCK0881
*                                                                     * ADCK0882
*     where :                                                         * ADCK0883
*                                                                     * ADCK0884
*       ICASE is  an  integer  input  variable  which  specifies  the * ADCK0885
*         considered function as follows :                            * ADCK0886
*                                                                     * ADCK0887
*         ICASE = 1  Bessel function of first kind, J.                * ADCK0888
*                                                                     * ADCK0889
*         ICASE = 2  The derivative of Bessel function of first kind. * ADCK0890
*                                                                     * ADCK0891
*         ICASE = 3  Bessel function of second kind, Y.               * ADCK0892
*                                                                     * ADCK0893
*         ICASE = 4  The derivative of Bessel function of second kind.* ADCK0894
*                                                                     * ADCK0895
*       X is the value of the independent variable.                   * ADCK0896
*                                                                     * ADCK0897
*                                                                     * ADCK0898
*     The COMMON-transferred variables used here are :                * ADCK0899
*                                                                     * ADCK0900
*       INF is an integer output variable set as follows :            * ADCK0901
*                                                                     * ADCK0902
*         INF = 1    The required answer is accurate.                 * ADCK0903
*                                                                     * ADCK0904
*         INF = 2    The evaluation of the function FNC has failed.   * ADCK0905
*                                                                     * ADCK0906
*       XNU is   an input  variable  which  defines  the order of the * ADCK0907
*         Bessel function.                                            * ADCK0908
*                                                                     * ADCK0909
*                                                                     * ADCK0910
*     Subprograms required :                                          * ADCK0911
*                                                                     * ADCK0912
*       RFSFNS-Supplied .... COULFG                                   * ADCK0913
*                            (by A.R. Barnett, 1982, CPC, 27, 147)    * ADCK0914
*                                                                     * ADCK0915
*                                                                     * ADCK0916
*                                 RFSFNS.  Version of 21-07-1995      * ADCK0917
*                                                                     * ADCK0918
*                                 Department of Mathematics           * ADCK0919
*                                 University of Patras                * ADCK0920
*                                 GR-261.10 Patras                    * ADCK0921
*                                 Greece                              * ADCK0922
*                                 Email : vrahatis@math.upatras.gr    * ADCK0923
*                                                                     * ADCK0924
*********************************************************************** ADCK0925
      IMPLICIT REAL*8(A-H, O-Z), INTEGER*8(I-N)                         ADCK0926
      PARAMETER (ZERO = 0.D0, ONE = 1.D0)                               ADCK0927
      INTEGER*8 INF                                                     ADCK0928
      COMMON / BLK2 /   GAMMA                                           ADCK0929
      DIMENSION FC(3),GC(3),FCP(3),GCP(3)                               ADCK0930
      REAL*8 XNU                                                        ADCK0931
      INF=1                                                             ADCK0932
*                                                                       ADCK0933
*     Set the starting values for COULFG.                               ADCK0934
*                                                                       ADCK0935
      KFN = 2                                                           ADCK0936
      ETA = ZERO                                                        ADCK0937
      XLMIN = XNU                                                       ADCK0938
      XLMAX = XLMIN                                                     ADCK0939
      IF (ICASE .EQ. 1) THEN                                            ADCK0940
         MODE = 3                                                       ADCK0941
      ELSE IF (ICASE .EQ. 2) THEN                                       ADCK0942
         MODE = 2                                                       ADCK0943
      ELSE                                                              ADCK0944
         MODE = 1                                                       ADCK0945
      ENDIF                                                             ADCK0946
*                                                                       ADCK0947
*     Evoke COULFG.                                                     ADCK0948
*                                                                       ADCK0949
      CALL COULFG(X,ETA,XLMIN,XLMAX,FC,GC,FCP,GCP,MODE,KFN,IFAIL)       ADCK0950
*                                                                       ADCK0951
*     Check for possible errors reported by COULFG.                     ADCK0952
*                                                                       ADCK0953
      IF (IFAIL .EQ. 0) THEN                                            ADCK0954
         INDEX = 1                                                      ADCK0955
      ELSE IF (IFAIL .GE. 1) THEN                                       ADCK0956
*                                                                       ADCK0957
*        In case COULFG fails recalculate the functions using also      ADCK0958
*        the order XNU+1.                                               ADCK0959
*                                                                       ADCK0960
         XLMAX = XLMAX+ONE                                              ADCK0961
         CALL COULFG(X,ETA,XLMIN,XLMAX,FC,GC,FCP,GCP,MODE,KFN,IFAIL)    ADCK0962
         IF (IFAIL .EQ. 0) THEN                                         ADCK0963
            INDEX = 2                                                   ADCK0964
         ELSE                                                           ADCK0965
            INF = 2                                                     ADCK0966
            RETURN                                                      ADCK0967
         END IF                                                         ADCK0968
      ELSE                                                              ADCK0969
         INF = 2                                                        ADCK0970
         RETURN                                                         ADCK0971
      END IF                                                            ADCK0972
*                                                                       ADCK0973
      IF (ICASE .EQ. 1) THEN                                            ADCK0974
         FNC = FC(INDEX)                                                ADCK0975
      ELSE IF (ICASE .EQ. 2) THEN                                       ADCK0976
         FNC = GC(INDEX)                                                ADCK0977
      ELSE IF (ICASE .EQ. 3) THEN                                       ADCK0978
         FNC = FCP(INDEX)                                               ADCK0979
      ELSE                                                              ADCK0980
         FNC = GCP(INDEX)                                               ADCK0981
      ENDIF                                                             ADCK0982
*                                                                       ADCK0983
      RETURN                                                            ADCK0984
*                                                                       ADCK0985
*     Last statement of the function FNC.                               ADCK0986
*                                                                       ADCK0987
      END                                                               ADCK0988
*---------------------------------------------------------------------* ADCK0989
                                                                        ADCK0990
      REAL*8 FUNCTION G(ICASE,X, XNU,INF)                               ADCK0991
*********************************************************************** ADCK0992
*                                                                     * ADCK0993
*     FUNCTION G                                                      * ADCK0994
*                                                                     * ADCK0995
*                                                                     * ADCK0996
*       The purpose of this subprogram is to compute the value of the * ADCK0997
*     integrand function G for the computation of the Kronecker-Picard* ADCK0998
*     integral.                                                       * ADCK0999
*                                                                     * ADCK1000
*                                                                     * ADCK1001
*     The function statement is :                                     * ADCK1002
*                                                                     * ADCK1003
*     REAL*8 FUNCTION G (ICASE, X)                                    * ADCK1004
*                                                                     * ADCK1005
*     where :                                                         * ADCK1006
*                                                                     * ADCK1007
*       ICASE is  an  integer  input  variable  which  specifies  the * ADCK1008
*         considered function as follows :                            * ADCK1009
*                                                                     * ADCK1010
*         ICASE = 1  Bessel function of first kind, J.                * ADCK1011
*                                                                     * ADCK1012
*         ICASE = 2  The derivative of Bessel function of first kind. * ADCK1013
*                                                                     * ADCK1014
*         ICASE = 3  Bessel function of second kind, Y.               * ADCK1015
*                                                                     * ADCK1016
*         ICASE = 4  The derivative of Bessel function of second kind.* ADCK1017
*                                                                     * ADCK1018
*       X is the value of the independent variable.                   * ADCK1019
*                                                                     * ADCK1020
*                                                                     * ADCK1021
*     The COMMON-transferred variables used here are :                * ADCK1022
*                                                                     * ADCK1023
*       INF is an integer output variable set as follows :            * ADCK1024
*                                                                     * ADCK1025
*         INF = 1    The required answer is accurate.                 * ADCK1026
*                                                                     * ADCK1027
*         INF = 2    The evaluation of the function G has failed.     * ADCK1028
*                                                                     * ADCK1029
*       XNU is an input variable which defines the order of the Bessel* ADCK1030
*         function.                                                   * ADCK1031
*                                                                     * ADCK1032
*       GAMMA  is  a  parameter  involved  in  the  Kronecker-Picard  * ADCK1033
*         integral. Its default value is taken equal to 1.D0 but the  * ADCK1034
*         user is able to change it in his own main program. If this  * ADCK1035
*         value is less than the machine precision EPSMCH, INTSUB     * ADCK1036
*         sets GAMMA equal to 1.D0.                                   * ADCK1037
*                                                                     * ADCK1038
*                                                                     * ADCK1039
*     Subprograms required :                                          * ADCK1040
*                                                                     * ADCK1041
*       RFSFNS-Supplied .... COULFG                                   * ADCK1042
*                            (by A.R. Barnett, 1982, CPC, 27, 147)    * ADCK1043
*                                                                     * ADCK1044
*                                                                     * ADCK1045
*                                 RFSFNS.  Version of 21-07-1995      * ADCK1046
*                                                                     * ADCK1047
*                                 Department of Mathematics           * ADCK1048
*                                 University of Patras                * ADCK1049
*                                 GR-261.10 Patras                    * ADCK1050
*                                 Greece                              * ADCK1051
*                                 Email : vrahatis@math.upatras.gr    * ADCK1052
*                                                                     * ADCK1053
*********************************************************************** ADCK1054
      IMPLICIT REAL*8(A-H, O-Z), INTEGER*8(I-N)                         ADCK1055
      PARAMETER (ZERO = 0.D0, ONE = 1.D0, TWO = 2.D0, THREE = 3.D0)     ADCK1056
      INTEGER*8  INF                                                    ADCK1057
      COMMON / BLK2 / GAMMA                                             ADCK1058
      DIMENSION FC(3),GC(3),FCP(3),GCP(3)                               ADCK1059
      REAL*8 XNU                                                        ADCK1060
      INF = 1                                                           ADCK1061
*                                                                       ADCK1062
*     Set the starting values for COULFG.                               ADCK1063
*                                                                       ADCK1064
      KFN = 2                                                           ADCK1065
      ETA = ZERO                                                        ADCK1066
      MODE = 1                                                          ADCK1067
      XLMIN = XNU                                                       ADCK1068
      XLMAX = XLMIN                                                     ADCK1069
*                                                                       ADCK1070
*     Evoke COULFG.                                                     ADCK1071
*                                                                       ADCK1072
      CALL COULFG(X,ETA,XLMIN,XLMAX,FC,GC,FCP,GCP,MODE,KFN,IFAIL)       ADCK1073
*                                                                       ADCK1074
*     Check for possible errors reported by COULFG.                     ADCK1075
*                                                                       ADCK1076
      IF (IFAIL .EQ. 0) THEN                                            ADCK1077
         INDEX = 1                                                      ADCK1078
      ELSE IF (IFAIL .GE. 1) THEN                                       ADCK1079
*                                                                       ADCK1080
*        In case COULFG fails recalculate the functions using also      ADCK1081
*        the order XNU+1.                                               ADCK1082
*                                                                       ADCK1083
         XLMAX = XLMAX+ONE                                              ADCK1084
         CALL COULFG(X,ETA,XLMIN,XLMAX,FC,GC,FCP,GCP,MODE,KFN,IFAIL)    ADCK1085
         IF (IFAIL .EQ. 0) THEN                                         ADCK1086
            INDEX = 2                                                   ADCK1087
         ELSE                                                           ADCK1088
            INF = 2                                                     ADCK1089
            RETURN                                                      ADCK1090
         END IF                                                         ADCK1091
      ELSE                                                              ADCK1092
         INF = 2                                                        ADCK1093
         RETURN                                                         ADCK1094
      END IF                                                            ADCK1095
*                                                                       ADCK1096
      IF (ICASE .EQ. 1 .OR. ICASE .EQ. 3) THEN                          ADCK1097
         F = FC(INDEX)                                                  ADCK1098
         FP = FCP(INDEX)                                                ADCK1099
      ELSE                                                              ADCK1100
         F = GC(INDEX)                                                  ADCK1101
         FP = GCP(INDEX)                                                ADCK1102
      ENDIF                                                             ADCK1103
*                                                                       ADCK1104
      COE1 = ONE/X                                                      ADCK1105
      COE2 = COE1*COE1                                                  ADCK1106
      COE3 = ONE-XNU*XNU*COE2                                           ADCK1107
      FPP = -COE1*FP-COE3*F                                             ADCK1108
      F2 = F*F                                                          ADCK1109
      FP2 = FP*FP                                                       ADCK1110
      FPP2 = FPP*FPP                                                    ADCK1111
*                                                                       ADCK1112
      IF (ICASE .EQ. 1 .OR. ICASE .EQ. 2) THEN                          ADCK1113
         G = (FP2-F*FPP)/(F2+GAMMA*GAMMA*FP2)                           ADCK1114
      ELSE                                                              ADCK1115
         FPPP = (TWO*COE2-COE3)*FP+COE1*(THREE*COE3-TWO)*F              ADCK1116
         G = (FPP2-FP*FPPP)/(FP2+GAMMA*GAMMA*FPP2)                      ADCK1117
      ENDIF                                                             ADCK1118
*                                                                       ADCK1119
      RETURN                                                            ADCK1120
*                                                                       ADCK1121
*     Last statement of the function G.                                 ADCK1122
*                                                                       ADCK1123
      END                                                               ADCK1124
*---------------------------------------------------------------------* ADCK1125
                                                                        ADCK1126
      SUBROUTINE COULFG(XX,ETA1,XLMIN,XLMAX, FC,GC,FCP,GCP,             ADCK1127
     *                  MODE1,KFN,IFAIL)                                ADCK1128
C                                                                       ADCK1129
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCADCK1130
C                                                                      CADCK1131
C  REVISED COULOMB WAVEFUNCTION PROGRAM USING STEED'S METHOD           CADCK1132
C                                                                      CADCK1133
C  A. R. BARNETT           MANCHESTER  MARCH   1981                    CADCK1134
C                                                                      CADCK1135
C  ORIGINAL PROGRAM 'RCWFN'      IN    CPC   8 (1974) 377-395          CADCK1136
C                 + 'RCWFN'      IN    CPC  11 (1976) 141-142          CADCK1137
C  FULL DESCRIPTION OF ALGORITHM IN    CPC  21 (1981) 297-314          CADCK1138
C  THIS VERSION WRITTEN UP       IN    CPC  27 (1982) 147-166          CADCK1139
*                                                                      *ADCK1140
*  Modified to be used with RFSFNS on July 21, 1995.                   *ADCK1141
*                                                                      *ADCK1142
C                                                                      CADCK1143
C  COULFG RETURNS F,G,F',G', FOR REAL XX.GT.0,REAL ETA1 (INCLUDING 0), CADCK1144
C   AND REAL LAMBDA(XLMIN) .GT. -1 FOR INTEGER-SPACED LAMBDA VALUES    CADCK1145
C   THUS GIVING POSITIVE-ENERGY SOLUTIONS TO THE COULOMB SCHRODINGER   CADCK1146
C   EQUATION,TO THE KLEIN-GORDON EQUATION AND TO SUITABLE FORMS OF     CADCK1147
C   THE DIRAC EQUATION ,ALSO SPHERICAL & CYLINDRICAL BESSEL EQUATIONS  CADCK1148
C                                                                      CADCK1149
C  FOR A RANGE OF LAMBDA VALUES (XLMAX - XLMIN) MUST BE AN INTEGER,    CADCK1150
C  STARTING ARRAY ELEMENT IS M1 = MAX0(IDINT(XLMIN+ACCUR),0) + 1       CADCK1151
C      SEE TEXT FOR MODIFICATIONS FOR INTEGER L-VALUES                 CADCK1152
C                                                                      CADCK1153
C  IF 'MODE' = 1  GET F,G,F',G'   FOR INTEGER-SPACED LAMBDA VALUES     CADCK1154
C            = 2      F,G      UNUSED ARRAYS MUST BE DIMENSIONED IN    CADCK1155
C            = 3      F               CALL TO AT LEAST LENGTH (1)      CADCK1156
C  IF 'KFN'  = 0 REAL        COULOMB FUNCTIONS ARE RETURNED            CADCK1157
C            = 1 SPHERICAL   BESSEL      >>     >>    >>               CADCK1158
C            = 2 CYLINDRICAL BESSEL      >>     >>    >>               CADCK1159
C  THE USE OF 'MODE' AND 'KFN' IS INDEPENDENT                          CADCK1160
C                                                                      CADCK1161
C  PRECISION:  RESULTS TO WITHIN 2-3 DECIMALS OF 'MACHINE ACCURACY'    CADCK1162
C   IN OSCILLATING REGION X .GE. ETA1 + SQRT(ETA1**2 + XLM(XLM+1))     CADCK1163
C   COULFG IS CODED FOR REAL*8 ON IBM OR EQUIVALENT  ACCUR = 10**-16   CADCK1164
C   USE AUTODBL + EXTENDED PRECISION ON HX COMPILER  ACCUR = 10**-33   CADCK1165
C   FOR MANTISSAS OF 56 & 112 BITS. FOR SINGLE PRECISION CDC (48 BITS) CADCK1166
C   REASSIGN DSQRT=SQRT ETC.  SEE TEXT FOR COMPLEX ARITHMETIC VERSION  CADCK1167
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCADCK1168
C                                                                       ADCK1169
      IMPLICIT REAL*8 (A-H,O-Z)                                         ADCK1170
      DIMENSION    FC(1),GC(1),FCP(1),GCP(1)                            ADCK1171
      LOGICAL      ETANEO,XLTURN                                        ADCK1172
      COMMON       /STEED/ PACCQ,NFP,NPQ,IEXP,M1                        ADCK1173
C***  COMMON BLOCK IS FOR INFORMATION ONLY.  NOT REQUIRED IN CODE       ADCK1174
C***  COULFG HAS CALLS TO: DSQRT,DABS,DMOD,IDINT,DSIGN,DFLOAT,DMINI     ADCK1175
*********************************************************************** ADCK1176
*     The following command has been added in the original code.        ADCK1177
*     COMMON block BLK3 is used to transfer the machine accuracy        ADCK1178
*     EPSMCH which has been computed within INTSUB.                     ADCK1179
*                                                                       ADCK1180
      COMMON / BLK3 / EPSMCH                                            ADCK1181
*********************************************************************** ADCK1182
      DATA ZERO,ONE,TWO,TEN2,ABORT / 0.0D0, 1.0D0, 2.0D0, 1.0D2, 2.0D4 /ADCK1183
      DATA HALF,TM30 / 0.5D0, 1.0D-30 /                                 ADCK1184
      DATA RT2DPI    / 0.79788456080286535587989211986876373D0/         ADCK1185
C *** THIS CONSTANT IS  DSQRT(TWO/PI):  USE Q0 FOR IBM REAL*16: D0 FOR  ADCK1186
C *** REAL*8 & CDC DOUBLE P:  E0 FOR CDC SINGLE P; AND TRUNCATE VALUE.  ADCK1187
C                                                                       ADCK1188
*********************************************************************** ADCK1189
*     The following command of the original code has been modified.     ADCK1190
*                       ACCUR = 1.0D-16                                 ADCK1191
*     The accuracy is set equal to the machine accuracy.                ADCK1192
*********************************************************************** ADCK1193
                        ACCUR = EPSMCH                                  ADCK1194
*********************************************************************** ADCK1195
C ***            CHANGE ACCUR TO SUIT MACHINE AND PRECISION REQUIRED    ADCK1196
      MODE  = 1                                                         ADCK1197
      IF(MODE1 .EQ. 2 .OR. MODE1 .EQ. 3 ) MODE = MODE1                  ADCK1198
      IFAIL = 0                                                         ADCK1199
      IEXP  = 1                                                         ADCK1200
      NPQ   = 0                                                         ADCK1201
      ETA   = ETA1                                                      ADCK1202
      GJWKB = ZERO                                                      ADCK1203
      PACCQ = ONE                                                       ADCK1204
      IF(KFN .NE. 0) ETA = ZERO                                         ADCK1205
                 ETANEO  = ETA .NE. ZERO                                ADCK1206
      ACC   = ACCUR                                                     ADCK1207
      ACC4  = ACC*TEN2*TEN2                                             ADCK1208
      ACCH  = DSQRT(ACC)                                                ADCK1209
C ***    TEST RANGE OF XX, EXIT IF.LE.DSQRT(ACCUR) OR IF NEGATIVE       ADCK1210
C                                                                       ADCK1211
      IF(XX .LE. ACCH)                          GO TO 100               ADCK1212
      X     = XX                                                        ADCK1213
      XLM   = XLMIN                                                     ADCK1214
      IF(KFN .EQ. 2)  XLM = XLM - HALF                                  ADCK1215
      IF(XLM .LE. -ONE .OR. XLMAX .LT. XLMIN)   GO TO 105               ADCK1216
      E2MM1 = ETA*ETA + XLM*XLM + XLM                                   ADCK1217
      XLTURN= X*(X - TWO*ETA) .LT. XLM*XLM + XLM                        ADCK1218
      DELL  = XLMAX - XLMIN + ACC                                       ADCK1219
      IF(DABS(DMOD(DELL,ONE)) .GT. ACC) WRITE(6,2040)XLMAX,XLMIN,DELL   ADCK1220
      LXTRA = IDINT(DELL)                                               ADCK1221
      XLL   = XLM + DFLOAT(LXTRA)                                       ADCK1222
C ***       LXTRA IS NUMBER OF ADDITIONAL LAMBDA VALUES TO BE COMPUTED  ADCK1223
C ***       XLL  IS MAX LAMBDA VALUE, OR 0.5 SMALLER FOR J,Y BESSELS    ADCK1224
C ***         DETERMINE STARTING ARRAY ELEMENT (M1) FROM XLMIN          ADCK1225
*********************************************************************** ADCK1226
*     The following command of the original code was replaced by the    ADCK1227
*     next one so that arrays of low dimension can be used.             ADCK1228
*                                                                       ADCK1229
*     M1    = MAX0(IDINT(XLMIN +ACC),0) + 1                             ADCK1230
*********************************************************************** ADCK1231
      M1    = MAX0(IDINT(DELL + ACC),0) + 1                             ADCK1232
*********************************************************************** ADCK1233
      L1    = M1 + LXTRA                                                ADCK1234
C                                                                       ADCK1235
C ***   EVALUATE CF1  =  F  \ FPRIME(XL,ETA,X)/F(XL,ETA,X)              ADCK1236
C                                                                       ADCK1237
      XI  = ONE/X                                                       ADCK1238
      FCL = ONE                                                         ADCK1239
      PK  = XLL + ONE                                                   ADCK1240
      PX  = PK  + ABORT                                                 ADCK1241
    2 EK  = ETA / PK                                                    ADCK1242
      F   = (EK + PK*XI)*FCL + (FCL - ONE)*XI                           ADCK1243
      PK1 = PK + ONE                                                    ADCK1244
C ***   TEST ENSURES B1 .NE. ZERO FOR NEGATIVE ETA; FIXUP IS EXACT.     ADCK1245
             IF(DABS(ETA*X + PK*PK1) .GT. ACC)  GO TO 3                 ADCK1246
             FCL  = (ONE + EK*EK)/(ONE + (ETA/PK1)**2)                  ADCK1247
             PK   =  TWO + PK                                           ADCK1248
      GO TO 2                                                           ADCK1249
    3 D   = ONE/((PK + PK1)*(XI + EK/PK1))                              ADCK1250
      DF  = -FCL*(ONE + EK*EK)*D                                        ADCK1251
            IF(FCL .NE. ONE )  FCL = -ONE                               ADCK1252
            IF(D   .LT. ZERO)  FCL = -FCL                               ADCK1253
      F   = F  + DF                                                     ADCK1254
C                                                                       ADCK1255
C ***   BEGIN CF1 LOOP ON PK = K = LAMBDA + 1                           ADCK1256
C                                                                       ADCK1257
      P     = ONE                                                       ADCK1258
    4 PK    = PK1                                                       ADCK1259
        PK1 = PK1 + ONE                                                 ADCK1260
        EK  = ETA / PK                                                  ADCK1261
        TK  = (PK + PK1)*(XI + EK/PK1)                                  ADCK1262
        D   =  TK - D*(ONE + EK*EK)                                     ADCK1263
              IF(DABS(D) .GT. ACCH)             GO TO 5                 ADCK1264
C             WRITE(*,1000) D,DF,ACCH,PK,EK,ETA,X                       ADCK1265
              P = P  +   ONE                                            ADCK1266
              IF( P .GT. TWO)                   GO TO 110               ADCK1267
*********************************************************************** ADCK1268
*     The following command of the original code was replaced by the    ADCK1269
*     next five commands to avoid interruption of the execution due     ADCK1270
*     to division by zero.                                              ADCK1271
*                                                                       ADCK1272
*   5 D     = ONE/D                                                     ADCK1273
*********************************************************************** ADCK1274
    5 IF (DABS(D) .LT. EPSMCH) THEN                                     ADCK1275
        IFAIL = 1                                                       ADCK1276
        RETURN                                                          ADCK1277
      ENDIF                                                             ADCK1278
      D     = ONE/D                                                     ADCK1279
*********************************************************************** ADCK1280
              IF (D .LT. ZERO) FCL = -FCL                               ADCK1281
        DF  = DF*(D*TK - ONE)                                           ADCK1282
        F   = F  + DF                                                   ADCK1283
              IF(PK .GT. PX)                    GO TO 110               ADCK1284
      IF(DABS(DF) .GE. DABS(F)*ACC)             GO TO 4                 ADCK1285
                  NFP = PK - XLL - 1                                    ADCK1286
      IF(LXTRA .EQ.0)                           GO TO 7                 ADCK1287
C                                                                       ADCK1288
C *** DOWNWARD RECURRENCE TO LAMBDA = XLM. ARRAY GC, IF PRESENT, STORES ADCK1289
C                                                                       ADCK1290
      FCL = FCL*TM30                                                    ADCK1291
      FPL = FCL*F                                                       ADCK1292
      IF(MODE .EQ. 1) FCP(L1) = FPL                                     ADCK1293
                      FC (L1) = FCL                                     ADCK1294
      XL  = XLL                                                         ADCK1295
      RL  = ONE                                                         ADCK1296
      EL  = ZERO                                                        ADCK1297
      DO 6  LP = 1,LXTRA                                                ADCK1298
         IF(ETANEO) EL = ETA/XL                                         ADCK1299
         IF(ETANEO) RL = DSQRT(ONE + EL*EL)                             ADCK1300
         SL    =  EL  + XL*XI                                           ADCK1301
         L     =  L1 - LP                                               ADCK1302
         FCL1  = (FCL *SL + FPL)/RL                                     ADCK1303
         FPL   =  FCL1*SL - FCL *RL                                     ADCK1304
         FCL   =  FCL1                                                  ADCK1305
         FC(L) =  FCL                                                   ADCK1306
         IF(MODE .EQ. 1) FCP(L) = FPL                                   ADCK1307
         IF(MODE .NE. 3 .AND. ETANEO) GC(L+1) = RL                      ADCK1308
    6 XL = XL - ONE                                                     ADCK1309
      IF(FCL .EQ. ZERO) FCL = ACC                                       ADCK1310
      F  = FPL/FCL                                                      ADCK1311
C ***    NOW WE HAVE REACHED LAMBDA = XLMIN = XLM                       ADCK1312
C ***    EVALUATE CF2 = P + I.Q AGAIN USING STEED'S ALGORITHM           ADCK1313
C ***    SEE TEXT FOR COMPACT COMPLEX CODE FOR SP CDC OR NON-ANSI IBM   ADCK1314
C                                                                       ADCK1315
    7 IF( XLTURN ) CALL JWKB(X,ETA,DMAX1(XLM,ZERO),FJWKB,GJWKB,IEXP)    ADCK1316
      IF( IEXP .GT. 1 .OR. GJWKB .GT. ONE/(ACCH*TEN2))  GO TO 9         ADCK1317
          XLTURN = .FALSE.                                              ADCK1318
      TA =  TWO*ABORT                                                   ADCK1319
      PK =  ZERO                                                        ADCK1320
      WI =  ETA + ETA                                                   ADCK1321
      P  =  ZERO                                                        ADCK1322
      Q  =  ONE - ETA*XI                                                ADCK1323
      AR = -E2MM1                                                       ADCK1324
      AI =  ETA                                                         ADCK1325
      BR =  TWO*(X - ETA)                                               ADCK1326
      BI =  TWO                                                         ADCK1327
      DR =  BR/(BR*BR + BI*BI)                                          ADCK1328
      DI = -BI/(BR*BR + BI*BI)                                          ADCK1329
      DP = -XI*(AR*DI + AI*DR)                                          ADCK1330
      DQ =  XI*(AR*DR - AI*DI)                                          ADCK1331
    8 P     = P  + DP                                                   ADCK1332
         Q  = Q  + DQ                                                   ADCK1333
         PK = PK + TWO                                                  ADCK1334
         AR = AR +PK                                                    ADCK1335
         AI = AI + WI                                                   ADCK1336
         BI = BI + TWO                                                  ADCK1337
         D  = AR*DR - AI*DI + BR                                        ADCK1338
         DI = AI*DR + AR*DI + BI                                        ADCK1339
         C  = ONE/(D*D + DI*DI)                                         ADCK1340
         DR =  C*D                                                      ADCK1341
         DI = -C*DI                                                     ADCK1342
         A  = BR*DR - BI*DI - ONE                                       ADCK1343
         B  = BI*DR + BR*DI                                             ADCK1344
         C  = DP*A  - DQ*B                                              ADCK1345
         DQ = DP*B  + DQ*A                                              ADCK1346
         DP = C                                                         ADCK1347
         IF(PK .GT. TA)                         GO TO 120               ADCK1348
      IF(DABS(DP)+DABS(DQ).GE.(DABS(P)+DABS(Q))*ACC)   GO TO 8          ADCK1349
                      NPQ   = PK/TWO                                    ADCK1350
                      PACCQ = HALF*ACC/DMIN1(DABS(Q),ONE)               ADCK1351
                      IF(DABS(P) .GT. DABS(Q)) PACCQ = PACCQ*DABS(P)    ADCK1352
C                                                                       ADCK1353
C *** SOLVE FOR FCM = F AT LAMBDA = XLM,THEN FIND NORM FACTOR W=W/FCM   ADCK1354
C                                                                       ADCK1355
      GAM = (F - P)/Q                                                   ADCK1356
            IF(Q .LE. ACC4*DABS(P))             GO TO 130               ADCK1357
      W   = ONE/DSQRT((F - P)*GAM + Q)                                  ADCK1358
            GO TO 10                                                    ADCK1359
C *** ARRIVE HERE IF G(XLM) .GT. 10**6 OR IEXP .GT. 70 &XLTURN = .TRUE. ADCK1360
    9 W   = FJWKB                                                       ADCK1361
      GAM = GJWKB*W                                                     ADCK1362
      P   = F                                                           ADCK1363
      Q   = ONE                                                         ADCK1364
C                                                                       ADCK1365
C *** NORMALISE FOR SPHERICAL OR CYLINDRICAL BESSEL FUNCTIONS           ADCK1366
C                                                                       ADCK1367
   10                     ALPHA = ZERO                                  ADCK1368
          IF(KFN  .EQ. 1) ALPHA = XI                                    ADCK1369
          IF(KFN  .EQ. 2) ALPHA = XI*HALF                               ADCK1370
                          BETA  = ONE                                   ADCK1371
          IF(KFN  .EQ. 1) BETA  = XI                                    ADCK1372
          IF(KFN  .EQ. 2) BETA  = DSQRT(XI)*RT2DPI                      ADCK1373
      FCM  = DSIGN(W,FCL)*BETA                                          ADCK1374
           FC(M1)  = FCM                                                ADCK1375
                      IF(MODE .EQ. 3)           GO TO 11                ADCK1376
           IF(.NOT. XLTURN)   GCL =  FCM*GAM                            ADCK1377
           IF(      XLTURN)   GCL =  GJWKB*BETA                         ADCK1378
           IF( KFN .NE. 0 )   GCL = -GCL                                ADCK1379
           GC(M1)  = GCL                                                ADCK1380
           GPL =  GCL*(P - Q/GAM) - ALPHA*GCL                           ADCK1381
                      IF(MODE .EQ. 2)           GO TO 11                ADCK1382
           GCP(M1) = GPL                                                ADCK1383
           FCP(M1) = FCM*(F - ALPHA)                                    ADCK1384
   11 IF(LXTRA .EQ. 0 ) RETURN                                          ADCK1385
C *** UPWARD RECURRENCE FROM GC(M1),GCP(M1)  STORED VALUE IS RL         ADCK1386
C *** RENORMALISE FC,FCP AT EACH LAMBDA AND CORRECT REGULAR DERIVATIVE  ADCK1387
C ***    XL   = XLM HERE  AND RL = ONE , EL = ZERO FOR BESSELS          ADCK1388
         W    = BETA*W/DABS(FCL)                                        ADCK1389
         MAXL = L1 - 1                                                  ADCK1390
      DO 12 L = M1,MAXL                                                 ADCK1391
                      IF(MODE .EQ. 3)           GO TO 12                ADCK1392
                      XL = XL + ONE                                     ADCK1393
         IF(ETANEO)   EL = ETA/XL                                       ADCK1394
         IF(ETANEO)   RL = GC(L+1)                                      ADCK1395
                      SL = EL + XL*XI                                   ADCK1396
         GCL1     = ((SL - ALPHA)*GCL - GPL)/RL                         ADCK1397
         GPL      =   RL*GCL -  (SL + ALPHA)*GCL1                       ADCK1398
         GCL      = GCL1                                                ADCK1399
         GC(L+1)  = GCL1                                                ADCK1400
                      IF(MODE .EQ. 2)           GO TO 12                ADCK1401
         GCP(L+1) = GPL                                                 ADCK1402
         FCP(L+1) = W*(FCP(L+1) - ALPHA*FC(L+1))                        ADCK1403
   12 FC(L+1)     = W* FC(L+1)                                          ADCK1404
      RETURN                                                            ADCK1405
1000  FORMAT(/' CF1 ACCURACY LOSS: D,DF,ACCH,K,ETA/K,ETA,X = ',1P7D9.2/)ADCK1406
C                                                                       ADCK1407
C ***    ERROR MESSAGES                                                 ADCK1408
C                                                                       ADCK1409
  100 IFAIL = -1                                                        ADCK1410
      WRITE(*,2000) XX,ACCH                                             ADCK1411
 2000 FORMAT(' FOR XX = ',1PD12.3,' TRY SMALL-X SOLUTIONS',             ADCK1412
     *' OR X NEGATIVE'/ ,' SQUARE ROOT ACCURACY PARAMETER =  ',D12.3/)  ADCK1413
      RETURN                                                            ADCK1414
  105 IFAIL = -2                                                        ADCK1415
      WRITE(*,2005) XLMAX,XLMIN,XLM                                     ADCK1416
 2005 FORMAT(/' PROBLEM WITH INPUT ORDER VALUES:XLMAX,XLMIN,XLM = ',    ADCK1417
     *1P3D15.6/)                                                        ADCK1418
      RETURN                                                            ADCK1419
  110 IFAIL =  1                                                        ADCK1420
      WRITE(*,2010) ABORT,F ,DF,PK,PX,ACC                               ADCK1421
 2010 FORMAT(' CF1 HAS FAILED TO CONVERGE AFTER ',F10.0,' ITERATIONS',/ ADCK1422
     *' F,DF,PK,PX,ACCUR =  ',1P5D12.3//)                               ADCK1423
      RETURN                                                            ADCK1424
  120 IFAIL =  2                                                        ADCK1425
      WRITE(*,2020) ABORT,P,Q,DP,DQ,ACC                                 ADCK1426
 2020 FORMAT(' CF2 HAS FAILED TO CONVERGE AFTER ',F7.0,' ITERATIONS',/  ADCK1427
     *' P,Q,DP,DQ,ACCUR =  ',1P4D17.7,D12.3//)                          ADCK1428
      RETURN                                                            ADCK1429
  130 IFAIL =  3                                                        ADCK1430
      WRITE(*,2030) P,Q,ACC,DELL,LXTRA,M1                               ADCK1431
 2030 FORMAT(' FINAL Q.LE.DABS(P)*ACC*10**4 , P,Q,ACC = ',1P3D12.3,4X,  ADCK1432
     *' DELL,LXTRA,M1 = ',D12.3,2I5 /)                                  ADCK1433
      RETURN                                                            ADCK1434
 2040 FORMAT('XLMAX - XLMIN = DELL NOT AN INTEGER ',1P3D20.10/)         ADCK1435
      END                                                               ADCK1436
C                                                                       ADCK1437
      SUBROUTINE JWKB(XX,ETA1,XL,FJWKB,GJWKB,IEXP)                      ADCK1438
      REAL*8          XX,ETA1,XL,FJWKB,GJWKB,DZERO                      ADCK1439
C *** COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS    FOR XL.GE. 0 ADCK1440
C *** AS MODIFIED BY BIEDENHARN ET AL. PHYS REV 97 (1955) 542-554       ADCK1441
C *** CALLS DMAX1,SQRT,ALOG,EXP,ATAN2,FLOAT,INT        BARNETT FEB 1981 ADCK1442
      DATA  ZERO,HALF,ONE,SIX,TEN/ 0.0E0, 0.5E0, 1.0E0, 6.0E0, 10.0E0 / ADCK1443
      DATA  DZERO, RL35, ALOGE   / 0.0D0, 35.0E0, 0.4342945E0 /         ADCK1444
      X     = XX                                                        ADCK1445
      ETA   = ETA1                                                      ADCK1446
      GH2   = X*(ETA + ETA - X)                                         ADCK1447
      XLL1  = DMAX1(XL*XL + XL,DZERO)                                   ADCK1448
      IF(GH2 + XLL1 .LE. ZERO) RETURN                                   ADCK1449
      HLL   = XLL1 + SIX/RL35                                           ADCK1450
      HL    = SQRT(HLL)                                                 ADCK1451
      SL    = ETA/HL + HL/X                                             ADCK1452
      RL2   = ONE + ETA*ETA/HLL                                         ADCK1453
      GH    = SQRT(GH2 + HLL)/X                                         ADCK1454
      PHI   = X*GH - HALF*( HL*ALOG((GH + SL)**2/RL2) - ALOG(GH) )      ADCK1455
          IF(ETA .NE. ZERO) PHI = PHI - ETA*ATAN2(X*GH,X - ETA)         ADCK1456
      PHI10 = -PHI*ALOGE                                                ADCK1457
      IEXP  =  INT(PHI10)                                               ADCK1458
      IF(IEXP .GT. 70) GJWKB = TEN**(PHI10 - FLOAT(IEXP))               ADCK1459
      IF(IEXP .LE. 70) GJWKB = EXP(-PHI)                                ADCK1460
      IF(IEXP .LE. 70) IEXP  = 0                                        ADCK1461
      FJWKB = HALF/(GH*GJWKB)                                           ADCK1462
      RETURN                                                            ADCK1463
      END                                                               ADCK1464
