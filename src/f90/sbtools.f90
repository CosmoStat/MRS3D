MODULE SBTools
	IMPLICIT NONE

	CONTAINS

	!===============================================================================
    !===============================================================================

    SUBROUTINE BRootFinder(qln, nnmax, nlmax)

	IMPLICIT NONE

	INTEGER :: nnmax, nlmax, l, nrr, maxrt, lastRoot, err
	PARAMETER(maxrt = 1000)
	REAL*8 :: rmax, order,A,B, A_before
	REAl*8, DIMENSION(0:nlmax,1:nnmax) :: qln
	REAL*8, DIMENSION(:), allocatable :: roots
	integer :: status


	!$OMP parallel &
	!$OMP shared(qln,rmax,nnmax,nlmax,status) &
	!$OMP private(roots,order,A,B,A_before,l,nrr,lastRoot,err)

	allocate(roots(1:maxrt),stat = status)

	!$OMP do schedule(dynamic,1)
	DO l = 0, nlmax

	    order = real(l) + 0.5
	    lastRoot = 1

		PRINT*, "ORDER l=",l

	    A=0.0
	    B=2000.0

		DO WHILE ( lastRoot .LT. nnmax)

			roots = 0.0
			A_before = A
		    CALL ROOTBESSJ( order, A, B, maxrt, nrr, err, roots )

			if(nrr .EQ. 0) then
				B = B + 1000.0
			else if(err .NE. 1) then
				print*,"err code :",err
				A = A_before
				B = (B - A)/2.0 + A
			else

				if(lastRoot + nrr .LE. nnmax) then
					qln(l,lastRoot:lastRoot + nrr - 1) = roots(1:nrr)
					lastRoot = lastRoot + nrr - 1
				else
					qln(l,lastRoot:nnmax) = roots(1: nnmax - lastRoot + 1)
					lastRoot = nnmax
				endif
				B = A + 2048.0
			endif

		ENDDO

	ENDDO
	!$OMP end do
	deallocate(roots)
	!$OMP end parallel

    END SUBROUTINE BRootFinder


    !===============================================================================
    !===============================================================================

    SUBROUTINE BJL(L,X,JL)


	    !!== MODIFIED SUBROUTINE FOR SPHERICAL BESSEL FUNCTIONS.                       ==!!
	    !!== CORRECTED THE SMALL BUGS IN PACKAGE CMBFAST&CAMB(for l=4,5, x~0.001-0.002)==!!
	    !!== CORRECTED THE SIGN OF J_L(X) FOR X<0 CASE                                 ==!!
	    !!== WORKS FASTER AND MORE ACCURATE FOR LOW L, X<<L, AND L<<X cases            ==!!
	    !!== zqhuang@astro.utoronto.ca                                                 ==!!
	    IMPLICIT NONE
	    INTEGER L
	    real*8 X,JL
	    real*8 AX,AX2
	    real,PARAMETER::LN2=0.6931471805599453094D0
	    real,PARAMETER::ONEMLN2=0.30685281944005469058277D0
	    real,PARAMETER::PID2=1.5707963267948966192313217D0
	    real,PARAMETER::PID4=0.78539816339744830961566084582D0
	    real,parameter::ROOTPI12 = 21.269446210866192327578D0
	    real,parameter::GAMMA1 =   2.6789385347077476336556D0 !/* Gamma function of 1/3 */
	    real,parameter::GAMMA2 =   1.3541179394264004169452D0 !/* Gamma function of 2/3 */
	    real,PARAMETER::PI=3.141592653589793238463D0
	    real*8 NU,NU2,BETA,BETA2,COSB
	    real*8 sx,sx2
	    real*8 cotb,cot3b,cot6b,secb,sec2b
	    real*8 trigarg,expterm,L3

	    IF(L.LT.0)THEN
		write(*,*) 'Can not evaluate Spherical Bessel Function with index l<0'
		STOP
	    ENDIF
	    AX=DABS(X)
	    AX2=AX**2
	    IF(L.LT.7)THEN
		IF(L.EQ.0)THEN
		    IF(AX.LT.1.D-1)THEN
			JL=1.D0-AX2/6.D0*(1.D0-AX2/20.D0)
		    ELSE
			JL=DSIN(AX)/AX
		    ENDIF

		ELSEIF(L.EQ.1)THEN
		    IF(AX.LT.2.D-1)THEN
			JL=AX/3.D0*(1.D0-AX2/10.D0*(1.D0-AX2/28.D0))
		    ELSE
			JL=(DSIN(AX)/AX-DCOS(AX))/AX
		    ENDIF
		ELSEIF(L.EQ.2)THEN
		    IF(AX.LT.3.D-1)THEN
			JL=AX2/15.D0*(1.D0-AX2/14.D0*(1.D0-AX2/36.D0))
		    ELSE
			JL=(-3.0D0*DCOS(AX)/AX-DSIN(AX)*(1.D0-3.D0/AX2))/AX
		    ENDIF
		ELSEIF(L.EQ.3)THEN
		    IF(AX.LT.4.D-1)THEN
			JL=AX*AX2/105.D0*(1.D0-AX2/18.D0*(1.D0-AX2/44.D0))
		    ELSE
			JL=(DCOS(AX)*(1.D0-15.D0/AX2)-DSIN(AX)*(6.D0-15.D0/AX2)/AX)/AX
		    ENDIF
		ELSEIF(L.EQ.4)THEN
		    IF(AX.LT.6.D-1)THEN
			JL=AX2**2/945.D0*(1.D0-AX2/22.D0*(1.D0-AX2/52.D0))
		    ELSE
			JL=(DSIN(AX)*(1.D0-(45.D0-105.D0/AX2)/AX2)+DCOS(AX)*(10.D0-105.D0/AX2)/AX)/AX
		    ENDIF
		ELSEIF(L.EQ.5)THEN
		    IF(AX.LT.1.D0)THEN
			JL=AX2**2*AX/10395.D0*(1.D0-AX2/26.D0*(1.D0-AX2/60.D0))
		    ELSE
			JL=(DSIN(AX)*(15.D0-(420.D0-945.D0/AX2)/AX2)/AX-DCOS(AX)*(1.D0-(105.D0-945.0d0/AX2)/AX2))/AX
		    ENDIF
		ELSE
		    IF(AX.LT.1.D0)THEN
			JL=AX2**3/135135.D0*(1.D0-AX2/30.D0*(1.D0-AX2/68.D0))
		    ELSE
			JL=(DSIN(AX)*(-1.D0+(210.D0-(4725.D0-10395.D0/AX2)/AX2)/AX2)+ &
			    DCOS(AX)*(-21.D0+(1260.D0-10395.D0/AX2)/AX2)/AX)/AX
		    ENDIF
		ENDIF
	    ELSE
		NU=0.5D0+L
		NU2=NU**2
		IF(AX.LT.1.D-40)THEN
		    JL=0.D0
		ELSEIF((AX2/L).LT.5.D-1)THEN
		    JL=DEXP(L*DLOG(AX/NU)-LN2+NU*ONEMLN2-(1.D0-(1.D0-3.5D0/NU2)/NU2/30.D0)/12.D0/NU) &
		       /NU*(1.D0-AX2/(4.D0*NU+4.D0)*(1.D0-AX2/(8.D0*NU+16.D0)*(1.D0-AX2/(12.D0*NU+36.D0))))
		ELSEIF((real(L)**2/AX).LT.5.D-1)THEN
		    BETA=AX-PID2*(L+1)
		    JL=(DCOS(BETA)*(1.D0-(NU2-0.25D0)*(NU2-2.25D0)/8.D0/AX2*(1.D0-(NU2-6.25)*(NU2-12.25D0)/48.D0/AX2)) &
		       -DSIN(BETA)*(NU2-0.25D0)/2.D0/AX* (1.D0-(NU2-2.25D0)*(NU2-6.25D0)/24.D0/AX2*(1.D0-(NU2-12.25)* &
			   (NU2-20.25)/80.D0/AX2)) )/AX
		ELSE
		    L3=NU**0.325
		    IF(AX .LT. NU-1.31*L3) then
			COSB=NU/AX
			SX = DSQRT(NU2-AX2)
			COTB=NU/SX
			SECB=AX/NU
			BETA=DLOG(COSB+SX/AX)
			COT3B=COTB**3
			COT6B=COT3B**2
			SEC2B=SECB**2
			EXPTERM=( (2.D0+3.D0*SEC2B)*COT3B/24.D0 &
			   - ( (4.D0+SEC2B)*SEC2B*COT6B/16.D0 &
			   + ((16.D0-(1512.D0+(3654.D0+375.D0*SEC2B)*SEC2B)*SEC2B)*COT3B/5760.D0 &
			   + (32.D0+(288.D0+(232.D0+13.D0*SEC2B)*SEC2B)*SEC2B)*SEC2B*COT6B/128.D0/NU)*COT6B/NU) &
			   /NU)/NU
			JL=DSQRT(COTB*COSB)/(2.D0*NU)*DEXP(-NU*BETA+NU/COTB-EXPTERM)

		    !          /**************** Region 2: x >> l ****************/

		    ELSEIF (AX .GT. NU+1.48*L3) then
			COSB=NU/AX
			SX=DSQRT(AX2-NU2)
			COTB=NU/SX
			SECB=AX/NU
			BETA=DACOS(COSB)
			COT3B=COTB**3
			COT6B=COT3B**2
			SEC2B=SECB**2
			TRIGARG=NU/COTB-NU*BETA-PID4 &
			       -((2.0+3.0*SEC2B)*COT3B/24.D0  &
			       +(16.D0-(1512.D0+(3654.D0+375.D0*SEC2B)*SEC2B)*SEC2B)*COT3B*COT6B/5760.D0/NU2)/NU
			EXPTERM=( (4.D0+sec2b)*sec2b*cot6b/16.D0 &
			       -(32.D0+(288.D0+(232.D0+13.D0*SEC2B)*SEC2B)*SEC2B)*SEC2B*COT6B**2/128.D0/NU2)/NU2
			JL=DSQRT(COTB*COSB)/NU*DEXP(-EXPTERM)*DCOS(TRIGARG)

		    !          /***************** Region 3: x near l ****************/

		    ELSE
			BETA=AX-NU
			BETA2=BETA**2
			SX=6.D0/AX
			SX2=SX**2
			SECB=SX**0.3333333333333333d0
			SEC2B=SECB**2
			JL=( GAMMA1*SECB + BETA*GAMMA2*SEC2B &
			      -(BETA2/18.D0-1.D0/45.D0)*BETA*SX*SECB*GAMMA1 &
			      -((BETA2-1.D0)*BETA2/36.D0+1.D0/420.D0)*SX*SEC2B*GAMMA2   &
			      +(((BETA2/1620.D0-7.D0/3240.D0)*BETA2+1.D0/648.D0)*BETA2-1.D0/8100.D0)*SX2*SECB*GAMMA1 &
			      +(((BETA2/4536.D0-1.D0/810.D0)*BETA2+19.D0/11340.D0)*BETA2-13.D0/28350.D0)*BETA*SX2*SEC2B*GAMMA2 &
			      -((((BETA2/349920.D0-1.D0/29160.D0)*BETA2+71.D0/583200.D0)*BETA2-121.D0/874800.D0)* &
			       BETA2+7939.D0/224532000.D0)*BETA*SX2*SX*SECB*GAMMA1)*DSQRT(SX)/ROOTPI12
		    ENDIF
		ENDIF
	    ENDIF
	    IF(X.LT.0.AND.MOD(L,2).NE.0)JL=-JL


    END SUBROUTINE BJL

    !===============================================================================
    !===============================================================================

END MODULE SBTools
