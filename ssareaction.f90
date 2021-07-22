! ==============================================================================
! Subroutine: SSAREACTION (CAP,ISEED,TOUT)
! 
! Purpose:
!
! Method:
!
! Arguments:
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            09/01/2013   Original Code
!
! Dependancies:
!
! Modules - SYSTEMVAR, CLASS_CAPSID
! Functions - RANDOM
! Subroutines -
!
! Author(s): Eric Dykeman
!
! ==============================================================================

      SUBROUTINE SSAREACTION (CAP,ISEED,TOUT)

        USE SystemVar, ONLY : psum,nsum,nrna,npro,time,&
                            & itot,tnext
        USE Class_Capsid

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        INTEGER, INTENT(INOUT) :: iseed
        TYPE (CAPSID_INT), INTENT(INOUT) :: cap(nrna)

        DOUBLE PRECISION, INTENT(IN) :: tout

        !=== VARIABLES ===!

        INTEGER :: i,j,k,n,n1,n2,nvirus

        DOUBLE PRECISION :: x,xp,r,ar,ac
        DOUBLE PRECISION :: tau,amax,random


        xp = DBLE(npro)

        !=== Capsid Reaction Probabilites ===!

        n = nsum / 2

        amax = psum(1,n) * xp + psum(0,n)


        !=== Compute Time Increment ===!

        tau = 1.0d10

        IF ( amax > 0.0d0 ) THEN

          r = RANDOM(iseed)

          tau = DLOG(1.0d0/r)
          tau = tau / amax

        ENDIF

        time = time + tau


        !=== Output Data ===!

!        IF ( MIN(time,tnext) > tout ) THEN

          !write(*,*)MIN(time,tnext),itot,npro

          !nvirus = 0

          !do i=1,nrna
          !if (cap(i)% nt >= 60 ) nvirus = nvirus + 1
          !enddo

          !open(unit=101,file='num_capsids.dat',access='append')

          !write(101,*)time,nvirus

          !close(unit=101)

!        ENDIF


        !=== Fire Reaction ===!

        r = RANDOM(iseed)
        amax = r * amax

        IF ( itot == 120000 ) tnext = 1.0d10 !!!MAX PROTEIN ALLOWED

        IF ( tnext < time .and. itot < 120000 ) THEN !!!MAX PROTEIN ALLOWED

          time = tnext

          !=== Add Coat Protein ===!

          npro = npro + 1
          itot = itot + 1

          !write(*,*)time,itot,npro
          !write(88,*)time,itot,npro

          !=== Compute Tnext ===!

          !=== 5min ===!

          ar = 4.00d2  !!!CHANGE RATE OF PROTEIN SYNTH

          r = RANDOM(iseed)

          tnext = DLOG(1.0d0/r)
          tnext = time + tnext / ar

        ELSE

          !=== Find Reaction to Fire ===!

          i = n

          DO WHILE ( MOD(n,2) == 0 )

            n = n / 2
            j = i - n

            r = psum(1,j) * xp + psum(0,j)

            IF ( r >= amax ) THEN

              i = j

            ELSE

              i = i + n
              amax = amax - r

            ENDIF

          ENDDO

          r = cap(i)% a(1) * xp + cap(i)% a(0)

          IF ( r >= amax ) THEN

            j = i

          ELSE

            j = i + 1
            amax = amax - r

          ENDIF

          !=== Fire Reaction ===!

          CALL CAPSID_FIRE (cap(j),amax)
          CALL CAPSID_REAC (cap(j))


          !=== Resum Probability Table ===!

          n = 1
          n1= 2
          n2= 4

          psum(0,i) = cap(i)% a(0) + cap(i+1)% a(0)
          psum(1,i) = cap(i)% a(1) + cap(i+1)% a(1)

          DO WHILE ( n1 < nsum )

            i = INT(i/n2) * n2 + n1

            j = i - n
            k = i + n

            psum(0,i) = psum(0,j) + psum(0,k)
            psum(1,i) = psum(1,j) + psum(1,k)

            n  = n1
            n1 = n2
            n2 = 2 * n2

          ENDDO

        ENDIF

        RETURN

      END SUBROUTINE SSAREACTION
