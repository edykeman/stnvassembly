! ==============================================================================
! Subroutine: READDATA (CAP)
! 
! Purpose: Reads in the list of "bond" interactions between capsid proteins
!          in the capsid, the RNA binding rates, and the hamiltonian path map.
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
! Functions -
! Subroutines -
!
! Author(s): Eric Dykeman
!
! ==============================================================================

      SUBROUTINE READDATA (CAP)

        USE SystemVar
        USE Class_Capsid

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        TYPE (CAPSID_INT), INTENT(INOUT) :: cap(nrna)

        !=== VARIABLES ===!

        INTEGER :: i,j,k,jj,kk
        CHARACTER (LEN=70) :: fmat

        DOUBLE PRECISION :: x,rcf,rcb
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ps,rf
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: psr

        DOUBLE PRECISION, PARAMETER :: cfac = 0.16605389210321896964d-8


        !=== Set Beta ===!

        beta = 1.0d0 / ( temp * gcons )


        !=== Set rkap1 and rkap2 ===!
        !=== convert from 1/M*s to 1/s ===!

        rkap1 = ratep1

        rkap2 = ratep2 * cfac / vol


        !=== Read in Parameter File ===!

        READ(1,*)ncp,nbond

        ALLOCATE (lbond(2,nbond),gbond(nbond))

        READ(1,*)(lbond(1,i),lbond(2,i),gbond(i),i=1,nbond)


        !=== Read In Hamiltonian Map ===!

        READ(2,*)nps,nnps

        ALLOCATE (maphp(nnps,nps))

        fmat = '(15I5)'

        DO i=1,nps

          READ(2,fmat)(maphp(j,i),j=1,nnps)
!	  PRINT*, (maphp(j,i),j=1,nnps)

        ENDDO

        !NOTE : for stnv - nps = 30 --> reset !
        nps = 30

        ALLOCATE (ps(nps),rf(nps),psr(2,nps))

        !=== Read in RNA PSs ===!

        fmat = '(30E11.4)'

        !ps(:) = 12.0d0!8.0d0
        !ps(1:5) = 12.00d0
        !rf(:) = 1.1d7
        DO i=1,nrna

          if ( i == 1 ) then
          READ(3,fmat)(rf(j),j=1,nps)
          READ(3,fmat)(ps(j),j=1,nps)
          !elseif ( i == 3001 ) then
          !ps(:) = 3.0d0
          endif

          DO j=1,nps

            rcf = rf(j)

            !=== convert from 1/M*s to 1/s ===!

            rcf = rcf * cfac / vol

            !=== Calculate rateB ===!

            x = ps(j) * beta
            x = DEXP(x)

            rcb = rf(j) / x

            psr(1,j) = rcf
            psr(2,j) = rcb

          ENDDO

          cap(i)% ps(:) = ps(:)
          cap(i)% psr(:,:) = psr(:,:)

        ENDDO


        !=== Form Neighbormap ===!

        ALLOCATE (nncp(ncp))

        nnmax = 0
        nncp(:) = 0

        DO i=1,nbond

          j = lbond(1,i)
          k = lbond(2,i)

          nncp(j) = nncp(j) + 1
          nncp(k) = nncp(k) + 1

        ENDDO

        DO i=1,ncp
        IF ( nncp(i) > nnmax ) THEN
          nnmax = nncp(i)
        ENDIF
        ENDDO

        ALLOCATE (mapnn(nnmax,ncp))
        ALLOCATE (gb(nnmax,ncp))

        nncp(:) = 0

        DO i=1,nbond

          j = lbond(1,i)
          k = lbond(2,i)

          nncp(j) = nncp(j) + 1
          nncp(k) = nncp(k) + 1

          jj = nncp(j)
          kk = nncp(k)

          mapnn(jj,j) = k
          mapnn(kk,k) = j

          gb(jj,j) = gbond(i)
          gb(kk,k) = gbond(i)

        ENDDO

        !=== Read in Autosteric Exclusion Map ===!

        ALLOCATE (mapauto(ncp))

        OPEN(unit=123,FILE='mut_exc.dat',Status='Old')

        DO i=1,ncp
        READ(123,*)j
        mapauto(i) = j
        ENDDO

        !=== Read in Two Fold Map ===!

        ALLOCATE (maptwo(2,ncp))

        OPEN(unit=124,FILE='two_fold.in',Status='Old')

        DO i=1,ncp
        READ(124,*)j,k
        maptwo(1,i) = j
        maptwo(2,i) = k
        ENDDO

        RETURN

      END SUBROUTINE READDATA
