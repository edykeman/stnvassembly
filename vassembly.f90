! ==============================================================================
! Program: VASSEMBLY
!
! Discription: Performs a "black box" integration of a system of
!              chemical reactions which describe virus assembly.
!
! NOTES:
!        FILE TREE
!
!        Unit = 1  --- Parameter File
!        Unit = 2  --- Hamiltonian Path File
!        Unit = 3  --- RNA Data File
!        Unit = 4  --- Option File
!        Unit = 5  --- Standard Out (NOT USED)
!        Unit = 6  --- Standard In (NOT USED)
!        Unit = 7  --- Output File
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
! Subroutines - READDATA
!
! Author(s): Eric Dykeman
!
! ==============================================================================

PROGRAM VASSEMBLY

	USE SystemVar
	USE Class_Capsid

	IMPLICIT NONE

	!=== Variable Declaration ===!

	TYPE(CAPSID_INT), DIMENSION(:), ALLOCATABLE :: cap

	INTEGER :: i,j,k,n,io,IHIST(60)
	INTEGER :: iseed,iargc,narg

	DOUBLE PRECISION :: x,r,dt,tout,random

	CHARACTER(LEN=256) :: optfile,parmfile,rnafile
	CHARACTER(LEN=256) :: hpfile,outfile,arg,fname_in,fname_out


	!=== WELCOME ===!

	!WRITE(*,*)'               WELCOME TO VASSEMBLY                 '
	!WRITE(*,*)' '
	!WRITE(*,*)'This program integrates a set of chemical reactions '
	!WRITE(*,*)'between capsid intermediates as a function of time.'
	!WRITE(*,*)' '

	!=== Default names ===!

	optfile = 'options.in'

	parmfile = 'capsid.parm'
	rnafile  = 'capsid.rna'
	hpfile   = 'capsid.hmap'
	outfile  = 'capsid.out'


        narg = IARGC ()

        DO i=1,narg,2

          CALL GETARG (i,arg)

          SELECT CASE (arg)

            CASE ('-p')
              CALL GETARG (i+1,optfile)
            CASE ('-o')
              CALL GETARG (i+1,fname_out)
            CASE ('-i')
              CALL GETARG (i+1,fname_in)
            CASE DEFAULT
              WRITE(*,*)arg,'Invalid Line Argument'
              STOP
          END SELECT

        ENDDO


	!=== SECTION 0 ===!

	!=== Read in Options ===!

	OPEN (UNIT = 4,FILE = optfile,STATUS='Unknown')

	READ(4,*)parmfile
	READ(4,*)hpfile
	READ(4,*)rnafile
	READ(4,*)outfile

	READ(4,*)nrna
	READ(4,*)npro
	READ(4,*)vol
	READ(4,*)temp
	READ(4,*)ratep1
	READ(4,*)ratep2
	READ(4,*)time
	READ(4,*)tfinal
	READ(4,*)iseed

	CLOSE (UNIT = 4)

	ALLOCATE (cap(nrna))

        parmfile = 'stnv.parm'
        hpfile = 'stnv.hmap'
        rnafile = fname_in
        outfile = fname_out

	!=== SECTION 1 ===!

	!=== Read In Parameters ===!

	OPEN (UNIT = 1,FILE = parmfile,STATUS='Unknown')
	OPEN (UNIT = 2,FILE =   hpfile,STATUS='Unknown')
	OPEN (UNIT = 3,FILE =  rnafile,STATUS='Unknown')

        !WRITE(*,*)'Reading in system data.'

	CALL READDATA (cap)

	CLOSE (UNIT = 1)
	CLOSE (UNIT = 2)
	CLOSE (UNIT = 3)


	!=== SECTION 2 ===!

	io = 1			!*
	dt = 1.0d-7
	tout = dt

	IF ( time >= dt ) THEN

		tout = DLOG(time) / DLOG(10.0d0)

		i = FLOOR(tout)

		dt = 10.0d0 ** i

		io = NINT(time/dt)

		io = io + 1
		tout = DBLE(io) * dt

	ENDIF

	!ECD-resetRAMP
	npro = 0 !180000
	time = 0.0d0

	io = 1			!*
	dt = 1.0d-7
	tout = dt

	x = 4.00d2!!!CHANGE RATE OF PROTEIN SYNTH

	r = RANDOM(iseed)

	tnext = DLOG(1.0d0/r)
	tnext = time + tnext / x
	!ECD-end

	DO i=1,nrna

		cap(i)% ipro(:) = 0
		cap(i)% irna(:) = 0

		cap(i)% i5 = 0
		cap(i)% i3 = 0
		cap(i)% nt = 0
		cap(i)% a(:) = 0.0d0

		CALL CAPSID_REAC (cap(i))

	ENDDO

	!ECD-START

	nsum = 2
	DO WHILE ( nsum < nrna )
		nsum = 2 * nsum
	ENDDO

	IF ( .not. ALLOCATED(psum) ) ALLOCATE (psum(0:1,nsum))

	psum(:,:) = 0.0d0

	DO i=1,nrna,2
		psum(0,i) = cap(i)% a(0) + cap(i+1)% a(0)
		psum(1,i) = cap(i)% a(1) + cap(i+1)% a(1)
	ENDDO

	k = 4
	n = 1
	DO WHILE ( k <= nsum )
		j = k / 2
		DO i=j,nsum,k
			psum(0,i) = psum(0,i-n) + psum(0,i+n)
			psum(1,i) = psum(1,i-n) + psum(1,i+n)
		ENDDO

		n = n * 2
		k = k * 2

	ENDDO

	!ECD-END

	!=== BEGIN STOCHASTIC SIMULATION ===!

	!WRITE(*,*)'Starting SSA ...'
	!write(*,*)time,tfinal,itot,npro
	! itot is the number of coat proteis synthesized
	! npro is the number of free proteins
	DO WHILE ( time < tfinal )

		CALL SSAREACTION (cap,iseed,tout)

		!=== Increment tout ===!

		IF ( time > tout ) THEN			!*
	
			tout = tout + dt

			io = io + 1

			IF ( io > 9 ) THEN
				io = 1
				dt = dt * 10.0d0
			ENDIF

		ENDIF					!*
	
	ENDDO

	!WRITE(*,*)'Stochastic Simulation complete.'

       !=== SECTION 4 ===!

        k = 0

        DO j=1,nrna
        IF ( cap(j)% nt == 60 ) k = k + 1
        ENDDO

        OPEN (UNIT=123,FILE=fname_out,STATUS='Unknown')
        WRITE(123,*)k
        CLOSE(UNIT=123)


END PROGRAM VASSEMBLY
