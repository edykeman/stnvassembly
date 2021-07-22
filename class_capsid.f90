! ==============================================================================
! Module: CLASS_CAPSID
! 
! Purpose: (MS2 VERSION)
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            09/01/2013   Original Code
!
! Contains:
!
! Modules -
! Functions -
! Subroutines -
!
! Author(s): Eric Dykeman
!
! ==============================================================================

      MODULE CLASS_CAPSID

        USE SystemVar, ONLY : mapnn,maphp,nncp,ncp,npro,ratep1,ratep2,&
                            & rkap1,rkap2,nps,nnps,nbond,gb,beta,mapauto

        IMPLICIT NONE

        PRIVATE

        PUBLIC :: CAPSID_REAC, CAPSID_FIRE

        TYPE, PUBLIC :: CAPSID_INT

          INTEGER :: ipro(60)
          INTEGER :: irna(30)
          INTEGER :: ireac(60)
          INTEGER :: iauto(60)

          INTEGER :: i5
          INTEGER :: i3
          INTEGER :: nt

          DOUBLE PRECISION :: a(0:1)
          DOUBLE PRECISION :: ps(30)
          DOUBLE PRECISION :: psr(2,30)

        END TYPE CAPSID_INT

        CONTAINS

        SUBROUTINE CAPSID_REAC (C)

          IMPLICIT NONE

          !=== ARGUMENTS ===!

          TYPE(CAPSID_INT), INTENT(INOUT) :: c

          !=== VARIABLES ===!

          LOGICAL :: ok5,ok3
          INTEGER :: i,j,k,m,i5,i3,nt
          INTEGER :: i0,i1,i2,is,nc,inum
          INTEGER :: num(ncp),low(ncp)
          INTEGER :: ls(2*nbond),lt(2*nbond)

          DOUBLE PRECISION :: x,dg,a(0:1)


          a(:) = 0.0d0

          i5 = c% i5
          i3 = c% i3
          nt = c% nt

          !=== Make Sure iauto is zeroed if nt == 0 ===!

          IF ( nt == 0 ) c% iauto(:) = 0

          is = 0

          IF ( nt /= 0 ) THEN

            i = 0
            is = 1
            ls(is) = 0
            lt(is) = 0

            DO WHILE ( ls(is) == 0 )

              i = i + 1
              j = c% ipro(i)

              IF ( j == 1 ) ls(is) = i

            ENDDO

          ENDIF

          ok5 = .false.
          IF ( i5 > 1 ) THEN
          IF ( c% irna(i5-1) == -1 ) ok5 = .true.
          ENDIF

          ok3 = .false.
          IF ( i3 < nps ) THEN
          IF ( c% irna(i3+1) == -1 ) ok3 = .true.
          ENDIF

          nc = 0
          inum = 0
          num(:) = 0
          low(:) = 0

          c% ireac(:) = c% ipro(:)

          !=== Depth First Search ===!

          DO WHILE ( is /= 0 )

            i1 = ls(is)
            i0 = lt(is)

            IF ( num(i1) == 0 ) THEN

              inum = inum + 1
              num(i1) = inum
              low(i1) = inum

              IF ( i0 /= 0 ) THEN
              IF ( num(i0) == 1 ) nc = nc + 1
              ENDIF

              !=== Add to stack ===!

              DO i=1,nncp(i1)

                i2 = mapnn(i,i1)
                j = c% ipro(i2)

                IF ( j == 0 ) c% ireac(i2) = 2

                IF ( j == 1 ) THEN
                IF ( num(i2) == 0 ) THEN

                  is = is + 1

                  ls(is) = i2
                  lt(is) = i1

                ELSEIF ( i2 /= i0 ) THEN

                  IF ( num(i2) < num(i1) ) THEN
                  low(i1) = MIN(low(i1),num(i2))
                  ENDIF

                ENDIF
                ENDIF

              ENDDO

            ELSE

              !=== Back Track LOW ===!

              IF ( i0 /= 0 ) THEN

                low(i0) = MIN(low(i0),low(i1))

                IF ( num(i0) == 1 ) THEN

                  IF ( nc > 1 ) c% ireac(i0) = 0

                ELSE

                  IF ( low(i1) >= num(i0) ) c% ireac(i0) = 0

                ENDIF

              ENDIF

              !=== Clear Stack ===!

              is = is - 1

            ENDIF

          ENDDO

          !=== Nucleation Check ===!

          IF ( i3-i5 == 1 ) THEN

            j = c% irna(i5)
            k = c% irna(i3)
            m = maphp(2,k)

            IF ( nt > 2 ) THEN
              c% ireac(j) = 0
              c% ireac(k) = 0
            ENDIF

          ELSEIF ( i3-i5 > 1 ) THEN

            DO i=i5+1,i3-1
            j = c% irna(i)
            c% ireac(j) = 0
            ENDDO

          ENDIF

          !=== Check for Autosteric CP Addition ===!

          DO i=1,ncp
          IF ( c% iauto(i) == 1 ) THEN
          IF ( c% ireac(i) == 2 ) c% ireac(i) = 3
          ENDIF
          ENDDO


          !=== COMPUTE REACTIONS ===!

          !=== RNA / CP Binding Reactions ===!

          DO i=1,nps

            !RNA i IS CP UNBOUND      irna(i) =  0
            !RNA i IS CP BOUND        irna(i) = -1
            !RNA i IN CAPSID COMPLEX  irna(i) >  0

            j = c% irna(i)

            IF ( j ==  0 ) a(1) = a(1) + c% psr(1,i)
            IF ( j == -1 ) a(0) = a(0) + c% psr(2,i)

          ENDDO


          !=== CAPSID Reactions ===!

          IF ( i5 == 0 .and. i3 == 0 ) THEN

            !=== Nucleation Reaction ===!

            DO i=1,nps-1

              IF ( c% irna(i)   /= -1 ) CYCLE
              IF ( c% irna(i+1) /= -1 ) CYCLE

              DO j=1,nnps

                !=== Only bring in neighbors which make a cp-cp bond ===!

                IF ( j == 1 .or. j == 3 ) a(0) = a(0) + rkap1

              ENDDO

            ENDDO

          ELSE

            !=== Elongation Reactions ===!

            DO i=1,ncp
            IF ( c% ireac(i) == 2 ) THEN

              DO j=1,nnps

                k = maphp(j,i)

                IF ( c% irna(i5) == k ) THEN
                IF ( ok5 ) a(0) = a(0) + rkap1
                ENDIF

                IF ( c% irna(i3) == k ) THEN
                IF ( ok3 ) a(0) = a(0) + rkap1
                ENDIF

              ENDDO

            ELSEIF ( c% ireac(i) == 3 ) THEN

              a(1) = a(1) + rkap2

            ELSEIF ( c% ireac(i) == 1 ) THEN

              dg = 0.0d0

              DO j=1,nncp(i)

                k = mapnn(j,i)
                m = c% ipro(k)

                IF ( m /= 0 ) dg = dg + gb(j,i)

              ENDDO

              x = dg * beta
              x = DEXP(x)

              IF ( c% iauto(i) /= 1 ) THEN
                x = ratep1 / x
              ELSE
                x = ratep2 / x
              ENDIF

              a(0) = a(0) + x

            ENDIF
            ENDDO

          ENDIF

          c% a(:) = a(:)

          RETURN

        END SUBROUTINE CAPSID_REAC

        SUBROUTINE CAPSID_FIRE (C,AMAX)

          IMPLICIT NONE

          !=== ARGUMENTS ===!

          TYPE(CAPSID_INT), INTENT(INOUT) :: c
          DOUBLE PRECISION, INTENT(IN) :: amax

          !=== VARIABLES ===!

          LOGICAL :: ok5,ok3
          INTEGER :: i,j,k,m,i5,i3,nt
          DOUBLE PRECISION :: x,xp,dg,atot


          atot = 0.0d0
          xp = DBLE(npro)

          i5 = c% i5
          i3 = c% i3
          nt = c% nt

          ok5 = .false.
          IF ( i5 > 1 ) THEN
          IF ( c% irna(i5-1) == -1 ) ok5 = .true.
          ENDIF

          ok3 = .false.
          IF ( i3 < nps ) THEN
          IF ( c% irna(i3+1) == -1 ) ok3 = .true.
          ENDIF


          !=== RNA / CP Binding Reactions ===!

          DO i=1,nps

            !RNA i IS CP UNBOUND      irna(i) =  0
            !RNA i IS CP BOUND        irna(i) = -1
            !RNA i IN CAPSID COMPLEX  irna(i) >  0

            j = c% irna(i)

            IF ( j ==  0 ) atot = atot + c% psr(1,i) * xp
            IF ( j == -1 ) atot = atot + c% psr(2,i)

            IF ( atot >= amax ) THEN

              IF ( j == 0 ) THEN
                npro = npro - 1
                c% irna(i) = -1
              ELSEIF ( j == -1 ) THEN
                npro = npro + 1
                c% irna(i) = 0
              ENDIF

              RETURN

            ENDIF

          ENDDO


          !=== CAPSID Reactions ===!

          IF ( i5 == 0 .and. i3 == 0 ) THEN

            !=== Nucleation Reaction ===!

            DO i=1,nps-1

              IF ( c% irna(i)   /= -1 ) CYCLE
              IF ( c% irna(i+1) /= -1 ) CYCLE

              DO j=1,nnps

                k = maphp(j,1)

                !=== Nucleate ===!

                IF ( j == 1 .or. j == 3 ) atot = atot + rkap1

                IF ( atot >= amax ) THEN

                  c% ipro(1) = 1 
                  c% irna(i) = 1

                  c% ipro(k) = 1
                  c% irna(i+1) = k

                  c% i5 = i
                  c% i3 = i + 1
                  c% nt = 2

                  !=== Autosteric marker ===!

                  m = mapauto(1)
                  c% iauto(m) = 1
                  m = mapauto(k)
                  c% iauto(m) = 1

                  RETURN

                ENDIF

              ENDDO

            ENDDO

          ELSE

            !=== Elongation Reactions ===!

            DO i=1,ncp
            IF ( c% ireac(i) == 2 ) THEN

              DO j=1,nnps

                k = maphp(j,i)

                m = c% irna(i5)

                IF ( m == k .and. ok5 ) THEN

                  atot = atot + rkap1

                  IF ( atot >= amax ) THEN

                    c% ipro(i) = 1
                    c% irna(i5-1) = i
                    c% i5 = c% i5 - 1
                    c% nt = c% nt + 1

                    !=== Autosteric marker ===!

                    m = mapauto(i)
                    c% iauto(m) = 1

                    RETURN

                  ENDIF

                ENDIF

                m = c% irna(i3)

                IF ( m == k .and. ok3 ) THEN

                  atot = atot + rkap1

                  IF ( atot >= amax ) THEN

                    c% ipro(i) = 1
                    c% irna(i3+1) = i
                    c% i3 = c% i3 + 1
                    c% nt = c% nt + 1

                    !=== Autosteric marker ===!

                    m = mapauto(i)
                    c% iauto(m) = 1

                    RETURN

                  ENDIF

                ENDIF

              ENDDO

            ELSEIF ( c% ireac(i) == 3 ) THEN

              atot = atot + rkap2 * xp

              IF ( atot >= amax ) THEN

                c% ipro(i) = 1
                c% nt = c% nt + 1
                npro = npro - 1

                RETURN

              ENDIF

            ELSEIF ( c% ireac(i) == 1 ) THEN

              dg = 0.0d0

              DO j=1,nncp(i)

                k = mapnn(j,i)
                m = c% ipro(k)

                IF ( m /= 0 ) dg = dg + gb(j,i)

              ENDDO

              x = dg * beta
              x = DEXP(x)

              IF ( c% iauto(i) /= 1 ) THEN
                x = ratep1 / x
              ELSE
                x = ratep2 / x
              ENDIF

              atot = atot + x

              IF ( atot >= amax ) THEN

                IF ( c% nt == 2 ) THEN

                  c% ipro(:) = 0
                  c% iauto(:) = 0
                  c% irna(i5) = -1
                  c% irna(i3) = -1
                  c% i5 = 0
                  c% i3 = 0
                  c% nt = 0

                ELSEIF ( c% irna(i5) == i ) THEN

                  c% ipro(i) = 0
                  c% irna(i5) = -1
                  c% i5 = c% i5 + 1
                  c% nt = c% nt - 1

                  !=== Remove Autosteric ===!

                  m = mapauto(i)
                  c% iauto(m) = 0

                ELSEIF ( c% irna(i3) == i ) THEN

                  c% ipro(i) = 0
                  c% irna(i3) = -1
                  c% i3 = c% i3 - 1
                  c% nt = c% nt - 1

                  !=== Remove Autosteric ===!

                  m = mapauto(i)
                  c% iauto(m) = 0

                ELSE

                  c% ipro(i) = 0
                  c% nt = c% nt - 1
                  npro = npro + 1

                ENDIF

                RETURN

              ENDIF

            ENDIF
            ENDDO

          ENDIF

          RETURN

        END SUBROUTINE CAPSID_FIRE

      END MODULE CLASS_CAPSID
