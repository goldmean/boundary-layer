MODULE FUNCTIONS

    USE CONSTANTS

    IMPLICIT NONE

    INTEGER, PARAMETER, PRIVATE :: VELOCITY_UNIT = 1, FRICTION_UNIT = 2, &
        PROFILES_UNIT = 3, RESIDUAL_UNIT = 4, MAXITER = 1E6

    DOUBLE PRECISION, PARAMETER, PRIVATE :: EPS = 1E-6, NU = MU / R0, &
        XSTEP = L / (NI - 1), YSTEP = H / (NJ - 1), A = 1.0 / U0 ** 2, &
        TSTEP = CFL / U0 * SQRT(XSTEP ** 2 + YSTEP ** 2)

    PRIVATE :: INITIAL_CONDITIONS, BOUNDARY_CONDITIONS, AUTOMODEL_VELOCITY_PROFILES
    PUBLIC :: CREATE_FILES, CREATE_MESH, NAVIER_STOKES_SOLVER, NUMERICAL_FRICTION, &
        BLASIUS_FRICTION

CONTAINS

    SUBROUTINE CREATE_MESH(XNODE, YNODE, XCELL, YCELL)
        IMPLICIT NONE

        INTEGER :: I, J
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(:, :) :: XNODE, YNODE
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(0:, 0:) :: XCELL, YCELL

        DO I = 1, NI
            DO J = 1, NJ
                XNODE(I, J) = (I - 1) * XSTEP
                YNODE(I, J) = (J - 1) * YSTEP
            END DO
        END DO

        XCELL(0, 1: NJ) = - XSTEP / 2.0
        YCELL(0, 1: NJ) = YNODE(1, 1: NJ) + YSTEP / 2.0

        YCELL(1: NI, 0) = - YSTEP / 2.0
        XCELL(1: NI, 0) = XNODE(1: NI, 1) + XSTEP / 2.0

        DO I = 1, NI
            DO J = 1, NJ
                XCELL(I, J) = XNODE(I, J) + XSTEP / 2.0
                YCELL(I, J) = YNODE(I, J) + YSTEP / 2.0
            END DO
        END DO

    END SUBROUTINE CREATE_MESH

    SUBROUTINE INITIAL_CONDITIONS(UCELL, VCELL, PCELL)
        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(OUT), DIMENSION(0:, 0:) :: UCELL, VCELL, PCELL

        UCELL(:, :) = U0
        VCELL(:, :) = V0
        PCELL(:, :) = P0

    END SUBROUTINE INITIAL_CONDITIONS

    SUBROUTINE BOUNDARY_CONDITIONS(UCELL, VCELL, PCELL)
        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(OUT), DIMENSION(0:, 0:) :: UCELL, VCELL, PCELL

        UCELL(:, 0) = - UCELL(:, 1)
        VCELL(:, 0) = - VCELL(:, 1)
        PCELL(:, 0) = PCELL(:, 1)

        UCELL(NI, :) = UCELL(NI - 1, :)
        VCELL(NI, :) = VCELL(NI - 1, :)
        PCELL(NI, :) = P0

        UCELL(0, :) = U0
        VCELL(0, :) = V0
        PCELL(0, :) = PCELL(1, :)

        VCELL(:, NJ) = VCELL(:, NJ - 1)

        WHERE (VCELL(:, NJ - 1) >= V0)
            UCELL(:, NJ) = UCELL(:, NJ - 1)
            PCELL(:, NJ) = P0
        ELSE WHERE
            UCELL(:, NJ) = U0
            PCELL(:, NJ) = PCELL(:, NJ - 1)
        END WHERE

    END SUBROUTINE BOUNDARY_CONDITIONS

    SUBROUTINE NAVIER_STOKES_SOLVER(UCELL, VCELL, PCELL)
        IMPLICIT NONE

        INTEGER :: I, J, ITERATION
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(0:, 0:) :: UCELL, VCELL, PCELL
        DOUBLE PRECISION :: UHATRIGHT, UHATLEFT, VHATUPPER, VHATBOUND, &
            URIGHT, VRIGHT, PRIGHT, ULEFT, VLEFT, PLEFT, UUPPER, VUPPER, PUPPER, &
            UBOUND, VBOUND, PBOUND

        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: RESU, RESV, RESP
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: MAXRESU, MAXRESV, MAXRESP
        INTEGER, ALLOCATABLE, DIMENSION(:) :: INDICES

        ALLOCATE(RESU(NI, NJ), RESV(NI, NJ), RESP(NI, NJ), MAXRESU(MAXITER), &
            MAXRESV(MAXITER), MAXRESP(MAXITER), INDICES(MAXITER))

        CALL INITIAL_CONDITIONS(UCELL, VCELL, PCELL)

        ITERATION = 0

        DO WHILE (ITERATION <= MAXITER)
            CALL BOUNDARY_CONDITIONS(UCELL, VCELL, PCELL)

            DO J = 1, NJ - 1
                DO I = 1, NI - 1
                    UHATRIGHT = 1.0 / 2.0 * (UCELL(I + 1, J) + UCELL(I, J))
                    UHATLEFT  = 1.0 / 2.0 * (UCELL(I - 1, J) + UCELL(I, J))
                    VHATUPPER = 1.0 / 2.0 * (VCELL(I, J + 1) + VCELL(I, J))
                    VHATBOUND = 1.0 / 2.0 * (VCELL(I, J - 1) + VCELL(I, J))

                    IF (UHATRIGHT >= 0.0) THEN
                        URIGHT = UCELL(I, J)
                        VRIGHT = VCELL(I, J)
                        PRIGHT = PCELL(I + 1, J)
                    ELSE
                        URIGHT = UCELL(I + 1, J)
                        VRIGHT = VCELL(I + 1, J)
                        PRIGHT = PCELL(I, J)
                    END IF

                    IF (UHATLEFT >= 0.0) THEN
                        ULEFT = UCELL(I - 1, J)
                        VLEFT = VCELL(I - 1, J)
                        PLEFT = PCELL(I, J)
                    ELSE
                        ULEFT = UCELL(I, J)
                        VLEFT = VCELL(I, J)
                        PLEFT = PCELL(I - 1, J)
                    END IF

                    IF (VHATUPPER >= 0.0) THEN
                        UUPPER = UCELL(I, J)
                        VUPPER = VCELL(I, J)
                        PUPPER = PCELL(I, J + 1)
                    ELSE
                        UUPPER = UCELL(I, J + 1)
                        VUPPER = VCELL(I, J + 1)
                        PUPPER = PCELL(I, J)
                    END IF

                    IF (VHATBOUND >= 0.0) THEN
                        UBOUND = UCELL(I, J - 1)
                        VBOUND = VCELL(I, J - 1)
                        PBOUND = PCELL(I, J)
                    ELSE
                        UBOUND = UCELL(I, J)
                        VBOUND = VCELL(I, J)
                        PBOUND = PCELL(I, J - 1)
                    END IF

                    RESU(I, J) = - ((UHATRIGHT * URIGHT - UHATLEFT * ULEFT) / XSTEP + &
                        (VHATUPPER * UUPPER - VHATBOUND * UBOUND) / YSTEP + &
                        (PRIGHT - PLEFT) / XSTEP - &
                        NU * (UCELL(I + 1, J) - 2.0 * UCELL(I, J) + UCELL(I - 1, J)) / &
                        XSTEP ** 2 - &
                        NU * (UCELL(I, J + 1) - 2.0 * UCELL(I, J) + UCELL(I, J - 1)) / &
                        YSTEP ** 2)

                    RESV(I, J) = - ((VHATUPPER * VUPPER - VHATBOUND * VBOUND) / YSTEP + &
                        (UHATRIGHT * VRIGHT - UHATLEFT * VLEFT) / XSTEP + &
                        (PUPPER - PBOUND) / YSTEP - &
                        NU * (VCELL(I + 1, J) - 2.0 * VCELL(I, J) + VCELL(I - 1, J)) / &
                        XSTEP ** 2 - &
                        NU * (VCELL(I, J + 1) - 2.0 * VCELL(I, J) + VCELL(I, J - 1)) / &
                        YSTEP ** 2)

                    IF (J == 1) VBOUND = 0.0

                    RESP(I, J) = - 1.0 / A * ((URIGHT - ULEFT) / XSTEP + &
                        (VUPPER - VBOUND) / YSTEP)

                    UCELL(I, J) = UCELL(I, J) + TSTEP * RESU(I, J)
                    VCELL(I, J) = VCELL(I, J) + TSTEP * RESV(I, J)
                    PCELL(I, J) = PCELL(I, J) + TSTEP * RESP(I, J)
                END DO
            END DO

            ITERATION = ITERATION + 1

            MAXRESU(ITERATION) = MAXVAL(ABS(RESU(1: NI - 1, 1: NJ - 1)))
            MAXRESV(ITERATION) = MAXVAL(ABS(RESV(1: NI - 1, 1: NJ - 1)))
            MAXRESP(ITERATION) = MAXVAL(ABS(RESP(1: NI - 1, 1: NJ - 1)))

            IF (MAXRESU(ITERATION) < EPS .AND. MAXRESV(ITERATION) < EPS .AND. &
                MAXRESP(ITERATION) < EPS) EXIT
        END DO

        FORALL (I = 1: ITERATION) INDICES(I) = I

        OPEN(UNIT=RESIDUAL_UNIT, FILE=RESIDUAL, ACTION="WRITE")

        WRITE(RESIDUAL_UNIT, *) 'VARIABLES = "ITERATION", "RESIDUAL"'
        WRITE(RESIDUAL_UNIT, *) "ZONE I = ", ITERATION, ", DATAPACKING = BLOCK"
        WRITE(RESIDUAL_UNIT, "(100I8)") INDICES(1: ITERATION)
        WRITE(RESIDUAL_UNIT, "(100E25.16)") MAXRESU(1: ITERATION)

        WRITE(RESIDUAL_UNIT, *) 'VARIABLES = "ITERATION", "RESIDUAL"'
        WRITE(RESIDUAL_UNIT, *) "ZONE I = ", ITERATION, ", DATAPACKING = BLOCK"
        WRITE(RESIDUAL_UNIT, "(100I8)") INDICES(1: ITERATION)
        WRITE(RESIDUAL_UNIT, "(100E25.16)") MAXRESV(1: ITERATION)

        WRITE(RESIDUAL_UNIT, *) 'VARIABLES = "ITERATION", "RESIDUAL"'
        WRITE(RESIDUAL_UNIT, *) "ZONE I = ", ITERATION, ", DATAPACKING = BLOCK"
        WRITE(RESIDUAL_UNIT, "(100I8)") INDICES(1: ITERATION)
        WRITE(RESIDUAL_UNIT, "(100E25.16)") MAXRESP(1: ITERATION)

        CLOSE(UNIT=RESIDUAL_UNIT)

        DEALLOCATE(RESU, RESV, RESP, MAXRESU, MAXRESV, MAXRESP, INDICES)

    END SUBROUTINE NAVIER_STOKES_SOLVER

    SUBROUTINE NUMERICAL_FRICTION(UCELL, CF)
        IMPLICIT NONE

        INTEGER :: I
        DOUBLE PRECISION, INTENT(IN), DIMENSION(0:, 0:) :: UCELL
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(:) :: CF

        DO I = 1, NI
            CF(I) = 2.0 * MU * (UCELL(I, 1) - UCELL(I, 0)) / &
                YSTEP / (R0 * U0 ** 2)  ! SECOND ORDER
        END DO

    END SUBROUTINE NUMERICAL_FRICTION

    SUBROUTINE BLASIUS_FRICTION(XNODE, RE, BL)
        IMPLICIT NONE

        INTEGER :: I
        DOUBLE PRECISION, INTENT(IN), DIMENSION(:, :) :: XNODE
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(:) :: RE, BL

        DO I = 1, NI
            RE(I) = U0 * XNODE(I, NJ) / NU
            BL(I) = 0.664 / SQRT(RE(I))
        END DO

    END SUBROUTINE BLASIUS_FRICTION

    SUBROUTINE AUTOMODEL_VELOCITY_PROFILES(YCELL, UCELL)
        IMPLICIT NONE

        INTEGER :: I, J, K
        DOUBLE PRECISION :: DELTA
        DOUBLE PRECISION, INTENT(IN), DIMENSION(0:, 0:) :: YCELL, UCELL
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: YAUTO, UAUTO

        ALLOCATE(YAUTO(NJ), UAUTO(NJ))

        DO K = 1, SIZE(SECTIONS)
            I = INT(SECTIONS(K) / XSTEP) + 1

            DO J = 1, NJ
                IF (UCELL(I, J) >= 0.99 * MAXVAL(UCELL(I, :))) THEN
                    DELTA = YCELL(I, J)
                    EXIT
                END IF
            END DO

            DO J = 1, NJ
                YAUTO(J) = YCELL(I, J) / DELTA
                UAUTO(J) = UCELL(I, J) / MAXVAL(UCELL(I, :))
            END DO

            WRITE(PROFILES_UNIT, *) 'VARIABLES = "U", "Y"'
            WRITE(PROFILES_UNIT, *) "ZONE I = ", NJ, ", DATAPACKING = BLOCK"
            WRITE(PROFILES_UNIT, "(100E25.16)") UAUTO
            WRITE(PROFILES_UNIT, "(100E25.16)") YAUTO
        END DO

        DEALLOCATE(YAUTO, UAUTO)

    END SUBROUTINE AUTOMODEL_VELOCITY_PROFILES

    SUBROUTINE CREATE_FILES(XNODE, YNODE, YCELL, UCELL, VCELL, PCELL, RE, CF, BL)
        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: RE, CF, BL
        DOUBLE PRECISION, INTENT(IN), DIMENSION(:, :) :: XNODE, YNODE
        DOUBLE PRECISION, INTENT(IN), DIMENSION(0:, 0:) :: YCELL, UCELL, VCELL, PCELL

        OPEN(UNIT=VELOCITY_UNIT, FILE=VELOCITY, ACTION="WRITE")
        OPEN(UNIT=FRICTION_UNIT, FILE=FRICTION, ACTION="WRITE")
        OPEN(UNIT=PROFILES_UNIT, FILE=PROFILES, ACTION="WRITE")

        WRITE(VELOCITY_UNIT, *) 'VARIABLES = "X", "Y", "U", "V", "P"'
        WRITE(VELOCITY_UNIT, *) "ZONE I = ", NI, ", J = ", NJ, ", DATAPACKING = BLOCK, &
            & VARLOCATION = ([3-20] = CELLCENTERED)"
        WRITE(VELOCITY_UNIT, "(100E25.16)") XNODE(1: NI, 1: NJ)
        WRITE(VELOCITY_UNIT, "(100E25.16)") YNODE(1: NI, 1: NJ)
        WRITE(VELOCITY_UNIT, "(100E25.16)") UCELL(1: NI - 1, 1: NJ - 1)
        WRITE(VELOCITY_UNIT, "(100E25.16)") VCELL(1: NI - 1, 1: NJ - 1)
        WRITE(VELOCITY_UNIT, "(100E25.16)") PCELL(1: NI - 1, 1: NJ - 1)

        WRITE(FRICTION_UNIT, *) 'VARIABLES = "RE", "CF"'
        WRITE(FRICTION_UNIT, *) "ZONE I = ", NI - 1, ", DATAPACKING = BLOCK"
        WRITE(FRICTION_UNIT, "(100E25.16)") RE(2: NI)
        WRITE(FRICTION_UNIT, "(100E25.16)") CF(2: NI)

        WRITE(FRICTION_UNIT, *) 'VARIABLES = "RE", "CF"'
        WRITE(FRICTION_UNIT, *) "ZONE I = ", NI - 1, ", DATAPACKING = BLOCK"
        WRITE(FRICTION_UNIT, "(100E25.16)") RE(2: NI)
        WRITE(FRICTION_UNIT, "(100E25.16)") BL(2: NI)

        CALL AUTOMODEL_VELOCITY_PROFILES(YCELL, UCELL)

        CLOSE(UNIT=VELOCITY_UNIT)
        CLOSE(UNIT=FRICTION_UNIT)
        CLOSE(UNIT=PROFILES_UNIT)

    END SUBROUTINE CREATE_FILES

END MODULE FUNCTIONS
