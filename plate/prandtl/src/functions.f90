MODULE FUNCTIONS

    USE CONSTANTS
    USE SYSTEM_SOLVERS, ONLY: SIMPLE_SOLVER

    IMPLICIT NONE

    INTEGER, PARAMETER, PRIVATE :: VELOCITY_UNIT = 1, FRICTION_UNIT = 2, &
        PROFILES_UNIT = 3

    DOUBLE PRECISION, PARAMETER, PRIVATE :: EPS = 1E-6, NU = MU / R0, &
        XSTEP = L / (NI - 1), YSTEP = H / (NJ - 1)

    PRIVATE :: INITIAL_CONDITIONS, BOUNDARY_CONDITIONS, AUTOMODEL_VELOCITY_PROFILES
    PUBLIC :: CREATE_FILES, CREATE_MESH, PRANDTL_SOLVER, NUMERICAL_FRICTION, &
        BLASIUS_FRICTION

CONTAINS

    SUBROUTINE CREATE_MESH(X, Y)
        IMPLICIT NONE

        INTEGER :: I, J
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(:, :) :: X, Y

        DO I = 1, NI
            DO J = 1, NJ
                X(I, J) = (I - 1) * XSTEP
                Y(I, J) = (J - 1) * YSTEP
            END DO
        END DO

    END SUBROUTINE CREATE_MESH

    SUBROUTINE INITIAL_CONDITIONS(U, V, P)
        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(OUT), DIMENSION(:, :) :: U, V, P

        U(1, :) = U0
        V(1, :) = V0
        P(:, :) = P0

    END SUBROUTINE INITIAL_CONDITIONS

    SUBROUTINE BOUNDARY_CONDITIONS(A, B, C, D, V)
        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(OUT), DIMENSION(:) :: A, B, C, D
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(:, :) :: V

        A(1) = 0.0
        B(1) = 1.0
        C(1) = 0.0
        D(1) = 0.0

        A(NJ) = 0.0
        B(NJ) = 1.0
        C(NJ) = 0.0
        D(NJ) = 1.0

        V(:, 1) = V0

    END SUBROUTINE BOUNDARY_CONDITIONS

    SUBROUTINE PRANDTL_SOLVER(U, V, P)
        IMPLICIT NONE

        INTEGER :: I, J
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(:, :) :: U, V, P
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: A, B, C, D, US, VS

        ALLOCATE(A(NJ), B(NJ), C(NJ), D(NJ), US(NJ), VS(NJ))

        CALL INITIAL_CONDITIONS(U, V, P)
        CALL BOUNDARY_CONDITIONS(A, B, C, D, V)

        DO I = 2, NI
            U(I, :) = U(I - 1, :)
            V(I, :) = V(I - 1, :)

            DO WHILE (.TRUE.)
                US(:) = U(I, :)
                VS(:) = V(I, :)

                DO J = 2, NJ - 1
                    A(J) = - V(I, J - 1) / (2.0 * YSTEP) - NU / YSTEP ** 2
                    B(J) = U(I, J) / XSTEP + 2.0 * NU / YSTEP ** 2
                    C(J) = V(I, J + 1) / (2.0 * YSTEP) - NU / YSTEP ** 2
                    D(J) = U(I - 1, J) ** 2 / XSTEP
                END DO

                CALL SIMPLE_SOLVER(B, C, A, D, NJ, U(I, :))

                DO J = 2, NJ
                    V(I, J) = V(I, J - 1) - YSTEP / (2.0 * XSTEP) * (U(I, J) - &
                        U(I - 1, J) + U(I, J - 1) - U(I - 1, J - 1))
                END DO

                IF (MAXVAL(ABS(U(I, :) - US(:))) / MAXVAL(ABS(U(I, :))) < EPS .AND. &
                    MAXVAL(ABS(V(I, :) - VS(:))) / MAXVAL(ABS(V(I, :))) < EPS) EXIT
            END DO
        END DO

        DEALLOCATE(A, B, C, D, US, VS)

    END SUBROUTINE PRANDTL_SOLVER

    SUBROUTINE NUMERICAL_FRICTION(U, CF)
        IMPLICIT NONE

        INTEGER :: I
        DOUBLE PRECISION, INTENT(IN), DIMENSION(:, :) :: U
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(:) :: CF

        DO I = 1, NI
            CF(I) = - 2.0 * MU * (3.0 * U(I, 1) - 4.0 * U(I, 2) + U(I, 3)) / &
                (2.0 * YSTEP) / (R0 * U0 ** 2)  ! SECOND ORDER
        END DO

    END SUBROUTINE NUMERICAL_FRICTION

    SUBROUTINE BLASIUS_FRICTION(X, RE, BL)
        IMPLICIT NONE

        INTEGER :: I
        DOUBLE PRECISION, INTENT(IN), DIMENSION(:, :) :: X
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(:) :: RE, BL

        DO I = 1, NI
            RE(I) = U0 * X(I, NJ) / NU
            BL(I) = 0.664 / SQRT(RE(I))
        END DO

    END SUBROUTINE BLASIUS_FRICTION

    SUBROUTINE AUTOMODEL_VELOCITY_PROFILES(Y, U)
        IMPLICIT NONE

        INTEGER :: I, J, K
        DOUBLE PRECISION :: DELTA
        DOUBLE PRECISION, INTENT(IN), DIMENSION(:, :) :: Y, U
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: YAUTO, UAUTO

        ALLOCATE(YAUTO(NJ), UAUTO(NJ))

        DO K = 1, SIZE(SECTIONS)
            I = INT(SECTIONS(K) / XSTEP) + 1

            DO J = 1, NJ
                IF (U(I, J) >= 0.99 * MAXVAL(U(I, :))) THEN
                    DELTA = Y(I, J)
                    EXIT
                END IF
            END DO

            DO J = 1, NJ
                YAUTO(J) = Y(I, J) / DELTA
                UAUTO(J) = U(I, J) / MAXVAL(U(I, :))
            END DO

            WRITE(PROFILES_UNIT, *) 'VARIABLES = "U", "Y"'
            WRITE(PROFILES_UNIT, *) "ZONE I = ", NJ, ", DATAPACKING = BLOCK"
            WRITE(PROFILES_UNIT, "(100E25.16)") UAUTO
            WRITE(PROFILES_UNIT, "(100E25.16)") YAUTO
        END DO

        DEALLOCATE(YAUTO, UAUTO)

    END SUBROUTINE AUTOMODEL_VELOCITY_PROFILES

    SUBROUTINE CREATE_FILES(X, Y, U, V, P, RE, CF, BL)
        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: RE, CF, BL
        DOUBLE PRECISION, INTENT(IN), DIMENSION(:, :) :: X, Y, U, V, P

        OPEN(UNIT=VELOCITY_UNIT, FILE=VELOCITY, ACTION="WRITE")
        OPEN(UNIT=FRICTION_UNIT, FILE=FRICTION, ACTION="WRITE")
        OPEN(UNIT=PROFILES_UNIT, FILE=PROFILES, ACTION="WRITE")

        WRITE(VELOCITY_UNIT, *) 'VARIABLES = "X", "Y", "U", "V", "P"'
        WRITE(VELOCITY_UNIT, *) "ZONE I = ", NI, ", J = ", NJ, ", DATAPACKING = BLOCK"
        WRITE(VELOCITY_UNIT, "(100E25.16)") X(1: NI, 1: NJ)
        WRITE(VELOCITY_UNIT, "(100E25.16)") Y(1: NI, 1: NJ)
        WRITE(VELOCITY_UNIT, "(100E25.16)") U(1: NI, 1: NJ)
        WRITE(VELOCITY_UNIT, "(100E25.16)") V(1: NI, 1: NJ)
        WRITE(VELOCITY_UNIT, "(100E25.16)") P(1: NI, 1: NJ)

        WRITE(FRICTION_UNIT, *) 'VARIABLES = "RE", "CF"'
        WRITE(FRICTION_UNIT, *) "ZONE I = ", NI - 1, ", DATAPACKING = BLOCK"
        WRITE(FRICTION_UNIT, "(100E25.16)") RE(2: NI)
        WRITE(FRICTION_UNIT, "(100E25.16)") CF(2: NI)

        WRITE(FRICTION_UNIT, *) 'VARIABLES = "RE", "CF"'
        WRITE(FRICTION_UNIT, *) "ZONE I = ", NI - 1, ", DATAPACKING = BLOCK"
        WRITE(FRICTION_UNIT, "(100E25.16)") RE(2: NI)
        WRITE(FRICTION_UNIT, "(100E25.16)") BL(2: NI)

        CALL AUTOMODEL_VELOCITY_PROFILES(Y, U)

        CLOSE(UNIT=VELOCITY_UNIT)
        CLOSE(UNIT=FRICTION_UNIT)
        CLOSE(UNIT=PROFILES_UNIT)

    END SUBROUTINE CREATE_FILES

END MODULE FUNCTIONS
