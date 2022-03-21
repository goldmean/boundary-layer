MODULE FUNCTIONS

    USE CONSTANTS
    USE MATRIX_OPERATIONS, ONLY: MATRIX_INVERSE

    IMPLICIT NONE

    INTEGER, PARAMETER, PRIVATE :: VELOCITY_UNIT = 1, FRICTION_UNIT = 2

    DOUBLE PRECISION, PARAMETER, PRIVATE :: EPS = 1E-6, NU = MU / R0, &
        XSTEP = L / (NI - 1), YSTEP = H / (NJ - 1)

    PRIVATE :: INITIAL_CONDITIONS, BOUNDARY_CONDITIONS
    PUBLIC :: CREATE_FILES, CREATE_MESH, SIMUNI_PRANDTL_SOLVER, &
        NUMERICAL_FRICTION, BLASIUS_FRICTION

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

        U(:, :) = U0
        V(:, :) = V0
        P(:, :) = P0

    END SUBROUTINE INITIAL_CONDITIONS

    SUBROUTINE BOUNDARY_CONDITIONS(A, F, D)
        IMPLICIT NONE

        INTEGER :: I, J
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(:) :: F, D
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(:, :) :: A

        FORALL (I = 1: NJ, J = 1: NJ) A(I, J) = 0.0

        A(1, 1) = 1.0
        A(1, 2) = - 1.0

        A(NJ, NJ - 1) = 0.0
        A(NJ, NJ) = 1.0

        FORALL (J = 1: NJ) F(J) = 0.0
        FORALL (J = 1: NJ) D(J) = 0.0

    END SUBROUTINE BOUNDARY_CONDITIONS

    SUBROUTINE SIMUNI_PRANDTL_SOLVER(U, V, P)
        IMPLICIT NONE

        INTEGER :: I, J
        DOUBLE PRECISION :: FINT, DINT, UINT
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(:, :) :: U, V, P
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: A, B
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: F, D, FTMP, DTMP, &
            US, VS, PS

        ALLOCATE(A(NJ, NJ), B(NJ, NJ), F(NJ), D(NJ), DTMP(NJ), FTMP(NJ), &
            US(NJ), VS(NJ), PS(NJ))

        CALL INITIAL_CONDITIONS(U, V, P)

        DO I = 2, NI
            U(I, :) = U(I - 1, :)
            V(I, :) = V(I - 1, :)
            P(I, :) = P(I - 1, :)

            DO WHILE (.TRUE.)
                US = U(I, :)
                VS = V(I, :)
                PS = P(I, :)

                CALL BOUNDARY_CONDITIONS(A, F, D)

                DO J = 2, NJ - 1
                    A(J, J) = U(I, J) / XSTEP + 2.0 * NU / YSTEP ** 2
                    A(J, J + 1) = V(I, J + 1) / (2.0 * YSTEP) - NU / YSTEP ** 2
                    A(J, J - 1) = - V(I, J - 1) / (2.0 * YSTEP) - NU / YSTEP ** 2
                END DO

                DO J = 2, NJ - 1
                    F(J) = 1 / XSTEP
                    D(J) = U(I - 1, J) ** 2 / XSTEP + P(I - 1, J) / XSTEP
                END DO

                CALL MATRIX_INVERSE(A, B, NJ)

                FTMP = MATMUL(B, F)
                DTMP = MATMUL(B, D)

                FINT = 2.0 * SUM(FTMP(2: NJ - 1)) + FTMP(1) + FTMP(NJ)
                DINT = 2.0 * SUM(DTMP(2: NJ - 1)) + DTMP(1) + DTMP(NJ)
                UINT = 2.0 * SUM(U(I - 1, 2: NJ - 1)) + U(I - 1, 1) + U(I - 1, NJ)

                P(I, :) = (DINT - UINT) / FINT

                DO J = 1, NJ
                    U(I, J) = DTMP(J) - FTMP(J) * P(I, J)
                END DO

                DO J = 2, NJ
                    V(I, J) = V(I, J - 1) - YSTEP / (2.0 * XSTEP) * (U(I, J) - &
                        U(I - 1, J) + U(I, J - 1) - U(I - 1, J - 1))
                END DO

                IF (MAXVAL(ABS(U(I, :) - US)) / MAXVAL(ABS(U(I, :))) < EPS .AND. &
                    MAXVAL(ABS(V(I, :) - VS)) / MAXVAL(ABS(V(I, :))) < EPS .AND. &
                    MAXVAL(ABS(P(I, :) - PS)) / MAXVAL(ABS(P(I, :))) < EPS) EXIT
            END DO
        END DO

        DEALLOCATE(A, B, F, D, FTMP, DTMP, US, VS, PS)

    END SUBROUTINE SIMUNI_PRANDTL_SOLVER

    SUBROUTINE NUMERICAL_FRICTION(U, CF)
        IMPLICIT NONE

        INTEGER :: I
        DOUBLE PRECISION, INTENT(IN), DIMENSION(:, :) :: U
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(:) :: CF

        DO I = 1, NI
            CF(I) = - 2.0 * MU * (3.0 * U(I, NJ) - 4.0 * U(I, NJ - 1) + U(I, NJ - 2)) / &
                (2.0 * YSTEP) / (R0 * MAXVAL(U(I, :)) ** 2)  ! SECOND ORDER
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

    SUBROUTINE CREATE_FILES(X, Y, U, V, P, RE, CF, BL)
        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: RE, CF, BL
        DOUBLE PRECISION, INTENT(IN), DIMENSION(:, :) :: X, Y, U, V, P

        OPEN(UNIT=VELOCITY_UNIT, FILE=VELOCITY, ACTION="WRITE")
        OPEN(UNIT=FRICTION_UNIT, FILE=FRICTION, ACTION="WRITE")

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

        CLOSE(UNIT=VELOCITY_UNIT)
        CLOSE(UNIT=FRICTION_UNIT)

    END SUBROUTINE CREATE_FILES

END MODULE FUNCTIONS
