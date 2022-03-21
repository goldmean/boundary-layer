MODULE SYSTEM_SOLVERS

    PUBLIC :: SUPPLE_SOLVER, SIMPLE_SOLVER

CONTAINS

    SUBROUTINE SUPPLE_SOLVER(A, B, C, D, N, X)
        IMPLICIT NONE

        INTEGER :: I, N
        DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: A, B, C, D
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(:) :: X
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: P, Q, R, S, T, W

        ALLOCATE(P(N), Q(N), R(N), S(N), T(N), W(N))

        P(2) = D(1) / A(1)
        Q(2) = - B(1) / A(1)
        R(2) = - C(1) / A(1)

        DO I = 3, N
            W(I) = A(I - 1) + Q(I - 1) * C(I - 1)
            P(I) = (D(I - 1) - P(I - 1) * C(I - 1)) / W(I)
            Q(I) = - B(I - 1) / W(I)
            R(I) = - R(I - 1) * C(I - 1) / W(I)
        END DO

        S(N - 1) = Q(N) + R(N)
        T(N - 1) = P(N)

        DO I = N - 2, 1, -1
            S(I) = Q(I + 1) * S(I + 1) + R(I + 1)
            T(I) = Q(I + 1) * T(I + 1) + P(I + 1)
        END DO

        X(N) = (D(N) - B(N) * T(1) - P(N) * C(N)) / &
            (A(N) + B(N) * S(1) + C(N) * (Q(N) + R(N)))

        DO I = 1, N - 1
            X(I) = S(I) * X(N) + T(I)
        END DO

        DEALLOCATE(P, Q, R, S, T)

    END SUBROUTINE SUPPLE_SOLVER

    SUBROUTINE SIMPLE_SOLVER(B, C, A, D, N, X)
        IMPLICIT NONE

        INTEGER :: I, N
        DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: A, B, C, D
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(:) :: X
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: P, Q, R

        ALLOCATE(P(N), Q(N), R(N))

        P(1) = B(1)
        Q(1) = - C(1) / P(1)
        R(1) = D(1) / P(1)

        DO I = 2, N - 1
            P(I) = B(I) + A(I) * Q(I - 1)
            Q(I) = - C(I) / P(I)
            R(I) = (D(I) - A(I) * R(I - 1)) / P(I)
        END DO

        P(N) = B(N) + A(N) * Q(N - 1)
        R(N) = (D(N) - A(N) * R(N - 1)) / P(N)

        X(N) = R(N)

        DO I = N - 1, 1, -1
            X(I) = Q(I) * X(I + 1) + R(I)
        END DO

        DEALLOCATE(P, Q, R)

    END SUBROUTINE SIMPLE_SOLVER

END MODULE SYSTEM_SOLVERS
