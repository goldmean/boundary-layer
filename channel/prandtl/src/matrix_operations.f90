MODULE MATRIX_OPERATIONS

    PUBLIC :: MATRIX_INVERSE

CONTAINS

    SUBROUTINE MATRIX_INVERSE(X, Y, N)
        IMPLICIT NONE

        INTEGER :: I, J, K, N
        DOUBLE PRECISION :: T
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: A, B, C
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: L, U
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(:, :) :: X, Y

        ALLOCATE(L(N, N), U(N, N), A(N), B(N), C(N))

        L = 0.0
        U = 0.0
        B = 0.0

        DO K = 1, N - 1
            DO I = K + 1, N
                T = X(I, K) / X(K, K)
                L(I, K) = T
                DO J = K + 1, N
                    X(I, J) = X(I, J) - T * X(K, J)
                END DO
            END DO
        END DO

        FORALL(I = 1: N) L(I, I) = 1.0

        DO J = 1, N
            DO I = 1, J
                U(I, J) = X(I, J)
            END DO
        END DO

        DO K = 1, N
            B(K) = 1.0
            A(1) = B(1)

            DO I = 2, N
                A(I) = B(I)
                DO J = 1, I - 1
                    A(I) = A(I) - L(I, J) * A(J)
                END DO
            END DO

            C(N) = A(N) / U(N, N)

            DO I = N - 1, 1, -1
                C(I) = A(I)
                DO J = N, I + 1, -1
                    C(I) = C(I) - U(I, J) * C(J)
                END DO
                C(I) = C(I) / U(I, I)
            END DO

            FORALL (I = 1: N) Y(I, K) = C(I)

            B(K) = 0.0
        END DO

        DEALLOCATE(L, U, A, B, C)

    END SUBROUTINE MATRIX_INVERSE

END MODULE MATRIX_OPERATIONS
