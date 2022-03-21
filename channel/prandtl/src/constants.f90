MODULE CONSTANTS
    IMPLICIT NONE

    INTEGER, PARAMETER, PUBLIC :: NI = 101, NJ = 41
    DOUBLE PRECISION, PARAMETER, PUBLIC :: L = 1.0, H = 0.2, U0 = 1.0, V0 = 0.0, &
        P0 = 0.0, MU = 20.0, R0 = 1E3

    CHARACTER(*), PARAMETER, PUBLIC :: VELOCITY = "out/velocity_components.plt", &
        FRICTION = "out/friction_coefficient.plt"

END MODULE CONSTANTS
