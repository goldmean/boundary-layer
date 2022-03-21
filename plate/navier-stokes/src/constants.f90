MODULE CONSTANTS
    IMPLICIT NONE

    INTEGER, PARAMETER, PUBLIC :: NI = 51, NJ = 41
    DOUBLE PRECISION, PARAMETER, PUBLIC :: L = 1.0, H = 0.2, U0 = 1.0, V0 = 0.0, &
        P0 = 0.0, MU = 1.0, R0 = 1E3, CFL = 1E-2

    CHARACTER(*), PARAMETER, PUBLIC :: VELOCITY = "out/velocity_components.plt", &
        FRICTION = "out/friction_coefficient.plt", &
        PROFILES = "out/velocity_profiles.plt", &
        RESIDUAL = "out/variable_residuals.plt"

    DOUBLE PRECISION, PARAMETER, PUBLIC, DIMENSION(*) :: SECTIONS = [0.2, 0.6, 0.8]

END MODULE CONSTANTS
