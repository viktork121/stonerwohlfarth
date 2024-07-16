MODULE sw_minimizer
    USE kinds
    IMPLICIT NONE

CONTAINS

    FUNCTION SW_energy(phi, gamma, H_ext) RESULT(energy)
        ! Energy of a Stoner-Wohlfarth particle. 
        ! Inputs: phi: angle of magnetization to easy axis,
        !         gamma: angle of external field to easy axis,
        !         H_ext: magnitude of external field

        REAL(kind=DBL), INTENT(IN) :: phi, gamma, H_ext
        REAL(kind=DBL) :: H_par, H_perp
        REAL(kind=DBL) :: energy

        H_par = H_ext*COS(gamma) ! Component of external field perpendicular to easy axis
        H_perp = H_ext*SIN(gamma) ! Component of external field parallel to easy axis

        energy = SIN(phi)**2 - (H_par * COS(phi) + H_perp * SIN(phi))
    END FUNCTION SW_energy

    FUNCTION SW_energy_grad(phi, gamma, H_ext) RESULT(grad)
        ! Helper function, calculates gradient of energy function
        ! with respect to phi using centered differences. 
        ! Inputs: same as above, new parameter h: grid spacing
        REAL(kind=DBL), INTENT(IN) :: phi, gamma, H_ext
        REAL(kind=DBL) :: grad, h = 1.0e-8 ! Due to the sensitivity of the minimization
                                           ! to the step size of the derivative, various
                                           ! numerical errors and inconsistensies can occur.
                                           ! Through various tests, h = 1.0e-8 gives the most
                                           ! accurate result.

        grad = (SW_energy(phi + h, gamma, H_ext) - SW_energy(phi - h, gamma, H_ext)) / (2.0 * h)
    END FUNCTION SW_energy_grad

    FUNCTION conjugate_gradient(gamma, H_ext, phi_start, tol, max_iter) RESULT(phi_min)
        ! Minimization using the Fletcher-Reeves algorithm of the conjugate gradient method, adapted to 1 dimension.
        ! For reference, see following link:
        ! https://indrag49.github.io/Numerical-Optimization/conjugate-gradient-methods-1.html#nonlinear-conjugate-gradient-algorithm
        ! New parameters:
        ! phi_start: initial guess for phi (or variable to be minimized by)
        ! tol: small number to check for convergence of method
        ! max_iter: number of iterations
        ! d: (conjugate) descent direction
        ! alpha: length of descent step
        ! chi : quotient of squared gradients

        REAL(kind=DBL), INTENT(IN) :: gamma, H_ext, phi_start, tol
        INTEGER, INTENT(IN) :: max_iter
        REAL(kind=DBL) :: phi_min, grad, new_grad, d, new_d, alpha = 0.01, chi
        INTEGER :: iter

        phi_min = phi_start
        grad = SW_energy_grad(phi_min, gamma, H_ext)
        d = -grad  ! Initial direction is steepest descent

        DO iter = 1, max_iter

            phi_min = phi_min + alpha * d
            new_grad = SW_energy_grad(phi_min, gamma, H_ext)
            
            ! After calculating gradient at new position, new (conjugate) descent
            ! direction is calculated via the squared gradients. 

            chi =  (new_grad**2)/(grad**2)
            new_d = -new_grad + chi * d
            grad = new_grad
            d = new_d

            IF (ABS(grad) < tol) EXIT  ! Convergence check
        END DO
    END FUNCTION conjugate_gradient

END module sw_minimizer
