PROGRAM sw_simulation
! Program runs a hysteresis simulation for a Stoner-Wohlfarth particle by minimizing its anisotropy energy.
! The energy functional and its parameters as well as the minimizer, are found in the sw_minimizer module.
    USE kinds
    USE sw_minimizer
    IMPLICIT NONE
    REAL(KIND=DBL) :: phi_min ! phi value for which minimal energy is obtained
    REAL(KIND=DBL) :: phi_0 ! initial guess for phi for minimizer
	REAL(KIND=DBL) :: tol ! small number to test convergence in minimizer
    REAL(KIND = DBL) :: PI = 3.14159265358979323846
    REAL(KIND = DBL), DIMENSION(200) :: H_range ! array for external magnetic field
    REAL(KIND = DBL), DIMENSION(3) :: gammas ! array for angle gamma
    REAL(KIND = DBL), DIMENSION(100) :: phi_range ! array for angle phi (energy profile plot)
    REAL(KIND = DBL), DIMENSION(8) :: H_values ! isolated values of H_ext (energy profile plot)
    REAL(KIND = DBL) :: energy_single_val 
    
    INTEGER :: i, j

    tol = 1.0e-5
    
    ! Fill array with 200 values from -2 to 2 and back to -2, linearly spaced
    DO i = 1, 100
        H_range(i) = -2.0 + (i-1)*0.04
    END DO
    DO i = 101, 200
        H_range(i) = 2.0 - (i-101)*0.04
    END DO
    

    gammas = [1D0*PI/180D0, 45D0*PI/180D0, 90D0*PI/180D0] ! gamma = {1 deg, 45 deg, 90 deg} 
	
	! Minimize energy and store values to plot hysteresis
	
    OPEN(unit = 1, file = "output.txt", status = "replace", action = "write")
    WRITE(3, *) "gamma	", "phi_min	", "H_ext"
    DO i = 1, 3
        phi_0 = 0.0D0
        DO j = 1, 200
            phi_min = conjugate_gradient(gammas(i), H_range(j), phi_0, tol, 2000)
            phi_0 = phi_min
            PRINT*, "Minimal angle for gamma =", gammas(i) * 180 / PI, "and H =", H_range(j), "is", phi_min * 180 / PI
            WRITE(3, *) gammas(i), phi_min, H_range(j)
        END DO
    END DO
    
    ! Calculate energy for several phi values to plot energy profile
    
    DO i = 1, 100
    	phi_range(i) = -6.0 + (i-1)*0.08
    END DO
    
    H_values = [0.25D0, 0.5D0, 0.75D0, 1D0, 1.25D0, 1.5D0, 1.75D0, 2D0]
    
    OPEN(unit = 2, file = "output_energies.txt", status = "replace", action = "write")
    WRITE(3, *) "phi	", "energy	",
    DO i = 0, 100
    
    DO j = 1, 8
    	DO i = 1, 100
    		energy_single_val = SW_energy(phi_range(i), 0D0, H_values(j))
    		WRITE(3, *) phi_range(i)*180/PI, energy_single_val
    	END DO
    END DO
    
    CLOSE(1)
    CLOSE(2)

END PROGRAM sw_simulation
