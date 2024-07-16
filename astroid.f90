PROGRAM sw_astroid
! Program calculates the Stoner Wohlfarth astroid (switching field curve)
! numerically by minimizing the magnetic energy, for a given set of material
! parameters, i.e. for the corresponding critical field.

    USE kinds
    USE sw_minimizer
    IMPLICIT NONE
    
    REAL(KIND = DBL), PARAMETER :: PI = 3.14159265358979323846
    REAL(KIND = DBL) :: Hc, H_par, H_perp
    REAL(KIND = DBL) :: phi, phi_0
    REAL(KIND = DBL) :: tol = 1E-05
    REAL(KIND = DBL), DIMENSION(361) :: gammas
    REAL(KIND = DBL), DIMENSION(100) :: H_values
    INTEGER :: i, j
    
    Hc = 2 ! Critical field - dimensionless units.
    ! If one were to use SI units, the critical
  	! field is calculated as Hc = 2*K/(mu_0*M_s),
  	! with K = anisotropy constant,
    ! M_s = saturation magnetization. Beware however,
    ! that in this case, the functions in the minimizer
    ! module must be modified to include the corresponding
    ! constants.
    
    DO i = 1, 100
    	H_values(i) = 0.001D0 + (i-1)*Hc/100.
    END DO
    
    DO i = 1, 361
  		gammas(i) = -91D0 + i ! gammas is the range [-90, 270]
  	END DO
  	    
  	OPEN(unit = 1, file = "output_astroid.txt", status = "replace", action = "write")
  	WRITE(1, *) "H_par	", "H_perp"
  	DO i = 1, 361
  		IF (gammas(i) < 91) THEN
  			phi = PI
  			DO j = 1, 100
  			  	H_par = H_values(j)*COS(gammas(i)*PI/180.)
  			  	H_perp = H_values(j)*SIN(gammas(i)*PI/180.)
  				phi = conjugate_gradient(gammas(i)*PI/180., H_values(j), phi, tol, 2000)
  				IF (COS(phi) > 0.) EXIT ! this condition essentially checks if the
  										! magnetization has a non-negative component
  										! in the direction of the ext. field, i.e.
  										! if the 'switching' has occured.
  			END DO
			WRITE(1, *) H_par, H_perp
  		ELSE
  			phi = 0
  			DO j = 1, 100
  			  	H_par = H_values(j)*COS(gammas(i)*PI/180.)
  			  	H_perp = H_values(j)*SIN(gammas(i)*PI/180.)
  				phi = conjugate_gradient(gammas(i)*PI/180., H_values(j), phi, tol, 2000)
  				IF (COS(phi) < 0.) EXIT ! the condition has to be inverted in this case
  									    ! since we are on the other "half" of the astroid
  			END DO
  			WRITE(1, *) H_par, H_perp
  		END IF 
  	END DO
  	
  	CLOSE(1)

END PROGRAM

