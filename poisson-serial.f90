PROGRAM main
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: u_exact, u0, u1, f
  DOUBLE PRECISION :: epsilon=1.e-3, resid
  INTEGER :: n=20, nit=100
!
  INTEGER :: i,j,it
  DOUBLE PRECISION :: temp
  DOUBLE PRECISION :: residue, error
!
  WRITE(*,*) 'Enter number of intervals n, nit, epsilon'
  READ(*,*) n, nit, epsilon
  ALLOCATE( f(n, n), &
       &    u_exact(0:n+1, 0:n+1), & 
       &    u0(0:n+1, 0:n+1), &
       &    u1(0:n+1, 0:n+1) )

!  Define RHS u_exact, f
  CALL RANDOM_NUMBER(u_exact)
  u_exact(0,:) = 0.0
  u_exact(n+1,:) = 0.0
  u_exact(:,0) = 0.0
  u_exact(:,n+1) = 0.0
  DO i=1,n
     DO j=1,n
        f(i,j) = u_exact(i-1,j) + u_exact(i,j+1) + &
             &   u_exact(i,j-1) + u_exact(i+1,j) - &
             &   4.0 * u_exact(i,j)
     END DO
  END DO
!
!  Jacobi iterations
  u0 = 0.0
  u1 = 0.0
  it=0
  DO
     CALL jacobi(u0, u1, f, n)
     CALL jacobi(u1, u0, f, n)
     it = it+2
     resid = residue(u0, f, n)
!!$!
!!$     WRITE(*,'(a,i5,a,1pe12.3)') 'it =', it, '  Residue =', resid
!
     IF( resid .LT. epsilon .OR. it .GT. nit) EXIT
  END DO

! Compute Errors
  error = 0.0
  DO i=1,n
     DO j=1,n
        temp = ABS(u_exact(i,j) - u0(i,j))
        error = MAX(error,temp)
     END DO
  END DO
  temp = error
  WRITE(*,'(a7,/(10f8.3))') 'u_exact ', u_exact(1:n,2)
  WRITE(*,'(a7,/(10f8.3))') 'u0 ', u0(1:n,2)
  WRITE(*,'(a,i12)')      '# of iterations =', it
  WRITE(*,'(a,1pe10.3)') 'Residue =', resid
  WRITE(*,'(a,1pe10.3)') 'Maximum error =', error
!
END PROGRAM main
                                                                  
SUBROUTINE jacobi(u0, u1, f, n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  DOUBLE PRECISION, INTENT(IN)  :: f(n,n)
  DOUBLE PRECISION, INTENT(IN)  :: u0(0:n+1,0:n+1)
  DOUBLE PRECISION, INTENT(OUT) :: u1(0:n+1,0:n+1)
!
  INTEGER :: i, j
!
  DO i = 1, n
     DO j = 1, n
        u1(i,j) = 0.25d0 * (u0(i-1,j) + u0(i,j+1) + &
             &   u0(i,j-1) + u0(i+1,j) - f(i,j))
     END DO
  END DO
END SUBROUTINE jacobi

FUNCTION residue(u, f, n)
  IMPLICIT NONE
  DOUBLE PRECISION :: residue
  INTEGER, INTENT(IN) :: n
  DOUBLE PRECISION, INTENT(IN) :: u(0:n+1,0:n+1)
  DOUBLE PRECISION, INTENT(IN) :: f(n,n)
!
  INTEGER :: i, j
!
  residue = 0.0
  DO i = 1, n
     DO j= 1, n
        residue=residue + ( u(i-1,j) + u(i,j+1) + &
             &   u(i,j-1) + u(i+1,j) - 4.0d0 * u(i,j) - f(i,j) )**2
     END DO
  END DO
  residue = SQRT(residue)
END FUNCTION residue
