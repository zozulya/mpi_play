PROGRAM main
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: u_exact, u0, u1, f
  DOUBLE PRECISION :: epsilon=1.e-3, resid
  INTEGER :: n=20, nit=100
!
  INTEGER :: s, e, nlocal         ! 1D Domain decomposition in x
  INTEGER :: me, npes, ierr ! Parallel environment
  INTEGER :: prev, next       ! Neighbors
!
  INTEGER :: i,j,it
  DOUBLE PRECISION :: temp, error
  DOUBLE PRECISION :: residue
!
!  Parallel environment
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, npes, ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)
  prev = me-1
  next = me+1
  IF( me .EQ. 0 )      prev = MPI_PROC_NULL
  IF( me .EQ. npes-1 ) next = MPI_PROC_NULL
!
!  Read and broadcast inputs
  IF( me.EQ.0 ) THEN
     WRITE(*,*) 'Enter number of points n, nit, epsilon'
     READ(*,*) n, nit, epsilon
  END IF
  CALL MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(nit, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(epsilon, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
!
!  1D Domain decomposition allong the 2nd dimension
  CALL dist1d(1, n, s, nlocal)
  e = s + nlocal - 1
  WRITE(*,'(a,i3.3,a,4i8)') 'PE', me, '  prev, next, s, e', &
       & prev, next, s, e
!
  ALLOCATE( f(1:n, s:e), &
       &    u_exact(0:n+1, s-1:e+1), & 
       &    u0(0:n+1, s-1:e+1), &     ! Ghost cells for solutions
       &    u1(0:n+1, s-1:e+1) )

!  Define u_exact and f
  CALL RANDOM_NUMBER(u_exact)
  u_exact(0,:) = 0.0
  u_exact(n+1,:) = 0.0
  IF( s.EQ.1 ) u_exact(:,s-1) = 0.0
  IF( e.EQ.n) u_exact(:,e+1) = 0.0
  CALL exchange(u_exact, s, e, n, prev, next)
  DO i=1,n
     DO j=s,e
        f(i,j) = u_exact(i-1,j) + u_exact(i,j+1) + &
             &   u_exact(i,j-1) + u_exact(i+1,j) - &
             &   4.0d0 * u_exact(i,j)
     END DO
  END DO
!
!  Jacobi iterations
  u0 = 0.0
  u1 = 0.0
  DO it=1,nit,2
     CALL jacobi(u0, u1, f, s, e, n)
         CALL exchange(u1, s, e, n, prev, next)
     CALL jacobi(u1, u0, f, s, e, n)
         CALL exchange(u0, s, e, n, prev, next)
     resid = residue(u0, f, s, e, n)
!
!!$     IF( me.EQ.0 ) THEN
!!$        WRITE(*,'(a,i5,a,1pe12.3)') 'it =', it+1, '  Residue =', resid
!!$     END IF
!
     IF( resid .LT. epsilon ) EXIT
  END DO

! Compute Errors
  IF( me.EQ.0 ) THEN
     WRITE(*,'(a/(10f8.3))') 'u_exact ', u_exact(1:n,s+1)
     WRITE(*,'(a/(10f8.3))') 'u0 ', u0(1:n,s+1)
   END IF
  error = 0.0
  DO i=1,n
     DO j=s,e
        temp = ABS(u_exact(i,j) - u0(i,j))
        error = MAX(error,temp)
     END DO
  END DO
  temp = error
  CALL MPI_ALLREDUCE(temp, error, 1, MPI_DOUBLE_PRECISION,&
          &   MPI_MAX, MPI_COMM_WORLD, ierr)
  IF( me.EQ.0 ) THEN
     WRITE(*,'(a,i12)')      '# of iterations =', it+1
     WRITE(*,'(a,1pe10.3)') 'Residue =', resid
     WRITE(*,'(a,1pe12.3)') 'Maximum error   =', error
  END IF
!
  CALL MPI_FINALIZE(ierr)
END PROGRAM main
                                                                  

SUBROUTINE dist1d(s0, ntot, s, nloc)
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, INTENT(in) :: s0, ntot
  INTEGER, INTENT(out) :: s, nloc
  INTEGER :: me, npes, ierr, naver, rem
!
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, npes, ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)
  naver = ntot/npes
  rem = MODULO(ntot,npes)
  s = s0 + MIN(rem,me) + me*naver
  nloc = naver
  IF( me.LT.rem ) nloc = nloc+1
!
END SUBROUTINE dist1d


SUBROUTINE jacobi(u0, u1, f, s, e, n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: s, e, n
  DOUBLE PRECISION, INTENT(IN)  :: f(1:n,s:e)
  DOUBLE PRECISION, INTENT(IN)  :: u0(0:n+1,s-1:e+1)
  DOUBLE PRECISION, INTENT(OUT) :: u1(0:n+1,s-1:e+1)
!
  INTEGER :: i, j
!
  DO i = 1, n
     DO j = s, e
        u1(i,j) = 0.25d0 * (u0(i-1,j) + u0(i,j+1) + &
             &   u0(i,j-1) + u0(i+1,j) - f(i,j))
     END DO
  END DO
END SUBROUTINE jacobi


FUNCTION residue(u, f, s, e, n)
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  DOUBLE PRECISION :: residue
  INTEGER, INTENT(IN) :: s, e, n
  DOUBLE PRECISION, INTENT(IN) :: u(0:n+1,s-1:e+1)
  DOUBLE PRECISION, INTENT(IN) :: f(1:n,s:e)
!
  DOUBLE PRECISION :: temp1, temp2
  INTEGER :: ierr
  INTEGER :: i, j
!
  temp1 = 0.0
  DO i = 1, n
     DO j= s, e
        temp1 = temp1 + ( u(i-1,j) + u(i,j+1) + &
             &   u(i,j-1) + u(i+1,j) - &
             &   4.0d0 * u(i,j) - f(i,j) )**2
     END DO
  END DO
  CALL MPI_ALLREDUCE(temp1, temp2, 1, MPI_DOUBLE_PRECISION,&
          &   MPI_SUM, MPI_COMM_WORLD, ierr)
  residue = SQRT(temp2)
END FUNCTION residue


SUBROUTINE exchange(u, s, e, n, prev, next)
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER, INTENT(IN) :: s, e, n, prev, next
  DOUBLE PRECISION, INTENT(INOUT) :: u(0:n+1,s-1:e+1)
!
  INTEGER :: stat(MPI_STATUS_SIZE, 4), req(4), ierr
!
  CALL MPI_ISEND(u(1,s),   n, MPI_DOUBLE_PRECISION, prev, &
       &0, MPI_COMM_WORLD, req(1), ierr)
  CALL MPI_ISEND(u(1,e),   n, MPI_DOUBLE_PRECISION, next, &
       &1, MPI_COMM_WORLD, req(2), ierr)
  CALL MPI_IRECV(u(1,e+1), n, MPI_DOUBLE_PRECISION, next, &
       &0, MPI_COMM_WORLD, req(3), ierr)
  CALL MPI_IRECV(u(1,s-1), n, MPI_DOUBLE_PRECISION, prev, &
       &1, MPI_COMM_WORLD, req(4), ierr)
  CALL MPI_WAITALL(4, req, stat, ierr)
END SUBROUTINE exchange
