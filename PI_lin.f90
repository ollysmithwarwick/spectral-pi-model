program PI_lin
! --------------------------------------------------------------------------!
! A program for solving the plasma interchange model using spectral methods !
! When compiling one must have FFTW3 loaded, include fftw3.f03 and link the fftw3 library !

  use, intrinsic :: iso_c_binding 
  implicit none
  include 'fftw3.f03'

! DECLARE VARIABLES !


  type(C_PTR) :: plan, plan2  ! plans for input to fftws - declared later
  complex(C_DOUBLE_COMPLEX), dimension(:), allocatable :: in, out                                                   ! memory allocated for input and outputs of ffts
  complex(C_DOUBLE_COMPLEX), dimension(:), allocatable :: n0, a_n0, n1, a_n1, phi0, a_phi0, phi1, a_phi1, E, b_plus ! 0 corresponds to slow varying background, 1 to fast varying part, 
                                                                                                                    ! a to fourier space, lack of a to real space, E to gradient of phi, and b_plus to n1+iphi1
  complex(C_DOUBLE_COMPLEX), dimension(:), allocatable :: dxxn1, a_dxxn1, dxxphi1, a_dxxphi1                        !Second derivatives of fast varying terms
  complex(C_DOUBLE_COMPLEX), dimension(:), allocatable :: new_n1, new_phi1                                          ! Temporary arrays to hold new n and phi     

  integer, parameter :: N_x = 101, N_t = 1000                ! N_x is number of points in x (and in fourier space) and N_t the number of timesteps to calculate
  integer, parameter :: dp = selected_real_kind(15,307)      ! Gives a double precision value to be used in variable definitions
  integer(8) :: j, k                                         ! Integers of do loops
  real(dp) :: f, f2                                          ! Helpful divisor terms
  real(dp), parameter :: pi = 4*atan(1.0_dp)                 
  real(dp), parameter :: L = 8*atan(1.0_dp)                  ! Length of periodic box
  real(dp), parameter :: dt = 0.001                          ! Length of timestep
  complex(dp), parameter :: i = (0.0_dp, 1.0_dp)             
  real(dp), parameter :: T_end = 1.0_dp

! Allocate sizes to all complex arrays !
  allocate(in(N_x))
  allocate(out(N_x))
  allocate(n0(N_x))
  allocate(a_n0(N_x))
  allocate(phi0(N_x))
  allocate(a_phi0(N_x))
  allocate(n1(N_x))
  allocate(a_n1(N_x))
  allocate(phi1(N_x))
  allocate(a_phi1(N_x))
  allocate(dxxn1(N_x))
  allocate(a_dxxn1(N_x))
  allocate(dxxphi1(N_x))
  allocate(a_dxxphi1(N_x))
  allocate(new_n1(N_x))
  allocate(new_phi1(N_x))
  allocate(E(N_x))
  allocate(b_plus(N_x))

! Declare plans - only difference is the direction of the transform
  plan = fftw_plan_dft_1d(N_x, in, out, FFTW_FORWARD, FFTW_ESTIMATE)
  plan2 = fftw_plan_dft_1d(N_x, in, out, FFTW_BACKWARD, FFTW_ESTIMATE)

! Initialise variables (update to read an input file later)
  phi0(:) = 0.0_dp
  n0(:) = 0.0_dp
  phi1(:) = 0.0_dp

! For now use a simple harmonic function of frequency 1
  do k = 1, N_x
     n1(k) = exp(i*2*(k-1)*L/N_x)
     b_plus(k) = phi1(k)+i*n1(k)
  end do

! Commonly used divisions taken out to reduce calculation time
  f = (1.0_dp/(L*L))
  f2 = sqrt(1.0_dp/N_x)

! Open output files for writing
  open(2, file = 'phi0.dat')
  open(3, file = 'phi1.dat')
  open(4, file = 'n0.dat')
  open(5, file = 'n1.dat')
  open(6, file = 'b_plus.dat')

! Write first lines of output files with inputs
  write(2,1) phi0
  write(3,1) phi1
  write(4,1) n0
  write(5,1) n1
  write(6,1) b_plus

! Transform into k space using fftw 
  in = phi1
  call fftw_execute_dft(plan, in, out)
  a_phi1 = out*f2

  in = n1
  call fftw_execute_dft(plan, in, out)
  a_n1 = out*f2
  
  write(*,*) a_n1

! Split based on whether N_x is even or odd - currently no significant difference but allows for adaptation later 
  if (mod(N_x,2) == 0) then
     ! Iterate over timesteps
     do j = 1, N_t
        ! Find the fourier space coefficents of second derivatives
        do k = 0, (N_x/2)
           a_dxxphi1(k+1) = -4 * pi * pi * (k) * (k) * f * a_phi1(k+1)
           a_dxxn1(k+1) = -4 * pi * pi * (k) * (k) * f * a_n1(k+1)
        end do
     
        do k = int(N_x/2+1), N_x-1
           a_dxxphi1(k+1) = -4 * pi * pi * (k-N_x) * (k-N_x) * f * a_phi1(k+1)
           a_dxxn1(k+1) = -4 * pi * pi * (k-N_x) * (k-N_x) * f * a_n1(k+1)
        end do
        
        ! Fourier transform to real space
        in = a_dxxphi1
        call fftw_execute_dft(plan2, in, out)
        dxxphi1 = out*f2
        
        in = a_dxxn1
        call fftw_execute_dft(plan2, in, out)
        dxxn1 = out*f2
        
        ! Apply the equations using a numerical method scheme
        do k = 1, N_x
           new_n1(k) = n1(k) + dt * ( -i * phi1(k) + dxxn1(k))
           new_phi1(k) = phi1(k) + dt * ( i * n1(k) + dxxphi1(k))
        end do
        
        ! Shift temporary data into permanent
        n1(:) = new_n1(:)
        phi1(:) = new_phi1(:)

        ! Transorm into F-space to give new fourier space terms (could move this to start of loop)
        in = phi1
        call fftw_execute_dft(plan, in, out)
        a_phi1 = out*f2

        in = n1
        call fftw_execute_dft(plan, in, out)
        a_n1 = out*f2

        ! Test the a_+ (called b_plus here) function described in the paper to see if it works as expected in this linear case
        do k = 1, N_x
           b_plus(k) = phi1(k)+i*n1(k)
        end do

        ! Write outputs
        write(2,1) phi0
        write(3,1) phi1
        write(4,1) n0
        write(5,1) n1
        write(6,1) b_plus
     end do
     
  else if (mod(N_x,2) == 1) then

     do j = 1, N_t

        do k = 0, (N_x/2)
           a_dxxphi1(k+1) = -4 * pi * pi * (k) * (k) * f * a_phi1(k+1)
           a_dxxn1(k+1) = -4 * pi * pi * (k) * (k) * f * a_n1(k+1)
        end do
     
        do k = int(N_x/2+1), N_x-1
           a_dxxphi1(k+1) = -4 * pi * pi * (k-N_x) * (k-N_x) * 1/L * a_phi1(k+1)
           a_dxxn1(k+1) = -4 * pi * pi * (k-N_x) * (k-N_x) * 1/L * a_n1(k+1)
        end do

        in = a_dxxphi1
        call fftw_execute_dft(plan2, in, out)
        dxxphi1 = out*f2
        
        in = a_dxxn1
        call fftw_execute_dft(plan2, in, out)
        dxxn1 = out*f2
 
        do k = 1, N_x
           new_n1(k) = n1(k) + dt * ( -i * phi1(k) + dxxn1(k))
           new_phi1(k) = phi1(k) + dt * ( i * n1(k) + dxxphi1(k))
        end do

        n1(:) = new_n1(:)
        phi1(:) = new_phi1(:)

        in = phi1
        call fftw_execute_dft(plan, in, out)
        a_phi1 = out*f2

        in = n1
        call fftw_execute_dft(plan, in, out)
        a_n1 = out*f2

        do k = 1, N_x
           b_plus(k) = phi1(k)+i*n1(k)
        end do

        write(2,1) phi0
        write(3,1) phi1
        write(4,1) n0
        write(5,1) n1
        write(6,1) b_plus
   
     end do
  end if

close(2)
close(3)
close(4)
close(5)

1 format(1000f12.5)

end program PI_lin
