module co_random_numbers
  use kinds
  implicit none
  private
  !  The following code, originally written for mpi by Waseem Kamleh
  !  and adapted to Coarray Fortran by Curtis Abell,
  !  creates independent streams of random numbers on each processor
  !  of a parallel Coarray Fortran program.
  !
  !  It employs a linear congruential generator that starts with a
  !  global integer seed (iseed) and creates image dependent seeds.
  !  Once each coarray has been independently seeded
  !  the Fortran intrinsic random_number generator can be used.
  !
  !  The process dependent seed avoids spurious correlations, and the
  !  global integer seed means that if you run again on the same
  !  number of processes you will generate the same sequence of random
  !  numbers. This provides the benefit of code repeatability, which
  !  time based seeds do not.

  public :: co_random_seed

contains

  subroutine co_random_seed(iseed, iImg)

    integer, intent(in) :: iseed, iImg

    !  local variables
    !
    integer :: N
    integer, dimension(:), allocatable :: seeder_array

    !  Find out the size of the seed array used in the Fortran 90
    !  random number generator
    !
    call random_seed(size=N)

    allocate ( seeder_array(N) )

    !  Initialize the seed array by setting each bit using a
    !  simple linear congruential generator (in bitrand_seeder)
    !
    !  Because this code executes identically in each process, we
    !  include an image dependent seed to ensure different random
    !  numbers are generated on each process.

    call bitrand_seeder(iseed+iImg-1, N, seeder_array)

    call random_seed( put=seeder_array(1:N) )

    deallocate( seeder_array )

  end subroutine co_random_seed


  !  Initialize each bit of the 32-bit integer seed array using a
  !  simple linear congruential generator.
  !
  subroutine bitrand_seeder(seed, N, seed_array)

    integer :: seed, N
    integer :: seed_array(N)

    !  local variables
    !
    integer :: i,j,k
    integer :: s,t,randx
    real(dp), parameter :: TWOTO31 = 2147483648.0   ! 2^{31}
    real(dp), parameter :: TWOTONEG31=(1.0/TWOTO31)
    integer, parameter :: A=1103515245, C=12345
    integer :: M
    ! data  M   /Z'7FFFFFFF'/
    M = Z'7FFFFFFF' ! Z interprets as hexadecimal

    !  Initialize the linear congruential generator RAND.
    !
    randx = seed

    do i = 1, N
       s = 0
       t = 1
       do j = 0, 30
          randx = iand((A * randx + C),M)
          if((randx * TWOTONEG31) .lt. 0.5) then
             s = s + t
          endif
          t = 2 * t
          !  Throw away a few random numbers just to foil any
          !  possible correlations that might affect RAND.
          do k=1,5
             randx = iand((A * randx + C),M)
          end do
       end do

       !  Throw away a few random numbers just to foil any
       !  possible correlations that might affect RAND.
       do k=1,13
          randx = iand((A * randx + C),M)
       end do
       seed_array(i) = s
    end do

  end subroutine bitrand_seeder

end module co_random_numbers
