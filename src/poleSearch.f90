program poleSearch
    use kinds
    use numFort
    use heft
    use SMatrix
    use bobyqa_module
    implicit none

    ! Data input
    integer :: nPointsData
    real(DP), dimension(:), allocatable :: E_data
    real(DP), dimension(:), allocatable :: etaData
    real(DP), dimension(:), allocatable :: etaDataErr
    real(DP), dimension(:), allocatable :: phaseShiftData
    real(DP), dimension(:), allocatable :: phaseShiftDataErr

    ! Scattering output
    integer :: nPoints = 121
    ! real(DP) :: E_init = 1.11_DP, E_fin = 1.8_DP
    real(DP), dimension(:), allocatable :: E_arr
    real(DP), dimension(:), allocatable :: eta_arr
    real(DP), dimension(:,:), allocatable :: phaseShift
    real(DP), dimension(:,:,:), allocatable :: crossSec

    ! Misc
    integer :: i, j, iii, jjj ! misc loops
    real(DP), dimension(2,2) :: id2 = reshape((/ 1,0,0,1 /) , shape(id2))

    real(DP), dimension(2) :: searchError = (/ 1e-6_DP, 1e-6_DP /)
    real(DP) :: searchMaxStepFactor = 100.0_DP

    integer :: zi, zj, zk
    ! integer, parameter :: nPoles = 2
    ! complex(DP), dimension(nPoles) :: E_poles
    complex(DP) :: E_pole_guess
    real(DP), dimension(2) :: E_pole2
    real(DP) :: poleInverse

    integer :: maxIterations

    ! bobyqa variables
    integer :: bq_npt
    real(DP) :: bq_rhobeg, bq_rhoend
    integer :: bq_iprint
    integer :: bq_maxfun
    logical :: useBOBYQA = .true.
    real(DP), dimension(2) :: pole_low, pole_high

    call initialiseHEFT()
    call printCurrentParameters(iParamChoice)
    call initialiseHEFTInfinite()


    if (IamRoot) then
        allocate( lowBound(n_ch), midBound(n_ch) )
        lowBound(:) = m_mes(:) + m_bar(:)
        absErr = 1.0E-6_DP
        relErr = 0.0_DP

        ! write(*,*) '--------------------TEST-------------------'
        ! ich = 1
        ! ibare = 1
        ! i_bare_pole = 1
        ! E_gl = 1.6_DP
        ! E_pole2(:) = (/ 1.2, -0.05 /)
        ! E_pole = cmplx(E_pole2(1), E_pole2(2))
        ! call calcSMatrixPole(2, E_pole2, poleInverse)
        ! write(*,*) poleInverse
        ! stop

        open(142, file='data/poleGridSearch.out', action='write')

        printFittingOutput = .false.
        do zk = 1, max(n_bare, 1)
           i_bare_pole = zk

           do zi = 1,5
              do zj = 1,5
                 E_pole_guess = cmplx(1.1_DP + (zi-1)*0.3_DP &
                     & , -0.01 - (zj-1)*0.02,DP)
                 E_pole_guess = cmplx(1.1_DP + (zi-1)*0.3_DP &
                     & , -0.05_DP * 2.0_DP**(zj-1.0_DP))

                 ! do zi = 1,3
                 !    do zj = 1,3
                 !       E_pole_guess = cmplx(1.2_DP + (zi-1)*0.2_DP &
                 !           & , -0.05_DP * 2.0_DP**(zj-1.0_DP))

                 ! do zi = 1, 2
                 !    do zj = 1, 1
                 !       if (zi.eq.1) then
                 !           ! E_pole_guess = cmplx(1.65065073840254_DP,-0.122720220940563_DP, DP)
                 !           E_pole_guess = cmplx(1.657_DP, -0.0555_DP, DP)
                 !           ! E_pole_guess = cmplx(1.500_DP, -0.050_DP, DP)
                 !       else
                 !           E_pole_guess = cmplx(1.658_DP, -0.056_DP, DP)
                 !       end if

                 E_pole2(:) = (/ real(E_pole_guess, DP), aimag(E_pole_guess) /)
                 write(*,*) 'Pole number:  ', zk, zi, zj
                 write(*,'(a,f8.6,a,f8.6,a)') 'Guess: ', E_pole2(1) &
                     & , ' - ', abs(E_pole2(2)), 'i'


                 ! maxIterations = 1000
                 ! call Minimize( calcSMatrixPole, 2, E_pole2, searchError &
                 !     & , searchMaxStepFactor, poleInverse, printFittingOutput &
                 !     & , maxIterations )
                 ! call Minimize( calcSMatrixPole, 2, E_pole2, searchError &
                 !     & , searchMaxStepFactor, poleInverse, printFittingOutput &
                 !     & , maxIterations )


                 ! initialise BOBYQA parameters
                 bq_npt = 2*2 + 1
                 bq_rhobeg = 1.0d-4
                 bq_rhoend = 1.0d-7
                 bq_iprint = 0
                 bq_maxfun = 100000
                 pole_low(:) = [0.1, -1.0]
                 pole_high(:) = [3.0, 0.01]

                 call bobyqa( 2, bq_npt, E_pole2 &
                     & , pole_low, pole_high, bq_rhobeg, bq_rhoend, bq_iprint &
                     & , bq_maxfun, calcSMatrixPole_bq )
                 poleInverse = chi2_pole_min

                 E_poles(zk) = cmplx(E_pole2(1), E_pole2(2))

                 ! write(*,*)
                 ! write(*,*) 'Search finished'
                 ! write(*,*) 'Pole Loc: ', E_pole_min
                 ! write(*,*) 'T-inverse: ', chi2_pole_min
                 ! stop

                 ! write(*,*) 'Pole loc:     ', E_poles(zk)
                 write(*,'(a,f8.6,a,f8.6,a)') 'Pole loc:     ', real(E_poles(zk)), &
                     & ' - ', abs(aimag(E_poles(zk))), 'i'
                 write(*,'(a,es10.3)') 'Pole inverse: ', poleInverse
                 write(142,'(3i4,3f19.14)') zk, zj, zi, E_pole2(1), E_pole2(2), poleInverse
              end do
           end do
        end do

        close(142)


        ! deallocate( E_arr, phaseShift, eta_arr, crossSec )
        deallocate( lowBound, midBound )
        call finaliseHEFT()
        close(101)
        close(102)
    end if
    !-----------------------------------------------------------------------!

end program poleSearch
