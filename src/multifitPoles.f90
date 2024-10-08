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
    integer :: i, j, ii, jj, iii, jjj ! misc loops
    character (len=8) :: dummyString
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

    integer :: io_error
    logical :: io_stop = .false.
    logical :: paramsFileExists
    integer :: nParamFits ! how many fits to loop over
    integer :: iFit, iFit_nImg
    real(DP), dimension(:), allocatable :: thisParam
    real(DP), dimension(:,:), allocatable :: paramArray
    character(len=128) :: multifit_params_file, multifit_params_file_full

    integer :: nPolesMultiFit
    real(DP), dimension(2) :: pole_read
    complex(DP), dimension(:), allocatable :: initialPoleGuess
    complex(DP), dimension(:), allocatable :: multifitPoles
    real(DP), dimension(:), allocatable :: multifitInverse
    character(len=128) :: fileName_initialPoles
    logical :: initialPolesFileExists = .false.

    call initialiseHEFT()
    call printCurrentParameters(iParamChoice)
    call initialiseHEFTInfinite()


    ! ----------------------Read in fits from file----------------------
    if (IamRoot) then
        call get_command_argument(1, multifit_params_file)
        ! check if there is an argument
        !
        if (multifit_params_file.eq.'') then
            write(*,*) 'Enter a param file as an argument'
            write(*,*) 'Stopping...'
            io_stop = .true.

        else
            ! if there is an argument, check it is valid
            !
            multifit_params_file_full = trim(multifit_params_file) // '.params'
            inquire(file=trim(multifit_params_file_full), exist=paramsFileExists)
            if (paramsFileExists) then
                ! invalid if called one of these
                !
                if (trim(multifit_params_file).eq.'allFits') then
                    write(*,*) 'allFits is not a valid input'
                    write(*,*) 'Stopping...'
                    io_stop = .true.
                else if (trim(multifit_params_file).eq.'singleFit') then
                    write(*,*) 'singleFit is not a valid input'
                    write(*,*) 'Stopping...'
                    io_stop = .true.
                end if
            else
                ! Invalid if the file doesn't even exist
                !
                write(*,*) 'That .params file does not exist'
                write(*,*) 'Stopping...'
                io_stop = .true.
            end if
        end if

        ! Check pole initial guesses file exists for this fit
        !
        fileName_initialPoles = 'data/poleSearch_fit' &
            & // trim(adjustl(int2str(iParamChoice)))//'.out'
        inquire(file=fileName_initialPoles, exist=initialPolesFileExists)
        if (.not. initialPolesFileExists) then
            write(*,*) 'File: ' // fileName_initialPoles &
                & // ' does not exist.'
            write(*,*) 'Please run poles.x for this fit number first.'
            io_stop = .true.
        end if

        ! Stop the program if no params file was given as an argument
        !
        if (io_stop) stop

        ! Read params from the file now we know its a valid input
        !
        write(*,*) 'Params file:  ' // trim(multifit_params_file_full)
        open(192, file=trim(multifit_params_file_full), action='read')
        read(192,*)
        read(192,*) dummyString, nParamFits

        allocate(paramArray(nParamFits,nParamsTotal))
        do ii = 1, nParamFits
           read(192,*) paramArray(ii,:)
        end do
        close(192)
    end if



    if (IamRoot) then
        ! integration settings
        allocate( lowBound(n_ch), midBound(n_ch) )
        lowBound(:) = m_mes(:) + m_bar(:)
        absErr = 1.0E-6_DP
        relErr = 0.0_DP

        ! minimisation settings
        bq_npt = 2*2 + 1 ! 2 because two params (real, imag)
        bq_rhobeg = 1.0d-4
        bq_rhoend = 1.0d-7
        bq_iprint = 0
        bq_maxfun = 100000
        pole_low(:) = [0.1, -1.0]
        pole_high(:) = [3.0, 0.001]

        ! Read in pole initial guesses from file
        open(143, file=fileName_initialPoles, action='read')
        read(143,*) dummyStr, nPolesMultiFit
        allocate(initialPoleGuess(nPolesMultiFit) &
            & , multifitPoles(nPolesMultiFit) &
            & , multifitInverse(nPolesMultiFit))

        do i = 1, nPolesMultiFit
           read(143,*) pole_read
           initialPoleGuess(i) = cmplx(pole_read(1), pole_read(2), DP)
        end do
        close(143)

        allocate(thisParam(nParamsTotal))
        open(142, file='data/poles_multifit_' &
            & //trim(multifit_params_file)//'.out', action='write')

        ! nPolesMultiFit = 4


        ! initialPoleGuess(1) = cmplx(1.21_DP, -0.04_DP, DP)
        ! initialPoleGuess(2) = cmplx(1.6_DP, -0.01_DP, DP)

        ! 1.20000000006397   0.00000000005408
        ! 1.21061604312481  -0.04906418450419

        ! initialPoleGuess(1) = cmplx(1.211_DP, -0.049_DP, DP)
        ! initialPoleGuess(2) = cmplx(1.2_DP, -0.01_DP, DP)


        ! ! Delta 1b2c
        ! initialPoleGuess(1) = cmplx(1.21_DP, -0.049_DP, DP)
        ! initialPoleGuess(2) = cmplx(1.43_DP, -0.21_DP, DP)
        ! initialPoleGuess(3) = cmplx(1.328_DP, -0.412_DP, DP)

        ! initialPoleGuess(1) = cmplx(1.208773_DP, -0.049467_DP, DP)
        ! initialPoleGuess(2) = cmplx(1.474522_DP, -0.140634_DP, DP)
        ! initialPoleGuess(3) = cmplx(1.321505_DP, -0.364162_DP, DP)



        ! N 2b3c
        ! initialPoleGuess(1) = cmplx(1.500_DP, -0.050_DP, DP)
        ! initialPoleGuess(2) = cmplx(1.658_DP, -0.056_DP, DP)
        ! initialPoleGuess(3) = cmplx(1.750_DP, -0.330_DP, DP)

        ! Delta 2b3c
        ! fit 4
        ! initialPoleGuess(1) = cmplx(1.207440_DP, -0.048613_DP, DP)
        ! initialPoleGuess(2) = cmplx(1.486456_DP, -0.131985_DP, DP)
        ! initialPoleGuess(3) = cmplx(1.900094_DP, -0.107395_DP, DP)
        ! initialPoleGuess(4) = cmplx(2.238719_DP, -0.239501_DP, DP)
        ! initialPoleGuess(5) = cmplx(1.327214_DP, -0.387872_DP, DP)
        ! initialPoleGuess(6) = cmplx(1.329421_DP, -0.433259_DP, DP)
        ! 1.025360 - 0.702534i
        ! 1.650208 - 0.620415i

        ! initialPoleGuess(1) = cmplx(1.208460_DP, -0.042714_DP, DP)
        ! initialPoleGuess(2) = cmplx(1.481200_DP, -0.091565_DP, DP)
        ! initialPoleGuess(3) = cmplx(2.302976_DP, -0.072803_DP, DP)
        ! initialPoleGuess(4) = cmplx(1.604583_DP, -0.236970_DP, DP)

        ! fit 9
        ! initialPoleGuess(1) = cmplx(1.210374_DP, -0.049242_DP, DP)
        ! initialPoleGuess(2) = cmplx(1.375035_DP, -0.216378_DP, DP)
        ! initialPoleGuess(3) = cmplx(1.880061_DP, -0.218347_DP, DP)
        ! initialPoleGuess(4) = cmplx(2.064616_DP, -0.341609_DP, DP)
        ! initialPoleGuess(5) = cmplx(1.299029_DP, -0.438981_DP, DP)
        ! initialPoleGuess(6) = cmplx(2.539204_DP, -0.618716_DP, DP)

        ! if (n_bare.gt.1) then
        !     initialPoleGuess(2) = cmplx(m_bare(2), -0.01_DP, DP)
        ! else
        !     initialPoleGuess(2) = cmplx(1.6, -0.01_DP, DP)
        ! end if
        multiFitPoles(:) = initialPoleGuess(:)

        printFittingOutput = .false.

        fitLoop : do iFit = 1, nParamFits

           thisParam(:) = paramArray(iFit,:)
           call unpackFitParams(thisParam)

           write(*,*) 'iFit: ', iFit

           do zk = 1, nPolesMultiFit
              i_bare_pole = zk
              E_pole2(:) = [ real(multiFitPoles(zk),DP), aimag(multiFitPoles(zk)) ]

              write(*,'(a,f8.6,a,f8.6,a)') '   Guess: ', E_pole2(1) &
                  & , ' - ', abs(E_pole2(2)), 'i'

              ! pole_low = E_pole2 - abs(E_pole2*0.1_DP)
              ! pole_high = E_pole2 + abs(E_pole2*0.1_DP)

              pole_low = E_pole2 - cmplx(0.01_DP, 0.005_DP, DP)
              pole_high = E_pole2 + cmplx(0.01_DP, 0.005_DP, DP)

              call bobyqa( 2, bq_npt, E_pole2 &
                  & , pole_low, pole_high, bq_rhobeg, bq_rhoend, bq_iprint &
                  & , bq_maxfun, calcSMatrixPole_bq )

              multifitInverse(zk) = chi2_pole_min
              multifitPoles(zk) = cmplx(E_pole2(1), E_pole2(2))

              ! write(*,*)
              ! write(*,*) 'Search finished'
              ! write(*,*) 'Pole Loc: ', E_pole_min
              ! write(*,*) 'T-inverse: ', chi2_pole_min
              ! stop

              if (abs(E_pole2(2)).lt.1.0d-4) E_pole2(2) = 0.0_DP

              ! write(*,*) 'Pole loc:     ', E_poles(zk)
              write(*,'(a,f8.6,a,f8.6,a,5x,es19.6)') '   Pole loc:     ', E_pole2(1), &
                  & ' - ', abs(E_pole2(2)), 'i', chi2_pole_min
              ! write(*,'(a,es10.3)') 'Pole inverse: ', poleInverse
              write(142,'(i6,3f19.14)') iFit &
                  & , E_pole2(1), E_pole2(2), chi2_pole_min
           end do
        end do fitLoop

        close(142)


        ! deallocate( E_arr, phaseShift, eta_arr, crossSec )
        deallocate( lowBound, midBound )
        call finaliseHEFT()
        close(101)
        close(102)
    end if
    !-----------------------------------------------------------------------!

end program poleSearch
