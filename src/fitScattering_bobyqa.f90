program fitInfiniteVol
    use kinds
    use numFort
    use heft
    use SMatrix
    use bobyqa_module
    use stdlib_string_type
    use stdlib_strings
    implicit none

    ! Data input
    integer :: nPointsData
    real(DP), dimension(:), allocatable :: E_data
    real(DP), dimension(:), allocatable :: etaData
    real(DP), dimension(:), allocatable :: etaDataErr
    real(DP), dimension(:), allocatable :: SrData
    real(DP), dimension(:), allocatable :: SrDataErr
    real(DP), dimension(:), allocatable :: phaseShiftData
    real(DP), dimension(:), allocatable :: phaseShiftDataErr

    ! Scattering output
    integer :: nPoints = 33
    ! real(DP) :: E_init = 1.11_DP, E_fin = 1.8_DP
    real(DP), dimension(:), allocatable :: E_arr
    real(DP), dimension(:), allocatable :: eta_arr
    real(DP), dimension(:), allocatable :: Sr_arr
    real(DP), dimension(:,:), allocatable :: phaseShift
    real(DP), dimension(:,:,:), allocatable :: crossSec

    ! Misc
    integer :: i, j, iii, jjj ! misc loops
    integer :: fitCounter = 0
    integer :: iMin
    integer :: maxIterations
    integer :: iRead = 42
    integer :: nMinimise = 5
    character(len=1) :: inputAppendFits
    character(len=32) :: allFitsFmt
    logical :: doInputAppendFits
    logical :: inputAppendExists
    logical :: enableFitting = .true. ! Disable this if testing non-fitting stuff

    ! Fitting stuff
    real(DP), parameter :: paramSearchMaxStepFactor = 1e3_DP
    real(DP) :: chi2, chi2_old
    real(DP) :: fit_t_start, fit_t_end, fit_t_mid
    logical :: limitFit = .false.
    logical :: fitPseudoData
    real(DP) :: pseudo_uncertainty
    real(DP) :: SrDataMin = 0.0005_DP
    type(string_type) :: fileName_scatData
    logical :: useRandomInit = .false.

    ! Fitting mins
    real(DP) :: chi2_min = HUGE(1.0_DP)
    real(DP), dimension(:), allocatable :: x_min

    ! bobyqa stuff (bq) (turn off limitFit)
    integer :: bq_npt
    real(DP) :: bq_rhobeg, bq_rhoend
    integer :: bq_iprint
    integer :: bq_maxfun

    call random_seed()
    call initialiseHEFT()
    call initialiseHEFTInfinite()
    call initialiseHEFTFitting()
    if (useRandomInit .and. IamRoot) then
        call scrambleActiveParams(activeFitParams, activeFitParamLow, activeFitParamHigh)
        write(*,*) 'Fit params randomised within bounds'
    end if
    call printCurrentParameters(iParamChoice)

    if (IamRoot) then
        if (verbose) write(*,*) 'Initialisations complete'
        allocate( lowBound(n_ch), midBound(n_ch) )
        allocate( x_min(nParams) )

        ! Stuff for phase continuity
        allocate( isClockwise(n_ch), openChPoint(n_ch), isOpenCh(n_ch) )
        allocate( old_SMat(n_ch,n_ch) )
        isOpenCh(:) = .false.

        ! -----------------------Read data for fitting----------------------
        fitPseudoData = .true.

        if (fitPseudoData) then
            fileName_scatData = 'dataInf_pseudo.in'
        else
            fileName_scatData = 'dataInf.in'
        end if
        write(*,*) 'Fitting to ' // fileName_scatData
        write(*,*)

        if (verbose) write(*,*) 'Reading ' // fileName_scatData
        open(104, file=char(fileName_scatData), action='read')
        if (verbose) write(*,'(a,i0)') ' Opened '//fileName_scatData//' as ID: ', 104
        read(104,*) nPoints
        ! Allocate the data arrays
        allocate( E_data(nPoints), phaseShiftData(nPoints) &
            & , phaseShiftDataErr(nPoints) &
            & , SrData(nPoints), SrDataErr(nPoints) &
            & , etaData(nPoints), etaDataErr(nPoints) )
        if (verbose) write(*,*) 'Data arrays allocated'
        do iii = 1,nPoints
           read(104,*) E_data(iii), phaseShiftData(iii) &
               & , phaseShiftDataErr(iii), SrData(iii), SrDataErr(iii)
        end do

        E_init = E_data(1)
        E_fin = E_data(nPoints)
        if (verbose) write(*,'(x,a,f6.3,a,f6.3)') 'Data read, E_init=', E_init &
            & , ', E_fin=', E_fin

        ! Convert from Sr to eta
        etaData(:) = sqrt(1.0_DP - SrData(:))

        if (fitPseudoData) then
            ! If using pseudo-data, set uncertainty to some percentage
            pseudo_uncertainty = 0.05_DP
            phaseShiftDataErr(:) = phaseShiftData(:) * pseudo_uncertainty
            SrDataErr(:) = SrData(:) * pseudo_uncertainty
            where (SrDataErr(:).lt.SrDataMin) SrDataErr(:) = SrDataMin
            etaDataErr(:) = etaData(:) * pseudo_uncertainty

            ! etaDataErr(:) = 0.5_DP / sqrt(1.0_DP - SrData(:)) * SrDataErr(:)
        else
            ! If using single-energy data, calculate error from Sr
            etaDataErr(:) = 0.5_DP / sqrt(1.0_DP - SrData(:)) * SrDataErr(:)
        end if

        close(104)
        ! ------------------------------------------------------------------


        ! -----------------------Initialise some stuff----------------------
        allocate( E_arr(nPoints), eta_arr(nPoints), Sr_arr(nPoints) &
            & , phaseShift(nPoints,n_ch), crossSec(nPoints,n_ch,n_ch) )
        if (verbose) write(*,*) 'Scattering observables allocated'
        if (E_data(1).gt.1000) E_data(:) = E_data(:)/1000.0_DP
        E_arr(:) = E_data(:)
        lowBound(:) = m_mes(:) + m_bar(:)
        absErr = 1.0E-6_DP ! absolute error of integrals
        relErr = 0.0_DP ! relative error of integrals
        ! ------------------------------------------------------------------


        ! ------------------------------------------------------------------
        call packFitParams()
        activeFitParamGuess(:) = activeFitParams(:)

        if (verbose) write(*,*)
        if (verbose) write(*,*) 'Fit params packed as:'
        if (verbose) write(*,*) activeFitParams(:)

        call cpu_time(fit_t_start)
        if (verbose) write(*,*) 'Begining minimisation process'

        ! ------------------------------------------------------------------
        ! ----------------------Do Fitting with bobyqa----------------------
        ! ------------------------------------------------------------------
        ! npt = number of interpolation conditions
        ! This must be > nParams, typically 2*nParams+1
        bq_npt = 2*nParams + 1

        ! rhobeg is the initial trust region
        ! Should be set to about 1/10 the greatest
        ! expected change to a variable
        bq_rhobeg = 1.0d-4

        ! bobyqa then shrinks the trust region periodically
        ! until it reaches this final value
        bq_rhoend = 1.0d-8

        ! How much printing to do
        ! 0 : no printing
        ! 1 : only print final results
        ! 2 : print updates as it runs
        bq_iprint = 1

        ! Maximum number of function evaluations
        bq_maxfun = 100000

        ! Run the bobyqa algorithm
        call bobyqa(nParams, bq_npt, activeFitParams &
            & , activeFitParamLow, activeFitParamHigh &
            & , bq_rhobeg, bq_rhoend, bq_iprint &
            & , bq_maxfun, parameterFit)

        ! calc final chi2
        chi2 =  sum( (phaseShift(:,ch_onshell) - phaseShiftData(:))**2 &
            & / phaseShiftDataErr(:)**2 )
        if (doEtaFitting .and. n_ch.gt.1) then
            chi2 = chi2 + sum( (eta_arr(:) - etaData(:))**2 &
                & / etaDataErr(:)**2 )
        end if

        ! printCurrentParameters prints m_bare_phys so gotta update it
        m_bare_phys = m_bare
        write(*,*)
        call printCurrentParameters()

        write(*,'(a,f0.2)') ' chi2 = ', chi2
        write(*,'(a,f0.2)') ' chi2/dof = ', chi2 / (size(phaseShiftData(:)) - nParams)
        write(*,*)
        call cpu_time(fit_t_end)
        ! ------------------------------------------------------------------


        ! ------------------------Write fits to file------------------------
        open(108, file='singleFit.params', action='write')
        write(108,*) 'nParam', nParamsTotal
        write(108,'(30f14.8)') fitParams
        close(108)
        ! ------------------------------------------------------------------


        ! -------------------------Print Data vs Fit------------------------
        write(*,'(17x,a)') 'Phase Shifts'
        write(*,'(15x,a4,6x,a4)') 'Data', ' Fit'
        do iii = 1, nPoints, nPoints/10
           write(*,'(13x,f6.2,5x,f6.2)') phaseShiftData(iii) &
               & , phaseShift(iii,ch_onshell)
        end do

        if (n_ch .gt. 1) then
            write(*,*)
            write(*,'(17x,a)') 'Inelasticity'
            write(*,'(15x,a4,6x,a4)') 'Data', ' Fit'
            do iii = 1, nPoints, nPoints/10
               write(*,'(13x,f6.2,5x,f6.2)') etaData(iii) &
                   & , eta_arr(iii)
            end do
        end if

        write(*,*)
        write(*,*) '                t = ' &
            & , timePrint(fit_t_end-fit_t_start)
        ! fitParams(fitParamIndex) = activeFitParams(:)
        fitParams(fitParamIndex) = x_min(:)
        call unpackFitParams(fitParams)
        ! ------------------------------------------------------------------

        ! -----------Potentially append the fit to the allFits file----------
        write(*,'(/,a)') ' Do you wish to append the fit to allFits.params? (y/n): '
        call request_y_or_n(doInputAppendFits)
        ! doInputAppendFits = .false.

        if (doInputAppendFits) then
            ! if (inputAppendFits.eq.'y') then
            !    check if the allFits file exists
            inquire(file='allFits.params', exist=inputAppendExists)
            open(108, file='allFits.params')

            if (inputAppendExists) then
                ! If it does exist, go to the end of the file
                !    to find the new fit number (iRead+1)
                read(108,*)
                read(108,*)
                read(108,*)
                read(108,*)
                do
                   read(108,*) iRead
                   if (iRead.eq.0) exit
                end do
                backspace(108)
                backspace(108)
                read(108,*) iRead
            else
                ! If it does not, set it up with the nParams at the top
                write(108,*) 'nParams', nParamsTotal
                write(108,*) 'iChoice', 1
                iRead = 0
            end if

            write(108,'(i5,'//trim(int2str(nParamsTotal)) &
                & //'f14.8, f14.2, 6x, a6, i0, a12, i0, a9)') &
                & iRead+1, fitParams, chi2, pot_choices//' | ' &
                & , nParams, '-param | Fit', iParamChoice, ' | bobyqa'
            write(108,'(i5,'//trim(int2str(nParamsTotal))//'f14.8, f14.2)') &
                & 0, fitParams*0.0_DP, 0.0_DP
            close(108)
        end if
        ! ------------------------------------------------------------------


        deallocate( x_min )
        deallocate( E_arr, phaseShift, eta_arr, Sr_arr, crossSec )
        deallocate( lowBound, midBound )
        deallocate( isClockwise, openChPoint, isOpenCh )
        deallocate( E_data, phaseShiftData, phaseShiftDataErr &
            & , etaData, etaDataErr )
        close(101)
        close(102)


    end if
    call finaliseHEFT()

    !-----------------------------------------------------------------------!

contains

    subroutine request_y_or_n(isYes)
        ! Subroutine to force a yes/no call and store the result in a logical
        implicit none
        logical, intent(out) :: isYes
        character(len=4) :: y_n
        do
           read(*,*) y_n
           if (y_n.eq.'y' .or. y_n.eq.'yes' .or. y_n.eq.'Y') then
               isYes = .true.
               exit
           else if (y_n.eq.'n' .or. y_n.eq.'no' .or. y_n.eq.'N') then
               isYes = .false.
               exit
           else
               write(*,'(a,\)') '  Please enter "y" or "n": '
           end if
        end do
    end subroutine request_y_or_n



    function colour_array(arr, colour, fmt, mask_in) result(out)
        implicit none
        real(DP), dimension(:), intent(in) :: arr
        character(len=*), intent(in) :: colour
        character(len=*), intent(in) :: fmt
        logical, dimension(size(arr)), intent(in), optional :: mask_in

        logical :: isMasked
        type(string_type) :: out
        character(len=7) :: col, reset
        integer :: i_arr

        ! Check if there is a mask
        isMasked = present(mask_in)

        ! Define the colours
        reset = achar(27)//'[0;00m'

        ! 1 bold, 2 dim, 4 underlines, 5 blink, 7 invert, 8 hidden
        select case (colour)
        case ('red')
            col = achar(27)//'[0;31m'
        case ('green')
            col = achar(27)//'[0;32m'
        case ('blue')
            col = achar(27)//'[0;34m'
        case ('cyan')
            col = achar(27)//'[0;36m'
        case ('magenta')
            col = achar(27)//'[0;35m'
        case ('yellow')
            col = achar(27)//'[0;33m'
        case ('black')
            col = achar(27)//'[0;30m'
        case ('grey')
            col = achar(27)//'[0;90m'
        case default
            col = reset
        end select

        ! Create the output string according to mask
        out = ''
        do i_arr = 1, size(arr)
           if (isMasked) then
               if (mask_in(i_arr)) then
                   out = out // col // to_string(arr(i_arr),fmt) &
                       & // reset
               else
                   out = out // to_string(arr(i_arr),fmt)
               end if

           else
               ! out = out // to_string(arr(i_arr),fmt)
               out = out // col // to_string(arr(i_arr),fmt) &
                   & // reset
           end if
        end do
    end function colour_array



    subroutine parameterFit(n, x, f)
        implicit none
        integer, intent(in) :: n
        real(DP), dimension(:), intent(in) :: x
        real(DP), intent(out) :: f ! chi^2
        ! Local variables
        logical :: fitVerb = .false.
        logical :: printMin = .true.

        fitCounter = fitCounter + 1
        if (fitVerb) write(*,*) 'x: ', x

        fitParams(fitParamIndex) = x(:)
        call unpackFitParams(fitParams)

        do i = 1, nPoints
           globalPoint = i
           call calcSMatrix(E_arr(i), phaseShift(i,:), eta_arr(i))
           if (fitVerb) write(*,*) i, E_arr(i), phaseShift(i,1)

           ! If the diff between this point and the previous points is more than 150 deg
           !    then just add 180 deg on. atanRangeFixer doesn't seem to work for >180 :(
           if (i.gt.1) then
               if (abs(phaseShift(i,ch_onshell) - phaseShift(i-1,ch_onshell)).gt.150.0_DP) then
                   phaseShift(i,ch_onshell) = phaseShift(i,ch_onshell) + 180.0_DP
               end if
           end if
        end do

        ! calculate chi^2
        f =  sum( (phaseShift(:,ch_onshell) - phaseShiftData(:))**2 &
            & / phaseShiftDataErr(:)**2 )
        if (doEtaFitting .and. n_ch.gt.1) then
            f = f + sum( (eta_arr(:) - etaData(:))**2 &
                & / etaDataErr(:)**2 )
        end if

        if (f.lt.chi2_min) then
            ! Print the new minimum found
            chi2_min = f
            x_min(:) = x(:)
            if (printMin) then
                ! colour_array makes the variables which arent being varied grey
                write(*,*) to_string(fitCounter,'i10') // '    x_min:' &
                    & // colour_array(fitParams, 'grey', 'f10.4', .not.isActiveParam)
                write(*,'(14x,a,f0.2,a)') 'chi2 = ', f, '   |   ' &
                    & //pot_choices//'   |   '//'fit'//to_string(iParamChoice) &
                    & //'   |   '//char(fileName_scatData)//'   |   bobyqa'
                write(*,*)
            end if
        end if

    end subroutine parameterFit




end program fitInfiniteVol
