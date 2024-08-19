module transFunctions
    use kinds
    implicit none

contains

    pure function atanTrans(x, xMin, xMax, xGuess, dim) result(y)
        implicit none
        integer, intent(in) :: dim
        real(DP), dimension(dim), intent(in) :: x, xMin, xMax, xGuess
        ! Locals
        real(DP), dimension(dim) :: y, yGuess
        real(DP), dimension(dim) :: alpha

        yGuess(:) = (xGuess(:) - xMin(:)) / (xMax(:) - xMin(:))
        alpha(:) = -log(yGuess(:)) / log(2.0_DP)
        y(:) = 1.0_DP/pi * atan(x(:)) + 0.5_DP
        y(:) = xMin(:) + (xMax(:) - xMin(:)) * y(:)**alpha(:)

    end function atanTrans

    pure function tanTrans(y, xMin, xMax, xGuess, dim) result(x)
        implicit none
        integer, intent(in) :: dim
        real(DP), dimension(dim), intent(in) :: y, xMin, xMax, xGuess
        ! Locals
        real(DP), dimension(dim) :: x, yGuess
        real(DP), dimension(dim) :: alpha

        yGuess(:) = (xGuess(:) - xMin(:)) / (xMax(:) - xMin(:))
        alpha(:) = -log(yGuess(:)) / log(2.0_DP)
        x(:) = tan(pi*((y(:) - xMin(:))/(xMax(:) &
            & - xMin(:)))**(1.0_DP/alpha(:)) &
            & - pi/2.0_DP)

    end function tanTrans

end module transFunctions


program fitInfiniteVol
    use kinds
    use numFort
    use heft
    use SMatrix
    use transFunctions
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
    logical :: limitFit = .true.
    logical :: fitPseudoData = .false.
    character(len=128) :: fileName_scatData

    ! Fitting mins
    real(DP) :: chi2_min = HUGE(1.0_DP)
    real(DP), dimension(:), allocatable :: x_min

    call initialiseHEFT()
    call printCurrentParameters(iParamChoice)
    call initialiseHEFTInfinite()
    call initialiseHEFTFitting()

    if (IamRoot) then
        if (verbose) write(*,*) 'Initialisations complete'
        allocate( lowBound(n_ch), midBound(n_ch) )
        allocate( x_min(nParams) )

        ! Stuff for phase continuity
        allocate( isClockwise(n_ch), openChPoint(n_ch), isOpenCh(n_ch) )
        allocate( old_SMat(n_ch,n_ch) )
        isOpenCh(:) = .false.

        ! -----------------------Read data for fitting----------------------
        if (fitPseudoData) then
            fileName_scatData = 'dataInf_pseudo.in'
        else
            fileName_scatData = 'dataInf.in'
        end if

        if (verbose) write(*,*) 'Reading ' // trim(fileName_scatData)
        open(104, file=trim(fileName_scatData), action='read')
        if (verbose) write(*,'(a,i0)') ' Opened dataInf.in as ID: ', 104
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
        etaDataErr(:) = 0.5_DP / sqrt(1.0_DP - SrData(:)) * SrDataErr(:)
        close(104)
        ! ------------------------------------------------------------------


        ! -----------------------Initialise some stuff----------------------
        allocate( E_arr(nPoints), eta_arr(nPoints), Sr_arr(nPoints) &
            & , phaseShift(nPoints,n_ch), crossSec(nPoints,n_ch,n_ch) )
        if (verbose) write(*,*) 'Scattering observables allocated'
        if (E_data(1).gt.1000) E_data(:) = E_data(:)/1000.0_DP
        E_arr(:) = E_data(:)
        lowBound(:) = m_mes(:) + m_bar(:)
        absErr = 1.0E-6_DP
        relErr = 0.0_DP
        ! ------------------------------------------------------------------


        ! --------------------------- ---------------------------------------
        ! -----------------------Do parameter fitting-----------------------
        ! ------------------------------------------------------------------
        call packFitParams()
        activeFitParamGuess(:) = activeFitParams(:)

        if (verbose) write(*,*)
        if (verbose) write(*,*) 'Fit params packed as:'
        if (verbose) write(*,*) activeFitParams(:)

        limitFit=.true.
        if (limitFit) then
            if (verbose) write(*,*) 'Fit limited to bounds using atan transformations'
            activeFitParams(:) = 0.0_DP
        end if

        ! ----------------------------Do Fitting----------------------------
        printFittingOutput = .false.
        call cpu_time(fit_t_start)
        if (verbose) write(*,*) 'Begining minimisation process'

        ! maxIterations = 100
        ! if (verbose) write(*,*) 'Max iterations:', maxIterations
        ! call Minimize( parameterFit, nParams, activeFitParams &
        !     & , activeFitParamErrors, paramSearchMaxStepFactor &
        !     & , chi2, printFittingOutput, maxIterations )

        ! call Minimize( parameterFit, nParams, activeFitParams &
        !     & , activeFitParamErrors, paramSearchMaxStepFactor &
        !     & , chi2, printFittingOutput )

        ! Finding that minf often stops early, so keep running it
        ! until the chi2 isnt being updated much anymore
        ! bit hacky but oh well
        if (enableFitting) then
            iMin = 1
            do
               maxIterations = iMin*500
               call Minimize( parameterFit, nParams, activeFitParams &
                   & , activeFitParamErrors, paramSearchMaxStepFactor &
                   & , chi2, printFittingOutput, maxIterations )
               fitCounter = 0

               write(*,*)
               if (iMin.gt.1) then
                   if (chi2/chi2_old.gt.0.99999_DP .or. iMin.eq.20) then
                       write(*,*) chi2/chi2_old
                       ! printCurrentParameters prints mbarephys for reasons
                       m_bare_phys = m_bare
                       call printCurrentParameters()
                       exit
                   end if
               end if
               chi2_old = chi2
               iMin = iMin + 1
            end do
        end if

        ! Transform the fit variables back to parameter space if needed
        if (limitFit) then
            activeFitParams(:) = atanTrans( activeFitParams, activeFitParamLow &
                & , activeFitParamHigh, activeFitParamGuess, nParams )
        end if
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
                & //'f14.8, f14.2, 6x, a3, a12, i0)') &
                & iRead+1, fitParams, chi2, pot_choices &
                & , '    From Fit', iParamChoice
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
           if (y_n.eq.'y' .or. y_n.eq.'yes') then
               isYes = .true.
               exit

           else if (y_n.eq.'n' .or. y_n.eq.'no') then
               isYes = .false.
               exit

           else
               write(*,'(a,\)') '  Please enter "y" or "n": '
           end if
        end do

    end subroutine request_y_or_n



    subroutine parameterFit(n, x, f)
        implicit none
        integer, intent(in) :: n
        real(DP), dimension(n), intent(in) :: x
        real(DP), intent(out) :: f ! chi^2
        ! Local variables
        real(DP), dimension(n) :: xTrans
        logical :: fitVerb = .false.

        fitCounter = fitCounter + 1

        ! Transform x in (-inf,inf) to (paramLow,paramHigh)
        if (limitFit) then
            xTrans(:) = atanTrans( x(:), activeFitParamLow &
                & , activeFitParamHigh, activeFitParamGuess, n )
        else
            xTrans(:) = x(:)
        end if

        if (fitVerb) write(*,*) 'x: ', x
        if (fitVerb) write(*,*) 'xTrans: ', xTrans

        fitParams(fitParamIndex) = xTrans(:)
        call unpackFitParams(fitParams)

        do i = 1, nPoints
           globalPoint = i
           call calcSMatrix(E_arr(i), phaseShift(i,:), eta_arr(i))
           if (fitVerb) write(*,*) i, E_arr(i), phaseShift(i,1)
        end do

        ! calculate chi^2
        f =  sum( (phaseShift(:,ch_onshell) - phaseShiftData(:))**2 &
            & / phaseShiftDataErr(:)**2 )
        if (doEtaFitting .and. n_ch.gt.1) then
            Sr_arr(:) = 1.0_DP - eta_arr(:)**2
            f = f + sum( (Sr_arr(:) - SrData(:))**2 &
                & / SrDataErr(:)**2 )
        end if

        if (f.lt.chi2_min) then
            chi2_min = f
            x_min(:) = xTrans(:)
            ! colour_array makes the variables which arent being varied grey
            write(*,'(a)') trim(int2str(fitCounter)) // '    x_min:  ' &
                & // colour_array(fitParams, 'grey', '(f10.4)', .not.isActiveParam)
            write(*,'(14x,a,f0.2,a)') 'chi2 = ', f, '   |   ' &
                & //pot_choices//'   |   '//'fit'//trim(int2str(iParamChoice)) &
                & //'   |   '//trim(fileName_scatData)//'   |   minf'
            write(*,*)
        end if

    end subroutine parameterFit



    function colour_array(arr, colour, fmt, mask_in) result(out)
        ! Converts an array of real numbers into a coloured string
        implicit none
        real(DP), dimension(:), intent(in) :: arr
        character(len=*), intent(in) :: colour
        character(len=*), intent(in) :: fmt
        logical, dimension(size(arr)), intent(in), optional :: mask_in

        logical :: isMasked
        character(len=:), allocatable :: out
        character(len=7) :: col, reset
        character(len=32) :: numString
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
               write(numString,fmt) arr(i_arr)
               if (mask_in(i_arr)) then
                   out = out // col // trim(adjustl(numString))  // '  ' // reset
               else
                   out = out // reset // trim(adjustl(numString)) // '  '
               end if

           else
               out = out // col // trim(adjustl(numString)) // '  ' // reset
           end if
        end do
    end function colour_array


end program fitInfiniteVol
