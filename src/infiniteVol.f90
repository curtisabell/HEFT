program infiniteVol
    use kinds
    use numFort
    use heft
    use SMatrix
    implicit none

    ! Data input
    integer :: nPointsData
    real(DP), dimension(:), allocatable :: E_data
    real(DP), dimension(:), allocatable :: etaData, SrData
    real(DP), dimension(:), allocatable :: etaDataErr, SrDataErr
    real(DP), dimension(:), allocatable :: phaseShiftData
    real(DP), dimension(:), allocatable :: phaseShiftDataErr

    ! Scattering output
    real(DP), dimension(:), allocatable :: E_arr
    real(DP), dimension(:), allocatable :: eta_arr, Sr_arr
    real(DP), dimension(:,:), allocatable :: phaseShift
    real(DP), dimension(:,:), allocatable :: Gamma_decay, Gamma_ratio
    real(DP), dimension(:,:,:), allocatable :: crossSec
    complex(DP), dimension(:), allocatable :: Tmatrix_onshell

    real(DP) :: E_img
    real(DP) :: eta_img, Sr_img

    integer :: i, j, iii, jjj, i_nImg, jImg ! misc loops
    real(DP) :: chi2
    real(DP), dimension(2,2) :: id2 = reshape((/ 1,0,0,1 /) , shape(id2))

    integer :: file_scattering, file_crossSec
    integer :: file_tmatrix
    character(len=128) :: fileName_scattering, fileName_crossSec, fileName_tmatrix

    logical :: usePseudoData = .false.
    real(DP) :: pseudo_uncertainty
    character(len=128) :: fileName_data

    call initialiseHEFT()
    call printCurrentParameters(iParamChoice)
    call initialiseHEFTInfinite()

    ! Set up file outputs
    if (IamRoot) then
        fileName_scattering = 'data/scattering_fit' &
            & // trim(adjustl(int2str(iParamChoice))) // '.out'
        fileName_crossSec = 'data/crossSection_fit' &
            & // trim(adjustl(int2str(iParamChoice))) // '.out'
        fileName_tmatrix = 'data/tmatrix_fit' &
            & // trim(adjustl(int2str(iParamChoice))) // '.out'

        if (usePseudoData) then
            fileName_data = 'dataInf_pseudo.in'
        else
            fileName_data = 'dataInf.in'
        end if
    end if


    if (IamRoot) then
        allocate( lowBound(n_ch) )

        ! Need these variables for phase shift continuity
        allocate( isClockwise(n_ch), openChPoint(n_ch), isOpenCh(n_ch) )
        allocate( old_SMat(n_ch,n_ch) )
        isOpenCh(:) = .false.

        ! ------------------------------------------------------------------
        ! -------------------------Read data points-------------------------
        ! ------------------------------------------------------------------
        if (useDataPoints) then
            open(142, file=fileName_data, action='read')
            read(142,*) nPoints_inf

            ! Allocate the arrays
            allocate( E_data(nPoints_inf), phaseShiftData(nPoints_inf)&
                & , phaseShiftDataErr(nPoints_inf) &
                & , etaData(nPoints_inf), etaDataErr(nPoints_inf) &
                & , SrData(nPoints_inf), SrDataErr(nPoints_inf) )

            do iii = 1,nPoints_inf
               read(142,*) E_data(iii), phaseShiftData(iii), phaseShiftDataErr(iii) &
                   & , SrData(iii), SrDataErr(iii)
            end do

            ! Convert from Sr to eta: S_r = 1 - eta**2
            etaData(:) = sqrt(1.0_DP - SrData(:))
            etaDataErr(:) = 0.5_DP / etaData(:) * SrDataErr(:)
            ! etaDataErr(:) = 0.5_DP / SrData(:) * SrDataErr(:)
            close(142)

            ! this is for when I artificially alter the inelasticity
            !    uncertainties for fittting
            if (usePseudoData) then
                pseudo_uncertainty = 0.01_DP
                phaseShiftDataErr(:) = phaseShiftData(:) * pseudo_uncertainty

                SrDataErr(:) = max(SrData(:) * pseudo_uncertainty, 0.01)
                etaDataErr(:) = etaData(:) * pseudo_uncertainty

                ! SrDataErr(:) = 2.0_DP*etaData * etaDataErr
            end if
        end if

        ! Initialisations
        allocate( E_arr(nPoints_inf) )
        allocate( eta_arr(nPoints_inf), Sr_arr(nPoints_inf) &
            & , phaseShift(nPoints_inf,n_ch) &
            & , crossSec(nPoints_inf,n_ch,n_ch) &
            & , Gamma_decay(nPoints_inf,n_ch) &
            & , Gamma_ratio(nPoints_inf,n_ch) &
            & , Tmatrix_onshell(nPoints_inf) )
        lowBound(:) = m_mes(:) + m_bar(:)
        absErr = 1.0E-7_DP
        relErr = 0.0_DP

        ! --------------------Construct the energy array--------------------
        if (useDataPoints) then
           ! Checks if the data is in MeV or GeV
           ! this wouldn't work for pi-pi scattering though since
           !    the pi-pi threshold is not >1 GeV
            if (E_data(1).gt.1000) then
                E_arr(:) = E_data(:)/1000.0_DP
            else
                E_arr(:) = E_data(:)
            end if
            E_init = E_arr(1)
            E_fin = E_arr(nPoints_inf)
        else
            E_arr(:) = (E_fin-E_init)/(nPoints_inf-1) &
                 & * [ (iii,iii=0,nPoints_inf-1) ] + E_init
        end if

        write(*,*) '=========================================================='
        write(*,*) '    E                      Phase'
        if (useDataPoints) then
            write(*,'(11x,2a14)') 'Data', 'Fit'
        end if

        ! ------------------------------------------------------------------
        ! -----------------Calculate scattering observables-----------------
        ! ------------------------------------------------------------------
        do i = 1, nPoints_inf
           globalPoint = i
           call calcSMatrix( E_arr(i), phaseShift(i,:), eta_arr(i) &
               & , crossSec(i,:,:), Gamma_decay(i,: ), Tmatrix_onshell(i) )

           ! If the diff between this phase and the previous phase is more than 150 deg
           !    then just add 180 deg on. atanRangeFixer doesn't seem to always work
           !    properly for more than 180 degrees though, should look into it
           if (i.gt.1) then
               where (abs(phaseShift(i,:) - phaseShift(i-1,:)).gt.150.0_DP)
                   phaseShift(i,:) = phaseShift(i,:) + 180.0_DP
               end where
           end if

           ! Convert decay rate to branching ratio
           Gamma_ratio(i,:) = Gamma_decay(i,:) / sum(Gamma_decay(i,:))

           ! Write data to terminal
           if (useDataPoints) then
               write(*,'(f8.4,11x,2f14.4)') E_arr(i), phaseShiftData(i), phaseShift(i,ch_onshell)
           else
               if (mod(i,nPoints_inf/10).eq.0) then
                   write(*,'(f8.3,11x,f14.4)') E_arr(i), phaseShift(i,ch_onshell)
               end if
           end if
        end do

        ! -------------------------Calculate chi**2-------------------------
        if (useDataPoints) then
            chi2 = sum( (phaseShift(:,ch_onshell) - phaseShiftData(:))**2 &
                 & / phaseShiftDataErr(:)**2 )
            ! if there is more than one channel, include inelasticity
            if (n_ch.gt.1) then
               chi2 = chi2 + sum( (eta_arr(:) - etaData(:))**2 &
                    & / etaDataErr(:)**2 )
            end if
        else
            chi2 = fit_chi2
        end if

        if (chi2.gt.10e4) then
            write(*,'(a,es9.2)') 'chi2 = ', chi2
        else
            write(*,'(a,f8.2)') 'chi2 = ', chi2

            ! TODO: store number of parameters used
            !       for each fit to be read in for this calculation
            ! write(*,'(a,f8.2)') 'chi2/dof = ', chi2/(size(phaseShiftData(:)) - 3.0_DP)
            ! write(*,'(a,f8.2)') 'chi2/dof = ', chi2/(33.0_DP-6.0_DP)
        end if
        write(*,*) '=========================================================='

        ! change the NaNs to -1s for easier processing
        where (crossSec(:,:,:).ne.crossSec(:,:,:)) crossSec(:,:,:) = -1.0_DP

        ! ------------------Write phase and cross-sec to file-----------------
        open( 101, file=fileName_scattering, action='write' )
        open( 102, file=fileName_crossSec, action='write' )
        open( 103, file=fileName_tmatrix, action='write' )

        write(101,'(8a14)') 'E (GeV)', ch_labels(:), 'Inelasticity'
        write(102,'(8a14)') 'E (GeV)', ch_labels(:), 'Cross Section'
        do i = 1, nPoints_inf
           ! write phase shift and inelasticity to file
           write(101,'(8f14.8)') E_arr(i), phaseShift(i,:), eta_arr(i)

           ! write cross section to file
           write(102,'(f14.8)', advance='no') E_arr(i)
           do ich = 1, n_ch
              if (ich.ne.1) write(102,'(f14.8)', advance='no') 0.0
              write(102,'(8f14.8)') crossSec(i,ich,:)
           end do

           ! write onshell t-matrix element to file
           write(103, '(3f12.6)') E_arr(i), real(Tmatrix_onshell(i)) &
               & , aimag(Tmatrix_onshell(i))
        end do

        close(101)
        close(102)

        deallocate( phaseShift, eta_arr, crossSec, Sr_arr )
        ! deallocate( Gamma_decay, crossSec )
        deallocate( lowBound, E_arr )
        deallocate( isClockwise, openChPoint, isOpenCh )
        if (useDataPoints) then
            deallocate( E_data, phaseShiftData, phaseShiftDataErr &
                & , etaData, etaDataErr, SrData, SrDataErr )
        end if


    end if
    call finaliseHEFT()

    !-----------------------------------------------------------------------!

end program infiniteVol
