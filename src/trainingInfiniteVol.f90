program trainingInfiniteVol
  use kinds
  use numFort
  use heft
  use SMatrix
  use co_random_numbers
  implicit none

  ! Data input
  integer :: nPointsData
  real(DP), dimension(:), allocatable :: E_data
  real(DP), dimension(:), allocatable :: etaData
  real(DP), dimension(:), allocatable :: etaDataErr
  real(DP), dimension(:), allocatable :: phaseShiftData
  real(DP), dimension(:), allocatable :: phaseShiftDataErr

  ! Scattering output
  ! integer :: nPoints = 121
  ! real(DP) :: E_init = 1.11_DP, E_fin = 1.8_DP
  real(DP), dimension(:), allocatable :: E_arr
  real(DP), dimension(:), allocatable :: eta_arr, Sr_arr
  real(DP), dimension(:,:), allocatable :: phaseShift
  real(DP), dimension(:,:,:), allocatable :: crossSec

  real(DP) :: E_img
  real(DP) :: eta_img, Sr_img

  ! Misc
  integer :: i, j, iii, jjj, i_nImg, jImg ! misc loops
  integer :: nProgressPrints = 10
  real(DP) :: chi2[*]
  real(DP), dimension(2,2) :: id2 = reshape((/ 1,0,0,1 /) , shape(id2))
  character(len=16) :: stringDump
  logical :: verboseTerminal = .true.
  real(DP) :: t_init, t_final

  ! training stuff
  integer :: n_MC_total, n_MC
  integer :: n_MC_img
  integer :: n_MC_run, n_runs
  integer :: iMC, iRun
  integer :: MC_counter = 0
  integer :: file_training_config = 127
  integer :: randomSeed = 42
  real(DP), dimension(2) :: m_bounds, g_bounds, Lam_bounds
  real(DP), dimension(2) :: v_bounds, Lam_v_bounds
  real(DP), dimension(:), allocatable :: rand_params[:]

  call initialiseHEFT()
  call initialiseHEFTInfinite()

  call co_random_seed(randomSeed, iImg)
  allocate(rand_params(nParamsTotal)[*])

  ! ------------Read in config file for training generation-----------
  if (IamRoot) then
     open(file_training_config, file='HEFTTraining.config', action='read')
     read(file_training_config,*)
     read(file_training_config,*) stringDump, m_bounds
     read(file_training_config,*) stringDump, g_bounds
     read(file_training_config,*) stringDump, Lam_bounds
     read(file_training_config,*) stringDump, v_bounds
     read(file_training_config,*) stringDump, Lam_v_bounds
     read(file_training_config,*)

     read(file_training_config,*)
     read(file_training_config,*) n_MC_total
     n_MC_total = (n_MC_total/nImg) * nImg ! make sure n_MC is divisible
     close(file_training_config)
     write(*,*) 'Number of total MC:', n_MC_total
     write(*,*)
  end if

  call co_broadcast(m_bounds, 1)
  call co_broadcast(g_bounds, 1)
  call co_broadcast(Lam_bounds, 1)
  call co_broadcast(v_bounds, 1)
  call co_broadcast(Lam_v_bounds, 1)
  call co_broadcast(n_MC_total, 1)
  n_MC = n_MC_total / nImg

  ! ----------------------Read in scattering data---------------------
  if (IamRoot) then
     open(142, file='dataInf_orig.in', action='read')
     read(142,*) nPoints_inf
  end if

  call co_broadcast(nPoints_inf, 1)
  ! Allocate the arrays
  allocate( E_data(nPoints_inf), phaseShiftData(nPoints_inf)&
       & , phaseShiftDataErr(nPoints_inf) &
       & , etaData(nPoints_inf), etaDataErr(nPoints_inf) )

  if (IamRoot) then
     do iii = 1,nPoints_inf
        read(142,*) E_data(iii) &
             & , phaseShiftData(iii), phaseShiftDataErr(iii) &
             & , etaData(iii), etaDataErr(iii)
     end do

     ! Convert from Sr to eta
     etaData(:) = sqrt(1.0_DP - etaData(:))
     etaDataErr(:) = 0.5_DP / etaData(:) * etaDataErr(:)
     close(142)
  end if

  call co_broadcast(E_data, 1)
  call co_broadcast(phaseShiftData, 1)
  call co_broadcast(phaseShiftDataErr, 1)
  call co_broadcast(etaData, 1)
  call co_broadcast(etaDataErr, 1)
  ! ------------------------------------------------------------------

  allocate( lowBound(n_ch), midBound(n_ch) )
  ! Stuff for phase continuity
  allocate( isClockwise(n_ch), openChPoint(n_ch), isOpenCh(n_ch) )
  allocate( old_SMat(n_ch,n_ch) )
  isOpenCh(:) = .false.

  ! Initialisations
  allocate( eta_arr(nPoints_inf) &
       & , phaseShift(nPoints_inf,n_ch) &
       & , crossSec(nPoints_inf,n_ch,n_ch) )
  allocate( E_arr(nPoints_inf) )
  lowBound(:) = m_mes(:) + m_bar(:)
  absErr = 1.0E-6_DP
  relErr = 0.0_DP

  ! --------------------Construct the energy array--------------------
  if (useDataPoints) then
     if (E_data(1).gt.1000) then
        E_arr(:) = E_data(:)/1000.0_DP
     else
        E_arr(:) = E_data(:)
     end if
     E_init = E_arr(1)
     E_fin = E_arr(nPoints_inf)
  else
     E_arr(:) = (E_fin-E_init)/(nPoints_inf-1) &
          & * (/ (iii,iii=0,nPoints_inf-1) /) + E_init
  end if

  ! --------------------File to save training data--------------------
  if (IamRoot) then
     open(101, file='/home/a1686392/TrainingData/'//trim(adjustl(projName)) &
          & //'_training_N'//trim(int2str(n_MC_total)) &
          & //'.data', action='write')
     write(101, *) '         HEFT ' // trim(int2str(n_bare)) // 'b'&
          & // trim(int2str(n_ch)) // 'c Training Data'
     call cpu_time(t_init)
  end if

  do iMC = 1, n_MC
     ! --------------------------Do Monte-Carlo--------------------------
     call random_number(rand_params)
     if (n_bare .gt. 0) then
        rand_params(1:bare_end) = (m_bounds(2)-m_bounds(1)) &
             & * rand_params(1:bare_end) &
             & + m_bounds(1)
        m_bare = rand_params(1:bare_end)

        rand_params(bare_end+1:g_end) = (g_bounds(2)-g_bounds(1)) &
             & * rand_params(bare_end+1:g_end) &
             & + g_bounds(1)
        gBare = vectorToMat(rand_params(bare_end+1:g_end) &
             & , n_ch, n_bare)

        rand_params(g_end+1:Lam_end) = (Lam_bounds(2)-Lam_bounds(1)) &
             & * rand_params(g_end+1:Lam_end) &
             & + Lam_bounds(1)
        Lambda = vectorToMat(rand_params(g_end+1 : Lam_end) &
             & ,  n_ch, n_bare)
     end if

     rand_params(Lam_end+1:v_end) = (v_bounds(2)-v_bounds(1)) &
          & * rand_params(Lam_end+1:v_end) &
          & + v_bounds(1)
     vCh = vectorToSymMat(rand_params(Lam_end+1 : v_end))

     rand_params(v_end+1 : Lam_v_end) &
          & = (Lam_v_bounds(2)-Lam_v_bounds(1)) &
          & * rand_params(v_end+1 : Lam_v_end) &
          & + Lam_v_bounds(1)
     Lambda_v = rand_params(v_end+1 : Lam_v_end)

     ! -----------------Calculate scattering observables-----------------
     do i = 1, nPoints_inf
        globalPoint = i
        call calcSMatrix(E_arr(i), phaseShift(i,:) &
             & , eta_arr(i), crossSec(i,:,:))
     end do

     chi2 = sum( (phaseShift(:,ch_onshell) - phaseShiftData(:))**2 &
          & / phaseShiftDataErr(:)**2 ) &
          & + sum( (eta_arr(:) - etaData(:))**2 / etaDataErr(:)**2 )

     ! ------------------Write the training data to file-----------------
     sync all
     if (IamRoot) then
        do jImg = 1, nImg
           write(101,'('//trim(int2str(nParamsTotal))//'f14.8,es22.8)') &
                & rand_params(:)[jImg], chi2[jImg]
        end do
     end if

     ! --------------------Write progress to terminal--------------------
     if (IamRoot .and. verboseTerminal) then
        write(*,'(5x,a1,i7,a1,i0,2x,a1,2x,i3,a4,a,\)') '|', iMC*nImg, '/' &
             & , n_MC*nImg, '|', int(100*real(iMC,DP)/real(n_MC,DP)) &
             & , '%  |', char(13)
        if ((iMC .eq. 1) .or. (mod(iMC, n_MC/nProgressPrints) .eq. 0) &
             & .or. (iMC .eq. n_MC)) then
           write(*,'(5x,a1,i7,a1,i0,2x,a1,2x,i3,a4)') '|', iMC*nImg, '/' &
                & , n_MC*nImg, '|', int(100*real(iMC,DP)/real(n_MC,DP)) &
                & , '%  |'
        end if
     end if
  end do


  ! --------------Write an info file about this training--------------
  if (IamRoot) then
     call cpu_time(t_final)
     open(103, file='/home/a1686392/TrainingData/'//trim(projName) &
          & //'_training_N'//trim(int2str(n_MC_total)) &
          & //'.info', action='write')
     write(103, *) '--------HEFT ' // trim(int2str(n_bare)) // 'b'&
          & // trim(int2str(n_ch)) // 'c Training Data--------'
     write(103,*) 'nImages: ' // trim(int2str(nImg))
     if (n_bare.gt.0) then
        write(103,*) 'Bare states: '
        do i = 1, n_bare
           write(103,*) '     ', trim(bare_labels(i))
        end do
     end if
     write(103,*) 'Scattering Channels: '
     do i = 1, n_ch
        write(103,*) '     ', trim(ch(1,i)%name) // '-' // trim(ch(2,i)%name)
     end do

     write(103,*)
     write(103,*) 'N_MC:  '//trim(int2str(n_MC_total))
     write(103,*) 'Parameter Bounds:'
     write(103,'(a8,2f12.6)') 'm'     , m_bounds
     write(103,'(a8,2f12.6)') 'g'     , g_bounds
     write(103,'(a8,2f12.6)') 'Lam'   , Lam_bounds
     write(103,'(a8,2f12.6)') 'v'     , v_bounds
     write(103,'(a8,2f12.6)') 'Lam_v' , Lam_v_bounds

     write(103,*)
     write(103,*) 'Time:  ', timePrint(t_final-t_init)
  end if

  ! -----------------------------Finalise-----------------------------
  if (IamRoot) then
     write(*,*) 'Time:  ', timePrint(t_final-t_init)
     write(*,*)
     close(101)
     close(103)
  end if
  deallocate( phaseShift, eta_arr, crossSec, rand_params)
  deallocate( lowBound, midBound, E_arr )
  deallocate( isClockwise, openChPoint, isOpenCh )
  deallocate( E_data, phaseShiftData, phaseShiftDataErr &
       & , etaData, etaDataErr )
  call finaliseHEFT()

  !-----------------------------------------------------------------------!

contains


  function vectoriseUpperTriangle(M) result(vec)
    implicit none
    ! Converts an upper triangular matrix into an array, column-wise
    real(DP), dimension(:,:) :: M
    real(DP), dimension((size(M,1)*(size(M,1)+1))/2) :: vec
    ! Local variables
    integer :: dimM
    integer :: iM, jM

    dimM = size(M,1)
    do jM = 1,dimM
       do iM = 1,jM
          vec((jM*(jM-1))/2 + iM) = M(iM,jM)
       end do
    end do

  end function vectoriseUpperTriangle

end program trainingInfiniteVol
