! In this file I have subroutines calcSMatrix() and calcSMatrixPole()
! calcSMatrix() is used for fitting the scattering observables, and therefore
!    only needs real valued energies/momenta
! calcSMatrixPole() is for pole searches, which require complex E and k
!
! In principle I could just do all the calculations in calcSMatrixPole() and
!    then have calcSMatrix() just call calcSMatrixPole(), converting the
!    results to real numbers. Doing  all of the calculations with
!    complex numbers  however has a sizable impact on the time
!    it took to fit Hamiltonian to the scattering observables.
module SMatrix
    use kinds
    use numFort
    use heft
    use, intrinsic :: ieee_arithmetic
    implicit none

    integer  :: ich, jch, ibare=1, jbare=1
    real(DP) :: absErr, relErr ! integral errors
    real(DP), dimension(:), allocatable :: lowBound ! integral bounds
    real(DP), dimension(:), allocatable :: midBound

    ! Pole search variables
    real(DP), parameter :: pole_phase = -70.0_DP / degrees
    real(DP), dimension(:), allocatable :: pole_phases ! todo: different pole phases per channel

    complex(DP) :: k_complex_theta = cmplx(cos(pole_phase), sin(pole_phase), DP)
    integer :: i_bare_pole

    ! phase continuity variables
    real(DP) :: phaseDiff
    real(DP) :: dtheta, thetaNew, thetaOld
    complex(DP) :: Sii_old, Sii
    complex(DP), dimension(:,:), allocatable :: old_SMat
    logical, dimension(:), allocatable :: isClockwise, isOpenCh
    integer, dimension(:), allocatable :: openChPoint
    integer :: globalPoint, atan_ch

    complex(DP) :: E_pole_min
    real(DP) :: chi2_pole_min = HUGE(1.0d0)

    integer :: testInt = 0

contains

    subroutine calcSMatrix(E_in, phase, eta, cSec, Gamma_ij, Tmat_out)
        implicit none
        real(DP), intent(in)                  :: E_in ! on-shell energy
        real(DP), dimension(:), intent(out)   :: phase ! phase shift
        real(DP), intent(out)                 :: eta ! inelasticity
        real(DP), dimension(:,:), intent(out), optional :: cSec ! cross section (onshell)
        real(DP), dimension(:), intent(out), optional :: Gamma_ij ! decay rates
        complex(DP), intent(out), optional :: Tmat_out ! tmatrix (onshell element)

        ! Local Variables
        integer  :: ii, jj, bi, bj
        real(DP), dimension(1,1) :: mat11

        complex(DP) :: sumTerm
        real(DP) :: brkpt
        ! complex(DP) :: Sigma_BB, SigmaI_BB ! Self energy and background self energy
        complex(DP), dimension(:,:), allocatable :: Sigma_BB, SigmaI_BB ! Self energy and background self energy
        complex(DP), dimension(:,:), allocatable :: g_f_iB
        complex(DP), dimension(:), allocatable :: rho_i
        real(DP), dimension(:), allocatable :: k_in
        real(DP), dimension(:), allocatable :: dkdE
        real(DP), dimension(:), allocatable :: fsqrd
        real(DP), dimension(:), allocatable :: Re_S, Im_S
        complex(DP), dimension(:), allocatable :: S_diag
        complex(DP), dimension(:,:), allocatable :: Ainv
        complex(DP), dimension(:,:), allocatable :: M_ii
        complex(DP), dimension(:,:), allocatable :: gbar_iB
        complex(DP), dimension(:,:), allocatable :: SMatrix
        complex(DP), dimension(:,:), allocatable :: TMat, ttMat, ttildeMat

        E_gl = E_in
        brkpt = E_gl * 4.0_DP

        allocate( k_in(n_ch), SMatrix(n_ch,n_ch), TMat(n_ch,n_ch) &
            &, ttMat(n_ch,n_ch), ttildeMat(n_ch,n_ch), M_ii(n_ch,n_ch) &
            &, rho_i(n_ch), dkdE(n_ch), fsqrd(n_ch), Re_S(n_ch), Im_S(n_ch) &
            &, S_diag(n_ch))
        allocate( Ainv(n_bare,n_bare), g_f_iB(n_ch,n_bare), gbar_iB(n_ch,n_bare) &
            &, Sigma_BB(n_bare,n_bare), SigmaI_BB(n_bare,n_bare) )

        ! Eq (2.54)
        k_in(:) = sqrt( E_gl**2/4.0_DP + (m_bar(:)**2 - m_mes(:)**2)**2 &
            & /4.0_DP/E_gl**2 - (m_mes(:)**2 + m_bar(:)**2)/2.0_DP )

        ! Principle value integrals need to have an integrand denominator
        !    of the form x - x_0, where x in the integral variable.
        ! As the integrals are all defined to integrate over k however,
        !    with integrand denominators of the form E - omega(E), we need
        !    to do a change of variables so that we integrate over E instead.
        dkdE(:) = sqrt(k_in(:)**2 + m_mes(:)**2)*sqrt(k_in(:)**2 + m_bar(:)**2) &
            & /k_in(:)/E_gl

        ! ----------------------Bare state Self-energy----------------------
        ! Here, and further on, I set the loop variables to a global variable.
        !    This allows the integrand functions to access them, so they can
        !    Access the correct hadron masses, couplings, etc.
        !
        ! Eq. (2.43)
        Sigma_BB(:,:) = cmplx(0.0_DP,0.0_DP,DP)
        do bi = 1,n_bare
           ibare = bi
           do bj = 1,n_bare
              jbare = bj
              do ii = 1,n_ch
                 ich = ii
                 if (E_gl .ge. lowBound(ich)) then
                     ! Need a principle value integral when the channel is open
                     Sigma_BB(ibare,jbare) = Sigma_BB(ibare,jbare) &
                         & + cmplx( integralPV(Sigma_integrandPV, E_gl &
                         & , lowBound(ich)+0.00001_DP, brkpt, absErr, relErr) &
                         & + integral(Sigma_integrand, brkpt &
                         & , Infty, absErr, relErr) &
                         & , -pi * dkdE(ich) * k_in(ich)**2 &
                         & * gBare(ich,ibare) * gBare(ich,jbare) &
                         & * g_i(k_in(ich),ich,ibare) &
                         & * g_i(k_in(ich),ich,jbare) , DP )

                 else
                     ! If the channel is closed, just do a normal integral
                     Sigma_BB(ibare,jbare) = Sigma_BB(ibare,jbare) &
                         & + integral(Sigma_integrand, lowBound(ich) &
                         & , Infty, absErr, relErr)
                 end if
              end do
           end do
        end do

        ! ----------------------------M integral----------------------------
        ! M needs to be a diagonal matrix for the matrix equation
        !    in Eq. (2.31) to be true
        !
        ! Eq. (2.30)
        M_ii(:,:) = 0.0_DP
        do ii = 1,n_ch
           ich = ii
           ! TODO: generalise this properly to allow for n_f-separable potentials
           mat11 = matmul(transpose(f_i(k_in(ich),ich)), f_i(k_in(ich),ich))
           fsqrd(ich) = mat11(1,1)

           if (E_gl.ge.lowBound(ich)) then
               M_ii(ich,ich) = &
                   & cmplx( integralPV(M_integrandPV, E_gl &
                   & , lowBound(ich)+0.00001_DP, brkpt, absErr, relErr) &
                   & + integral(M_integrand, brkpt, Infty, absErr, relErr) &
                   & , -pi * dkdE(ich) * k_in(ich)**2 * fsqrd(ich), DP )
           else
               M_ii(ich,ich) = integral(M_integrand, lowBound(ich) &
                   & , Infty, absErr, relErr)
           end if
        end do

        ! ------------------------Background t matrix-----------------------
        !
        ! Eq. (2.31)
        ! I do this in two steps purely for clarity
        ttildeMat(:,:) = id(:,:) - matmul(vCh(:,:), M_ii(:,:))
        ttildeMat(:,:) = matmul(inv(ttildeMat(:,:)), vCh(:,:))

        ! Eq. (2.32)
        ! TODO: should be able to do this with matrix multiplication
        do jj = 1, n_ch
           do ii = 1, n_ch
              ttMat(ii,jj) = sqrt(fsqrd(ii)) * ttildeMat(ii,jj) * sqrt(fsqrd(jj))
           end do
        end do

        ! ---------------------------g_f integral---------------------------
        !
        ! Eq. (2.39)
        ! TODO: this needs to be generalised for n_f-separable potentials,
        !       g_f_iB should be a rank-3 object.
        g_f_iB(:,:) = cmplx(0.0_DP,0.0_DP,DP)
        do bi = 1, n_bare
           ibare = bi
           do ii = 1, n_ch
              ich = ii
              if (E_gl .ge. lowBound(ich)) then
                  mat11(:,:) = f_i(k_in(ich),ich)
                  g_f_iB(ich,ibare) = cmplx( integralPV(g_f_integrandPV, E_gl &
                      & , lowBound(ii)+0.00001_DP, brkpt, absErr, relErr) &
                      & + integral(g_f_integrand, brkpt, Infty, absErr, relErr) &
                      & , -pi * dkdE(ich) * k_in(ich)**2 * mat11(1,1) &
                      & * gBare(ich,ibare) * g_i(k_in(ich),ich,ibare), DP )
              else
                  g_f_iB(ich,ibare) = integral(g_f_integrand, lowBound(ich) &
                      & , Infty, absErr, relErr)
              end if
           end do
        end do

        ! ----------------------Background self energy----------------------
        !
        ! Eq. (2.46)
        SigmaI_BB(:,:) = cmplx(0.0_DP,0.0_DP,DP)
        do bi = 1,n_bare
           ibare = bi
           do bj = 1,n_bare
              jbare = bj
              do jj = 1,n_ch
                 jch = jj
                 do ii = 1,n_ch
                    ich = ii
                    SigmaI_BB(ibare,jbare) = SigmaI_BB(ibare,jbare) + g_f_iB(ich,ibare) &
                        & * ttildeMat(ich,jch) * g_f_iB(jch,jbare)
                 end do
              end do
           end do
        end do

        ! ----------------------Calculate Bare T-matrix---------------------
        !
        ! Eq. (2.40)
        ! Used to denote it \bar{G}, but I've changed it to \mathcal{G}
        !    to better differentiate with the finite-volume \bar{G}.
        ! TODO: change name of variable to represent this
        do bi = 1,n_bare
           ibare = bi
           do ii = 1,n_ch
              sumTerm = cmplx(0.0_DP,0.0_DP,DP)
              ich = ii
              do jj = 1,n_ch
                 sumTerm = sumTerm + ttildeMat(ii,jj) * g_f_iB(jj,ibare)
              end do
              ! TODO: generalise for n_f-separable potentials,
              !       should have a sum over i_f (n_f functions)
              mat11(:,:) = f_i(k_in(ich),ich)
              gbar_iB(ii,ibare) = gBare(ich,ibare) * g_i(k_in(ich),ich,ibare) &
                  & + mat11(1,1) * sumTerm
           end do
        end do

        ! Eq. (2.42)
        Ainv(:,:) = 0.0_DP
        do bi = 1, n_bare
           ibare = bi
           Ainv(ibare,ibare) = E_gl - m_bare(ibare)
        end do
        Ainv(:,:) = Ainv(:,:) - Sigma_BB(:,:) - SigmaI_BB(:,:)

        ! ----------------Add bare contributions to T matrix----------------
        !
        ! Eq. (2.34)
        TMat(:,:) = ttMat(:,:)
        if (n_bare.gt.0) then ! only add the bare component if there is a bare state
            TMat(:,:) = TMat(:,:) + matmul( matmul(gBar_iB(:,:),inv(Ainv(:,:))) &
                & , transpose(gBar_iB(:,:)) )
        end if

        ! -------------------Convert T-matrix to S-matrix-------------------
        !
        ! Eq. (2.55)
        ! Calculate density of states
        ! Wu et al: Finite-volume Hamiltonian method for
        !     coupled channel interactions in lattice QCD
        rho_i(:) = pi * k_in(:) / E_gl &
            & * sqrt(k_in(:)**2 + m_mes(:)**2)*sqrt(k_in(:)**2 + m_bar(:)**2)

        do jj = 1,n_ch
           do ii = 1,n_ch
              ! Eq. (2.56)
              SMatrix(ii,jj) = id(ii,jj) - 2.0_DP*cmplx(0.0_DP,1.0_DP,DP) &
                  & * sqrt(rho_i(ii)*rho_i(jj)) * TMat(ii,jj)
           end do

           ! Calculate decay rates using Fermi's golden rule,
           ! TODO: Double check which one was correct using our parametrisation
           !       of the S-matrix, think it's the second one
           if (present(Gamma_ij)) then
               ! Using T = i(1-S)
               ! Gamma_ij(jj) = 2.0_DP*pi * rho_i(jj) &
               !     & * abs(cmplx(0.0_DP,1.0_DP,DP)*(1.0_DP - SMatrix(ch_onshell,jj)))**2

               ! Using T = i(1-S)/2
               Gamma_ij(jj) = 2.0_DP*pi * rho_i(jj) &
                   & * abs(cmplx(0.0_DP,0.5_DP,DP)*(1.0_DP - SMatrix(ch_onshell,jj)))**2

               ! Using T matrix
               ! Gamma_ij(jj) = 2.0_DP*pi * abs(TMat(ch_onshell,jj))**2 * rho_i(jj)
           end if
        end do

        ! -----------------Calculate onshell Tmatrix element----------------
        !
        if (present(Tmat_out)) then
            ! Using T = i(1-S)  (S = 1 + iT)
            ! seems like I need to divide by 2 to match P33 WI08 solution for some reason
            Tmat_out = cmplx(0.0_DP, 0.5_DP, DP) &
                & * (1.0_DP - SMatrix(ch_onshell,ch_onshell))

            ! Using T = i(S-1)  (S = 1 - iT)
            ! Tmat_out = cmplx(0.0_DP,1.0_DP,DP) * (SMatrix(ch_onshell,ch_onshell) - 1.0_DP)

            ! Using calculated T-matrix
            ! Tmat_out = TMat(ch_onshell,ch_onshell)

            ! Output S-matrix for testing
            ! Tmat_out = SMatrix(ch_onshell,ch_onshell)
        end if

        ! -----------------Calculate Scattering Observables-----------------
        !
        ! There's an issue here with extracting the phase shifts. The atan/atan2
        !    functions both have a range of spread over 2pi (todo: check ranges)
        ! The phase shift however isn't limited to that range, and therefore to
        !    have a continuous phase shift when it leaves the range we need to
        !    "fix" the phase shift to the correct value.
        !    This turns out to be easier said than done automatically,
        !    and I've got it working well enough using the atanRangeFixer function.
        !
        ! TODO: Must be a better way than what I've done, especially because this method
        !       prevents parallelisation since you need to calculate
        !       the observables sequentially this way.
        !
        ! Eq. (2.58)
        do ii = 1,n_ch
           atan_ch = ii
           S_diag(ii) = SMatrix(ii,ii)
           Re_S(ii) = real(S_diag(ii))
           Im_S(ii) = aimag(S_diag(ii))
           phase(ii) = 0.5_DP*atan(real(Im_S(ii),DP)/real(Re_S(ii),DP)) * degrees

           ! Check phase shift contnuity for open channels
           if (isOpenCh(ii)) then
               Sii = SMatrix(ii,ii)
               Sii_old = old_SMat(ii,ii)
               ! This doesn't seem to always work when it goes above 180 deg :(
               phase(ii) = phase(ii) &
                   & + atanRangeFixer(Sii_old,Sii)*degrees/2.0_DP
               old_SMat(ii,ii) = Sii
           end if

           ! If the channel wasn't already open, or it was i=1, check if the channel
           !     needs to be opened for phase shift continuity
           if (.not.ieee_is_nan(phase(ii)) .and. .not.isOpenCh(ii)) then
               isOpenCh(ii) = .true.
               openChPoint(ii) = globalPoint
               old_SMat(ii,ii) = SMatrix(ii,ii)
           end if
        end do

        ! Eq. (2.58)
        ! Calculate inelasticity
        if (n_ch.gt.1) then
            eta = sqrt(Re_S(ch_onshell)**2 + Im_S(ch_onshell)**2)
        end if

        ! Eq. (2.53)
        ! Caluclate cross section
        if (present(cSec)) then
            do jj = 1,n_ch
               do ii = 1,n_ch
                  cSec(ii,jj) = 4.0_DP*pi**3/E_gl**2 * k_in(ii)/k_in(jj) &
                      & * sqrt(k_in(ii)**2 + m_mes(ii)**2)*sqrt(k_in(ii)**2 + m_bar(ii)**2) &
                      & * sqrt(k_in(jj)**2 + m_mes(jj)**2)*sqrt(k_in(jj)**2 + m_bar(jj)**2) &
                      & * abs(TMat(ii,jj))**2 * GeV2_mb
               end do
            end do
        end if

        deallocate( k_in, SMatrix,  TMat, ttMat, ttildeMat, M_ii &
            & , rho_i, dkdE, fsqrd, Re_S, Im_S )

        if (n_bare .gt. 0) then
            deallocate( Ainv, g_f_iB, gbar_iB, Sigma_BB, SigmaI_BB )
        end if

    end subroutine calcSMatrix


    subroutine calcSMatrixPole(n, x, f)
        implicit none
        integer, intent(in) :: n
        real(DP), dimension(n), intent(in) :: x
        real(DP), intent(out) :: f

        real(DP)                 :: E_in

        ! Local Variables
        integer  :: ii, jj, bi, bj
        complex(DP), dimension(1,1) :: mat11

        complex(DP) :: sumTerm
        complex(DP) :: detAinv
        ! complex(DP) :: Sigma_BB, SigmaI_BB ! Self energy and background self energy
        real(DP), dimension(:), allocatable      :: Re_S, Im_S
        complex(DP), dimension(:,:), allocatable :: Sigma_BB, SigmaI_BB ! Self energy and background self energy
        complex(DP), dimension(:), allocatable   :: rho_i
        complex(DP), dimension(:), allocatable   :: k_in
        complex(DP), dimension(:), allocatable   :: dkdE
        complex(DP), dimension(:), allocatable   :: fsqrd
        complex(DP), dimension(:), allocatable   :: S_diag
        complex(DP), dimension(:,:), allocatable :: g_f_iB
        complex(DP), dimension(:,:), allocatable :: Ainv
        complex(DP), dimension(:,:), allocatable :: M_ii
        complex(DP), dimension(:,:), allocatable :: gbar_iB
        complex(DP), dimension(:,:), allocatable :: SMatrix
        complex(DP), dimension(:,:), allocatable :: TMat, ttMat
        complex(DP), dimension(:,:), allocatable :: ttildeMat, ttildeMat_inv

        complex(DP) :: mu_1, mu_2, S11, S22, S12, S21

        real(DP) :: Re_E, Im_E
        real(DP) :: tempMin

        logical :: verboseMin = .false.

        Re_E = x(1)
        Im_E = x(2)
        E_pole = cmplx(Re_E, Im_E, DP)

        ! if (isEnergyDepPotential) vCh(:,:) = C_ij_E()

        allocate( k_in(n_ch), SMatrix(n_ch,n_ch), TMat(n_ch,n_ch) &
            & , ttMat(n_ch,n_ch), ttildeMat(n_ch,n_ch), M_ii(n_ch,n_ch) &
            & , rho_i(n_ch), dkdE(n_ch), fsqrd(n_ch), Re_S(n_ch), Im_S(n_ch) &
            & , Ainv(n_bare,n_bare), g_f_iB(n_ch,n_bare), gbar_iB(n_ch,n_bare) &
            & , S_diag(n_ch), Sigma_BB(n_bare,n_bare), SigmaI_BB(n_bare,n_bare) &
            & , ttildeMat_inv(n_ch,n_ch))

        k_in(:) = sqrt( E_pole**2/4.0_DP + (m_bar(:)**2 - m_mes(:)**2)**2 &
            & /4.0_DP/E_pole**2 - (m_mes(:)**2 + m_bar(:)**2)/2.0_DP )
        dkdE(:) = sqrt(k_in(:)**2 + m_mes(:)**2)*sqrt(k_in(:)**2 + m_bar(:)**2) &
            & /k_in(:)/E_pole

        ! ---------------------------Self Energy---------------------------
        Sigma_BB(:,:) = 0.0_DP
        if (n_bare.gt.0) then
            do bj = 1,n_bare
               jbare = bj
               do bi = 1,n_bare
                  ibare = bi
                  do ii = 1,n_ch
                     ich = ii
                     Sigma_BB(ibare,jbare) = Sigma_BB(ibare,jbare) + cmplx( &
                         &   integral(Sigma_k_integrand_real, 0.0_DP, Infty, absErr, relErr) &
                         & , integral(Sigma_k_integrand_imag, 0.0_DP, Infty, absErr, relErr), DP)
                  end do
               end do
            end do
        end if


        ! ----------------------------M integral----------------------------
        M_ii(:,:) = 0.0_DP
        do ii = 1,n_ch
           ich = ii
           M_ii(ich,ich) = cmplx( &
               &   integral(M_k_integrand_real, 0.0_DP, Infty, absErr, relErr) &
               & , integral(M_k_integrand_imag, 0.0_DP, Infty, absErr, relErr), DP)
        end do

        ttildeMat_inv(:,:) = id(:,:) - matmul(vCh(:,:), M_ii(:,:))
        ttildeMat(:,:) = matmul(inv(ttildeMat_inv(:,:)), vCh(:,:))

        do jj = 1,n_ch
           do ii = 1,n_ch
              ttMat(ii,jj) = sqrt(fsqrd(ii)) * ttildeMat(ii,jj) * sqrt(fsqrd(jj))
           end do
        end do

        ! ---------------------------g_f integral---------------------------
        g_f_iB(:,:) = cmplx(0.0_DP,0.0_DP,DP)
        do bi = 1,n_bare
           ibare = bi
           do ii = 1,n_ch
              ich = ii
              g_f_iB(ich,ibare) = cmplx( &
                  &   integral(g_f_k_integrand_real, 0.0_DP, Infty, absErr, relErr) &
                  & , integral(g_f_k_integrand_imag, 0.0_DP, Infty, absErr, relErr),DP)
           end do
        end do

        ! ----------------------Background self energy----------------------
        SigmaI_BB(:,:) = cmplx(0.0_DP,0.0_DP,DP)
        do bi = 1,n_bare
           ibare = bi
           do bj = 1,n_bare
              jbare = bj
              do jj = 1,n_ch
                 jch = jj
                 do ii = 1,n_ch
                    ich = ii
                    SigmaI_BB(ibare,jbare) = SigmaI_BB(ibare,jbare) + g_f_iB(ich,ibare) &
                        & * ttildeMat(ich,jch) * g_f_iB(jch,jbare)
                 end do
              end do
           end do
        end do


        ! ------------------------Calculate S-matrix------------------------
        do bi = 1,n_bare
           ibare = bi
           do ii = 1,n_ch
              sumTerm = cmplx(0.0_DP,0.0_DP,DP)
              ich = ii
              ! Still needs to be generalised for n_f seperable functions
              do jj = 1,n_ch
                 sumTerm = sumTerm + ttildeMat(ii,jj) * g_f_iB(jj,ibare)
              end do
              mat11(:,:) = f_i(k_in(ich),ich)
              gbar_iB(ii,ibare) = gBare(ich,ibare) * g_i(k_in(ich),ich,ibare) &
                  & + mat11(1,1) * sumTerm
           end do
        end do

        Ainv(:,:) = 0.0_DP
        do bi = 1, n_bare
           ibare = bi
           Ainv(ibare,ibare) = E_pole - m_bare(ibare)
        end do
        Ainv(:,:) = Ainv(:,:) - Sigma_BB(:,:) - SigmaI_BB(:,:)
        detAinv = det(Ainv(:,:))

        TMat(:,:) = ttMat(:,:)

        if (n_bare.gt.0) then
            TMat(:,:) = TMat(:,:) + matmul( matmul(gBar_iB(:,:),inv(Ainv(:,:))) &
                & , transpose(gBar_iB(:,:)) )
        end if

        rho_i(:) = sqrt(k_in(:)**2 + m_mes(:)**2)*sqrt(k_in(:)**2 + m_bar(:)**2) &
            & / E_pole * k_in(:)

        do jj = 1,n_ch
           do ii = 1,n_ch
              SMatrix(ii,jj) = id(ii,jj) - 2.0_DP*cmplx(0.0_DP,1.0_DP,DP)*pi &
                  & * sqrt(rho_i(ii)*rho_i(jj)) * TMat(ii,jj)
           end do
        end do


        f = HUGE(1.0_DP)
        if (n_bare.eq.0) then
            f = abs(det(ttildeMat_inv(:,:)))**2
        else
            f = abs(detAinv)**2
        end if



        if (f.lt.chi2_pole_min) then
            chi2_pole_min = f
            E_pole_min = E_pole
            if (verboseMin) write(*,*) E_pole_min, chi2_pole_min
            ! stop
        end if

        ! Some testing variables for the algebraic poles from
        !    the multi bare state paper
        !
        ! if (testInt.eq.0) then
        !     S11 = Sigma_BB(1,1) + SigmaI_BB(1,1)
        !     S22 = Sigma_BB(2,2) + SigmaI_BB(2,2)
        !     S12 = Sigma_BB(1,2) + SigmaI_BB(1,2)
        !     S21 = Sigma_BB(2,1) + SigmaI_BB(2,1)
        !     mu_1 = m_bare(1) + m_bare(2) + S11 + S22
        !     mu_2 = m_bare(2)*S11 + m_bare(1)*S22 + m_bare(1)*m_bare(2) &
        !         & + S11*S22 - S12*S21

        !     write(*,*)
        !     write(*,*) 'detA(E)num', Ainv(1,1)*Ainv(2,2) - Ainv(1,2)*Ainv(2,1)
        !     write(*,*) 'detA(E)num', detAinv
        !     write(*,*) 'detA(E)alg', E_pole**2 - E_pole*mu_1 + mu_2
        !     write(*,*)
        !     write(*,*) 'sqrt   ', sqrt(mu_1**2 - 4.0_DP*mu_2)
        !     write(*,*) 'm plus ', 0.5_DP*( sqrt(mu_1**2 - 4.0_DP*mu_2) + mu_1)
        !     write(*,*) 'm plus ', 0.5_DP*( sqrt(mu_1**2 - 4.0_DP*mu_2) + mu_1)
        !     write(*,*) 'm minus', 0.5_DP*(-sqrt(mu_1**2 - 4.0_DP*mu_2) + mu_1)
        !     write(*,*) 'pole1  ', cmplx(1.5_DP, -0.05_DP, DP)
        !     write(*,*) 'pole2  ', cmplx(1.657_DP, -0.056_DP, DP)
        !     write(*,*)
        !     stop

        !     write(*,*)
        !     write(*,*) 'E_pole_guess = ', E_pole
        !     write(*,*) 'Sigma_BB', Sigma_BB
        !     write(*,*) 'g_f', g_f_iB
        !     write(*,*) 'M_ii', M_ii
        !     write(*,*) 'tt_tilde', ttildeMat
        !     write(*,*) 'Sigma_I', SigmaI_BB
        !     write(*,*) 'Tmat', TMat
        !     write(*,*) 'Smat', SMatrix
        !     write(*,*) 'T-inverse', f
        ! end if
        ! testInt = testInt+1

        deallocate( k_in, SMatrix,  TMat, ttMat, ttildeMat, M_ii &
            & , rho_i, dkdE, fsqrd, Re_S, Im_S, Ainv, g_f_iB, gbar_iB &
            & , S_diag, Sigma_BB, SigmaI_BB, ttildeMat_inv )

    end subroutine calcSMatrixPole


    subroutine calcSMatrixPole_bq(n, x, f)
        ! The new minimisation routine using BOBYQA uses a slightly
        !    different function structure to minfun. This just converts
        !    to the BOBYQA format
        implicit none
        integer, intent(in) :: n
        real(DP), dimension(:), intent(in) :: x
        real(DP), intent(out) :: f ! chi^2

        call calcSMatrixPole(n, x, f)

        if (f.lt.chi2_pole_min) then
            chi2_pole_min = f
        end if
    end subroutine calcSMatrixPole_bq



    ! -------------------------Sigma integrands-------------------------
    ! Eq. (2.43)

    function Sigma_integrandPV(E)
        implicit none
        real(DP), intent(in) :: E
        real(DP) :: Sigma_integrandPV

        real(DP) :: q, dqdE

        q = sqrt( E**2/4.0_DP + (m_bar(ich)**2 - m_mes(ich)**2)**2 &
            & /4.0_DP/E**2 - (m_mes(ich)**2 + m_bar(ich)**2)/2.0_DP )

        dqdE = sqrt(q**2 + m_mes(ich)**2)*sqrt(q**2 + m_bar(ich)**2) / q / E
        Sigma_integrandPV = -q**2 * dqdE &
            & * gBare(ich,ibare) * gBare(ich,jbare) &
            & * g_i(q,ich,ibare) * g_i(q,ich,jbare)

    end function Sigma_integrandPV

    function Sigma_integrand(E)
        implicit none
        real(DP), intent(in) :: E
        real(DP) :: Sigma_integrand

        Sigma_integrand = Sigma_integrandPV(E) / (E - E_gl)

    end function Sigma_integrand


    function Sigma_k_integrand_real(k_in)
        implicit none
        real(DP), intent(in) :: k_in
        real(DP) :: Sigma_k_integrand_real
        complex(DP) :: Sigma_k_integrand_complex

        complex(DP) :: q, E

        q = k_in * k_complex_theta
        E = sqrt(q**2 + m_mes(ich)**2) + sqrt(q**2 + m_bar(ich)**2)
        Sigma_k_integrand_complex = -q**2 * gBare(ich,ibare) * gBare(ich,jbare) &
            & * g_i(q,ich,ibare) * g_i(q,ich,jbare) &
            & / (E - E_pole) * k_complex_theta
        Sigma_k_integrand_real = real(Sigma_k_integrand_complex)

    end function Sigma_k_integrand_real

    function Sigma_k_integrand_imag(k_in)
        implicit none
        real(DP), intent(in) :: k_in
        real(DP) :: Sigma_k_integrand_imag
        complex(DP) :: Sigma_k_integrand_complex

        complex(DP) :: q, E

        q = k_in * k_complex_theta
        E = sqrt(q**2 + m_mes(ich)**2) + sqrt(q**2 + m_bar(ich)**2)
        Sigma_k_integrand_complex = -q**2 * gBare(ich,ibare) * gBare(ich,jbare) &
            & * g_i(q,ich,ibare) * g_i(q,ich,jbare) &
            & / (E - E_pole) * k_complex_theta
        Sigma_k_integrand_imag = aimag(Sigma_k_integrand_complex)

    end function Sigma_k_integrand_imag

    ! --------------------------M_i integrands--------------------------
    ! Eq. (2.30)

    function M_integrandPV(E)
        implicit none
        real(DP), intent(in) :: E
        real(DP) :: M_integrandPV
        ! Local variables
        real(DP) :: q, dqdE
        real(DP), dimension(1,1) :: f_i11

        q = sqrt( E**2/4.0_DP + (m_bar(ich)**2 - m_mes(ich)**2)**2 &
            & /4.0_DP/E**2 - (m_mes(ich)**2 + m_bar(ich)**2)/2.0_DP )
        dqdE = sqrt(q**2 + m_mes(ich)**2)*sqrt(q**2 + m_bar(ich)**2) / q / E
        f_i11 = f_i(q,ich)
        M_integrandPV = -q**2 * f_i11(1,1)**2 * dqdE

    end function M_integrandPV

    function M_integrand(E)
        implicit none
        real(DP), intent(in) :: E
        real(DP) :: M_integrand

        M_integrand = M_integrandPV(E) / (E - E_gl)

    end function M_integrand

    function M_k_integrand_real(k_in)
        implicit none
        real(DP), intent(in) :: k_in
        real(DP) :: M_k_integrand_real
        ! Local variables
        complex(DP) :: M_k_integrand_cmplx
        complex(DP) :: q, dqdE, E
        complex(DP), dimension(1,1) :: f_i11

        q = k_in * k_complex_theta
        E = sqrt(q**2 + m_mes(ich)**2) + sqrt(q**2 + m_bar(ich)**2)
        f_i11 = f_i(q,ich)
        M_k_integrand_cmplx = -q**2 * f_i11(1,1)**2 / (E - E_pole) &
            & * k_complex_theta
        M_k_integrand_real = real(M_k_integrand_cmplx,DP)

    end function M_k_integrand_real

    function M_k_integrand_imag(k_in)
        implicit none
        real(DP), intent(in) :: k_in
        real(DP) :: M_k_integrand_imag
        ! Local variables
        complex(DP) :: M_k_integrand_cmplx
        complex(DP) :: q, E
        complex(DP), dimension(1,1) :: f_i11

        q = k_in * k_complex_theta
        E = sqrt(q**2 + m_mes(ich)**2) + sqrt(q**2 + m_bar(ich)**2)
        f_i11 = f_i(q,ich)
        M_k_integrand_cmplx = -q**2 * f_i11(1,1)**2 / (E - E_pole) &
            & * k_complex_theta
        M_k_integrand_imag = aimag(M_k_integrand_cmplx)

    end function M_k_integrand_imag


    ! --------------------------g_f integrands--------------------------
    ! Eq. (2.39)

    function g_f_integrandPV(E)
        implicit none
        real(DP), intent(in) :: E
        real(DP) :: g_f_integrandPV
        ! Local variables
        real(DP) :: q, dqdE
        real(DP), dimension(1,1) :: f_i11

        q = sqrt( E**2/4.0_DP + (m_bar(ich)**2 - m_mes(ich)**2)**2 &
            & /4.0_DP/E**2 - (m_mes(ich)**2 + m_bar(ich)**2)/2.0_DP )
        dqdE = sqrt(q**2 + m_mes(ich)**2)*sqrt(q**2 + m_bar(ich)**2) / q / E
        f_i11 = f_i(q,ich)
        g_f_integrandPV = -q**2 * f_i11(1,1) &
            & * gBare(ich,ibare) * g_i(q,ich,ibare) * dqdE

    end function g_f_integrandPV

    function g_f_integrand(E)
        implicit none
        real(DP), intent(in) :: E
        real(DP) :: g_f_integrand

        g_f_integrand = g_f_integrandPV(E) / (E - E_gl)

    end function g_f_integrand

    function g_f_k_integrand_real(k_in)
        implicit none
        real(DP), intent(in) :: k_in
        real(DP) :: g_f_k_integrand_real
        ! Local variables
        complex(DP) :: q, dqdE, E
        complex(DP) :: g_f_k_integrand_cmplx
        complex(DP), dimension(1,1) :: f_i11

        q = k_in * k_complex_theta
        E = sqrt(q**2 + m_mes(ich)**2) + sqrt(q**2 + m_bar(ich)**2)
        f_i11 = f_i(q,ich)
        g_f_k_integrand_cmplx = -1.0_DP * q**2 * f_i11(1,1) &
            & * gBare(ich,ibare) * g_i(q,ich,ibare) / (E - E_pole) &
            & * k_complex_theta
        g_f_k_integrand_real = real(g_f_k_integrand_cmplx,DP)

    end function g_f_k_integrand_real

    function g_f_k_integrand_imag(k_in)
        implicit none
        real(DP), intent(in) :: k_in
        real(DP) :: g_f_k_integrand_imag
        ! Local variables
        complex(DP) :: q, dqdE, E
        complex(DP) :: g_f_k_integrand_cmplx
        complex(DP), dimension(1,1) :: f_i11

        q = k_in * k_complex_theta
        E = sqrt(q**2 + m_mes(ich)**2) + sqrt(q**2 + m_bar(ich)**2)
        f_i11 = f_i(q,ich)
        g_f_k_integrand_cmplx = -1.0_DP * q**2 * f_i11(1,1) &
            & * gBare(ich,ibare) * g_i(q,ich,ibare) / (E - E_pole) &
            & * k_complex_theta
        g_f_k_integrand_imag = aimag(g_f_k_integrand_cmplx)

    end function g_f_k_integrand_imag


    function atanRangeFixer(z1, z2) result(fix)
        ! No recollection of how this works unfortunateky
        ! Must've come to me in a dream
        implicit none
        complex(DP), intent(in) :: z1, z2
        real(DP) :: fix

        integer :: Q1, Q2
        real(DP) :: theta1, theta2, dtheta

        theta1 = atan2(aimag(z1),real(z1))
        if (theta1.lt.0.0_DP) theta1 = theta1 + 2.0_DP*pi
        theta2 = atan2(aimag(z2),real(z2))
        if (theta2.lt.0.0_DP) theta2 = theta2 + 2.0_DP*pi
        dtheta = theta2 - theta1

        if (globalPoint .eq. openChPoint(atan_ch)+1) then
            if (dtheta.ge.0.0_DP) then
                isClockwise(atan_ch) = .false.
            else
                isClockwise(atan_ch) = .true.
            end if
        else

            if (abs(dtheta).gt.pi) then
                if (dtheta.ge.0.0_DP) then
                    isClockwise(atan_ch) = .true.
                else
                    isClockwise(atan_ch) = .false.
                end if
            end if
        end if

        Q1 = floor(theta1/(pi/2.0_DP) + 1.0_DP)
        Q2 = floor(theta2/(pi/2.0_DP) + 1.0_DP)
        fix = 0.0_DP

        if (Q2.eq.1) then
            if (isClockwise(atan_ch)) fix = -2.0_DP*pi
        else if (Q2.eq.2) then
            if (isClockwise(atan_ch)) then
                fix = -pi
            else
                fix = pi
            end if
        else if (Q2.eq.3) then
            if (isClockwise(atan_ch)) then
                fix = -pi
            else
                fix = pi
            end if
        else if (Q2.eq.4) then
            if (.not.isClockwise(atan_ch)) then
                fix = 2.0_DP*pi
            end if
        end if

    end function atanRangeFixer


end module SMatrix
