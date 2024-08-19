module hamiltonian
  use kinds
  use heft
  implicit none

contains

    subroutine generateHamiltonian(H, omega, k, L)
        implicit none
        real(DP), dimension(:,:), intent(out) :: H
        real(DP), dimension(:), intent(out)   :: omega
        real(DP), dimension(:,:) :: k
        real(DP), intent(in) :: L

        ! ith, jth elements of H
        integer :: ii, jj
        integer :: ch_i, ch_j ! loop over channels
        integer :: bare_i, bare_j ! loop over bare states
        integer :: ik, jk ! loop over momentum
        real(DP) :: k_ik, k_jk ! store momentum in these for quick access
        integer :: N_H
        real(DP) :: g_i_temp, V_ij_temp

        N_H = size(H,1)
        H(:,:) = 0.0_DP
        omega(:) = 0.0_DP

        ! Set the first diagonal(s) to bare state mass(es)
        if (n_bare.gt.0) then
            omega(:n_bare) = m_bare(:)
            do ii = 1,n_bare
               H(ii,ii) = omega(ii)
            end do
        end if

        do jj = 1,N_H
           ! channel index in j'th position
           ch_j = mod(jj-1-n_bare,n_ch) + 1
           ! Momentum index in j'th position
           jk = (jj-1-n_bare)/n_ch + 1
           k_jk = k(jk, ch_j)
           ! bare index in j'th pos
           if (jj.le.n_bare) bare_j = jj

           do ii = 1,N_H
              ! Current channel in i'th position
              ch_i = mod(ii-1-n_bare,n_ch) + 1
              ! Momentum index in i'th position
              ik = (ii-1-n_bare)/n_ch + 1
              k_ik = k(ik, ch_i)

              ! Non-interacting Hamiltonian diagonals (H_0)
              if (ii.eq.jj .and. ii.gt.n_bare) then
                  omega(ii) = sqrt(k_ik**2 + m_mes(ch_i)**2) &
                      & + sqrt(k_ik**2 + m_bar(ch_i)**2)
                  H(ii,jj) = omega(ii)
              end if

              ! Bare -> channel couplings e.g. Delta->N+pi
              if (jj.le.n_bare .and. ii.gt.n_bare .and. n_bare.gt.0) then
                  g_i_temp = gBare(ch_i,bare_j) &
                      & * g_i(k_ik, ch_i, bare_j)
                  ! Finite-volume factors: Eq. (2.75)
                  g_i_temp = g_i_temp * sqrt(C_3packed(ik)/4.0_DP/pi) &
                      & * (2.0_DP*pi*hbar_c/L)**(1.5_DP)
                  H(ii,jj) = g_i_temp
                  H(jj,ii) = g_i_temp
              end if

              ! Channel -> channel couplings e.g. N+pi->Delta+pi
              if (ii.gt.n_bare .and. jj.gt.n_bare) then
                  V_ij_temp = V_ij(k_ik, ch_i, k_jk, ch_j)
                  ! Eq. (2.76)
                  V_ij_temp = sqrt(real(C_3packed(ik)*C_3packed(jk),DP)) &
                      & /4.0_DP/pi * (2.0_DP*pi*hbar_c/L)**3 &
                      & * V_ij_temp
                  H(ii,jj) = H(ii,jj) + V_ij_temp
              end if

           end do
        end do

    end subroutine generateHamiltonian

end module hamiltonian
