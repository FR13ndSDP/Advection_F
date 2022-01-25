module face_velocity_module
  implicit none
  private
  public :: get_face_velocity
contains

  ! pro_lo is the physics location
  subroutine get_face_velocity(time, &
                               vx, vxlo, vxhi, &
                               vy, vylo, vyhi, &
                               dx, prob_lo)
    use amrex_base_module

    real(amrex_real), intent(in) :: time
    integer, intent(in) :: vxlo(2), vxhi(2), vylo(2), vyhi(2)
    real(amrex_real), intent(out) :: vx(vxlo(1):vxhi(1), vxlo(2):vxhi(2))
    real(amrex_real), intent(out) :: vy(vylo(1):vyhi(1), vylo(2):vyhi(2))
    real(amrex_real), intent(in) :: dx(2), prob_lo(2)

    integer :: i, j, plo(2), phi(2)
    real(amrex_real) :: x, y
    real(amrex_real), pointer, contiguous :: psi(:, :)
    real(amrex_real), parameter :: M_PI = 3.1415926536
    real(amrex_real) :: qtrdx(2)
    qtrdx = 0.25d0/dx

    plo(1) = min(vxlo(1) - 1, vylo(1) - 1)
    plo(2) = min(vxlo(2) - 1, vylo(2) - 1)
    phi(1) = max(vxhi(1), vyhi(1) + 1)
    phi(2) = max(vxhi(2) + 1, vyhi(2))

    call amrex_allocate(psi, plo, phi)

    do j = plo(2), phi(2)
      y = (dble(j) + 0.5d0)*dx(2) + prob_lo(2)
      do i = plo(1), phi(1)
        x = (dble(i) + 0.5d0)*dx(1) + prob_lo(1)
        ! streamfunction psi
        psi(i, j) = sin(M_PI*x)**2*sin(M_PI*y)**2*cos(M_PI*time/2.d0)*(1.d0/M_PI)
      end do
    end do

    do j = vxlo(2), vxhi(2)
      do i = vxlo(1), vxhi(1)
        ! x velocity
        vx(i, j) = -((psi(i, j + 1) + psi(i - 1, j + 1)) - (psi(i, j - 1) + psi(i - 1, j - 1))) * qtrdx(2)
      end do
    end do

    do j = vylo(2), vyhi(2)
      do i = vylo(1), vyhi(1)
        ! y velocity
        vy(i, j) = ((psi(i + 1, j) + psi(i + 1, j - 1)) - (psi(i - 1, j) + psi(i - 1, j - 1))) * qtrdx(1)
      end do
    end do

    call amrex_deallocate(psi)
  end subroutine get_face_velocity

end module face_velocity_module
