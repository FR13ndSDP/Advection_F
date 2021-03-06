module advect_module
  use iso_fortran_env, only : stderr => error_unit
  use amrex_base_module
  implicit none
  private

  public :: advect

contains

  subroutine advect(lo, hi, &
                    uin, ui_lo, ui_hi, &
                    uout, uo_lo, uo_hi, &
                    vx, vx_lo, vx_hi, &
                    vy, vy_lo, vy_hi, &
                    flxx, fx_lo, fx_hi, &
                    flxy, fy_lo, fy_hi, &
                    dx, dt)

    use compute_flux_module

    integer, intent(in) :: lo(2), hi(2), ui_lo(2), ui_hi(2), uo_lo(2), uo_hi(2), &
                           vx_lo(2), vx_hi(2), vy_lo(2), vy_hi(2), &
                           fx_lo(2), fx_hi(2), fy_lo(2), fy_hi(2)
    real(amrex_real), intent(in) :: dx(2), dt
    real(amrex_real), intent(in) :: uin(ui_lo(1):ui_hi(1), ui_lo(2):ui_hi(2))
    real(amrex_real), intent(inout) :: uout(uo_lo(1):uo_hi(1), uo_lo(2):uo_hi(2))
    real(amrex_real), intent(in) :: vx(vx_lo(1):vx_hi(1), vx_lo(2):vx_hi(2))
    real(amrex_real), intent(in) :: vy(vy_lo(1):vy_hi(1), vy_lo(2):vy_hi(2))
    real(amrex_real), intent(out) :: flxx(fx_lo(1):fx_hi(1), fx_lo(2):fx_hi(2))
    real(amrex_real), intent(out) :: flxy(fy_lo(1):fy_hi(1), fy_lo(2):fy_hi(2))

    integer :: i, j
    integer :: glo(2), ghi(2)
    real(amrex_real) :: dtdx(2), umax, vmax

    real(amrex_real), dimension(:, :), pointer, contiguous :: phix_1d, phiy_1d, phix, phiy, slope
    dtdx = dt/dx
    glo = lo - 1
    ghi = hi + 1

    call amrex_allocate(phix_1d, glo, ghi)
    call amrex_allocate(phiy_1d, glo, ghi)
    call amrex_allocate(phix, glo, ghi)
    call amrex_allocate(phiy, glo, ghi)
    call amrex_allocate(slope, glo, ghi)

    ! check for CFL condition
    umax = maxval(vx)
    vmax = maxval(vy)
    if (umax*dt >= dx(1) .or. vmax*dt >= dx(2)) then
      write(stderr, *) "umax = ", umax, ", vmax = ", vmax, ", dt = ", dt, ", dx = ", dx
      call amrex_error("CFL violation")
    end if

    call compute_flux(lo, hi, dt, dx, &
                      uin, ui_lo, ui_hi, &
                      vx, vx_lo, vx_hi, &
                      vy, vy_lo, vy_hi, &
                      flxx, fx_lo, fx_hi, &
                      flxy, fy_lo, fy_hi, &
                      phix_1d, phiy_1d, phix, phiy, slope, glo, ghi)

    do j = lo(2), hi(2)
      do i = lo(1), hi(1)
        uout(i, j) = uin(i, j) &
                     + (flxx(i, j) - flxx(i + 1, j))*dtdx(1) &
                     + (flxy(i, j) - flxy(i, j + 1))*dtdx(2)
      end do
    end do

    do j = lo(2), hi(2)
      do i = lo(1), hi(1)+1
        flxx(i,j) = flxx(i,j) * (dt * dx(2))
      end do
    end do

    do j = lo(2), hi(2) + 1
      do i = lo(1), hi(1)
        flxy(i,j) = flxy(i,j) * (dt * dx(1))
      end do
    end do

    call amrex_deallocate(phix_1d)
    call amrex_deallocate(phiy_1d)
    call amrex_deallocate(phix)
    call amrex_deallocate(phiy)
    call amrex_deallocate(slope)

  end subroutine advect

end module advect_module
