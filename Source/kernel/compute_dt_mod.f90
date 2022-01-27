module compute_dt_module
  use amrex_amr_module
  use amrex_fort_module, only: rt => amrex_real
  implicit none
  private
  public :: compute_dt

contains

  subroutine compute_dt()
    use my_amr_module, only: t_new, dt, stop_time, nsubsteps

    integer :: lev, nlevs, n_factor
    real(rt) :: dt_0, eps
    real(rt), allocatable :: dt_tmp(:)
    real(rt), parameter :: change_max = 1.1_rt

    nlevs = amrex_get_numlevels()
    allocate (dt_tmp(0:nlevs - 1))
    do lev = 0, nlevs - 1
      dt_tmp(lev) = est_timestep(lev, t_new(lev))
    end do
    call amrex_parallel_reduce_min(dt_tmp, nlevs)

    dt_0 = dt_tmp(0)
    n_factor = 1
    do lev = 0, nlevs - 1
      dt_tmp(lev) = min(dt_tmp(lev), change_max*dt(lev))
      n_factor = n_factor*nsubsteps(lev)
      dt_0 = min(dt_0, n_factor*dt_tmp(lev))
    end do

    eps = 1.e-3_rt*dt_0
    if (t_new(0) + dt_0 > stop_time - eps) then
      dt_0 = stop_time - t_new(0)
    end if

    dt(0) = dt_0
    do lev = 1, nlevs - 1
      dt(lev) = dt(lev - 1)/nsubsteps(lev)
    end do
  end subroutine compute_dt

  function est_timestep(lev, time) result(dt)
    use my_amr_module, only: phi_new, cfl
    use face_velocity_module, only: get_face_velocity

    real(rt) :: dt
    integer, intent(in) :: lev
    real(rt), intent(in) :: time

    real(rt) :: dt_est, umax
    type(amrex_box) :: bx
    type(amrex_fab) :: u
    type(amrex_mfiter) :: mfi
    real(rt), contiguous, pointer :: p(:, :, :, :)

    dt_est = huge(1._rt)

    !$omp parallel reduction(min:dt_est) private(umax, bx, u, mfi, p)
    call u%reset_omp_private()
    call amrex_mfiter_build(mfi, phi_new(lev), tiling=.true.)
    do while (mfi%next())
      bx = mfi%nodaltilebox()
      call u%resize(bx, amrex_spacedim)
      p => u%dataptr()

      call get_face_velocity(time, &
                             p(:, :, :, 1), bx%lo, bx%hi, &
                             p(:, :, :, 2), bx%lo, bx%hi, &
                             amrex_geom(lev)%dx, amrex_problo)

      umax = u%norminf(1, 1)
      if (umax > 1.e-100_rt) then
        dt_est = min(dt_est, amrex_geom(lev)%dx(1)/umax)
      end if

      umax = u%norminf(2, 1)
      if (umax > 1.e-100_rt) then
        dt_est = min(dt_est, amrex_geom(lev)%dx(2)/umax)
      end if

    end do
    call amrex_mfiter_destroy(mfi)
    call amrex_fab_destroy(u)
    !$omp end parallel

    dt = dt_est*cfl

  end function est_timestep

end module compute_dt_module
