! wrapper of amrex functions
module fillpatch_module
  use iso_c_binding
  use amrex_amr_module
  use bc_module, only: lo_bc, hi_bc
  use amrex_fort_module, only: rt => amrex_real
  implicit none

  private

  public :: fillpatch, fillcoarsepatch

contains
  ! Fill phi with data from phi_old and phi_new of current level and one level below.
  ! use interpolater conservative linear cell interp.
  ! Fill boundary as well as ghost cells
  subroutine fillpatch(lev, time, phi)
    use amr_data_module, only: t_old, t_new, phi_old, phi_new
    integer, intent(in) :: lev
    real(rt), intent(in) :: time
    type(amrex_multifab), intent(inout) :: phi
    !specify the components number
    integer, parameter :: src_comp = 1, dst_comp = 1, num_comp = 1

    if (lev .eq. 0) then
      call amrex_fillpatch(phi, t_old(lev), phi_old(lev), &
                           t_new(lev), phi_new(lev), &
                           amrex_geom(lev), fill_physbc, &
                           time, src_comp, dst_comp, num_comp)
    else
      call amrex_fillpatch(phi, t_old(lev - 1), phi_old(lev - 1), &
                           t_new(lev - 1), phi_new(lev - 1), &
                           amrex_geom(lev - 1), fill_physbc, &
                           t_old(lev), phi_old(lev), &
                           t_new(lev), phi_new(lev), &
                           amrex_geom(lev), fill_physbc, &
                           time, src_comp, dst_comp, num_comp, &
                           amrex_ref_ratio(lev - 1), amrex_interp_cell_cons, &
                           lo_bc, hi_bc)
    end if
  end subroutine fillpatch

  ! Fill patch only with data from one level below
  subroutine fillcoarsepatch(lev, time, phi)
    use amr_data_module, only: t_old, t_new, phi_old, phi_new
    integer, intent(in) :: lev
    real(rt), intent(in) :: time
    type(amrex_multifab), intent(inout) :: phi
    !specify the components number
    integer, parameter :: src_comp = 1, dst_comp = 1, num_comp = 1

    call amrex_fillcoarsepatch(phi, t_old(lev - 1), phi_old(lev - 1), &
                               t_new(lev - 1), phi_new(lev - 1), &
                               amrex_geom(lev - 1), fill_physbc, &
                               amrex_geom(lev), fill_physbc, &
                               time, src_comp, dst_comp, num_comp, &
                               amrex_ref_ratio(lev - 1), amrex_interp_cell_cons, &
                               lo_bc, hi_bc)
  end subroutine fillcoarsepatch

  ! How to fill boudary  cells
  subroutine fill_physbc(pmf, scomp, ncomp, time, pgeom) bind(c)
    use amrex_geometry_module, only: amrex_is_all_periodic
    use amrex_filcc_module, only: amrex_filcc

    type(c_ptr), value :: pmf, pgeom
    integer(c_int), value :: scomp, ncomp
    real(rt), value :: time

    type(amrex_geometry) :: geom
    type(amrex_multifab) :: mf
    type(amrex_mfiter) :: mfi
    real(rt), contiguous, pointer :: p(:, :, :, :)
    integer :: p_lo(4), p_hi(4)

    if (.not. amrex_is_all_periodic()) then
      geom = pgeom
      mf = pmf

      !$omp parallel private(mfi, p, p_lo, p_hi)
      call amrex_mfiter_build(mfi, mf, tiling=.false.)
      do while (mfi%next())
        p => mf%dataptr(mfi)
        ! part of this box is outside of the domain
        if (.not. geom%domain%contains(p)) then
          p_lo = lbound(p)
          p_hi = ubound(p)
          call amrex_filcc(p, p_lo, p_hi, &
                           geom%domain%lo, geom%domain%hi, &
                           geom%dx, &
                           geom%get_physical_location(p_lo), &
                           lo_bc, hi_bc)
        end if
      end do
      !$omp end parallel
    end if
  end subroutine fill_physbc
end module fillpatch_module
