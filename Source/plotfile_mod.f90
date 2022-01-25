module plotfile_module

  use amrex_amr_module
  use amr_data_module, only: phi_new, t_new
  use my_amr_module, only: plot_file, stepno

  implicit none

  private
  public :: write_plotfile

contains

  subroutine write_plotfile()
    integer :: nlevs
    character(len=127) :: name
    character(len=16) :: current_step
    type(amrex_string) :: varname(1)

    write (current_step, fmt='(i5.5)') stepno(0)
    name = trim(plot_file)//current_step

    nlevs = amrex_get_numlevels()

    call amrex_string_build(varname(1), "phi")
    call amrex_write_plotfile(name, nlevs, phi_new, varname, amrex_geom, &
                              t_new(0), stepno, amrex_ref_ratio)

  end subroutine write_plotfile
end module plotfile_module
