module initdata_module
  use amrex_amr_module
  use my_amr_module
  use plotfile_module, only: write_plotfile
  use averagedown_module, only: averagedown

  implicit none

  private

  public :: initdata

contains

  subroutine initdata()
    call amrex_init_from_scratch(0.0_amrex_real)
    call averagedown()

    if (plot_int .gt. 0) call write_plotfile
  end subroutine initdata

end module initdata_module
