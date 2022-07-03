module bc_module
  use amrex_base_module
  implicit none
  ! boundary conditions
  ! second parameter is number of components
  ! default to periodic boundary, see amrex_bc_types_module for bc list.
  integer, save :: lo_bc(amrex_spacedim, 1) = amrex_bc_ext_dir
  integer, save :: hi_bc(amrex_spacedim, 1) = amrex_bc_ext_dir
end module bc_module
