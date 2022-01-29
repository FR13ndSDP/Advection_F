module tagging_module

  use iso_c_binding
  use amrex_fort_module, only: rt => amrex_real

  implicit none
  private

  public :: tag_phi_error, tag_phi_grad

contains

  subroutine tag_phi_error(lo, hi, phi, philo, phihi, tag, taglo, taghi, phierr, settag, cleartag)
    integer, intent(in) :: lo(3), hi(3), philo(4), phihi(4), taglo(4), taghi(4)
    real(rt), intent(in) :: phi(philo(1):phihi(1), philo(2):philo(2), philo(3):phihi(3))
    character(kind=c_char), intent(inout) :: tag(taglo(1):taghi(1), taglo(2):taghi(2), taglo(3):taghi(3))

    real(rt), intent(in) :: phierr
    character(kind=c_char), intent(in) :: settag, cleartag

    integer :: i, j, k

    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          if (phi(i, j, k) >= phierr) then
            tag(i, j, k) = settag
          else
            tag(i, j, k) = cleartag
          end if
        end do
      end do
    end do
  end subroutine tag_phi_error

  subroutine tag_phi_grad(lo, hi, phi, philo, phihi, tag, taglo, taghi, phierr, settag, cleartag)
    integer, intent(in) :: lo(3), hi(3), philo(4), phihi(4), taglo(4), taghi(4)
    real(rt), intent(in) :: phi(philo(1):phihi(1), philo(2):philo(2), philo(3):phihi(3))
    character(kind=c_char), intent(inout) :: tag(taglo(1):taghi(1), taglo(2):taghi(2), taglo(3):taghi(3))

    real(rt), intent(in) :: phierr
    character(kind=c_char), intent(in) :: settag, cleartag

    integer :: i, j, k
    real(rt) :: ax, ay, az

    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          ax = abs(phi(i + 1, j, k) - phi(i, j, k))
          ay = abs(phi(i, j + 1, k) - phi(i, j, k))
          az = abs(phi(i, j, k + 1) - phi(i, j, k))

          ax = max(ax, abs(phi(i, j, k) - phi(i - 1, j, k)))
          ay = max(ay, abs(phi(i, j, k) - phi(i, j - 1, k)))
          az = max(az, abs(phi(i, j, k) - phi(i, j, k - 1)))

          if (max(ax, ay, az) >= phierr) then
            tag(i, j, k) = settag
          else
            tag(i, j, k) = cleartag
          end if
        end do
      end do
    end do
  end subroutine tag_phi_grad
end module tagging_module
