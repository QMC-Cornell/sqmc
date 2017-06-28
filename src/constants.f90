module constants
  use types, only: rk
  implicit none
  save
  private
  public :: pi

  real(kind = rk), parameter ::            pi = 4.0_rk * atan(1.0_rk)

end module constants
