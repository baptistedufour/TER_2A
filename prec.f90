!> @author Y. Coudière
!> @date Time-stamp: <2014-05-07 22:04:55 yves>
!> @brief Module de définition de la précision des réels.
!>
!> @details Précision des nombres à vrigule flottante et différentes
!>          constantes numériques universelles (pi, etc).
module prec
  implicit none

  public

  !> paramètre pour déclarer les réels tous avec la même précision
  !> ailleurs dans le code.
  integer, parameter :: pr = kind(1.D0) ! selected_real_kind(15)
  !> le paramètre de la précision de Matlab.
  integer, parameter :: matlab_pr = kind(1.D0) ! selected_real_kind(15)
  !> écart entre 1. et le nombre à virgule flottante suivant
  real(pr), parameter :: eps = epsilon(1.0_pr)

  !> Pi
  real(pr), parameter :: pi = 3.141592653589793238462643383279502884_pr

end module prec
