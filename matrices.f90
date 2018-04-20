!> @brief Paramètres communs à tous les types de matrices
!> @date Time-stamp: <2014-04-04 20:11:54 yves>
!
!> @author Y. Coudière
!> - Code démarré vers 2006
!> @date mai 2013 mise à plat
module matrices

  use prec
  implicit none

  private

  integer, parameter, public :: mat_add = 0 ! paramètres de construction
  integer, parameter, public :: mat_set = 1 ! d'une matrice

end module matrices

