!> @brief Quelques algorithmes de tri des entiers et des réels.
!> @date Time-stamp: <2017-10-26 03:08:16 yves>
!> @author Yves.Coudiere@univ-nantes.fr
!> @details Tri par insertion et shell
!> @date Septembre 2015 mise à plat
!> - changement des noms des procédures (notation avec des _)
!> - documentation doxygen
module sort

  use prec
  implicit none
  
  private

  !> @brief Tri par insertion, ordre décroissant
  interface sort_insertion
     module procedure sort_insertion_i,sort_insertion_f
  end interface

  !> @brief Tri par shell, ordre décroissant. Plus rapide que par
  !> insertion sur les tableaux de taille moyenne à grande, moins
  !> rapide que quicksort sur les grand tableaux.
  interface sort_shell
     module procedure sort_shell_i,sort_shell_f
  end interface

  public :: sort_insertion, sort_shell

contains

  !> @brief Tri par insertion (wikipedia : le plus rapide pour moins
  !> de 15 elts), ordre décroissant. Le tableau ne change pas, mais on
  !> renvoie les indices à suivre pour aller dans l'ordre décroissant.
  subroutine sort_insertion_i(t,new_i)

    implicit none
    integer, dimension(:),intent(in)   :: t
    integer, dimension(:), intent(out) :: new_i

    integer :: n,i,j,z

    n = size(t)
    i = 0
    new_i = (/ (i, i=1,n) /)

    do i=2,n
       z = t(i)                    ! Valeur à insérer
       j = i-1
       ! Pas de boucle do car on ne connait pas la valeur de j en sortie
       do while (j>=1)
          if (z>t(new_i(j))) then 
             new_i(j+1)=new_i(j)   ! Decalage
             j = j-1
          else
             exit
          end if
       end do
       new_i(j+1)=i                ! Insertion
    end do

  end subroutine sort_insertion_i


  !> @brief Tri 'ShellSort' (wikipedia : le plus rapide pour moins de
  !> 15 elts), ordre décroissant. En pratique le plus rapide sur des
  !> tableaux de quelques centaines ? Le tableau ne change pas, mais
  !> on renvoie les indices à suivre pour aller dans l'ordre
  !> croissant.
  subroutine sort_shell_i(t,new_i)

    implicit none
    integer, dimension(:),intent(in)   :: t
    integer, dimension(:), intent(out) :: new_i

    integer :: n,i,j,z,gap,i_z

    n = size(t)
    i = 0
    new_i = (/ (i, i=1,n) /)

    ! Recherche du Gap optimal, defini par U(n) = 3U(n-1)+1 et U(n)<n
    gap = 0
    do while (gap<n)
       gap = 3*gap+1
    end do

    do while (gap>2)

       gap = gap/3                          ! Calcul du gap
!print *,'GAP=',gap
       ! Tri par insertion

       do i=1+gap,n
!print "('i=',I3)",i
          i_z = new_i(i)                    ! Indice de la valeur a inserer
          z = t(i_z)                        ! Valeur à insérer
!print "('valeur a inserer = ',I3)",z
          j = i-gap
          do while (j>=1)
!print "('test avec j=',I3)",j
             if (z>t(new_i(j))) then 
                new_i(j+gap)=new_i(j)       ! Decalage
                j = j-gap
             else
                exit
             end if
          end do
          new_i(j+gap)=i_z                  ! Insertion
!print "('i=',I3,' j=',I3,' new_i=',20I3)",i,j,new_i
       end do
!print "(6X,' FIN  new_i=',20I3)",new_i
    end do

  end subroutine sort_shell_i
  

  !> @brief Tri par insertion (wikipedia : le plus rapide pour moins
  !> de 15 elts), ordre décroissant. Le tableau ne change pas, mais on
  !> renvoie les indices à suivre pour aller dans l'ordre croissant.
  subroutine sort_insertion_f(t,new_i)

    implicit none
    real(pr), dimension(:),intent(in)   :: t
    integer, dimension(:), intent(out) :: new_i

    integer :: n,i,j
    real(pr) :: z

    n = size(t)
    i = 0
    new_i = (/ (i, i=1,n) /)

    do i=2,n
       z = t(i)                    ! Valeur à insérer
       j = i-1
       ! Pas de boucle do car on ne connait pas la valeur de j en sortie
       do while (j>=1)
          if (z>t(new_i(j))) then 
             new_i(j+1)=new_i(j)   ! Decalage
             j = j-1
          else
             exit
          end if
       end do
       new_i(j+1)=i                ! Insertion
    end do

  end subroutine sort_insertion_f


  !> @brief Tri 'ShellSort' (wikipedia : le plus rapide pour moins de
  !> 15 elts), ordre décroissant. En pratique le plus rapide sur des
  !> tableaux de quelques centaines ? Le tableau ne change pas, mais
  !> on renvoie les indices à suivre pour aller dans l'ordre
  !> croissant.
  subroutine sort_shell_f(t,new_i)

    implicit none
    real(pr), dimension(:),intent(in)   :: t
    integer, dimension(:), intent(out) :: new_i

    integer :: n,i,j,gap,i_z
    real(pr) :: z

    n = size(t)
    i = 0
    new_i = (/ (i, i=1,n) /)

    ! Recherche du Gap optimal, U(n) = 3U(n-1)+1
    gap = 0
    do while (gap<n)
       gap = 3*gap+1
    end do

    do while (gap>2)

       gap = gap/3

       do i=1+gap,n
          i_z = new_i(i)               ! Indice de la valeur a inserer
          z = t(i_z)                  ! Valeur à insérer
          j = i-gap
          do while (j>0)
             if (z>t(new_i(j))) then
                new_i(j+gap)=new_i(j) ! Decalage
                j = j-gap
             else
                exit
             end if
          end do
          new_i(j+gap)=i_z
       end do

    end do

  end subroutine sort_shell_f
  
end module sort
