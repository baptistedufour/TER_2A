program test_csc

  use prec
  use mat_list
  use mat_csc
!  use iso_c_binding
  use umfpack_calls
  
  ! We want to test how to factorize and solve a linear system in the csc
  ! format with umfpack.

  ! The matrix:
  ! A = [4 2 0  1]
  !     [1 4 1  0]
  !     [0 0 1 -1]
  !     [0 0 3  6]
  !
  ! It CSC arrays
  !
  ! col = 1   3   5     8      11
  ! row = 1 2 1 2 2 3 4 1  3 4
  ! val = 4 1 2 4 1 1 3 1 -1 6

  integer :: n = 4
  integer, dimension(5), target :: col = [1,3,5,8,11]-1
  integer, dimension(10), target :: row = [1,2,1,2,2,3,4,1,3,4]-1
  real(pr), dimension(10), target :: val = &
       & [4._pr,1._pr,2._pr,4._pr,1._pr,1._pr,3._pr,1._pr,-1._pr,6._pr]
  type(c_ptr) :: ptr_symbolic, ptr_numeric

  real(pr), dimension(4), target :: x,b = [7._pr, 6._pr, 0._pr, 9._pr]

  type(m_list) :: M_l
  type(csc) :: M_csc
  
  integer :: status

  print"('Using the mat_csc module')"
  ! A = [4 2 0  1]
  !     [1 4 1  0]
  !     [0 0 1 -1]
  !     [0 0 3  6]
  call zeros(M_l, n,n)
  call mat_add_value(M_l, 1,1, 4._pr)
  call mat_add_value(M_l, 1,2, 2._pr)
  call mat_add_value(M_l, 1,4, 1._pr)
  call mat_add_value(M_l, 2,1, 1._pr)
  call mat_add_value(M_l, 2,2, 4._pr)
  call mat_add_value(M_l, 2,3, 1._pr)
  call mat_add_value(M_l, 3,3, 1._pr)
  call mat_add_value(M_l, 3,4,-1._pr)
  call mat_add_value(M_l, 4,3, 3._pr)
  call mat_add_value(M_l, 4,4, 6._pr)
  call print(M_l)
  M_csc = M_l
  call print(M_csc)
  call factorize(M_csc)

  print "('Using explicit calls to umfpack')"
  status = umfpack_di_symbolic(n,n, c_loc(col), c_loc(row), c_loc(val), &
       ptr_symbolic, c_null_ptr,c_null_ptr)
  print *,'di_symbolic : ',status
  status = umfpack_di_numeric(c_loc(col), c_loc(row), c_loc(val), &
       ptr_symbolic, ptr_numeric, c_null_ptr,c_null_ptr)
  print *,'di_numeric : ',status
  call umfpack_di_free_symbolic(ptr_symbolic)

  status = umfpack_di_solve(UMFPACK_A, c_loc(col), c_loc(row), c_loc(val), &
       & c_loc(x), c_loc(b), ptr_numeric, c_null_ptr,c_null_ptr)
  print *,'di_solve : ',status
  call umfpack_di_free_numeric(ptr_numeric)

  print "(6E12.5)", x

end program test_csc
