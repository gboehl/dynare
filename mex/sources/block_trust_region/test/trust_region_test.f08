program trust_region_test
  use trust_region
  use iso_fortran_env
  implicit none (type, external)

  real(real64), dimension(2), parameter :: rosenbrock_guess = [ -10**4, 1 ]
  real(real64), dimension(2), parameter :: rosenbrock_solution = [ 1, 1 ]
  real(real64), dimension(2) :: rosenbrock_x

  real(real64), dimension(4), parameter :: powell1_guess = [ 3, -1, 0, 1 ]
  real(real64), dimension(4), parameter :: powell1_solution = 0
  real(real64), dimension(4) :: powell1_x

  integer :: info

  rosenbrock_x = rosenbrock_guess
  call trust_region_solve(rosenbrock_x, rosenbrock, info)
  if (info /= 1 .or. maxval(abs(rosenbrock_x - rosenbrock_solution)) > 1e-11_real64) &
       error stop "Failed to solve rosenbrock"

  powell1_x = powell1_guess
  call trust_region_solve(powell1_x, powell1, info, tolf = 1e-8_real64)
  if (info /= 1 .or. maxval(abs(rosenbrock_x - rosenbrock_solution)) > 4e-5_real64) &
       error stop "Failed to solve powell1"
contains
  subroutine rosenbrock(x, fvec, fjac)
    real(real64), dimension(:), intent(in) :: x
    real(real64), dimension(size(x)), intent(out) :: fvec
    real(real64), dimension(size(x),size(x)), intent(out), optional :: fjac

    fvec(1) = 1 - x(1)
    fvec(2) = 10*(x(2)-x(1)**2)

    if (present(fjac)) then
       fjac(1,1) = -1
       fjac(1,2) = 0
       fjac(2,1) = -20*x(1)
       fjac(2,2) = 10
    end if
  end subroutine rosenbrock

  subroutine powell1(x, fvec, fjac)
    real(real64), dimension(:), intent(in) :: x
    real(real64), dimension(size(x)), intent(out) :: fvec
    real(real64), dimension(size(x),size(x)), intent(out), optional :: fjac

    fvec(1) = x(1)+10*x(2)
    fvec(2) = sqrt(5._real64)*(x(3)-x(4))
    fvec(3) = (x(2)-2*x(3))**2
    fvec(4) = sqrt(10._real64)*(x(1)-x(4))**2

    if (present(fjac)) then
       fjac(1,1) = 1
       fjac(1,2) = 10
       fjac(2,3) = sqrt(5._real64)
       fjac(2,4) = -fjac(2,3)
       fjac(3,2) = 2*(x(2)-2*x(3))
       fjac(3,3) = -2*fjac(3,2)
       fjac(4,1) = 2*sqrt(10._real64)*(x(1)-x(4))
       fjac(4,4) = -fjac(4,1)
    end if
  end subroutine powell1
end program trust_region_test
