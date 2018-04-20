MODULE fonction

  use prec

  implicit none

contains
 function get_lambda(xi,yi) result(lambda)

    real(pr):: E,nu
    real(pr), intent(in)::xi, yi
    real(pr)::lambda

    !pour un lambda constant
    E=62000000000._pr
    nu=0.3_pr

    !lambda=(E*nu)/((1_pr+nu)*(1_pr - 2*nu))
    lambda=1
    !Si on veut avoir un lmabda qui change selon x on peut faire une fonction ici

  end function get_lambda



  function get_mu(xi,yi) result(mu)

    real(pr):: E,nu
    real(pr), intent(in)::xi, yi
    real(pr)::mu

    !pour un lambda constant
    E=62000000000._pr
    nu=0.3_pr
    mu=1
    !mu=E/(2_pr*(1_pr+nu))
    !Si on veut avoir un mu qui change selon x on peut faire une fonction ici

  end function get_mu


END MODULE fonction
