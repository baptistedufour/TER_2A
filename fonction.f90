MODULE fonction

  use prec

  implicit none

contains
 function get_lambda(xi,yi) result(lambda)

    real(pr):: E1,nu, E2
    real(pr), intent(in)::xi, yi
    real(pr)::lambda

    !pour un lambda constant
    E1=74000._pr
    E2=210000._pr
    nu=0.3_pr

    !lambda=(E1*nu)/((1_pr+nu)*(1_pr - 2*nu))
    !lambda=1
    !Si on veut avoir un lmabda qui change selon x on peut faire une fonction ici

    if (yi<-65) then
        lambda=(E2*nu)/((1_pr+nu)*(1_pr - 2*nu))
      else if (yi>-45) then
        lambda=(E1*nu)/((1_pr+nu)*(1_pr - 2*nu))
      else if ((yi>-65) .AND. (yi<-45)) then
        lambda=((E2*nu)/((1_pr+nu)*(1_pr - 2*nu)))*(yi+45)/(-20)+((E1*nu)/((1_pr+nu)*(1_pr - 2*nu)))*(yi+65)/20
    end if

  end function get_lambda

  function get_mu(xi,yi) result(mu)

    real(pr):: E1, E2,nu
    real(pr), intent(in)::xi, yi
    real(pr)::mu

    !pour un lambda constant
    E1=74000._pr
    E2=210000._pr
    nu=0.3_pr
    !mu=1
    !mu=E1/(2_pr*(1_pr+nu))
    !Si on veut avoir un mu qui change selon x on peut faire une fonction ici

    if (yi<-65) then
        mu=E2/(2_pr*(1_pr+nu))
    else if (yi>-45) then
        mu=E1/(2_pr*(1_pr+nu))
    else if ((yi>-65) .AND. (yi<-45)) then
        mu=(E2/(2_pr*(1_pr+nu)))*(yi+45)/(-20)+(E1/(2_pr*(1_pr+nu)))*(yi+65)/20
    end if

  end function get_mu

END MODULE fonction
