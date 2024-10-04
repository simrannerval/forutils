module sincosintegrals

    implicit none
    
    contains
    
    subroutine cisia ( x, ci, si )

    !*****************************************************************************80
    !
    !! CISIA computes cosine Ci(x) and sine integrals Si(x).
    !
    !  Licensing:
    !
    !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
    !    they give permission to incorporate this routine into a user program 
    !    provided that the copyright is acknowledged.
    !
    !  Modified:
    !
    !    03 July 2012
    !
    !  Author:
    !
    !    Shanjie Zhang, Jianming Jin
    !
    !  Reference:
    !
    !    Shanjie Zhang, Jianming Jin,
    !    Computation of Special Functions,
    !    Wiley, 1996,
    !    ISBN: 0-471-11963-6,
    !    LC: QA351.C45.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the argument of Ci(x) and Si(x).
    !
    !    Output, real ( kind = 8 ) CI, SI, the values of Ci(x) and Si(x).
    !
    implicit none

    real ( kind = 8 ) bj(101)
    real ( kind = 8 ) ci
    real ( kind = 8 ) el
    real ( kind = 8 ) eps
    integer ( kind = 4 ) k
    integer ( kind = 4 ) m
    real ( kind = 8 ) p2
    real ( kind = 8 ) si
    real ( kind = 8 ) x
    real ( kind = 8 ) x2
    real ( kind = 8 ) xa
    real ( kind = 8 ) xa0
    real ( kind = 8 ) xa1
    real ( kind = 8 ) xcs
    real ( kind = 8 ) xf
    real ( kind = 8 ) xg
    real ( kind = 8 ) xg1
    real ( kind = 8 ) xg2
    real ( kind = 8 ) xr
    real ( kind = 8 ) xs
    real ( kind = 8 ) xss

    p2 = 1.570796326794897D+00
    el = 0.5772156649015329D+00
    eps = 1.0D-15
    x2 = x * x

    if ( x == 0.0D+00 ) then

        ci = -1.0D+300
        si = 0.0D+00

    else if ( x <= 16.0D+00 ) then

        xr = -0.25D+00 * x2
        ci = el + log ( x ) + xr
        do k = 2, 40
        xr = -0.5D+00 * xr * ( k - 1 ) / ( k * k * ( 2 * k - 1 ) ) * x2
        ci = ci + xr
        if ( abs ( xr ) < abs ( ci ) * eps ) then
            exit
        end if
        end do

        xr = x
        si = x
        do k = 1, 40
        xr = -0.5D+00 * xr * ( 2 * k - 1 ) / k / ( 4 * k * k + 4 * k + 1 ) * x2
        si = si + xr
        if ( abs ( xr ) < abs ( si ) * eps ) then
            return
        end if
        end do

    else if ( x <= 32.0D+00 ) then

        m = int ( 47.2D+00 + 0.82D+00 * x )
        xa1 = 0.0D+00
        xa0 = 1.0D-100
        do k = m, 1, -1
        xa = 4.0D+00 * k * xa0 / x - xa1
        bj(k) = xa
        xa1 = xa0
        xa0 = xa
        end do
        xs = bj(1)
        do k = 3, m, 2
        xs = xs + 2.0D+00 * bj(k)
        end do
        bj(1) = bj(1) / xs
        do k = 2, m
        bj(k) = bj(k) / xs
        end do
        xr = 1.0D+00
        xg1 = bj(1)
        do k = 2, m
        xr = 0.25D+00 * xr * ( 2.0D+00 * k - 3.0D+00 ) **2 &
            / ( ( k - 1.0D+00 ) * ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) * x
        xg1 = xg1 + bj(k) * xr
        end do

        xr = 1.0D+00
        xg2 = bj(1)
        do k = 2, m
        xr = 0.25D+00 * xr * ( 2.0D+00 * k - 5.0D+00 )**2 &
            / ( ( k-1.0D+00 ) * ( 2.0D+00 * k - 3.0D+00 ) ** 2 ) * x
        xg2 = xg2 + bj(k) * xr
        end do

        xcs = cos ( x / 2.0D+00 )
        xss = sin ( x / 2.0D+00 )
        ci = el + log ( x ) - x * xss * xg1 + 2.0 * xcs * xg2 - 2.0 * xcs * xcs
        si = x * xcs * xg1 + 2.0 * xss * xg2 - sin ( x )

    else

        xr = 1.0D+00
        xf = 1.0D+00
        do k = 1, 9
        xr = -2.0D+00 * xr * k * ( 2 * k - 1 ) / x2
        xf = xf + xr
        end do
        xr = 1.0D+00 / x
        xg = xr
        do k = 1, 8
        xr = -2.0D+00 * xr * ( 2 * k + 1 ) * k / x2
        xg = xg + xr
        end do
        ci = xf * sin ( x ) / x - xg * cos ( x ) / x
        si = p2 - xf * cos ( x ) / x - xg * sin ( x ) / x

    end if

    return
    end




    subroutine cisib ( x, ci, si )

    !*****************************************************************************80
    !
    !! CISIB computes cosine and sine integrals.
    !
    !  Licensing:
    !
    !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
    !    they give permission to incorporate this routine into a user program 
    !    provided that the copyright is acknowledged.
    !
    !  Modified:
    !
    !    20 March 2012
    !
    !  Author:
    !
    !    Shanjie Zhang, Jianming Jin
    !
    !  Reference:
    !
    !    Shanjie Zhang, Jianming Jin,
    !    Computation of Special Functions,
    !    Wiley, 1996,
    !    ISBN: 0-471-11963-6,
    !    LC: QA351.C45.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the argument of Ci(x) and Si(x).
    !
    !    Output, real ( kind = 8 ) CI, SI, the values of Ci(x) and Si(x).
    !
    implicit none

    real ( kind = 8 ) ci
    real ( kind = 8 ) fx
    real ( kind = 8 ) gx
    real ( kind = 8 ) si
    real ( kind = 8 ) x
    real ( kind = 8 ) x2

    x2 = x * x

    if ( x == 0.0D+00 ) then

        ci = -1.0D+300
        si = 0.0D+00

    else if ( x <= 1.0D+00 ) then

        ci = (((( -3.0D-08        * x2 &
                + 3.10D-06     ) * x2 &
                - 2.3148D-04   ) * x2 &
                + 1.041667D-02 ) * x2 &
                - 0.25D+00     ) * x2 + 0.577215665D+00 + log ( x )

        si = (((( 3.1D-07        * x2 &
                - 2.834D-05    ) * x2 &
                + 1.66667D-03  ) * x2 &
                - 5.555556D-02 ) * x2 + 1.0D+00 ) * x

    else

        fx = (((( x2              &
        + 38.027264D+00  ) * x2 &
        + 265.187033D+00 ) * x2 &
        + 335.67732D+00  ) * x2 &
        + 38.102495D+00  ) /    &
        (((( x2                 &
        + 40.021433D+00  ) * x2 &
        + 322.624911D+00 ) * x2 &
        + 570.23628D+00  ) * x2 &
        + 157.105423D+00 )

        gx = (((( x2               &
        + 42.242855D+00  ) * x2  &
        + 302.757865D+00 ) * x2  &
        + 352.018498D+00 ) * x2  &
        + 21.821899D+00 ) /      &
        (((( x2                  &
        + 48.196927D+00   ) * x2 &
        + 482.485984D+00  ) * x2 &
        + 1114.978885D+00 ) * x2 &
        + 449.690326D+00  ) / x

        ci = fx * sin ( x ) / x - gx * cos ( x ) / x

        si = 1.570796327D+00 - fx * cos ( x ) / x - gx * sin ( x ) / x

    end if

    return
    end

end module sincosintegrals