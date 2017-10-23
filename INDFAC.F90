
real(8) function itp(z, tbp, r, rh, r1, rp)
real(8) z, r, r1, rh, rp, tbp
real(8) mu, mup, muh, mu1, u, up, uh, u1
real(8) s, rer, ni, a, b, y, y1, yh, yp
real(8) f1, f2, f12, dit, it, rrr

if (r /= rp) then
    rrr= rp*tbp
	mu = r / rrr
    mup = rp / rrr
    muh = rh / rrr
    mu1 = r1 / rrr

    u = dSQRt(1. + mu**2)
    up = dSQRt(1. + mup**2)
    uh = dSQRt(1. + muh**2)
    u1 = dSQRt(1. + mu1**2)

    s = SGN(r - rp)
    rer = 1. - rp / r
    ni = ((((u - 1.) / mu) * (mup / (up - 1.)))**(z * s)) * dEXP(z * dABS(u - up))
    f1 = (s / (2. * z * mup)) * dsqrt(up / u)
    f2 = 1. / (ni - 1.) - (s / (24. * z)) * ((9. * mup**2 + 2.) / up**3 + (3. * mu**2 - 2.) / u**3) * dLOG(1. + (1. / (ni - 1.)))
    f12 = f1 * f2
    IF (r > rp) THEN
        it = z * rer * (1. + 2 * mup * z * f12)
      ELSE
        it = 2 * (z**2) * rer * mup * f12
    END IF
else
	it= dsin(datan(tbp))
end if

!    y = yfun(u)
!    yp = yfun(up)

!	if (rh < 1d-3) then
!		b= 0.
!	else
!		yh = yfun(uh)
!	    b = 1. / (1. - dEXP(-z * (y + yp - 2 * yh)))
!	end if
!   if (r1 > 1000.) then
!		a= 0.
!	else
!		y1 = yfun(u1)
!		a = 1. / (1. - dEXP(-z * (2 * y1 - y - yp)))
!	end if

!    dit = -z * rer * dsqrt(up / u) * (a - b)

    itp = it !+ dit

end function


real(8) function iap(z, tbp, r, rh, r1, rp)
real(8) z, r, r1, rh, rp, tbp
real(8) mu, it, itp

    mu = r / (rp*tbp)
	it = itp(z, tbp, r, rh, r1, rp)
    iap = mu * (it + z * (rp / r - 1.))

end function


real(8) function sgn(x)
    IF (x < 0.) THEN
		sgn= -1.
	ELSE IF (x == 0.) THEN
		sgn= 0.
	ELSE
		sgn= 1.
    END IF
end function


real(8) FUNCTION yfun (u)
real(8) u
   yfun = u - .5 * LOG((u + 1.) / (u - 1.))
END FUNCTION
