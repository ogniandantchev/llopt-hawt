program lloptwt

use msimslmd

real(8) pi, V, n, z, R, r1, rh, omega, rho, ni
real(8) D, dh, tif, c2tk2, sum, li, lt, J
real(8), external::	itp, iap, jp, Gfx, Gdx
real(8), allocatable:: s(:,:), a(:), RHS(:)
real(8), allocatable:: xii(:,:), xkk(:,:)
!real(8), allocatable:: dua1(:), dut1(:), dvr(:)
real(8), allocatable:: cl(:), cd(:), c(:)
real(8) drh, e1, cp, ct, kt, kt1
real(8) alfa, beta
integer(4) nn, k
integer(4) M, No
!logical(1) tb
real(8), external:: ktpif
real(8)	xdata(8), fdata(8)
real(8)	break(8), cscoef(4,8)
real(8), allocatable :: xa(:), ya(:), xb(:), yb(:)

real(8), parameter:: B1=0.0733, B2=0.0174
real(8)	deltah, deltat
real(8), allocatable:: deltac(:), delta(:), phi(:)
real(8) chh, ch1, ch2, cht, chr1, chr2
real(8)	xnaca(17), fcnacaa08b005(17), ftnaca66(17)
data xnaca /0., .0125, .025, .05, .075, .1, .2, .3, .4, &
	.45, .5, .6, .7, .8, .9, .95, 1./
data fcnacaa08b005 /0., .0686, .142, .282, .389, .475, .725, &
	.881, .97, .992, 1., .971, .877, .69, .352, .168, 0./
data ftnaca66 /0., .231, .306, .419, .508, .584, .8, .927, &
	.99, 1., .992, .931, .807, .622, .375, .229, .006/
integer(4), parameter:: maxft= 10

common /sp1/ break, cscoef
common /cho/ chh, ch1, ch2, cht, drh, chr1, chr2

pi = 4.*atan(1.)

V = 10.
n = 3.25
z = 2.
D = 7.
!T = ..
!read*, KT
!KT= 0.33
KT1= 0.
rho= 1.2
ni= 15d-6
R = D/2.
r1 = 1000000.
e1= 0.001
J= V/(n*D)

read*, P
!P = 1500
cp = 8.*P/(pi*rho*D**2*V**3)
tmp = dsqrt(4.*cp/3.)
ct = -2.* tmp * dsinh(dasinh((-2.*cp**2)/tmp**3)/3.)
KT = pi*ct*J**2/8.
print*, cp, ct, KT

drh = 0.2
rh= drh*R
dh= rh * 2.

omega = 2. * pi * n

chh= 0.07	!14
ch1= 0.09	!16
ch2= 0.09	!16
cht= 0.07	!14
chr1= 0.4
chr2= 0.7
deltah=	0.13
deltat=	0.04
xs= 0.
ts= 0.

alfa= 0.5
beta= 1.5

lt= V/(omega*R)
li= 0.944*lt-0.065		!???
!li = .26*lt**2 + .812*lt + .35*KT + .03

np= 40
M= 7
No= M+1

allocate(s(1:M+1,0:M), a(0:M), RHS(1:M+1))
allocate( xii(0:2,1:No), xkk(0:2,1:No))
allocate( cl(1:No), cd(1:No), c(1:No))
allocate( deltac(1:No), delta(1:No), phi(1:No))

do k= 1, No
	xkk(0,k)= dcos(pi*(2.*k-1.)/(2.*No+1.))
	xkk(1,k)= x2dr(xkk(0,k), drh)
	xkk(2,k)= R*xkk(1,k)
!	print*, xkk(0,k), xkk(1,k), xkk(2,k)
end do
print*, ' '
do i= 1, No
	xii(0,i)= dcos(2*pi*i/(2*No+1.))
	xii(1,i)= x2dr(xii(0,i), drh)
	xii(2,i)= R*xii(1,i)
!	print*, xii(0,i), xii(1,i), xii(2,i)
end do
do i= 1, No
	delta(i)= (1.-log10(9.*(xii(1,i)-drh)/(1.-drh)+1.))* &
		(deltah-deltat)+deltat
	c(i)= chord(xii(1,i))
end do

do while (abs(KT-KT1) > e1)

  do i= 1, No
	  do nn= 0, M
		sum= 0.
		  do k= 1, No
	  		c2tk2= dcos(pi*(2.*k-1.)/(2.*(2.*No+1.)))**2
		  	tif= itp(z, li/xkk(1,k), xii(1,i), 1d-4, 1d4, xkk(1,k))
		  	sum= sum + c2tk2*jp(nn, -.5,.5,xkk(0,k))*tif/(xkk(0,k)-xii(0,i))
		  end do
		s(i,nn)= (1/(1.-drh))*(4.*pi/(2.*No+1.))*sum
	  end do
	  s(i,0)= s(i,0) + pi*z/(2.*xii(1,i))
!	RHS(i)= 2.*((lt-li)/lt)*(xii(1,i)*li/(xii(1,i)**2+li**2))
    be = datan(lt/xii(1,i))
    bei= datan(li/xii(1,i))
	  RHS(i)= dsin(be-bei)*dsin(bei)/dsin(be)
  end do

call dlsarg(M+1, s, M+1, RHS, 1, a)

do i= 1, No
  be = datan(lt/xii(1,i))
  bei= datan(li/xii(1,i))
  u1=	dsin(be-bei)*dsin(bei)/dsin(be)
  v1= 1 - dtan(bei)*(pi*xii(1,i)/J+u1)
	xdata(No-i+1)= xii(1,i)
	tempg= Gfx(xii(0,i), M, alfa, beta, a)
!	tempvrm= pi*xii(1,i)/J + ((lt-li)/lt)*(xii(1,i)*li/(xii(1,i)**2+li**2))
!	tempvrm= pi*xii(1,i)/J + dsin(be-bei)*dsin(bei)/dsin(be)
	cbi= dcos(bei)
  vri= dsqrt((1-v1)**2+(pi*xii(1,i)/J)**2)
	rns= vri*c(i)*R/ni
	cd(i)= .05808*(1. + 2.3*delta(i))/rns**(.1458)
	cl(i)= 2.*pi*tempg/(c(i)*vri)
	fdata(No-i+1)= tempg*vri*cbi*(1.-(cd(i)/cl(i))*li/xii(1,i))
end do

  CALL dCSINT (No, XDATA, FDATA, BREAK, CSCOEF)
  CALL dQDAGS (ktpif, drh, 1., 0., 0.0001, RES1, ERR1)

	KT1 = -z*J**2.*pi*RES1/2.

  print*, KT1, li

  li= li*KT/KT1
end do

CALL dWRRRN ('a', 1, M+1, a, 1, 0)
print*, lt, li

do i= 1, No
	deltac(i)= B1 * cl(i)
	phi(i)= datan(li/xii(1,i)) - B2*cl(i)
!	delta(i)= - delta(i)
end do


allocate( xa(33), ya(33), xb(np+1), yb(np+1))
ttt = 8.*datan(1.)/dfloat(np)

open(1, FILE='t1.prg')
do k= 1, No
	do i= 17, 33
		xa(i)= xnaca(i-16)
		ya(i)= 2.*deltac(k)*fcnacaa08b005(i-16) + &
			delta(k)*ftnaca66(i-16)
	end do
	do i= 1, 17
		xa(18-i)= - xnaca(i)
		ya(18-i)= 2.*deltac(k)*fcnacaa08b005(i) - &
			delta(k)*ftnaca66(i)
	end do
	do i=0, (np / 2) -1
		xb(i+1) = (1+dcos(dfloat(i)*ttt))/2.
		xb(np-i+1) = - (1+dcos(dfloat(i)*ttt))/2.
	end do
	xb(np/2+1) = 0.
	CALL dCSIEZ (33, xa, ya, np+1, xb, yb)
	do i=1, np+1
		xb(i) = 1. - 2.*abs(xb(i))
	end do

	write(1,113)
	write(1,114)
	rr= xii(1,k)
	cc= c(k)*R/2.
	do jj= 1, np+1
		x= xs + (xb(jj)*dsin(phi(k))+yb(jj)*dcos(phi(k)))*cc
		tita= ts + (c(k)/(2.*rr))* &
			(xb(jj)*dcos(phi(k))-yb(jj)*dsin(phi(k)))
		y= rr*R*dcos(tita)
		z= rr*R*dsin(tita)
		write(1,111) x, y, z
	end do
	write(1,115)
	write(1,116)
end do
write(1,112)
close(1)

!deallocate( )

open(1, FILE='1.plt')
do k= 1, No
	do i= 17, 33
		xa(i)= xnaca(i-16)
		ya(i)= 2.*deltac(k)*fcnacaa08b005(i-16) + &
			delta(k)*ftnaca66(i-16)
	end do
	do i= 1, 17
		xa(18-i)= - xnaca(i)
		ya(18-i)= 2.*deltac(k)*fcnacaa08b005(i) - &
			delta(k)*ftnaca66(i)
	end do
	do i=0, (np / 2) -1
		xb(i+1) = (1+dcos(dfloat(i)*ttt))/2.
		xb(np-i+1) = - (1+dcos(dfloat(i)*ttt))/2.
	end do
	xb(np/2+1) = 0.
	CALL dCSIEZ (33, xa, ya, np+1, xb, yb)
	do i=1, np+1
		xb(i) = 1. - 2.*abs(xb(i))
	end do

	rr= xii(1,k)
	cc= c(k)*R/2.
	do jj= 1, np+1
		x= (xb(jj)*dsin(phi(k))+yb(jj)*dcos(phi(k)))*cc
		y= (xb(jj)*dcos(phi(k))-yb(jj)*dsin(phi(k)))*cc
		write(1,*) x, y
	end do
end do
close(1)

open(1, FILE='delta')
do i= 1, No
	write(1, *) xii(1, i), cl(i), cd(i), xii(2,i), c(i), deltac(i), &
		delta(i), phi(i)*180./pi
end do
close(1)

open(1, FILE='prof')
do i= 1, No
	do k= 1, 17
		write(1,*) xnaca(k)*c(i)*R, &
	 -(2.*deltac(i)*fcnacaa08b005(k) + delta(i)*ftnaca66(k))*c(i)*R, &
	 -(2.*deltac(i)*fcnacaa08b005(k) - delta(i)*ftnaca66(k))*c(i)*R
	end do
end do
close(1)

open(1, FILE = 'g')
do i= 1, 50
	write(1,*) (dfloat(i)*(2/51.)-1.), -Gfx((dfloat(i)*(2./51.)-1.), M, alfa, beta, a), -Gdx((dfloat(i)*(2./51.)-1.), M, alfa, beta, a)
end do
close(1)

111 format('K : ', F12.6, ',', F12.6, ',', F12.6)
112 format('E : end')
113 format('K : $ /cr m3 sp')
114 format('K : KEY')
115 format('K : DON')
116 format('K : N')
end program lloptwt


real(8)	function chord(r)
real(8) chh, ch1, ch2, cht, drh, chr1, chr2
common /cho/ chh, ch1, ch2, cht, drh, chr1, chr2
	if (r > chr2) then
		chord= (r-chr2)*((cht-ch2)/(1.-chr2)) + ch2
	else if (r > chr1) then
		chord= (r-chr1)*((ch2-ch1)/(chr2-chr1)) + ch1
	else
		chord= (r-drh)*((ch1-chh)/(chr1-drh)) + chh
	end if
end function

real(8) function ktpif(dr)
real(8)	dr, break(8), cscoef(4,8)
common /sp1/ break, cscoef
	ktpif= dCSVAL(dr, 7, BREAK, CSCOEF)
end function

real(8) FUNCTION bn(p, n)
real(8) p, a, f
integer(4) n
   IF (n == 0.) THEN
        bn = 1.
   ELSE
        a = 1.
        do i= 0, n-1
            a = a * (p - i)
        end do
        f = 1.
        do i= 1, n
            f = f * i
        end do
        bn = a / f
   END IF
END FUNCTION

real(8) FUNCTION jp(n, a, b, x)
integer(4) n, m
real(8) a, b, x, s
   s = 0.
   do m= 0, n
      s = s + bn(n+a, m)*bn(n+b, n-m)*((x-1.)**(n-m))*(x+1.)**m
   end do
   jp = s / (2**n)
END FUNCTION

real(8)	function x2dr(x, drh)
real(8)	x, drh
	x2dr= (x*(1.-drh)+1.+drh)/2.
end function

real(8)	function dr2x(dr, drh)
real(8)	dr, drh
	dr2x= (2.*dr-1.-drh)/(1.-drh)
end function

real(8) function Gfx(x, m, a, b, coefs)
integer(4) m, i
real(8), external:: jp
real(8) sum, h, x, a, b, coefs(0:m)
	sum= 0.
	do i= 1, m
		sum= sum + coefs(i)*jp(i-1, a, b, x)/dfloat(i)
	end do
	h= -coefs(0)*(.5*(b-a+1.)*dacos(x)+(b-a)*dsqrt(1.-x**2.)+.5*(b-a-1.)*x*dsqrt(1-x**2.))
	Gfx= h - .5*((1.-x)**a)*((1.+x)**b)*sum
end function

real(8) function Gdx(x, m, a, b, coefs)
integer(4) m, i
real(8), external:: jp
real(8) sum, x, a, b, coefs(0:m)
	sum= 0.
	do i= 0, m
		sum= sum + coefs(i)*jp(i, a-1., b-1., x)
	end do
	Gdx= ((1.-x)**(a-1.))*((1.+x)**(b-1.))*sum
end function
