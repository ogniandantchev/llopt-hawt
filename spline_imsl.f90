!==============================================================================
! MODULE spline_imsl
!
! Drop-in replacement for the IMSL CSINT / CSVAL / CSDER routine family.
! All routines use the IMSL pp-form convention:
!
!   On interval [BREAK(i), BREAK(i+1)], with dx = x - BREAK(i):
!     S(x) = CSCOEF(1,i)
!           + CSCOEF(2,i)*dx
!           + CSCOEF(3,i)*dx**2
!           + CSCOEF(4,i)*dx**3
!
! Public routines
! ---------------
!   CSINT  — fit not-a-knot cubic spline, return pp-form (BREAK, CSCOEF)
!   CSVAL  — evaluate spline (or any derivative) at a single point
!   CSDER  — differentiate pp-form: return pp-form of the k-th derivative
!
! Usage example
! -------------
!   USE spline_imsl
!   REAL(8) :: x(n), f(n), brk(n), coef(4,n)
!   CALL CSINT(n, x, f, brk, coef)
!   val = CSVAL(xq, n, brk, coef)          ! evaluate S(xq)
!   d1  = CSVAL(xq, n, brk, coef, IDERIV=1)! evaluate S'(xq)
!   CALL CSDER(1, n, brk, coef, dbrk, dcoef, nd)   ! first-derivative pp-form
!
! Verified against scipy.interpolate.CubicSpline(bc_type='not-a-knot').
!==============================================================================

MODULE spline_imsl
  IMPLICIT NONE
  PRIVATE

  ! Floating-point kind — change to REAL(4) or KIND=selected_real_kind(15)
  ! throughout if a different precision is needed.
  INTEGER, PARAMETER :: WP = KIND(1.0D0)   ! working precision (double)

  PUBLIC :: WP
  PUBLIC :: CSINT, CSVAL, CSDER, CSIEZ, CSITG

  !----------------------------------------------------------------------------
  INTERFACE CSINT
    MODULE PROCEDURE csint_impl
  END INTERFACE CSINT

  !----------------------------------------------------------------------------
  INTERFACE CSVAL
    MODULE PROCEDURE csval_impl
  END INTERFACE CSVAL

  !----------------------------------------------------------------------------
  INTERFACE CSDER
    MODULE PROCEDURE csder_impl
  END INTERFACE CSDER

  !----------------------------------------------------------------------------
  INTERFACE CSIEZ
    MODULE PROCEDURE csiez_impl
  END INTERFACE CSIEZ

  !----------------------------------------------------------------------------
  INTERFACE CSITG
    MODULE PROCEDURE csitg_impl
  END INTERFACE CSITG

CONTAINS

  !============================================================================
  ! CSINT — Compute cubic spline with not-a-knot end conditions.
  !
  ! CALL CSINT(NDATA, XDATA, FDATA, BREAK, CSCOEF)
  !
  !   NDATA   (in)  INTEGER   Number of data points (>= 2).
  !   XDATA   (in)  REAL(WP)  Abscissas, shape (NDATA).  Strictly increasing.
  !   FDATA   (in)  REAL(WP)  Ordinates, shape (NDATA).
  !   BREAK   (out) REAL(WP)  Breakpoints, shape (NDATA).  Copy of XDATA.
  !   CSCOEF  (out) REAL(WP)  Coefficients, shape (4, NDATA).
  !                           Column j unused for j = NDATA (zeroed sentinel).
  !============================================================================
  SUBROUTINE csint_impl(NDATA, XDATA, FDATA, BREAK, CSCOEF)
    INTEGER,     INTENT(IN)  :: NDATA
    REAL(WP),    INTENT(IN)  :: XDATA(NDATA), FDATA(NDATA)
    REAL(WP),    INTENT(OUT) :: BREAK(NDATA), CSCOEF(4,NDATA)

    INTEGER  :: i, n
    REAL(WP) :: h(NDATA-1)   ! h(i) = x_{i+1} - x_i
    REAL(WP) :: s(NDATA-1)   ! s(i) = (f_{i+1}-f_i)/h(i)
    REAL(WP) :: bd(NDATA)    ! main diagonal (overwritten during Thomas sweep)
    REAL(WP) :: ud(NDATA)    ! upper diagonal
    REAL(WP) :: r(NDATA)     ! right-hand side
    REAL(WP) :: m(NDATA)     ! second derivatives S''(x_i)
    REAL(WP) :: alpha, beta, piv

    n = NDATA

    IF (n < 2) ERROR STOP 'CSINT: NDATA must be >= 2'

    BREAK = XDATA

    ! ------------------------------------------------------------------
    ! n = 2: linear, stored as a degenerate cubic
    ! ------------------------------------------------------------------
    IF (n == 2) THEN
      IF (XDATA(2) <= XDATA(1)) &
        ERROR STOP 'CSINT: XDATA must be strictly increasing'
      CSCOEF(1,1) = FDATA(1)
      CSCOEF(2,1) = (FDATA(2) - FDATA(1)) / (XDATA(2) - XDATA(1))
      CSCOEF(3,1) = 0.0_WP
      CSCOEF(4,1) = 0.0_WP
      CSCOEF(1,2) = FDATA(2)
      CSCOEF(2:4,2) = 0.0_WP
      RETURN
    END IF

    ! ------------------------------------------------------------------
    ! Interval widths and first divided differences
    ! ------------------------------------------------------------------
    DO i = 1, n-1
      h(i) = XDATA(i+1) - XDATA(i)
      IF (h(i) <= 0.0_WP) THEN
        WRITE(*,'(A,I0)') 'CSINT: XDATA not strictly increasing at index ', i
        ERROR STOP
      END IF
      s(i) = (FDATA(i+1) - FDATA(i)) / h(i)
    END DO

    ! ------------------------------------------------------------------
    ! n = 3: NAK forces constant second derivative m = S''
    !   Derived by setting m(1)=m(2)=m(3) and solving the single
    !   interior continuity equation:
    !     m*(3*h(1) + 3*h(2)) = 6*(s(2) - s(1))
    ! ------------------------------------------------------------------
    IF (n == 3) THEN
      m(1) = 2.0_WP * (s(2) - s(1)) / (h(1) + h(2))
      m(2) = m(1)
      m(3) = m(1)

    ELSE
      ! ----------------------------------------------------------------
      ! General case n >= 4.
      !
      ! Solve a reduced (n-2)x(n-2) tridiagonal system for m(2)..m(n-1)
      ! after eliminating m(1) and m(n) via the NAK conditions:
      !
      !   NAK left:   m(1) = [(h1+h2)*m(2) - h1*m(3)] / h2
      !     alpha = h(1)/h(2)
      !     => row 2: bd = (h1+h2)*(2+alpha),  ud = h2 - alpha*h1
      !
      !   NAK right:  m(n) = [(h_{n-2}+h_{n-1})*m(n-1) - h_{n-1}*m(n-2)] / h_{n-2}
      !     beta = h(n-1)/h(n-2)
      !     => row n-1: sub = h_{n-2} - beta*h_{n-1},
      !                 bd  = (h_{n-2}+h_{n-1})*(2+beta)
      ! ----------------------------------------------------------------
      alpha = h(1)   / h(2)
      beta  = h(n-1) / h(n-2)

      ! Row 2 (NAK-modified left endpoint)
      bd(2) = (h(1) + h(2)) * (2.0_WP + alpha)
      ud(2) = h(2) - alpha*h(1)
      r(2)  = 6.0_WP * (s(2) - s(1))

      ! Standard interior rows 3 .. n-2
      DO i = 3, n-2
        bd(i) = 2.0_WP * (h(i-1) + h(i))
        ud(i) = h(i)
        r(i)  = 6.0_WP * (s(i) - s(i-1))
      END DO

      ! Row n-1 (NAK-modified right endpoint)
      bd(n-1) = (h(n-2) + h(n-1)) * (2.0_WP + beta)
      ud(n-1) = 0.0_WP
      r(n-1)  = 6.0_WP * (s(n-1) - s(n-2))

      ! Thomas algorithm — forward sweep (rows 3 .. n-1)
      ! The sub-diagonal is h(i-1) for interior rows and
      ! (h(n-2) - beta*h(n-1)) for the last row.
      DO i = 3, n-1
        IF (i == n-1) THEN
          piv = (h(n-2) - beta*h(n-1)) / bd(i-1)
        ELSE
          piv = h(i-1) / bd(i-1)
        END IF
        bd(i) = bd(i) - piv*ud(i-1)
        r(i)  = r(i)  - piv*r(i-1)
      END DO

      ! Back substitution
      m(n-1) = r(n-1) / bd(n-1)
      DO i = n-2, 2, -1
        m(i) = (r(i) - ud(i)*m(i+1)) / bd(i)
      END DO

      ! Recover end second derivatives from NAK expressions
      m(1) = ((h(1)+h(2))*m(2) - h(1)*m(3)) / h(2)
      m(n) = ((h(n-2)+h(n-1))*m(n-1) - h(n-1)*m(n-2)) / h(n-2)

    END IF

    ! ------------------------------------------------------------------
    ! Convert second derivatives to pp-form coefficients.
    !
    ! On [x_i, x_{i+1}], dx = x - x_i:
    !   C(1,i) = f_i
    !   C(2,i) = s_i - h_i*(2*m_i + m_{i+1})/6
    !   C(3,i) = m_i/2
    !   C(4,i) = (m_{i+1} - m_i)/(6*h_i)
    ! ------------------------------------------------------------------
    DO i = 1, n-1
      CSCOEF(1,i) = FDATA(i)
      CSCOEF(2,i) = s(i) - h(i)*(2.0_WP*m(i) + m(i+1))/6.0_WP
      CSCOEF(3,i) = m(i)/2.0_WP
      CSCOEF(4,i) = (m(i+1) - m(i))/(6.0_WP*h(i))
    END DO

    ! Zeroed sentinel column (IMSL: evaluator only uses cols 1..NDATA-1)
    CSCOEF(1,n)   = FDATA(n)
    CSCOEF(2:4,n) = 0.0_WP

  END SUBROUTINE csint_impl


  !============================================================================
  ! CSVAL — Evaluate a cubic spline (or one of its derivatives) at a point.
  !
  ! result = CSVAL(X, NDATA, BREAK, CSCOEF [, IDERIV])
  !
  !   X       (in) REAL(WP)   Evaluation point.
  !   NDATA   (in) INTEGER    Number of breakpoints (same as NDATA in CSINT).
  !   BREAK   (in) REAL(WP)   Breakpoints, shape (NDATA).
  !   CSCOEF  (in) REAL(WP)   Coefficients, shape (4, NDATA).
  !   IDERIV  (in) INTEGER    Optional.  Order of derivative to evaluate:
  !                             0 = S(x)   [default]
  !                             1 = S'(x)
  !                             2 = S''(x)
  !                             3 = S'''(x)
  !                            >3 = 0 (spline is piecewise cubic)
  !
  ! X is clamped to [BREAK(1), BREAK(NDATA)] (extrapolation uses the
  ! end cubic pieces, consistent with IMSL behaviour).
  !============================================================================
  FUNCTION csval_impl(X, NDATA, BREAK, CSCOEF, IDERIV) RESULT(val)
    REAL(WP),    INTENT(IN)           :: X
    INTEGER,     INTENT(IN)           :: NDATA
    REAL(WP),    INTENT(IN)           :: BREAK(NDATA), CSCOEF(4,NDATA)
    INTEGER,     INTENT(IN), OPTIONAL :: IDERIV
    REAL(WP)                          :: val

    INTEGER  :: i, id
    REAL(WP) :: dx

    id = 0
    IF (PRESENT(IDERIV)) id = IDERIV

    ! Locate the interval containing X (binary search)
    i = csint_locate(X, NDATA, BREAK)

    dx = X - BREAK(i)

    SELECT CASE (id)
    CASE (0)
      val = CSCOEF(1,i) + dx*(CSCOEF(2,i) + dx*(CSCOEF(3,i) + dx*CSCOEF(4,i)))
    CASE (1)
      val = CSCOEF(2,i) + dx*(2.0_WP*CSCOEF(3,i) + dx*3.0_WP*CSCOEF(4,i))
    CASE (2)
      val = 2.0_WP*CSCOEF(3,i) + dx*6.0_WP*CSCOEF(4,i)
    CASE (3)
      val = 6.0_WP*CSCOEF(4,i)
    CASE DEFAULT
      val = 0.0_WP
    END SELECT

  END FUNCTION csval_impl


  !============================================================================
  ! CSDER — Differentiate a pp-form spline k times.
  !
  ! Returns the pp-form of the k-th derivative, which is a piecewise
  ! polynomial of degree (3-k) on the same breakpoint sequence.
  !
  ! CALL CSDER(K, NDATA, BREAK, CSCOEF, DBREAK, DCSCOEF, NDBREAK)
  !
  !   K        (in)  INTEGER   Order of differentiation (1, 2, or 3).
  !   NDATA    (in)  INTEGER   Number of breakpoints.
  !   BREAK    (in)  REAL(WP)  Breakpoints, shape (NDATA).
  !   CSCOEF   (in)  REAL(WP)  Coefficients, shape (4, NDATA).
  !   DBREAK   (out) REAL(WP)  Breakpoints of derivative, shape (NDATA).
  !                            Equal to BREAK (same knot sequence).
  !   DCSCOEF  (out) REAL(WP)  Coefficients of derivative, shape (4, NDATA).
  !                            Rows (K+1)..4 of each column are zero.
  !   NDBREAK  (out) INTEGER   Number of breakpoints (= NDATA; provided for
  !                            convenience when chaining CSDER calls).
  !============================================================================
  SUBROUTINE csder_impl(K, NDATA, BREAK, CSCOEF, DBREAK, DCSCOEF, NDBREAK)
    INTEGER,     INTENT(IN)  :: K, NDATA
    REAL(WP),    INTENT(IN)  :: BREAK(NDATA), CSCOEF(4,NDATA)
    REAL(WP),    INTENT(OUT) :: DBREAK(NDATA), DCSCOEF(4,NDATA)
    INTEGER,     INTENT(OUT) :: NDBREAK

    INTEGER  :: i, d
    REAL(WP) :: c(4)

    IF (K < 1 .OR. K > 3) &
      ERROR STOP 'CSDER: K must be 1, 2, or 3'

    DBREAK  = BREAK
    NDBREAK = NDATA
    DCSCOEF = 0.0_WP

    DO i = 1, NDATA-1
      c = CSCOEF(:,i)
      ! Differentiate d times in place:  c = d/dx of current polynomial
      DO d = 1, K
        c(1) = c(2)
        c(2) = 2.0_WP*c(3)
        c(3) = 3.0_WP*c(4)
        c(4) = 0.0_WP
      END DO
      DCSCOEF(:,i) = c
    END DO

    ! Sentinel column: value at final knot recoverable from column NDATA-1
    DCSCOEF(1,NDATA) = CSVAL(BREAK(NDATA), NDATA, BREAK, CSCOEF, IDERIV=K)
    DCSCOEF(2:4,NDATA) = 0.0_WP

  END SUBROUTINE csder_impl


  !============================================================================
  ! CSIEZ — Convenience: fit a cubic spline and evaluate at specified points.
  !
  ! Drop-in replacement for IMSL dCSIEZ.
  !
  ! CALL CSIEZ(NDATA, XDATA, FDATA, NOUT, XOUT, YOUT)
  !
  !   NDATA  (in)  INTEGER   Number of data points (>= 2).
  !   XDATA  (in)  REAL(WP)  Abscissas, shape (NDATA).  Strictly increasing.
  !   FDATA  (in)  REAL(WP)  Ordinates, shape (NDATA).
  !   NOUT   (in)  INTEGER   Number of evaluation points.
  !   XOUT   (in)  REAL(WP)  Points at which to evaluate, shape (NOUT).
  !   YOUT   (out) REAL(WP)  Spline values at XOUT, shape (NOUT).
  !============================================================================
  SUBROUTINE csiez_impl(NDATA, XDATA, FDATA, NOUT, XOUT, YOUT)
    INTEGER,  INTENT(IN)  :: NDATA, NOUT
    REAL(WP), INTENT(IN)  :: XDATA(NDATA), FDATA(NDATA)
    REAL(WP), INTENT(IN)  :: XOUT(NOUT)
    REAL(WP), INTENT(OUT) :: YOUT(NOUT)

    REAL(WP) :: BRK(NDATA), COEF(4, NDATA)
    INTEGER  :: j

    CALL CSINT(NDATA, XDATA, FDATA, BRK, COEF)
    DO j = 1, NOUT
      YOUT(j) = CSVAL(XOUT(j), NDATA, BRK, COEF)
    END DO

  END SUBROUTINE csiez_impl


  !============================================================================
  ! CSITG — Analytically integrate a cubic spline over [A, B].
  !
  ! result = CSITG(A, B, NDATA, BREAK, CSCOEF)
  !
  ! Since S(x) is a piecewise cubic, the integral is computed exactly
  ! by summing the integrals of each polynomial piece that overlaps [A, B].
  !============================================================================
  FUNCTION csitg_impl(A, B, NDATA, BREAK, CSCOEF) RESULT(integral)
    REAL(WP), INTENT(IN) :: A, B
    INTEGER,  INTENT(IN) :: NDATA
    REAL(WP), INTENT(IN) :: BREAK(NDATA), CSCOEF(4, NDATA)
    REAL(WP)             :: integral

    INTEGER  :: i, i_lo, i_hi
    REAL(WP) :: x_lo, x_hi, dx_lo, dx_hi

    integral = 0.0_WP

    ! Clamp bounds to the spline domain
    x_lo = MAX(A, BREAK(1))
    x_hi = MIN(B, BREAK(NDATA))
    IF (x_lo >= x_hi) RETURN

    ! Find intervals containing x_lo and x_hi
    i_lo = csint_locate(x_lo, NDATA, BREAK)
    i_hi = csint_locate(x_hi, NDATA, BREAK)

    ! Sum the exact integral of each polynomial piece
    DO i = i_lo, i_hi
      IF (i == i_lo) THEN
        dx_lo = x_lo - BREAK(i)
      ELSE
        dx_lo = 0.0_WP
      END IF

      IF (i == i_hi) THEN
        dx_hi = x_hi - BREAK(i)
      ELSE
        dx_hi = BREAK(i+1) - BREAK(i)
      END IF

      ! Integral of C1 + C2*dx + C3*dx^2 + C4*dx^3  from dx_lo to dx_hi
      integral = integral + &
        CSCOEF(1,i) * (dx_hi - dx_lo) + &
        CSCOEF(2,i) * (dx_hi**2 - dx_lo**2) / 2.0_WP + &
        CSCOEF(3,i) * (dx_hi**3 - dx_lo**3) / 3.0_WP + &
        CSCOEF(4,i) * (dx_hi**4 - dx_lo**4) / 4.0_WP
    END DO

  END FUNCTION csitg_impl


  !============================================================================
  ! csint_locate (private) — binary search for the pp interval containing X.
  !
  ! Returns index i such that BREAK(i) <= X < BREAK(i+1), with edge cases:
  !   X <  BREAK(1)    => i = 1        (extrapolate left)
  !   X >= BREAK(NDATA) => i = NDATA-1 (extrapolate right / right endpoint)
  !============================================================================
  PURE FUNCTION csint_locate(X, N, BREAK) RESULT(i)
    REAL(WP), INTENT(IN) :: X
    INTEGER,  INTENT(IN) :: N
    REAL(WP), INTENT(IN) :: BREAK(N)
    INTEGER              :: i, lo, hi, mid

    IF (X <= BREAK(1)) THEN
      i = 1; RETURN
    END IF
    IF (X >= BREAK(N)) THEN
      i = N-1; RETURN
    END IF

    lo = 1; hi = N
    DO WHILE (hi - lo > 1)
      mid = (lo + hi) / 2
      IF (X >= BREAK(mid)) THEN
        lo = mid
      ELSE
        hi = mid
      END IF
    END DO
    i = lo

  END FUNCTION csint_locate

END MODULE spline_imsl
