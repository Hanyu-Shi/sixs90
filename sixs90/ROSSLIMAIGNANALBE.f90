subroutine rlmaignanalbe(p1, p2, p3, brdfalb)
  implicit none
  real(8) :: p1, p2, p3, brdfalb
  brdfalb = p1 + p2*0.11627450 - p3*1.6798152
  if (brdfalb<0) then
    write (6, *) 'warning white sky albedo <0  , reset to 0.'
    brdfalb = 0.
  end if
  return
end subroutine rlmaignanalbe
