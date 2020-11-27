subroutine average_subpixel_rays (npixels,nrays, weighted_stokes, nstokes, &
                              pixel_index, observables)

  implicit none
  integer npixels, nrays, nstokes, pixel_index(nrays)
!f2py intent(in) pixel_index, npixels, nrays
  real weighted_stokes(nstokes, nrays)
!f2py intent(in) stokes
  real observables(nstokes, npixels)
!f2py intent(out) observables

  integer i,j, iray
  double precision temp(nstokes)

  iray=1
  do i=1,npixels
    temp = 0.0D0
    do while (pixel_index(iray) + 1 .EQ. i)
      temp(:) = temp(:) + weighted_stokes(:,iray)
      iray = iray + 1
    enddo
    observables(:,i) = temp(:)
  enddo
  return
end subroutine average_subpixel_rays

subroutine nearest_binary(size, reference, val, index_low, &
  index_high)

  implicit none
  integer size, index_low, index_high, upper, lower
!f2py intent(out) index_low, index_high
  double precision reference(size), val
!f2py intent(in) reference, val
  logical go
  integer index_test

  upper = size
  lower = 1
  go = .true.

  do while (go)
    index_test = (upper + lower) / 2 !rounds down.
    if (reference(index_test) .eq. val) then
      index_low = index_test
      index_high = index_test + 1
      go = .false.
    else if (reference(index_test) > val) then
      upper = index_test
    elseif (reference(index_test) < val) then
      lower = index_test
    endif

    if (upper - lower .eq. 1) then
      index_low = lower
      index_high = upper
      go = .false.
    endif
    enddo

  return
end subroutine nearest_binary

subroutine inner_product(nx1,ny1,nz1, &
  nx2,ny2,nz2, values1, values2, xgrid1,ygrid1,zgrid1, &
  xgrid2,ygrid2,zgrid2, reference, other)

  double precision xgrid1(nx1),ygrid1(ny1),zgrid1(nz1)
!f2py intent(in) xgrid1, ygrid1, zgrid1
  double precision xgrid2(nx2),ygrid2(ny2),zgrid2(nz2)
!f2py intent(in) xgrid2, ygrid2, zgrid2
  double precision values1(nx1,ny1,nz1), values2(nx2,ny2,nz2)
!f2py intent(in) values1, values2
  double precision ref_vol(nx1-1,ny1-1,nz1-1), ref_val(nx1-1,ny1-1,nz1-1)
!f2py intent(out) reference
  double precision other_vol(nx1-1,ny1-1,nz1-1), other_val(nx1-1,ny1-1,nz1-1)
!f2py intent(out) other

  double precision dx1,dy1,dz1,dx2,dy2,dz2
  double precision min_x,min_y,min_z,max_x,max_y,max_z
  integer ic,jc,kc, i,j,k
  integer index_low(3), index_high(3), temp_low, temp_high
  integer gridptr(3,8), t, index_low1(3), index_high1(3)
  double precision f(8), location(3)

  other = 0.0D0
  reference = 0.0D0
  do kc=1,nz1 - 1
    do jc=1,ny1 - 1
      do ic=1,nx1 - 1

        dx1 = xgrid1(ic+1) - xgrid1(ic)
        dy1 = ygrid1(jc+1) - ygrid1(jc)
        dz1 = zgrid1(kc+1) - zgrid1(kc)
        location(1) = xgrid1(ic) + dx1/2
        location(2) = ygrid1(jc) + dy1/2
        location(3) = zgrid1(kc) + dz1/2
        index_low1(1) = ic
        index_low1(2) = jc
        index_low1(3) = kc
        index_high1(1) = ic + 1
        index_high1(2) = jc + 1
        index_high1(3) = kc + 1
        ! call util_locate_point(nx1,ny1,nz1,xgrid1,ygrid1,zgrid1, &
        !     location, gridptr)
        call util_get_interp_kernel2(index_low1, index_high1, location(1), location(2), &
          location(3), f,nx1,ny1,nz1,xgrid1,ygrid1,zgrid1)

        ref_val(ic,jc,kc) = ref_val(ic,jc,kc) +  dx1*dy1*dz1*f(1)* &
                                  values1(ic,jc,kc)
        ref_val(ic,jc,kc) = ref_val(ic,jc,kc) +  dx1*dy1*dz1*f(2)* &
                                  values1(ic+1,jc,kc)
        ref_val(ic,jc,kc) = ref_val(ic,jc,kc) +  dx1*dy1*dz1*f(3)* &
                                  values1(ic+1,jc+1,kc)
        ref_val(ic,jc,kc) = ref_val(ic,jc,kc) +  dx1*dy1*dz1*f(4)* &
                                  values1(ic,jc+1,kc)
        ref_val(ic,jc,kc) = ref_val(ic,jc,kc) +  dx1*dy1*dz1*f(5)* &
                                  values1(ic,jc,kc+1)
        ref_val(ic,jc,kc) = ref_val(ic,jc,kc) +  dx1*dy1*dz1*f(6)* &
                                  values1(ic+1,jc,kc+1)
        ref_val(ic,jc,kc) = ref_val(ic,jc,kc) +  dx1*dy1*dz1*f(7)* &
                                  values1(ic+1,jc+1,kc+1)
        ref_val(ic,jc,kc) = ref_val(ic,jc,kc) +  dx1*dy1*dz1*f(8)* &
                                  values1(ic,jc+1,kc+1)


        ref_vol(ic,jc,kc) = dx1*dy1*dz1

        !for further speed up, replace these searches with an increment scheme
        !as index_low_new = index_high_old when ic is incremented (i think.)
        call nearest_binary(nx2, xgrid2, xgrid1(ic), index_low(1), &
          temp_high)
        call nearest_binary(nx2, xgrid2, xgrid1(ic+1), temp_low, &
          index_high(1))
        call nearest_binary(ny2, ygrid2, ygrid1(jc), index_low(2), &
          temp_high)
        call nearest_binary(ny2, ygrid2, ygrid1(jc+1), temp_low, &
          index_high(2))
        call nearest_binary(nz2, zgrid2, zgrid1(kc), index_low(3), &
          temp_high)
        call nearest_binary(nz2, zgrid2, zgrid1(kc+1), temp_low, &
          index_high(3))

        do k=index_low(3),index_high(3) - 1
            do j=index_low(2),index_high(2) - 1
              do i=index_low(1),index_high(1) - 1

                min_x = min(max(xgrid2(i), xgrid1(ic)),xgrid1(ic+1))
                min_y = min(max(ygrid2(j), ygrid1(jc)),ygrid1(jc+1))
                min_z = min(max(zgrid2(k), zgrid1(kc)),zgrid1(kc+1))

                max_x = max(min(xgrid2(i+1), xgrid1(ic+1)),xgrid1(ic))
                max_y = max(min(ygrid2(j+1), ygrid1(jc+1)),ygrid1(jc))
                max_z = max(min(zgrid2(k+1), zgrid1(kc+1)),zgrid1(kc))
                dx2 = max_x - min_x
                dy2 = max_y - min_y
                dz2 = max_z - min_z

                location(1) = 0.5*(min_x + max_x)
                location(2) = 0.5*(min_y + max_y)
                location(3) = 0.5*(min_z + max_z)
                ! call util_locate_point(nx2,ny2,nz2,xgrid2,ygrid2,zgrid2, &
                !     location, gridptr)
                call util_get_interp_kernel2(index_low, index_high, location(1), location(2), &
                  location(3), f,nx2,ny2,nz2,xgrid2,ygrid2,zgrid2)
                !
                other_val(ic,jc,kc) = other_val(ic,jc,kc) +  dx2*dy2*dz2*f(1)* &
                                          values2(i,j,k)
                other_val(ic,jc,kc) = other_val(ic,jc,kc) +  dx2*dy2*dz2*f(2)* &
                                          values2(i+1,j,k)
                other_val(ic,jc,kc) = other_val(ic,jc,kc) +  dx2*dy2*dz2*f(3)* &
                                          values2(i+1,j+1,k)
                other_val(ic,jc,kc) = other_val(ic,jc,kc) +  dx2*dy2*dz2*f(4)* &
                                          values2(i,j+1,k)
                other_val(ic,jc,kc) = other_val(ic,jc,kc) +  dx2*dy2*dz2*f(5)* &
                                          values2(i,j,k+1)
                other_val(ic,jc,kc) = other_val(ic,jc,kc) +  dx2*dy2*dz2*f(6)* &
                                          values2(i+1,j,k+1)
                other_val(ic,jc,kc) = other_val(ic,jc,kc) +  dx2*dy2*dz2*f(7)* &
                                          values2(i+1,j+1,k+1)
                other_val(ic,jc,kc) = other_val(ic,jc,kc) +  dx2*dy2*dz2*f(8)* &
                                          values2(i,j+1,k+1)

                other_vol(ic,jc,kc)  = other_vol(ic,jc,kc) + dx2*dy2*dz2

              enddo
            enddo
          enddo


      enddo
    enddo
  enddo
end subroutine inner_product


subroutine util_integrate_rays(nrays, nx, ny, nz, xgrid, ygrid, zgrid, &
  start_positions, end_positions, weights, paths)

  implicit none
  integer nrays, nx, ny, nz
  double precision xgrid(nx), ygrid(ny), zgrid(nz)
!f2py intent(in) xgrid, ygrid, zgrid
  double precision start_positions(3, nrays), end_positions(3, nrays)
!f2py intent(in) start_positions, end_positions
  real weights(nx,ny,nz)
!f2py intent(in) weights
  real paths(nrays)
!f2py intent(out) paths

  integer n
  real tau


  paths = -999.0
  do n=1,nrays
    tau = 0.0
    call util_integrate_ray(start_positions(:,n), end_positions(:,n),&
                            tau, weights, xgrid, ygrid, zgrid,     &
                            nx,ny,nz)
    paths(n) = tau
  enddo
end subroutine util_integrate_rays

subroutine util_integrate_ray(start_position, end_position,       &
  tau, weights, xgrid, ygrid, zgrid, nx,ny,nz)

  implicit none
  integer nrays, nx, ny, nz
  double precision xgrid(nx), ygrid(ny), zgrid(nz)
!f2py intent(in) xgrid, ygrid, zgrid
  double precision start_position(3), end_position(3)
!f2py intent(in) start_position, end_position
  real weights(nx,ny,nz), tau, weight1,weight2
!f2py intent(in) weights
!f2py intent(out) tau
  double precision cx, cy, cz, cxinv, cyinv, czinv, eps, r
  double precision sox, soy, soz, so,xn,yn,zn,xe,ye,ze,f(8)
  integer gridptr(3,8), ngrid, maxgrid,i,j
  integer gridptr_end(3,8)
  logical done, a,b,c, final_cell

  tau = 0.0
  eps = 1.0E-5
  ngrid = 0
  maxgrid = 50 * max(nx,ny,nz)

  r = sqrt(sum((end_position - start_position)**2))
  cx = (end_position(1) - start_position(1))/r
  cy = (end_position(2) - start_position(2))/r
  cz = (end_position(3) - start_position(3))/r

  if (abs(cx) .gt. 1.0e-6) then
    cxinv = 1.0D0/cx
  else
    cx = 0.0
    cxinv = 1.0e6
  endif
  if (abs(cy) .gt. 1.0e-6) then
    cyinv = 1.0D0/cy
  else
    cy = 0.0
    cyinv = 1.0e6
  endif
  if (abs(cz) .gt. 1.0e-6) then
    czinv = 1.0D0/cz
  else
    cz = 0.0
    czinv = 1.0e6
  endif

  xe = start_position(1)
  ye = start_position(2)
  ze = start_position(3)

  if (xe .le. xgrid(1)) then
    xe = xgrid(1) + eps
  else if (xe .ge.xgrid(nx)) then
    xe = xgrid(nx) - eps
  endif
  if (ye .le. ygrid(1)) then
    ye = ygrid(1) + eps
  else if (ye .ge.ygrid(ny)) then
    ye = ygrid(nx) - eps
  endif
  if (ze .le. zgrid(1)) then
    ze = zgrid(1) + eps
  else if (ze .ge.zgrid(nz)) then
    ze = zgrid(nz) - eps
  endif

  call util_locate_point(nx, ny, nz, xgrid, ygrid, zgrid,   &
    start_position, gridptr)

  call util_locate_point(nx, ny, nz, xgrid, ygrid, zgrid,   &
    end_position, gridptr_end)

  done = .false.
  do while (.not. done)
    call util_get_interp_kernel(gridptr,xe,ye, ze,f,&
      nx,ny,nz, &
      xgrid, ygrid, zgrid)

    weight1 = 0.0
    do i=1,8
      weight1 = weight1 + f(i)*weights(gridptr(1,i), &
      gridptr(2,i), gridptr(3,i))
    enddo

    if (cx .ge. 1.0e-6) then
      sox = (xgrid(gridptr(1,2)) - xe)*cxinv
    elseif (cx .le. -1.0e-6) then
      sox = (xgrid(gridptr(1,1)) - xe)*cxinv
    else
      sox = 1.0e20
    endif

    if (cy .ge. 1.0e-6) then
      soy = (ygrid(gridptr(2,3)) - ye)*cyinv
    elseif (cx .le. -1.0e-6) then
      soy = (ygrid(gridptr(2,2)) - ye)*cyinv
    else
      soy = 1.0e20
    endif
    if (cz .ge. 1.0e-6) then
      soz = (zgrid(gridptr(3,5)) - ze)*czinv
    elseif (cx .le. -1.0e-6) then
      soz = (zgrid(gridptr(3,4)) - ze)*czinv
    else
      soz = 1.0e20
    endif

    so = min(sox, soy, soz)

    if (so .eq. sox) then
      a = .true.
      b = .false.
      c = .false.
    elseif (so .eq. soy) then
      a = .false.
      b = .true.
      c = .false.
    elseif (so .eq. soz) then
      a = .false.
      b = .false.
      c = .true.
    endif
    xn = xe + so*cx
    yn = ye + so*cy
    zn = ze + so*cz

    final_cell = .true.
    do j=1,8
      do i=1,3
        if (.not. gridptr(i,j) .eq. gridptr_end(i,j)) then
          final_cell = .false.
        endif
      enddo
    enddo

    if (final_cell) then
      xn = end_position(1)
      yn = end_position(2)
      zn = end_position(3)
      done = .true.
    endif

    if (xn > xgrid(nx) + eps) then
      xn = xgrid(nx)
      stop 'went beyond grid +x'
    endif
    if (yn > ygrid(ny) + eps) then
      yn = ygrid(ny)
      stop 'went beyond grid +y'
    endif
    if (zn > zgrid(nz) + eps) then
      zn = zgrid(nz)
      stop 'went beyond grid +z'
    endif
    if (xn < xgrid(1) - eps) then
      xn = xgrid(1)
      stop 'went beyond grid -x'
    endif
    if (yn < ygrid(1)- eps) then
      yn = ygrid(1)
      stop 'went beyond grid -y'
    endif
    if (zn < zgrid(1)- eps) then
      zn = zgrid(1)
      stop 'went beyond grid -z'
    endif

    call util_get_interp_kernel(gridptr,xn,yn, zn,f,&
      nx,ny,nz, &
      xgrid, ygrid, zgrid)

    weight2 = 0.0
    do i=1,8
      weight2 = weight2 + f(i)*weights(gridptr(1,i), &
      gridptr(2,i), gridptr(3,i))
    enddo

    tau = tau + so*0.5*(weight1+weight2)

    if (.not. done) then
      call util_next_cell(a,b,c,cx,cy,cz,gridptr)
    endif

    xe = xn
    ye = yn
    ze = zn
    ngrid = ngrid + 1
    if (ngrid .ge. maxgrid) then
      done = .true.
    endif

  enddo

end subroutine util_integrate_ray

subroutine util_next_cell(a,b,c,cx,cy,cz,gridptr)

  implicit none
  logical a,b,c
!f2py intent(in) a,b,c
  double precision cx, cy, cz
!f2py intent(in) cx, cy, cz
  integer gridptr(3,8)
!f2py intent(in, out) gridptr

  integer i

  if (a) then
    if (cx .ge. 1e-6) then
      do i=1,8
        gridptr(1,i) = gridptr(1,i) + 1
      enddo
    elseif (cx .le. -1e-6) then
      do i=1,8
        gridptr(1,i) = gridptr(1,i) - 1
      enddo
    endif
  endif
  if (b) then
    if (cy .ge. 1e-6) then
      do i=1,8
        gridptr(2,i) = gridptr(2,i) + 1
      enddo
    elseif (cy .le. -1e-6) then
      do i=1,8
        gridptr(2,i) = gridptr(2,i) - 1
      enddo
    endif
  endif
  if (c) then
    if (cz .ge. 1e-6) then
      do i=1,8
        gridptr(3,i) = gridptr(3,i) + 1
      enddo
    elseif (cz .le. -1e-6) then
      do i=1,8
        gridptr(3,i) = gridptr(3,i) - 1
      enddo
    endif
  endif
  return

end subroutine util_next_cell


subroutine util_locate_point(nx, ny, nz, xgrid, ygrid, zgrid, &
  location, gridptr)
  implicit none

  integer nx, ny, nz, gridptr(3,8)
!f2py intent(out) gridptr
  double precision xgrid(nx), ygrid(ny), zgrid(nz)
!f2py intent(in) xgrid, ygrid, zgrid
  double precision location(3)
!f2py intent(in) location

  integer i,j
  gridptr = 1

  call nearest_binary(nx, xgrid, location(1),i,j)
  gridptr(1,1) = i
  gridptr(1,2) = j
  gridptr(1,3) = j
  gridptr(1,4) = i
  gridptr(1,5) = i
  gridptr(1,6) = j
  gridptr(1,7) = j
  gridptr(1,8) = i

  call nearest_binary(ny, ygrid, location(2),i,j)
  gridptr(2,1) = i
  gridptr(2,2) = i
  gridptr(2,3) = j
  gridptr(2,4) = j
  gridptr(2,5) = i
  gridptr(2,6) = j
  gridptr(2,7) = j
  gridptr(2,8) = i

  call nearest_binary(nz, zgrid, location(3),i,j)
  gridptr(3,1) = i
  gridptr(3,2) = i
  gridptr(3,3) = i
  gridptr(3,4) = i
  gridptr(3,5) = j
  gridptr(3,6) = j
  gridptr(3,7) = j
  gridptr(3,8) = j

  ! do i=1,nx
  !   if ((location(1) .ge. xgrid(i)) .and.  &
  !         (location(1) .le. xgrid(i+1))) then
  !     gridptr(1,1) = i
  !     gridptr(1,2) = i + 1
  !     gridptr(1,3) = i + 1
  !     gridptr(1,4) = i
  !     gridptr(1,5) = i
  !     gridptr(1,6) = i + 1
  !     gridptr(1,7) = i + 1
  !     gridptr(1,8) = i
  !   endif
  ! enddo
  !
  ! do i=1,ny
  !   if ((location(2) .ge. ygrid(i)) .and.  &
  !         (location(2) .le. ygrid(i+1))) then
  !     gridptr(2,1) = i
  !     gridptr(2,2) = i
  !     gridptr(2,3) = i + 1
  !     gridptr(2,4) = i + 1
  !     gridptr(2,5) = i
  !     gridptr(2,6) = i + 1
  !     gridptr(2,7) = i + 1
  !     gridptr(2,8) = i
  !   endif
  ! enddo
  !
  ! do i=1,nz
  !   if ((location(3) .ge. zgrid(i)) .and.  &
  !         (location(3) .le. zgrid(i+1))) then
  !     gridptr(3,1) = i
  !     gridptr(3,2) = i
  !     gridptr(3,3) = i
  !     gridptr(3,4) = i
  !     gridptr(3,5) = i + 1
  !     gridptr(3,6) = i + 1
  !     gridptr(3,7) = i + 1
  !     gridptr(3,8) = i + 1
  !   endif
  ! enddo
  return
end subroutine util_locate_point

subroutine util_get_interp_kernel2(index_low,index_high,xn,yn, zn,f,nx,ny,nz, &
  xgrid, ygrid, zgrid)
  implicit none

  integer nx,ny,nz,index_low(3), index_high(3)
!f2py intent(in) index_low, index_high
  double precision xn, yn, zn, f(8)
!f2py intent(in) xn, yn, zn
!f2py intent(out) f
  double precision xgrid(nx), ygrid(ny), zgrid(nz)
!f2py intent(in) xgrid, ygrid, zgrid
  double precision u,v,w, delx,dely, delz
  double precision invdelx, invdely, invdelz

  delx = xgrid(index_high(1)) - xgrid(index_low(1))
  dely = ygrid(index_high(2)) - xgrid(index_low(2))
  delz = zgrid(index_high(3)) - xgrid(index_low(3))

  if (delx .le. 0.0) then
    invdelx = 1.0
  else
    invdelx = 1.0D0/delx
  endif
  if (dely .le. 0.0) then
    invdely = 1.0
  else
    invdely = 1.0D0/dely
  endif
  if (delz .le. 0.0) then
    invdelz = 1.0
  else
    invdelz = 1.0D0/delz
  endif

  u = (xn - xgrid(index_low(1)))*invdelx
  v = (yn - ygrid(index_low(2)))*invdely
  w = (zn - zgrid(index_low(3)))*invdely

  f(1) = (1-W)*(1-V)*(1-U)
  f(2) = (1-W)*(1-V)*U
  f(3) = (1-W)*V*(1-U)
  f(4) = (1-W)*V*U
  f(5) = W*(1-V)*(1-U)
  f(6) = W*(1-V)*U
  f(7) = W*V*(1-U)
  f(8) = W*V*U

  return
end subroutine util_get_interp_kernel2




subroutine util_get_interp_kernel(gridptr,xn,yn, zn,f,nx,ny,nz, &
  xgrid, ygrid, zgrid)
  implicit none

  integer nx,ny,nz, gridptr(3,8)
!f2py intent(in) gridptr
  double precision xn, yn, zn, f(8)
!f2py intent(in) xn, yn, zn
!f2py intent(out) f
  double precision xgrid(nx), ygrid(ny), zgrid(nz)
!f2py intent(in) xgrid, ygrid, zgrid
  double precision u,v,w, delx,dely, delz
  double precision invdelx, invdely, invdelz

  delx = xgrid(gridptr(1,3)) - xgrid(gridptr(1,2))
  dely = ygrid(gridptr(2,2)) - xgrid(gridptr(2,1))
  delz = zgrid(gridptr(3,5)) - xgrid(gridptr(3,4))

  if (delx .le. 0.0) then
    invdelx = 1.0
  else
    invdelx = 1.0D0/delx
  endif
  if (dely .le. 0.0) then
    invdely = 1.0
  else
    invdely = 1.0D0/dely
  endif
  if (delz .le. 0.0) then
    invdelz = 1.0
  else
    invdelz = 1.0D0/delz
  endif

  u = (xn - xgrid(gridptr(1,2)))*invdelx
  v = (yn - ygrid(gridptr(2,1)))*invdely
  w = (zn - zgrid(gridptr(3,1)))*invdely

  f(1) = (1-W)*(1-V)*(1-U)
  f(2) = (1-W)*(1-V)*U
  f(3) = (1-W)*V*(1-U)
  f(4) = (1-W)*V*U
  f(5) = W*(1-V)*(1-U)
  f(6) = W*(1-V)*U
  f(7) = W*V*(1-U)
  f(8) = W*V*U

  return
end subroutine util_get_interp_kernel
