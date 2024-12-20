! -----  Init 3D ----
! Generate initial data for opencfd-sc, Ver 1.3
        implicit none
        integer:: nx, ny, nz, j, k, m, N_SPEC,equalpos
        real*8, allocatable, dimension(:, :):: d, u, v, w, T
        real*8, allocatable:: di(:, :, :)
        real*8:: tmp, d1, u1, v1, w1, T1, di1(100)
        character(len=64) :: string,keyString,numberString

        open(12,file='opencfd-comb2.in')
        read(12,*)
        keyString=' '
        do while(trim(keyString).ne.'N_SPEC')
          read(12,'(A)')string
          equalpos = index(string, "=")
          keyString = adjustl(string(1:equalPos-1))
          numberString = adjustl(string(equalPos+1:64))

          if(trim(keyString)=='nx_global') then
            read(numberString, *) nx
          elseif(trim(keyString)=='ny_global') then
            read(numberString, *) ny
          elseif(trim(keyString)=='nz_global') then
            read(numberString, *) nz
          elseif(trim(keyString)=='N_SPEC') then
            read(numberString, *) N_SPEC
          endif
        enddo
        print*,' nx    :',nx
        print*,' ny    :',ny
        print*,' nz    :',nz
        print*,' N_SPEC:',N_SPEC

        close(12)
        print*,' >> opencfd-comb2.in'

        allocate (d(nx, ny), u(nx, ny), v(nx, ny), w(nx, ny), T(nx, ny))
        allocate (di(nx, ny, N_SPEC))
!--------------------

        open (99, file="flow1d-inlet-comb.dat")
        read (99, *)
        do j = 1, ny
           read (99, *) tmp, d1, u1, v1, w1, T1, (di1(m), m=1, N_SPEC)
           d(:, j) = d1
           u(:, j) = u1
           v(:, j) = v1
           w(:, j) = w1
           T(:, j) = T1
           do m = 1, N_SPEC
              di(:, j, m) = di1(m)
           end do
        end do
        close (99)

        open (99, file="opencfd-comb.dat", form="unformatted")
        write (99) 0, 0.d0
        call write3d1(99, nx, ny, nz, d)
        call write3d1(99, nx, ny, nz, u)
        call write3d1(99, nx, ny, nz, v)
        call write3d1(99, nx, ny, nz, w)
        call write3d1(99, nx, ny, nz, T)
        do m = 1, N_SPEC
           call write3d1(99, nx, ny, nz, di(1, 1, m))
        end do
        close (99)
        print*,' << opencfd-comb.dat'
        deallocate (d, u, v, w, T, di)
     end

     subroutine write3d1(no, nx, ny, nz, u)
        implicit none
        integer:: no, nx, ny, nz, k
        real*8:: u(nx, ny)
        do k = 1, nz
           write (no) u
        end do
     end

