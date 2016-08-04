program diff
! test diffusion in a 3d lattice

implicit none

integer, parameter:: b8 = selected_real_kind(14)

integer :: i, j, nt, ntTotal, run, runTotal
integer :: prtclTotal
real(b8) :: xmin, xmax, ymin, ymax, zmin, zmax
real(b8) :: r
real(b8) :: dr(3), rsim(3,2)
real(b8), allocatable :: prtclArray(:,:)

integer :: cellTotal
real(b8) :: meanCount, varCount
integer,  allocatable :: countArray(:)
real(b8), allocatable :: cellArray(:,:,:), cellPolar(:,:), timeCount(:)


! open/create file
open(unit=10, file='prtclLocation.dat', action='write', status='replace')
write(10,*) ' '
close(10)

! set total number of runs
runTotal = 200
ntTotal  = 3000
! set total possible number of particles in system
prtclTotal = 1000
! initialize simulation space size
do i = 1, 3
    rsim(i,1) = 0.0_b8
    rsim(i,2) = 1.0_b8
end do
! initialize particle movement step size
dr(1) = (rsim(1,2) - rsim(1,1)) / 50.0_b8
dr(2) = (rsim(2,2) - rsim(2,1)) / 50.0_b8
dr(3) = (rsim(3,2) - rsim(3,1)) / 50.0_b8
! set total number of cell in system
cellTotal = 1
! allocate memory
allocate( prtclArray( prtclTotal, 4))
allocate( cellArray( cellTotal, 3, 2))
allocate( cellPolar( cellTotal, 3))
allocate( countArray( cellTotal))
allocate( timeCount( ntTotal))

write(*,*) 'prtclTotal =', prtclTotal
write(*,*) '  Ninitial =', prtclTotal/2
do i = 1, 3
    write(*,*) 'dr =', dr(1)
enddo
write(*,*)

call init_random_seed()

do run = 1, runTotal
    write(*,*) ' run', run

    timeCount(:)     = 0.0_b8
    countArray(:)    = 0
    cellPolar(:,:)   = 0.0_b8
    cellArray(:,:,:) = 0.0_b8
    ! initialize cell position
    call initOneCell( cellTotal, rsim, cellArray)
    ! initialize particle positions
    nt = 1
    prtclArray(:,:) = 0.0_b8
    do i = 1, (prtclTotal/2)
        prtclArray(i,4) = 1.0_b8
        do j = 1, 3
            call random_number(r) ! in cpm code use ran1() function
            prtclArray(i,j) = r *(rsim(j,2) - rsim(j,1)) + rsim(j,1)
        enddo
    enddo
    ! call wrtPrtclLocation( prtclTotal, nt, prtclArray)

    ! move particles for some number of timesteps
    do nt = 2, ntTotal+3000
    ! do while( nt < 40)
        ! nt = nt + 1
        ! update particle location and check boundary conditions
        call prtclUpdate( prtclTotal, dr, rsim, prtclArray)
        ! write(*,*) '  move  '
        do i = 1, prtclTotal
            if ( prtclArray(i,4) == 1.0_b8 ) then
                ! write(*,*) 'i=', i, prtclArray(i,1:3)
            end if
        enddo
        ! add flux of particles
        call prtclFlux( prtclTotal, dr, rsim, prtclArray)

        ! if ( mod(nt,3000) == 0 ) then
        !     call wrtPrtclLocation( prtclTotal, nt, prtclArray)
        ! end if
        ! INSTANTANEOUS
        ! if ( nt == ntTotal ) then
        !     call cellCount( cellTotal, prtclTotal, cellArray, prtclArray, countArray)
        !     write(100,*) countArray(1), run
        !     exit
        ! end if
        ! LONG TIME: count the particles within a cell
        ! if ( nt >= 3000 .AND. mod(nt,50) == 0 ) then
        if ( nt > 3000 ) then
        ! if ( nt >= 300 .AND. nt < 301 ) then
            call cellCount( cellTotal, prtclTotal, cellArray, prtclArray, countArray)
            timeCount(nt-3000) = float(countArray(1))
            ! write(100+run,*) countArray(1), nt-3000
            ! call cellpolarMW( cellTotal, prtclTotal, cellArray, prtclArray, cellPolar)
            ! call wrtPlrTotal( run, cellTotal, cellPolar, nt-3000)
        end if
    enddo

    ! calculate and output time averaged molecule count
    meanCount = sum(timeCount) / float(ntTotal)
    write(100,*) meanCount, run

enddo

contains
    ! gradient is assumed to be in the x-direction (1)
    ! add flux of particles at top x-boundary
    subroutine prtclFlux( N, dr, rsim, prtclArray)
        implicit none
        integer,  intent(in)    :: N
        real(b8), intent(in)    :: dr(:), rsim(:,:)
        real(b8), intent(inout) :: prtclArray(:,:)
        real(b8) :: a, d, g
        integer  :: i, j, k, nJ
        ! add flux to top boundary
        nJ = N / 100 ! nJ is the number of particles added
        if ( nJ == 0 ) then
            nJ = 1
        end if
        do i = 1, nJ
            j = 1
            do while( prtclArray(j,4) == 1.0_b8 )
                j = j + 1
                if ( j > N ) then
                    exit
                end if
            enddo
            if ( j > N ) then
                exit
            end if
            prtclArray(j,4) = 1.0_b8
            call random_number(r)
            prtclArray(j,1) = rsim(1,2) - r*dr(1)
            do k = 2, 3
                call random_number(r) ! cpm code used ran1() function
                prtclArray(j,k) = r *(rsim(k,2) - rsim(k,1)) + rsim(k,1)
            enddo
        enddo
    end subroutine prtclFlux


    ! move particles and check boundary conditions
    subroutine prtclUpdate( N, dr, rsim, prtclArray)
        implicit none
        integer,  intent(in)    :: N
        real(b8), intent(in)    :: dr(:), rsim(:,:)
        real(b8), intent(inout) :: prtclArray(:,:)
        integer  :: i, dim
        real(b8) :: r

        do i = 1, N
            ! check whether array index corresponds to a particle
            if ( prtclArray(i,4) == 1.0_b8 ) then
                ! move the particle by distance dr(j) in random direction
                call random_number(r)
                if ( r < (1.0_b8/6.0_b8) ) then
                    dim = 1
                    prtclArray(i,1) = prtclArray(i,1) - dr(1)
                elseif( r < (2.0_b8/6.0_b8) )then
                    dim = 1
                    prtclArray(i,1) = prtclArray(i,1) + dr(1)
                elseif( r < (3.0_b8/6.0_b8) )then
                    dim = 2
                    prtclArray(i,2) = prtclArray(i,2) - dr(2)
                elseif( r < (4.0_b8/6.0_b8) )then
                    dim = 2
                    prtclArray(i,2) = prtclArray(i,2) + dr(2)
                elseif( r < (5.0_b8/6.0_b8) )then
                    dim = 3
                    prtclArray(i,3) = prtclArray(i,3) - dr(3)
                else
                    dim = 3
                    prtclArray(i,3) = prtclArray(i,3) + dr(1)
                endif

                ! check boundary conditions
                if ( dim /= 1 ) then
                ! if ( dim > 0 ) then
                    ! reflective boundaries perpendicular to gradient
                    do j = 2, 3
                        if( prtclArray(i,j) < rsim(j,1) )then
                            ! write(*,*) 'too small', prtclArray(i,j), i, j
                            prtclArray(i,j) = rsim(j,2) - (rsim(j,1) - prtclArray(i,j))
                            exit
                        elseif( prtclArray(i,j) > rsim(j,2) )then
                            ! write(*,*) 'too big', prtclArray(i,j), i, j
                            prtclArray(i,j) = (prtclArray(i,j) - rsim(j,2)) + rsim(j,1)
                            exit
                        endif
                    enddo
                else
                    ! absorbing boundary parallel to gradient
                    if( (prtclArray(i,dim)<rsim(dim,1)) .OR. (prtclArray(i,dim)>rsim(dim,2)) )then
                        ! write(*,*) 'absorbed', prtclArray(i,dim), i
                        prtclArray(i,4) = 0.0_b8
                    endif
                endif
            endif
        enddo
    end subroutine prtclUpdate


    ! output particle location data
    subroutine wrtPrtclLocation( N, nt, prtclArray)
        implicit none
        integer,  intent(in) :: N, nt
        real(b8), intent(in) :: prtclArray(:,:)
        integer :: i
        open(unit=11, file='prtclLocation.dat', action='write', status='old')
        do i = 1, N
            if ( prtclArray(i,4) == 1.0_b8 ) then
                write(11,*) prtclArray(i,1:3), nt
            end if
        enddo
    end subroutine wrtPrtclLocation


    ! calculate MW cell polarization
    subroutine cellpolarMW( cellTotal, prtclTotal, cellArray, prtclArray, cellPolar)
        implicit none
        integer,  intent(in)  :: cellTotal, prtclTotal
        real(b8), intent(in)  :: cellArray(:,:,:), prtclArray(:,:)
        real(b8), intent(out) :: cellPolar(:,:)
        real(b8) :: center(3), check, q(3), qtot(3)
        integer :: i, j, k
        cellPolar(:,:) = 0.0_b8
        qtot = 0.0_b8
        do i = 1, cellTotal
            do j = 1, 3
                center(j) = (cellArray(i,j,2) - cellArray(i,j,1)) / (2.0_b8)
            enddo
            qtot = 0.0_b8
            do j = 1, prtclTotal
                q = 0.0_b8
                if ( prtclArray(j,4) == 1.0_b8 ) then
                    check = (prtclArray(j,1)-cellArray(i,1,1))*(prtclArray(j,1)-cellArray(i,1,2))
                    if ( check < 0.0 ) then
                        check = (prtclArray(j,2)-cellArray(i,2,1))*(prtclArray(j,2)-cellArray(i,2,2))
                        if ( check < 0.0 ) then
                            check = (prtclArray(j,3)-cellArray(i,3,1))*(prtclArray(j,3)-cellArray(i,3,2))
                            if ( check < 0.0 ) then
                                do k = 1, 3
                                    q(k) = prtclArray(j,k) - center(k)
                                enddo
                                q = q / sqrt(dot_product(q,q))
                                ! q = q
                            endif
                        endif
                    endif
                endif
                qtot = qtot + q
            enddo
            ! cellPolar(i,:) = qtot / sqrt(dot_product(qtot,qtot))
            cellPolar(i,:) = qtot
        enddo

    end subroutine cellpolarMW


    ! count the number of particels within cells volume
    subroutine cellCount( cellTotal, prtclTotal, cellArray, prtclArray, countArray)
        implicit none
        integer,  intent(in)  :: cellTotal, prtclTotal
        real(b8), intent(in)  :: cellArray(:,:,:), prtclArray(:,:)
        integer,  intent(out) :: countArray(:)
        integer  :: i, j
        real(b8) :: check
        do i = 1, cellTotal
            countArray(i) = 0
            do j = 1, prtclTotal
                if ( prtclArray(j,4) == 1.0_b8 ) then
                    check = (prtclArray(j,1)-cellArray(i,1,1))*(prtclArray(j,1)-cellArray(i,1,2))
                    if ( check < 0.0 ) then
                        ! write(*,*) 'passed x-check:', prtclArray(j,1:3)
                        check = (prtclArray(j,2)-cellArray(i,2,1))*(prtclArray(j,2)-cellArray(i,2,2))
                        if ( check < 0.0 ) then
                            ! write(*,*) 'passed y-check:', prtclArray(j,1:3)
                            check = (prtclArray(j,3)-cellArray(i,3,1))*(prtclArray(j,3)-cellArray(i,3,2))
                            if ( check < 0.0 ) then
                                ! write(*,*) 'passed z-check:', prtclArray(j,1:3)
                                countArray(i) = countArray(i) + 1
                            endif
                        endif
                    endif
                endif
            enddo
        enddo
    end subroutine cellCount


    subroutine initOneCell( N, rsim, cellArray)
        implicit none
        integer,  intent(in)  :: N
        real(b8), intent(in)  :: rsim(:,:)
        real(b8), intent(out) :: cellArray(:,:,:)
        real(b8) :: a, center(3)
        integer :: i
        if ( N > 1 ) then
            write(*,*) 'error: cellTotal > 1'
            write(*,*) 'only 1 cell will be initialized'
            write(*,*)
        end if
        ! calculate the center of the cell
        do i = 1, 3
            center(i) = (rsim(i,2) - rsim(i,1)) / 2.0_b8
        enddo
        ! set the cell length
        a = (rsim(1,2) - rsim(1,1)) / 2.0_b8
        do i = 1, 3
            ! cellArray(1,i,1) = center(i) - ( a / 2.0_b8 )
            ! cellArray(1,i,2) = center(i) + ( a / 2.0_b8 )
            cellArray(1,i,1) = center(i) - ( 0.250_b8 )
            cellArray(1,i,2) = center(i) + ( 0.250_b8 )
        enddo
    end subroutine initOneCell


    ! output total cluster polarization to fort.141
    subroutine wrtPlrTotal( nRun, cellTotal, cellPolar, nt)
        implicit none
        integer,  intent(in) :: cellTotal, nRun, nt
        real(b8), intent(in), dimension(:,:) :: cellPolar
        real(b8) :: px, py, pz
        integer  :: i

        px = 0.0_b8
        py = 0.0_b8
        pz = 0.0_b8
        do i = 1, cellTotal
            px = px + cellPolar(i,1)
            py = py + cellPolar(i,2)
            pz = pz + cellPolar(i,3)
        enddo
        write(100+nRun,"(E16.8)", advance="no") px
        write(100+nRun,"(E17.8)", advance="no") py
        write(100+nRun,"(E17.8)", advance="no") pz
        write(100+nRun,"(I10)", advance="no")   nt
        write(100+nRun,*) ''
    end subroutine wrtPlrTotal
end program


! initialize RANDOM_SEED
subroutine init_random_seed()
    integer :: values(1:8), k
    integer, dimension(:), allocatable :: seed
    real(8) :: r

    call date_and_time(values=values)

    call random_seed(size=k)
    allocate(seed(1:k))
    seed(:) = values(8)
    call random_seed(put=seed)
end subroutine init_random_seed


! returns random number between 0 - 1
! real(b8) function ran1()
!     implicit none
!     real(b8) x
!     call random_number(x) ! built in fortran 90 random number function
!     ran1=x
! end function ran1
