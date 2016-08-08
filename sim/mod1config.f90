module sysconfig
! module for config of cell cluster

! b8 will be used to define reals with 14 digits
integer,  parameter :: b8 = selected_real_kind(14)
! set radius of the cell
real(b8), parameter :: rCell = 0.200_b8

contains

    ! initialize cell positions to form a cluster
    ! set rsim space such that the cluster is located at the center
    subroutine itlClusterSys( cellTotal, cellArray, rsim)
        implicit none
        integer,  intent(in)  :: cellTotal
        real(b8), intent(out) :: rsim(:,:)
        real(b8), intent(out) :: cellArray(:,:,:)
        real(b8), allocatable :: center(:,:), dr(:,:), test(:)
        real(b8) :: clusterLength, r
        integer :: cellCheck, i, ii, j, k, n

        allocate( dr(6,3) )
        allocate( test(3) )
        allocate( center(cellTotal,3) )
        center(:,:)      = 0.0_b8
        cellArray(:,:,:) = 0.0_b8
        ! set dr, the displacement vector between cell centers
        dr(:,:) = 0.0_b8
        do i = 1, 3
            dr(i,i) = rCell * 2.0_b8
        enddo
        do i = 4, 6
            dr(i,i-3) = -rCell * 2.0_b8
        enddo
        ! set center and cellArray for cell 1
        do i = 1, 3
            center(1,i) = 0.0_b8
        enddo
        do i = 1, 3
            cellArray(1,i,1) = center(1,i) - rCell
            cellArray(1,i,2) = center(1,i) + rCell
        enddo
        if ( cellTotal > 1 ) then
            ! set center and cellArray for cells 2 to cellTotal
            n = 1
            i = 2
            do while ( i <= cellTotal .AND. n <= cellTotal )
                call random_number(r)
                j = 1 + floor( r * 6.0 )
                do k = 1, 6
                    ! set test location for cell center i
                    test(1) = center(n,1) + dr(j,1)
                    test(2) = center(n,2) + dr(j,2)
                    test(3) = center(n,3) + dr(j,3)
                    cellCheck = 0
                    do ii = 1, i-1
                        ! check that no other cell center occupy test location
                        if ( center(ii,1) == test(1) .AND. center(ii,2) == test(2) .AND. center(ii,3) == test(3) ) then
                            cellCheck = 1
                        end if
                    enddo
                    if ( cellCheck == 0 ) then
                        center(i,:) = test
                        do ii = 1, 3
                            cellArray(i,ii,1) = center(i,ii) - rCell
                            cellArray(i,ii,2) = center(i,ii) + rCell
                        enddo
                        i = i + 1
                    end if
                    if ( i > cellTotal ) then
                        exit
                    else
                        j = j + 1
                        if( j > 6 )then
                            j = 1
                        endif
                    endif
                enddo
                n = n + 1
            enddo
        end if
        ! cluster is now centered at [ 0.0, 0.0, 0.0]
        ! get cluster size, set system size, shift cluster to be centerd in the system
        do i = 1, 3
            clusterLength = maxval( cellArray(:,i,2)) - minval( cellArray(:,i,1))
            rsim(i,1) = 0.0_b8
            rsim(i,2) = clusterLength + (rCell * 4.0_b8)
            do j = 1, cellTotal
                cellArray(j,i,1) = cellArray(j,i,1) + ((rsim(i,2) - rsim(i,1)) / 2.0_b8)
                cellArray(j,i,2) = cellArray(j,i,2) + ((rsim(i,2) - rsim(i,1)) / 2.0_b8)
            enddo
        enddo
        deallocate( dr )
        deallocate( test )
        deallocate( center )
    end subroutine itlClusterSys


    ! output cell locations and system size
    subroutine wrtOutClusterSys( cellTotal, cellArray, rsim)
        implicit none
        integer,  intent(in) :: cellTotal
        real(b8), intent(in) :: cellArray(:,:,:), rsim(:,:)
        integer :: i
        write(*,*)
        do i = 1, cellTotal
            write(*,*) 'cell', i
            write(*,*) '   x', cellArray(i,1,:), cellArray(i,1,2)-cellArray(i,1,1)
            write(*,*) '   y', cellArray(i,2,:), cellArray(i,2,2)-cellArray(i,2,1)
            write(*,*) '   z', cellArray(i,3,:), cellArray(i,3,2)-cellArray(i,3,1)
        enddo
        write(*,*)
        write(*,*) 'rsim x', rsim(1,:)
        write(*,*) 'rsim y', rsim(2,:)
        write(*,*) 'rsim z', rsim(3,:)
    end subroutine wrtOutClusterSys


    ! initialize positions of cells to form a cluster
    subroutine itlCellCluster( cellTotal, cellArray, rsim)
        implicit none
        integer,  intent(in)  :: cellTotal
        real(b8), intent(out) :: cellArray(:,:,:)
        real(b8), intent(in)  :: rsim(:,:)
        ! real(b8) :: rCell, test(3), center(cellTotal,3)
        real(b8) :: r, rCell
        real(b8), allocatable :: center(:,:), dr(:,:), test(:)
        integer :: cellCheck, i, ii, j, k, n

        allocate( dr(6,3) )
        allocate( test(3) )
        allocate( center(cellTotal,3) )
        center(:,:)      = 0.0_b8
        cellArray(:,:,:) = 0.0_b8
        ! set radius of the cell
        rCell = 0.200_b8
        ! set dr
        dr(:,:) = 0.0_b8
        do i = 1, 3
            dr(i,i) = 2.0 * rCell
        enddo
        do i = 4, 6
            dr(i,i-3) = -2.0 * rCell
        enddo
        ! set center and cellArray for cell 1
        do i = 1, 3
            center(1,i) = (rsim(i,2) - rsim(i,1)) / 2.0_b8
        enddo
        do i = 1, 3
            cellArray(1,i,1) = center(1,i) - rCell
            cellArray(1,i,2) = center(1,i) + rCell
        enddo
        if ( cellTotal == 1 ) then
            return
        end if
        ! set center and cellArray for cells 2 to cellTotal
        n = 1
        i = 2
        do while ( i <= cellTotal .AND. n <= cellTotal )
            call random_number(r)
            j = 1 + floor( r * 6.0 )
            do k = 1, 6
                ! set test location for cell center i
                test(1) = center(n,1) + dr(j,1)
                test(2) = center(n,2) + dr(j,2)
                test(3) = center(n,3) + dr(j,3)
                cellCheck = 0
                do ii = 1, i-1
                    ! check that no other cell center occupy test location
                    if ( center(ii,1) == test(1) .AND. center(ii,2) == test(2) .AND. center(ii,3) == test(3) ) then
                        cellCheck = 1
                    end if
                enddo
                if ( cellCheck == 0 ) then
                    center(i,:) = test
                    do ii = 1, 3
                        cellArray(i,ii,1) = test(ii) - rCell
                        cellArray(i,ii,2) = test(ii) + rCell
                    enddo
                    i = i + 1
                end if
                if ( i > cellTotal ) then
                    exit
                else
                    j = j + 1
                    if( j > 6 )then
                        j = 1
                    endif
                endif
            enddo
            n = n + 1
        enddo
        deallocate( dr )
        deallocate( test )
        deallocate( center )
    end subroutine itlCellCluster


    ! make a list of cells which are on the perimeter
    ! edgeList(i) = 0 means that cell i is not an edge cell
    ! edgeList(i) = 1 means that cell i is an edge cell
    subroutine clusterEdgeList( cellTotal, cellArray, rsim, edgeList)
        implicit none
        integer,  intent(in)  :: cellTotal
        integer,  intent(out) :: edgeList(:)
        real(b8), intent(in)  :: cellArray(:,:,:), rsim(:,:)
        real(b8), allocatable :: center(:,:), dr(:,:)
        real(b8) :: dcell
        integer  :: i, j, k, n1, n2, check

        allocate( dr(6,3) )
        allocate( center(cellTotal,3) )
        ! set dr
        dr(:,:) = 0.0_b8
        do i = 1, 3
            dr(i,i) = 2.0 * rCell
        enddo
        do i = 4, 6
            dr(i,i-3) = -2.0 * rCell
        enddo
        ! set cell center from cellArray
        center(:,:)      = 0.0_b8
        do i = 1, cellTotal
            do j = 1, 3
                center(i,j) = cellArray(i,j,1) + (cellArray(i,j,2) - cellArray(i,j,1)) / 2.0_b8
            enddo
            write(*,*) 'cell', i, 'center =', center(i,:)
        enddo
        edgeList(:) = 1
        do n1 = 1, cellTotal
            check = 0
            do n2 = 1, cellTotal
                if ( n1 == n2 ) then
                    cycle
                end if
                dcell = 0.0_b8
                do j = 1, 3
                    dcell = dcell + (center(n1,j) - center(n2,j))**2
                enddo
                dcell = dsqrt(dcell)
                if ( dcell >= 2.0*rCell-0.000001 .AND. dcell <= 2.0*rCell+0.000001 ) then
                    check = check + 1
                    write(*,*) 'cell', n1, 'check =', check, 'cell 2 =', n2
                end if
                if ( check > 5 ) then
                    edgeList(n1) = 0
                    exit
                end if
            enddo
        enddo
        deallocate( dr )
        deallocate( center )
    end subroutine clusterEdgeList


    ! initialize 1 cell. rsim must be set beforehand.
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


end module
