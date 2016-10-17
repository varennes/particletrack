module sysconfig
! module for config of cell cluster
use parameters

!!!  SIMULATION PARAMETERS  [start] !!!
integer,  parameter ::   geoTotal = 1      ! total number of cluster geometries to iterate through
integer,  parameter ::   runTotal = 1      ! total number of runs
integer,  parameter ::  cellTotal = 10     ! total number of cells in system
integer,  parameter :: prtclTotal = 10000  ! total possible number of particles in system
!!!  SIMULATION PARAMETERS  [end]   !!!

contains

    subroutine getProbTimeScale( ntItl, dtReal, p, q)
        implicit none
        real(b8), intent(out) :: dtReal, p, q
        integer,  intent(out) :: ntItl
        ! calculate time step in real units
        dtReal = min( bReal*bReal / dReal / 6.0_b8, 1.0_b8 / kReal)
        ! calculate probailities of diffusion and production events
        p = 6.0_b8 * dReal * dtReal / bReal**2
        q = kReal * dtReal
        ! calculate the number of timesteps needed for equilibration
        ntItl = 10 * ceiling( max( lReal**2 / (dReal * dtReal), (kReal * dtReal)**(-1)) )
    end subroutine getProbTimeScale


    subroutine getSysLengthScales( dr, rsim)
        implicit none
        real(b8) , intent(out) :: dr(3), rsim(3,2)

        dr(:) = bReal
        rsim(:,1) = 0.0_b8
        rsim(1,2) =  lReal
        rsim(2,2) = syReal
        rsim(3,2) = szReal
    end subroutine getSysLengthScales


    ! initialize 3D cluster of cells using nearest-neighbor (NN) method
    subroutine itl3DClusterNN( cellArray, rsim)
        implicit none
        real(b8), intent(out) :: cellArray(:,:,:)
        real(b8), intent(in)  :: rsim(:,:)
        real(b8), allocatable :: dList(:), dTest(:)
        integer,  allocatable :: nnList(:,:), nnTest(:,:)
        real(b8) :: a, nnd, nndmax, x0, y0, z0
        integer :: center( cellTotal, 3), test(3)
        integer :: i, itest, inn, j, k, k1, k2, k3, kmax, n, centerCheck

        allocate(  dList(400))
        allocate(  dTest(400))
        allocate( nnList(400,3))
        allocate( nnTest(400,3))
        cellArray(:,:,:) = 0.0_b8
        ! set center and cellArray for cell 1
        x0 = (rsim(1,2) - rsim(1,1)) / 2.0_b8
        y0 = (rsim(2,2) - rsim(2,1)) / 2.0_b8
        z0 = (rsim(3,2) - rsim(3,1)) / 2.0_b8
        do i = 1, 3
            cellArray(1,i,1) = (rsim(i,2) - rsim(i,1)) / 2.0_b8 - rReal
            cellArray(1,i,2) = (rsim(i,2) - rsim(i,1)) / 2.0_b8 + rReal
        enddo
        center(:,:) = 0
        ! write(*,*) center(1,1:3)
        if ( cellTotal == 1 ) then
            return
        end if
        n = 2
        kmax = 1
        dList = 0.0_b8
        nnList = 0
        do while ( n <= cellTotal )
            ! check for free NN sites
            ! write(*,*) '  kmax =', kmax
            k = 1
            test = 0
            nnList = 0
            dList  = 0.0_b8
            do k1 = kmax, -kmax, -1
                do k2 = kmax, -kmax, -1
                    do k3 = kmax, -kmax, -1
                        do j = 1, n-1
                            centerCheck = 0
                            test(1:3) = [center(1,1) + k1, center(1,2) + k2, center(1,3) + k3]
                            if ( center(j,1) == test(1) .AND. center(j,2) == test(2) .AND. center(j,3) == test(3) ) then
                                centerCheck = 1
                                exit
                            endif
                        enddo
                        if ( centerCheck == 0 ) then
                            nnList(k,:) = test
                            dList(k) = 0.0
                            do j = 1, 3
                                dList(k) = dList(k) + float(nnList(k,j)-center(1,j))**2
                            enddo
                            dList(k) = sqrt( dList(k) )
                            k = k + 1
                        endif
                    enddo
                enddo
            enddo
            inn = k - 1
            nndmax = float(ceiling(minval( dList(1:inn)) + 0.00001))
            ! write(*,*)
            ! do i = 1, inn
            !     write(*,*) i, nnList(i,:), dList(i)
            ! enddo
            ! write(*,*)
            do while ( inn > 0 .AND. n <= cellTotal )
                nnd = minval( dList(1:inn))
                if ( nnd == maxval( dList(1:inn)) .OR. nnd >= nndmax ) then
                    exit
                end if
                ! write(*,*) 'nnd =', nnd, '  inn =', inn
                 dTest =  dList
                nnTest = nnList
                do i = 1, inn
                    if ( dList(i) == nnd .AND. n <= cellTotal) then
                        center(n,1:3) = nnList(i,1:3)
                        n = n + 1
                        ! remove nnList and dList entries that have been added to center
                         dTest(i)   = 0.0
                        nnTest(i,:) = [ 0, 0, 0]
                    end if
                enddo
                j = inn
                do i = inn, 1, -1
                    if ( dTest(i) == 0.0 ) then
                        if ( i /= inn ) then
                             dList(i:inn-1)   =  dList(i+1:inn)
                            nnList(i:inn-1,:) = nnList(i+1:inn,:)
                        end if
                        j = j - 1
                    end if
                enddo
                inn = j
                ! write(*,*)
                ! do i = 1, inn
                !     write(*,*) i, nnList(i,:), dList(i)
                ! enddo
                ! write(*,*)
            enddo
            kmax = kmax + 1
        enddo
        itest = cellTotal
        do i = 2, itest
            cellArray(i,1,1:2) = [ x0 + (2.0*float(center(i,1))-1.0) * rReal, x0 + (2.0*float(center(i,1))+1.0) * rReal]
            cellArray(i,2,1:2) = [ y0 + (2.0*float(center(i,2))-1.0) * rReal, y0 + (2.0*float(center(i,2))+1.0) * rReal]
            cellArray(i,3,1:2) = [ z0 + (2.0*float(center(i,3))-1.0) * rReal, z0 + (2.0*float(center(i,3))+1.0) * rReal]
        enddo
        if ( minval(cellArray(1:cellTotal,1:3,1)) < 0.0 .OR. maxval(cellArray(1:cellTotal,1:3,2)) > rsim(1,2) ) then
            write(*,*) '           INITIALIZATION ERROR             '
            write(*,*) 'ERROR: SYSTEM SIZE TOO SMALL FOR CELL CLUSTER'
        end if

        deallocate( dList)
        deallocate( dTest)
        deallocate(nnList)
        deallocate(nnTest)
    end subroutine itl3DClusterNN


    ! initialize 2D cluster of cells using nearest-neighbor (NN) method
    subroutine itl2DClusterNN( cellArray, rsim)
        implicit none
        real(b8), intent(out) :: cellArray(:,:,:)
        real(b8), intent(in)  :: rsim(:,:)
        real(b8), allocatable :: dList(:), dTest(:)
        integer,  allocatable :: nnList(:,:)
        real(b8) :: a, hCell, nnd, nndmax, x0, y0
        integer :: center( cellTotal, 3), test(2)
        integer :: i, itest, inn, j, k, k1, k2, kmax, n, centerCheck

        allocate(  dList(200))
        allocate(  dTest(200))
        allocate( nnList(200,3))
        cellArray(:,:,:) = 0.0_b8
        ! set cell height
        hCell = 0.10_b8
        ! set center and cellArray for cell 1
        x0 = (rsim(1,2) - rsim(1,1)) / 2.0_b8
        y0 = (rsim(2,2) - rsim(2,1)) / 2.0_b8
        cellArray(1,3,1) = ((rsim(3,2) - rsim(3,1)) / 2.0_b8) - hCell
        cellArray(1,3,2) = ((rsim(3,2) - rsim(3,1)) / 2.0_b8) + hCell
        do i = 1, 2
            cellArray(1,i,1) = (rsim(i,2) - rsim(i,1)) / 2.0_b8 - rReal
            cellArray(1,i,2) = (rsim(i,2) - rsim(i,1)) / 2.0_b8 + rReal
        enddo
        center(:,:) = 0
        ! write(*,*) center(1,1:2)
        if ( cellTotal == 1 ) then
            return
        end if
        n = 2
        kmax = 1
        dList = 0.0_b8
        nnList = 0
        do while ( n <= cellTotal )
            ! check for free NN sites
            ! write(*,*) '  kmax =', kmax
            k = 1
            test = 0
            nnList = 0
            dList  = 0.0_b8
            do k1 = kmax, -kmax, -1
                do k2 = kmax, -kmax, -1
                    do j = 1, n-1
                        centerCheck = 0
                        test = [center(1,1) + k1, center(1,2) + k2]
                        if ( center(j,1) == test(1) .AND. center(j,2) == test(2) ) then
                            centerCheck = 1
                            exit
                        end if
                    enddo
                    if ( centerCheck == 0 ) then
                        nnList(k,:) = test
                        dList(k) = sqrt( float(nnList(k,1)-center(1,1))**2 + float(nnList(k,2)-center(1,2))**2)
                        k = k + 1
                    end if
                enddo
            enddo
            inn = k - 1
            nndmax = float(ceiling(minval( dList(1:inn)) + 0.00001))
            do while ( inn > 0 .AND. n <= cellTotal )
                nnd = minval( dList(1:inn))
                if ( nnd == maxval( dList(1:inn)) .OR. nnd >= nndmax ) then
                    exit
                end if
                ! write(*,*) 'nnd =', nnd, '  inn =', inn
                dTest =  dList
                do i = 1, inn
                    if ( dList(i) == nnd .AND. n <= cellTotal) then
                        center(n,:) = nnList(i,:)
                        n = n + 1
                        ! remove nnList and dList entries that have been added to center
                        dTest(i)   = 0.0
                    end if
                enddo
                j = inn
                do i = inn, 1, -1
                    if ( dTest(i) == 0.0 ) then
                        if ( i /= inn ) then
                             dList(i:inn-1)   =  dList(i+1:inn)
                            nnList(i:inn-1,:) = nnList(i+1:inn,:)
                        end if
                        j = j - 1
                    end if
                enddo
                inn = j
            enddo
            kmax = kmax + 1
        enddo
        itest = cellTotal
        do i = 2, itest
            cellArray(i,1,1:2) = [ x0 + (2.0*float(center(i,1))-1.0) * rReal, x0 + (2.0*float(center(i,1))+1.0) * rReal]
            cellArray(i,2,1:2) = [ y0 + (2.0*float(center(i,2))-1.0) * rReal, y0 + (2.0*float(center(i,2))+1.0) * rReal]
            cellArray(i,3,1:2) = cellArray(1,3,1:2)
        enddo
        if ( minval(cellArray(1:cellTotal,1:3,1)) < 0.0 .OR. maxval(cellArray(1:cellTotal,1:3,2)) > rsim(1,2) ) then
            write(*,*) '           INITIALIZATION ERROR             '
            write(*,*) 'ERROR: SYSTEM SIZE TOO SMALL FOR CELL CLUSTER'
        end if

        deallocate(  dList)
        deallocate(  dTest)
        deallocate( nnList)
    end subroutine itl2DClusterNN


    ! initialize positions of cells to form a cluster
    ! cells are randomly placed along all possible sides of the cluster
    subroutine itl3DRandom( cellTotal, cellArray, rsim)
        implicit none
        integer,  intent(in)  :: cellTotal
        real(b8), intent(out) :: cellArray(:,:,:)
        real(b8), intent(in)  :: rsim(:,:)
        real(b8), allocatable :: center(:,:), dr(:,:), test(:)
        real(b8) :: r
        integer :: cellCheck, i, ii, j, k, n

        allocate( dr(6,3) )
        allocate( test(3) )
        allocate( center(cellTotal,3) )
        center(:,:)      = 0.0_b8
        cellArray(:,:,:) = 0.0_b8
        ! set dr
        dr(:,:) = 0.0_b8
        do i = 1, 3
            dr(i,i) = 2.0 * rReal
        enddo
        do i = 4, 6
            dr(i,i-3) = -2.0 * rReal
        enddo
        ! set center and cellArray for cell 1
        do i = 1, 3
            center(1,i) = (rsim(i,2) - rsim(i,1)) / 2.0_b8
        enddo
        do i = 1, 3
            cellArray(1,i,1) = center(1,i) - rReal
            cellArray(1,i,2) = center(1,i) + rReal
        enddo

        if ( cellTotal == 1 ) then
            return
        elseif ( cellTotal == 2 ) then
            center(2,1) = center(1,1) + (2.0 * rReal)
            center(2,2:3) = center(1,2:3)
            do i = 1, 3
                cellArray(2,i,1) = center(2,i) - rReal
                cellArray(2,i,2) = center(2,i) + rReal
            enddo
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
                        cellArray(i,ii,1) = test(ii) - rReal
                        cellArray(i,ii,2) = test(ii) + rReal
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
        if ( minval(cellArray(1:cellTotal,1:3,1)) < 0.0 .OR. maxval(cellArray(1:cellTotal,1:3,2)) > rsim(1,2) ) then
            write(*,*) '           INITIALIZATION ERROR             '
            write(*,*) 'ERROR: SYSTEM SIZE TOO SMALL FOR CELL CLUSTER'
        end if
        deallocate( dr )
        deallocate( test )
        deallocate( center )
    end subroutine itl3DRandom


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
            dr(i,i) = 2.0 * rReal
        enddo
        do i = 4, 6
            dr(i,i-3) = -2.0 * rReal
        enddo
        ! set cell center from cellArray
        center(:,:)      = 0.0_b8
        do i = 1, cellTotal
            do j = 1, 3
                center(i,j) = cellArray(i,j,1) + (cellArray(i,j,2) - cellArray(i,j,1)) / 2.0_b8
            enddo
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
                if ( dcell >= 2.0*rReal-0.000001 .AND. dcell <= 2.0*rReal+0.000001 ) then
                    check = check + 1
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


    ! calculate center of mass of cell cluster, neccasary for EC polarization
    subroutine clusterCenter( cellArray, clstrCOM)
        implicit none
        real(b8), intent(in)  :: cellArray(:,:,:)
        real(b8), intent(out) :: clstrCOM(3)
        integer :: i, j
        clstrCOM = 0.0_b8
        do i = 1, cellTotal
            do j = 1, 3
                clstrCOM(j) = clstrCOM(j) + cellArray(i,j,1) + (cellArray(i,j,2) - cellArray(i,j,1)) / (2.0_b8)
            enddo
        enddo
        clstrCOM(:) = clstrCOM(:) / float(cellTotal)
    end subroutine clusterCenter


    ! read in cell locations from file and copy them to cellArray
    subroutine readCellLocation( cellArray)
        implicit none
        real(b8), intent(out) :: cellArray(:,:,:)
        integer  :: i, j, k
        open( unit=11, file='cellLocation.dat', status="old", action="read")

        cellArray = 0.0_b8
        do i = 1, cellTotal
            read( 11, *) j, k, cellArray(i,1,2), cellArray(i,1,1)
            read( 11, *) j, k, cellArray(i,2,2), cellArray(i,2,1)
            read( 11, *) j, k, cellArray(i,3,2), cellArray(i,3,1)
            ! write(*,*) 'i', i, cellArray(i,:,:)
        enddo
        close(11)
    end subroutine readCellLocation


    subroutine wrtCellLocation( cellArray)
        implicit none
        real(b8), intent(in) :: cellArray(:,:,:)
        real(b8) :: center
        integer  :: i, j

        do i = 1, cellTotal
            do j = 1, 3
                center = cellArray(i,j,1) + (cellArray(i,j,2) - cellArray(i,j,1)) / (2.0_b8)
                write(210,'(I4)', advance='no') i
                write(210,'(I4)', advance='no') j
                write(210,'(F7.3,F7.3)', advance='no') cellArray(i,j,2), cellArray(i,j,1)
                write(210,'(F7.3)', advance='no') center
                write(210,*) ''
            enddo
        enddo
    end subroutine wrtCellLocation


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

end module
