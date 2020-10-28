!
! Command line compilation via: gfortran -Ofast -fdefault-real-8 -fdefault-integer-8 -ffree-line-length-none -o SU2.x SU2.f90
!
!
! Global variable declarations
!
    module globals
    integer :: particle_num,bond_num
    integer :: max_bond_num,max_iteration
    real :: lx,ly,lz,zero
    integer, allocatable, dimension(:) :: id,cid,nn,local_label
    integer, allocatable, dimension(:,:) :: nl,nb,bond_list
    real, allocatable, dimension(:) :: x,y,z
    end module globals
!
! SU(2) main program (Version 1.0)
!
    program SU2
    use globals
    implicit none

    integer :: i_particle,j_particle,i_bond,i_nn,ii_nn,jj_nn
    integer :: id_val,cid_val,nn_val,cn_num
    integer :: label(0:30,3),i_label,label_population(0:31)
    integer, allocatable, dimension(:) :: nl_val,cnl,cn_hist
    real :: x_val,y_val,z_val,xlo,xhi,ylo,yhi,zlo,zhi
    real :: average_q,rad_t1,rad_t2,radius(2)
    character*300 :: input1,input2,input3,input4,voro_command
    logical :: found_label

    zero=0.0
!
! Maximum number of iterations of modified Voronoi tessellation
!
    max_iteration=5
!
! Maximum number of bonds per atom (n.b. no bounds check)
!
    max_bond_num=45
    allocate (nl_val(max_bond_num),cnl(max_bond_num),cn_hist(max_bond_num))
!
! Radii of atom of type 1 and type 2 (currently set for Lennard-Jones Wahntroem potential)
!
    rad_t1=0.55
    rad_t2=0.45
    radius(:)=(/rad_t1,rad_t2/)
!
! Read in lammps configuration from the terminal. Assumes an orthorombic simulation cell and only two atom types
!
! n.b. the lammps header format might differ from the following assumed format
!
    print *,
    print *, 'Reading in atomic coordinates file'
    do i_particle=1,2
        read (5,*)
    end do
    read (5,*) particle_num
    allocate (id(particle_num),cid(particle_num),x(particle_num),y(particle_num),z(particle_num))
    print *, 'particle number....................',particle_num
    do i_particle=1,2
        read (5,*)
    end do
    read (5,*) xlo,xhi
    read (5,*) ylo,yhi
    read (5,*) zlo,zhi
    lx=xhi-xlo
    ly=yhi-ylo
    lz=zhi-zlo
    print *, 'lx, ly, lz.........................',lx,ly,lz
    do i_particle=1,13
        read (5,*)
    end do
    do i_particle=1,particle_num
        read (5,*) id_val,cid_val,x_val,y_val,z_val
        id(id_val)=id_val
        cid(id_val)=cid_val
        x(id_val)=x_val
        y(id_val)=y_val
        z(id_val)=z_val
    end do
!
! Generate configuration file for voro++
!
    open (unit=10,file='positions.dat')
    do i_particle=1,particle_num
        if (cid(i_particle)==1) then
            write (10,*) id(i_particle),x(i_particle),y(i_particle),z(i_particle),rad_t1
        else
            write (10,*) id(i_particle),x(i_particle),y(i_particle),z(i_particle),rad_t2
        end if
    end do
    close (10)
!
! Generate system command to call voro++ (n.b. system dependent path, assumes configuration has full periodic boundary
! conditions)
!
    write (input1,*) zero
    write (input2,*) lx
    write (input3,*) ly
    write (input4,*) lz
    voro_command='~/Software/voro++-0.4.6/src/voro++ -r -p -c "%i %x %y %z %s %n" '//trim(adjustl(input1))//' '// &
        trim(adjustl(input2))//' '//trim(adjustl(input1))//' '//trim(adjustl(input3))//' '//trim(adjustl(input1))// &
        ' '//trim(adjustl(input4))//' positions.dat'
    print *,
    print *, 'calling voro++'
    call execute_command_line(voro_command)
!
! Read in relevant data from voro++ generated file
!
    allocate (nn(particle_num),nl(particle_num,max_bond_num),nb(particle_num,max_bond_num))
    open (unit=10,file='positions.dat.vol')
    do i_particle=1,particle_num
        read (10,*) id_val,x_val,y_val,z_val,nn_val,(nl_val(i_nn),i_nn=1,nn_val)
        nn(id_val)=nn_val
        nl(id_val,1:nn_val)=nl_val(1:nn_val)
    end do
    close (10)
!
! Perform bond order analysis using standard Voronoi tessellation
!
    print *,
    print *, 'Bond order analysis using standard Voronoi tessellation'
!
! Calculate the number of common neighbours of atoms of a particular bond (prior to removing triplets)
!
    cn_hist(:)=0.0
    bond_num=0
    do i_particle=1,particle_num
        do i_nn=1,nn(i_particle)
            j_particle=nl(i_particle,i_nn)
            cn_num=0
            do ii_nn=1,nn(i_particle)
                do jj_nn=1,nn(j_particle)
                    if (nl(j_particle,jj_nn)==nl(i_particle,ii_nn)) cn_num=cn_num+1
                end do
            end do
            nb(i_particle,i_nn)=cn_num
            if (j_particle<i_particle) cycle
            cn_hist(cn_num)=cn_hist(cn_num)+1
            bond_num=bond_num+1
        end do
    end do
    print *, 'number of bonds .............',bond_num
!
! Calculate average coordination (prior to removing triplets)
!
    average_q=0.0
    do cn_num=1,max_bond_num
        average_q=average_q+float(cn_hist(cn_num)*cn_num)
    end do
    average_q=average_q/sum(cn_hist(:))
    print *, 'average_q..........................',average_q
    print *, 'average coordination (Eueler) .....',12.0/(6.0-average_q)
    print *, 'average coordination...............',float(sum(nn(1:500)))/500.0
!
! Calculate histogram of common neighbours of atoms of a particular bond
!
    print *,
    print *, 'histogram of bond order'
    do cn_num=1,max_bond_num
        if (cn_hist(cn_num)/=0) print *, cn_num,cn_hist(cn_num),100.0*float(cn_hist(cn_num))/float(bond_num)
    end do
!
! Perform modified Voronoi tesselation
!
    print *,
    print *, 'calling modRV'
    call modRV()
!
! Perform bond order analysis using modified Voronoi tessellation
!
    print *,
    print *, 'Bond order analysis using modified Voronoi tessellation'
!
! Calculate the number of common neighbours of atoms of a particular bond with triplets removed
!
    cn_hist(:)=0.0
    bond_num=0
    do i_particle=1,particle_num
        do i_nn=1,nn(i_particle)
            j_particle=nl(i_particle,i_nn)
            cn_num=0
            do ii_nn=1,nn(i_particle)
                do jj_nn=1,nn(j_particle)
                    if (nl(j_particle,jj_nn)==nl(i_particle,ii_nn)) cn_num=cn_num+1
                end do
            end do
            nb(i_particle,i_nn)=cn_num
            if (j_particle<i_particle) cycle
            cn_hist(cn_num)=cn_hist(cn_num)+1
            bond_num=bond_num+1
        end do
    end do
    print *, 'number of bonds .............',bond_num
!
! Calculate average coordination
!
    average_q=0.0
    do cn_num=1,max_bond_num
        average_q=average_q+float(cn_hist(cn_num)*cn_num)
    end do
    average_q=average_q/sum(cn_hist(:))
    print *, 'average_q..........................',average_q
    print *, 'average coordination (Eueler) .....',12.0/(6.0-average_q)
    print *, 'average coordination...............', float(sum(nn(1:500)))/500.0
!
! Calculate histogram of common neighbours of atoms of a particular bond
!
    print *,
    print *, 'histogram of bond order'
    do cn_num=1,max_bond_num
        if (cn_hist(cn_num)/=0) print *, cn_num,cn_hist(cn_num),100.0*float(cn_hist(cn_num))/float(bond_num)
    end do
!
! Enumerate bond list
!
    allocate (bond_list(bond_num,3))
    i_bond=0
    do i_particle=1,particle_num
        do i_nn=1,nn(i_particle)
            j_particle=nl(i_particle,i_nn)
            if (j_particle<i_particle) cycle
            i_bond=i_bond+1
            bond_list(i_bond,:)=(/nb(i_particle,i_nn),id(i_particle),id(j_particle)/)
        end do
    end do
!
! Define a list of local environments
!
    label(0,:)=(/0,12,0/)
    label(1,:)=(/0,12,2/)
    label(2,:)=(/0,12,3/)
    label(3,:)=(/0,12,4/)
    label(4,:)=(/0,12,5/)
    label(5,:)=(/0,12,6/)
    label(6,:)=(/1,10,2/)
    label(7,:)=(/1,10,3/)
    label(8,:)=(/1,10,4/)
    label(9,:)=(/1,10,5/)
    label(10,:)=(/1,10,6/)
    label(11,:)=(/1,10,7/)
    label(12,:)=(/2,8,0/)
    label(13,:)=(/2,8,1/)
    label(14,:)=(/2,8,2/)
    label(15,:)=(/2,8,3/)
    label(16,:)=(/2,8,4/)
    label(17,:)=(/2,8,5/)
    label(18,:)=(/2,8,6/)
    label(19,:)=(/2,8,7/)
    label(20,:)=(/2,8,8/)
    label(21,:)=(/3,6,0/)
    label(22,:)=(/3,6,1/)
    label(23,:)=(/3,6,2/)
    label(24,:)=(/3,6,3/)
    label(25,:)=(/3,6,4/)
    label(26,:)=(/3,6,5/)
    label(27,:)=(/3,6,6/)
    label(28,:)=(/3,6,7/)
    label(29,:)=(/3,6,8/)
    label(30,:)=(/3,6,9/)
!
! Determine local SU(2) topology of those atomic environments which satisfy Eulers theorem
!
    allocate(local_label(particle_num))
    label_population(:)=0
    do i_particle=1,particle_num
        cn_hist(:)=0
        do i_nn=1,nn(i_particle)
            cn_hist(nb(i_particle,i_nn))=cn_hist(nb(i_particle,i_nn))+1
        end do
        if (cn_hist(4)+cn_hist(5)+cn_hist(6)==nn(i_particle).and.12-cn_hist(4)+cn_hist(6)==nn(i_particle)) then
            found_label=.false.
            do i_label=0,30
                if (cn_hist(4)==label(i_label,1).and.cn_hist(5)==label(i_label,2).and.cn_hist(6)==label(i_label,3)) then
                    local_label(i_particle)=i_label
                    label_population(i_label)=label_population(i_label)+1
                    found_label=.true.
                    exit
                end if
            end do
            if (.not.found_label) then
                local_label(i_particle)=31
                label_population(31)=label_population(31)+1
            end if
        else
            local_label(i_particle)=31
            label_population(31)=label_population(31)+1
        end if
    end do
    close (10)
!
! Output SU(2) histogram
!
    print *,
    print *, 'SU2 histogram'
    do i_label=0,30
        print '(a2,i2,a1,i2,a1,i2,a1,1x,e16.10)', ' (',label(i_label,1),',',label(i_label,2),',',label(i_label,3),')',float(label_population(i_label))/float(particle_num)
    end do
    print *, ' other',float(label_population(31))/float(particle_num)
!
! Output local environments (via molecular id) and bonds to a lammps file
!
    open (unit=10,file='su2.lammps')
    write (10,'(a)') 'LAMMPS data file created by make_bond_polyhedra.f90'
    write (10,'(a)')
    write (10,'(i6,a)') particle_num,' atoms'
    write (10,'(i6,a)') bond_num,' bonds'
    write (10,'(a)') '2 atom types'
    write (10,'(a)') '8 bond types'
    write (10,'(a)')
    write (10,'(f10.7,1x,f10.7,a)') xlo,xhi,' xlo xhi'
    write (10,'(f10.7,1x,f10.7,a)') ylo,yhi,' ylo yhi'
    write (10,'(f10.7,1x,f10.7,a)') zlo,zhi,' zlo zhi'
    write (10,'(a)')
    write (10,'(a)') 'Masses'
    write (10,'(a)')
    write (10,'(a)') '1 2'
    write (10,'(a)') '2 1'
    write (10,'(a)')
    write (10,'(a)') 'Atoms'
    write (10,'(a)')
    do i_particle=1,particle_num
        write (10,*) id(i_particle),local_label(i_particle),cid(i_particle),x(i_particle),y(i_particle),z(i_particle)
    end do
    write (10,'(a)')
    write (10,'(a)') 'Bonds'
    write (10,'(a)')
    do i_bond=1,bond_num
           write (10,*) i_bond,bond_list(i_bond,:)
    end do
    close (10)

    end
!
!*****************************************************************************************
! Exlude voro++ bonds which pass through a triplet of neighbouring atoms
!
! This is done by identifying all possible neighbour triplets and using
! the tetrahedral construction to determine if a bond intersects the plane
! defined by the triplet positions.
!*****************************************************************************************
!
    subroutine modRV()
    use globals
    implicit none

    interface
        function tetrahedral_volume(id1,id2,id3,id4)
        integer :: id1,id2,id3,id4
        real :: tetrahedral_volume
        end function
    end interface

    logical :: ij_member,ik_member,jk_member
    integer :: bonds_to_exclude_num,i_particle,j_particle,iteration
    integer :: i_nn,ii_nn,iii_nn,j_nn,jj_nn,jjj_nn
    integer :: cnl_num,cnl(45),new_nl(45)
    integer :: i_triplet,j_triplet,k_triplet
    real :: tv_abde,tv_bcde,tv_cade

    iteration=0
    do 
        bonds_to_exclude_num=0
        do i_particle=1,particle_num
            do i_nn=1,nn(i_particle)
                j_particle=nl(i_particle,i_nn)
                if (j_particle<i_particle) cycle
                !
                ! Construct list of neighbours common to nearest neighbour particles 'i_particle' and 'j_particle'
                !
                cnl_num=0
                do ii_nn=1,nn(i_particle)
                    do jj_nn=1,nn(j_particle)
                        if (nl(j_particle,jj_nn)==nl(i_particle,ii_nn)) then
                            cnl_num=cnl_num+1
                            cnl(cnl_num)=nl(j_particle,jj_nn)
                        end if
                    end do
                end do
                !
                ! Cycle through all triplets of these common neighbours
                !
triplet:        do i_triplet=1,cnl_num
                    do j_triplet=i_triplet+1,cnl_num
                        do k_triplet=j_triplet+1,cnl_num
                            !
                            ! Check if a triplet of these common neighbours (to particles 'i_particle' and 'j_particle')
                            ! are common neighbours to each other
                            ! 
                            ij_member=.false.
                            do ii_nn=1,nn(cnl(i_triplet))
                                if (nl(cnl(i_triplet),ii_nn)==cnl(j_triplet)) ij_member=.true.
                            end do
                            ik_member=.false.
                            do ii_nn=1,nn(cnl(i_triplet))
                                if (nl(cnl(i_triplet),ii_nn)==cnl(k_triplet)) ik_member=.true.
                            end do
                            jk_member=.false.
                            do jj_nn=1,nn(cnl(j_triplet))
                                if (nl(cnl(j_triplet),jj_nn)==cnl(k_triplet)) jk_member=.true.
                            end do
                            !
                            ! If the triple are common neighbours to themselves, then check if the plane that they define
                            ! is intersected by the bond between particle's 'i_particle' and 'j_particle'
                            !
                            if (ij_member.and.ik_member.and.jk_member) then
                                tv_abde=tetrahedral_volume(cnl(i_triplet),cnl(j_triplet),i_particle,j_particle)
                                tv_bcde=tetrahedral_volume(cnl(j_triplet),cnl(k_triplet),i_particle,j_particle)
                                tv_cade=tetrahedral_volume(cnl(k_triplet),cnl(i_triplet),i_particle,j_particle)
                                !
                                ! If the intersection occurs omit the nearest-neighbour bond between particles 'i_particle' and 'j_particle'
                                !
                                if ((tv_abde<0.0.and.tv_bcde<0.0.and.tv_cade<0.0).or. &
                                    (tv_abde>0.0.and.tv_bcde>0.0.and.tv_cade>0.0)) then
                                    bonds_to_exclude_num=bonds_to_exclude_num+1
                                    ! modify the neighbour lists of the atoms involved in the central bond
                                    iii_nn=0
                                    do ii_nn=1,nn(i_particle)
                                        if (nl(i_particle,ii_nn)/=j_particle) then
                                            iii_nn=iii_nn+1
                                            new_nl(iii_nn)=nl(i_particle,ii_nn)
                                        end if
                                    end do
                                    nl(i_particle,1:iii_nn)=new_nl(1:iii_nn)
                                    nn(i_particle)=iii_nn
                                    jjj_nn=0
                                    do jj_nn=1,nn(j_particle)
                                        if (nl(j_particle,jj_nn)/=i_particle) then
                                            jjj_nn=jjj_nn+1
                                            new_nl(jjj_nn)=nl(j_particle,jj_nn)
                                        end if
                                    end do
                                    nl(j_particle,1:jjj_nn)=new_nl(1:jjj_nn)
                                    nn(j_particle)=jjj_nn
                                    exit triplet
                                end if
                            end if
                        end do
                    end do
                end do triplet
            end do
        end do
        iteration=iteration+1
        print *, 'iteration..........................',iteration
        print *, 'number of excluded bonds...........',bonds_to_exclude_num
        if (bonds_to_exclude_num==0.or.iteration>max_iteration) exit
    end do

    end subroutine
!
! Function to calculate volume of tetrahedron defined by four points (delivered as vectors)
!
    function tetrahedral_volume(id1,id2,id3,id4)
    use globals
    implicit none

    integer :: id1,id2,id3,id4
    real :: dx,dy,dz
    real :: veca(3),vecb(3),vecc(3),vecd(3),tetrahedral_volume
    real :: vecba(3),vecca(3),vecda(3),cp_vecca_vecda(3)

    veca(:)=0.0

    dx=x(id2)-x(id1)
    if (dx> lx/2.0) dx=dx-lx
    if (dx<-lx/2.0) dx=dx+lx
    dy=y(id2)-y(id1)
    if (dy> ly/2.0) dy=dy-ly
    if (dy<-ly/2.0) dy=dy+ly
    dz=z(id2)-z(id1)
    if (dz> lz/2.0) dz=dz-lz
    if (dz<-lz/2.0) dz=dz+lz
    vecb(:)=(/dx,dy,dz/)

    dx=x(id3)-x(id1)
    if (dx> lx/2.0) dx=dx-lx
    if (dx<-lx/2.0) dx=dx+lx
    dy=y(id3)-y(id1)
    if (dy> ly/2.0) dy=dy-ly
    if (dy<-ly/2.0) dy=dy+ly
    dz=z(id3)-z(id1)
    if (dz> lz/2.0) dz=dz-lz
    if (dz<-lz/2.0) dz=dz+lz
    vecc(:)=(/dx,dy,dz/)

    dx=x(id3)-x(id1)
    if (dx> lx/2.0) dx=dx-lx
    if (dx<-lx/2.0) dx=dx+lx
    dy=y(id3)-y(id1)
    if (dy> ly/2.0) dy=dy-ly
    if (dy<-ly/2.0) dy=dy+ly
    dz=z(id3)-z(id1)
    if (dz> lz/2.0) dz=dz-lz
    if (dz<-lz/2.0) dz=dz+lz
    vecc(:)=(/dx,dy,dz/)

    dx=x(id4)-x(id1)
    if (dx> lx/2.0) dx=dx-lx
    if (dx<-lx/2.0) dx=dx+lx
    dy=y(id4)-y(id1)
    if (dy> ly/2.0) dy=dy-ly
    if (dy<-ly/2.0) dy=dy+ly
    dz=z(id4)-z(id1)
    if (dz> lz/2.0) dz=dz-lz
    if (dz<-lz/2.0) dz=dz+lz
    vecd(:)=(/dx,dy,dz/)

    vecba(:)=vecb(:)-veca(:)
    vecca(:)=vecc(:)-veca(:)
    vecda(:)=vecd(:)-veca(:)

    cp_vecca_vecda(1)=vecca(2)*vecda(3)-vecca(3)*vecda(2)
    cp_vecca_vecda(2)=-(vecca(1)*vecda(3)-vecca(3)*vecda(1))
    cp_vecca_vecda(3)=vecca(1)*vecda(2)-vecca(2)*vecda(1)

    tetrahedral_volume=dot_product(vecba(:),cp_vecca_vecda(:))/6.0

    end function tetrahedral_volume

