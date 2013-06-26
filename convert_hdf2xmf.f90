program hdf2xml
  implicit none
  integer, parameter :: pr = 8
  integer :: nx, ny, nz, nscal,nvec,ivec,ndt,idt, io_error=0,iscal
  real (kind=pr) :: time, xl, yl, zl
  character (len=5) :: basefilename
  character prefixes_scalar(100)*10 ! maximal 100 prefixes à 10 chars length
  character prefixes_vector(100)*10 ! maximal 100 prefixes à 10 chars length
  character timesteps(10000)*10 ! maximal 10000 time steps à 10 chars length
  LOGICAL :: file_exists, scalars_exist, vectors_exist
  
  write (*,*) "-------------------------"
  write (*,*) " Fortran XMF generator"
  write (*,*) "-------------------------"
  
  !---------------------------------------------------------------------
  ! read in LIST OF scalar prefixes.  note it is mandatory to follow
  ! file naming convention.  for a file mask_00010.h5 "mask" is in
  ! this file
  ! ---------------------------------------------------------------------
  INQUIRE(FILE='prefixes_scalar.in', EXIST=scalars_exist) 
  if(scalars_exist) then
     io_error=0
     open (15, file='prefixes_scalar.in', action='read', status='old' )
     nscal = 1
     write (*,'(A)',advance='no') "scalars = "
     do while (io_error==0)
        read (15,'(A)',iostat=io_error) prefixes_scalar(nscal)
        if (io_error ==0) then ! pay a little attention here, else we go one too far
           write (*,'(A)',advance='no') prefixes_scalar(nscal)
           nscal = nscal+1
        endif
     enddo
     nscal = nscal-1
     write(*,*) " "
     close (15)
  endif
  ! now we know how many scalars we have and what their names are

  !---------------------------------------------------------------------
  ! read in LIST OF vector prefixes.  note vectors have 3 components
  ! and one basename. this file contains the basename only i.e. you
  ! save usx_00010.h5 usy_00010.h5 usz_00010.h5 the basename is "us"
  ! note it is mandatory to follow file naming convention.
  ! ---------------------------------------------------------------------

  INQUIRE(FILE='prefixes_vector.in', EXIST=vectors_exist)
  if(vectors_exist) then
     io_error=0
     open (15, file='prefixes_vector.in', action='read', status='old' )
     nvec = 1
     write (*,'(A)',advance='no') "vectors = "
     do while (io_error==0)
        read (15,'(A)',iostat=io_error) prefixes_vector(nvec)
        if (io_error ==0) then ! pay a little attention here, else we go one too far
           write (*,'(A)',advance='no') prefixes_vector(nvec)
           nvec = nvec+1
        endif
     enddo
     nvec = nvec-1
     write(*,*) " "
     close (15)
  endif
  ! now we know how many vectors we have and what their names are

  if((vectors_exist .eqv. .false.) .and.  (scalars_exist .eqv. .false.)) then
     write(*,*) "Neither scalars nor vectors found; unable to proceed."
     stop
  endif

  !--------------------------------------------------
  ! Read in LIST OF different time steps, generated externally
  ! file looks like
  ! 00010
  ! 00015
  ! 00020 etc
  !--------------------------------------------------
  
  INQUIRE(FILE='timesteps.in', EXIST=file_exists) 
  if(file_exists) then
     io_error=0
     open (15, file='timesteps.in', action='read', status='old' )
     ndt = 1
     write (*,'(A)',advance='no') "time steps= "
     do while (io_error==0)
        read (15,'(A)',iostat=io_error) timesteps(ndt)
        if (io_error ==0) then ! pay a little attention here, else we
                               ! go one too far
           write (*,'(A)',advance='no') timesteps(ndt)
           ndt = ndt+1
        endif
     enddo
     ndt = ndt-1
     write(*,*) " "
     close (15)
  else
     write(*,*) "timesteps.in does not exist; unable to proceed."
     stop
  endif
  ! Now we know how many vectors we have and what their names are

  !--------------------------------------
  ! Create the main *.xmf file
  !--------------------------------------
  open (14, file='ALL.xmf', status='replace')

  ! Read the resolution so that we can call BeginFile.
  if(scalars_exist) then
     call Fetch_attributes (trim(adjustl(prefixes_scalar(1)))//"_"//trim(adjustl(timesteps(1)))//".h5",&
          trim(adjustl(prefixes_scalar(1))), nx,ny,nz,xl,yl,zl,time)
  endif  

  if(vectors_exist) then
     call Fetch_attributes(trim(adjustl(prefixes_vector(1)))//"x_"//trim(adjustl(timesteps(1)))//".h5",&
          trim(adjustl(prefixes_vector(1)))//"x", nx,ny,nz,xl,yl,zl,time)
  endif

  ! File header, requires resolution
  call BeginFile(nx,ny,nz)

  ! begin loops over time steps
  do idt = 1, ndt  
     ! Get the time from any of these files
     if(scalars_exist) then
        call Fetch_attributes (trim(adjustl(prefixes_scalar(1)))//"_"//trim(adjustl(timesteps(idt)))//".h5",&
             trim(adjustl(prefixes_scalar(1))), nx,ny,nz,xl,yl,zl,time)  
     endif

     if(vectors_exist) then
        ! FIXME: this is probably wrong
        call Fetch_attributes (trim(adjustl(prefixes_vector(1)))//"x_"//trim(adjustl(timesteps(idt)))//".h5",&
             trim(adjustl(prefixes_vector(1)))//"x", nx,ny,nz,xl,yl,zl,time)  
     endif

     ! Time step header
     call BeginTimeStep(nx,ny,nz,xl,yl,zl,time)  

     ! Scalars
     if(scalars_exist) then
        do iscal = 1, nscal
           call WriteScalar(trim(adjustl(timesteps(idt))),&
                trim(adjustl(prefixes_scalar(iscal))),&
                nx,ny,nz,xl,yl ,zl,time)
        end do
     endif

     ! Vectors
     if(vectors_exist) then
        do ivec = 1, nvec  
           call WriteVector(trim(adjustl(timesteps(idt))),&
                trim(adjustl(prefixes_vector(ivec))),&
                nx,ny,nz,xl,yl,zl,time)
        enddo
     endif

     ! Time step footer
     call EndTimeStep()  
  enddo

  ! file footer
  call EndFile()  

  close (14)
end program hdf2xml


subroutine BeginFile(nx,ny,nz)
  implicit none
  integer, parameter :: pr = 8
  integer, intent (in) :: nx, ny, nz
  character (len=128) :: nxyz_str 

  ! note odd indices. paraview requirement
  write (nxyz_str,'(3(i4.4,1x))') nz, ny, nx 
  write (*,'("Resolution: ",3(i4,1x) )') nx, ny, nz

  write (14,'(A)') '<?xml version="1.0" ?>'
  write (14,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" ['
  write (14,'(A)') '<!ENTITY nxnynz "'//trim(adjustl(nxyz_str))//'">'
  write (14,'(A)') ']>'
  write (14,'(A)') '<Xdmf Version="2.0">  '
  write (14,'(A)') '<Domain>  '
  write (14,'(A)') '<Grid Name="Box" GridType="Collection" CollectionType="Temporal">'
end subroutine BeginFile


subroutine EndFile()
  implicit none
  write (14,'(A)') '</Grid>'
  write (14,'(A)') '</Domain>'
  write (14,'(A)') '</Xdmf>  '
end subroutine EndFile


subroutine BeginTimeStep ( nx, ny, nz, xl, yl ,zl, time )
  implicit none
  integer, parameter :: pr = 8
  integer, intent (in) :: nx, ny, nz
  real (kind=pr), intent(in) :: xl,yl,zl, time
  character (len=128) :: time_string

  write (time_string,'(es15.8)') time

  write (14,'(A)') '<!-- beginning time step -->    '
  write (14,'(A)') '<Grid Name="FLUSI_cartesian_grid" GridType="Uniform">'
  write (14,'(A)') '    <Time Value="'//trim(adjustl(time_string))//'" />'
  write (14,'(A)') '    <Topology TopologyType="3DCoRectMesh" Dimensions="&nxnynz;" />'
  write (14,'(A)') ' '
  write (14,'(A)') '    <Geometry GeometryType="Origin_DxDyDz">'
  write (14,'(A)') '    <DataItem Dimensions="3" NumberType="Float" Format="XML">'
  write (14,'(A)') '    0 0 0'
  write (14,'(A)') '    </DataItem>'
  write (14,'(A)') '    <DataItem Dimensions="3" NumberType="Float" Format="XML">'
  write (14,'(4x,3(es15.8,1x))') xl / dble (nx), yl/dble(ny), zl/dble(nz)
  write (14,'(A)') '    </DataItem>'
  write (14,'(A)') '    </Geometry>'  
end subroutine BeginTimeStep


subroutine EndTimeStep ( )
  implicit none
  write (14,*) '</Grid>'
end subroutine EndTimeStep


subroutine WriteVector ( basefilename, prefix, nx, ny, nz, xl, yl ,zl, time )
  implicit none
  integer, parameter :: pr = 8
  integer, intent (in) :: nx, ny, nz
  real (kind=pr), intent(in) :: xl,yl,zl, time
  character (len=*), intent (in) :: basefilename, prefix

  !   write (*,*) "vector: ", prefix, basefilename

  write (14,'(A)') '    <!--Vector-->'
  write (14,'(A)') '    <Attribute Name="'//prefix//'" AttributeType="Vector" Center="Node">'
  write (14,'(A)') '    <DataItem ItemType="Function" Function="JOIN($0, $1, $2)" Dimensions="&nxnynz; 3" NumberType="Float">'
  write (14,'(A)') '        <DataItem Dimensions="&nxnynz;" NumberType="Float" Format="HDF">'
  write (14,'(A)') '        '//prefix//'x_'//basefilename//'.h5:/'//prefix//'x'
  write (14,'(A)') '        </DataItem>'
  write (14,'(A)') '        '
  write (14,'(A)') '        <DataItem Dimensions="&nxnynz;" NumberType="Float" Format="HDF">'
  write (14,'(A)') '        '//prefix//'y_'//basefilename//'.h5:/'//prefix//'y'
  write (14,'(A)') '        </DataItem>'
  write (14,'(A)') '        '
  write (14,'(A)') '        <DataItem Dimensions="&nxnynz;" NumberType="Float" Format="HDF">'
  write (14,'(A)') '        '//prefix//'z_'//basefilename//'.h5:/'//prefix//'z'
  write (14,'(A)') '        </DataItem>     '
  write (14,'(A)') '    </DataItem>'
  write (14,'(A)') '    </Attribute>    '
end subroutine WriteVector


subroutine WriteScalar ( basefilename, prefix, nx, ny, nz, xl, yl ,zl, time )
  implicit none
  integer, parameter :: pr = 8
  integer, intent (in) :: nx, ny, nz
  real (kind=pr), intent(in) :: xl,yl,zl, time
  character (len=*), intent (in) :: basefilename, prefix

  !   write (*,*) "writing scalar ", basefilename, " ", prefix

  write (14,'(A)') '    <!--Scalar-->'
  write (14,'(A)') '    <Attribute Name="'//prefix//'" AttributeType="Scalar" Center="Node">'
  write (14,'(A)') '    <DataItem Dimensions="&nxnynz;" NumberType="Float" Format="HDF">'
  write (14,'(A)') '    '//prefix//'_'//basefilename//'.h5:/'//prefix
  write (14,'(A)') '    </DataItem>'
  write (14,'(A)') '    </Attribute>'  
end subroutine WriteScalar


subroutine Fetch_attributes( filename, dsetname,  nx, ny, nz, xl, yl ,zl, time )
  !----------------------------------------------------
  ! this routine fetches the resolution, the domain size and the time form a *.h5 file
  !----------------------------------------------------
  ! filename: a *.h5 file to read from.
  ! dsetname: a dataset inside the file. in our system, we only have one per file
  !           and this matches the prefix: mask_00010.h5  --> dsetname = "mask"
  ! note:
  !           the file must contain the dataset
  !           but especially the attributes "nxyz", "time", "domain_size"
  !----------------------------------------------------
  use hdf5
  implicit none
  integer, parameter :: pr = 8
  integer, intent (out) :: nx, ny, nz
  real (kind=pr), intent(out) :: xl,yl,zl, time

  character(len=*) :: filename  ! file name
  character(len=*) :: dsetname  ! dataset name
  character(len=4) :: aname     ! attribute name
  character(len=11) :: aname2

  integer(hid_t) :: file_id       ! file identifier
  integer(hid_t) :: dset_id       ! dataset identifier
  integer(hid_t) :: attr_id       ! attribute identifier
  integer(hid_t) :: aspace_id     ! attribute dataspace identifier
  integer(hid_t) :: atype_id      ! attribute dataspace identifier
  integer(hsize_t), dimension(1) :: adims = (/1/) ! attribute dimension
  integer     ::   arank = 1                      ! attribure rank
  integer(size_t) :: attrlen    ! length of the attribute string

  real (kind=pr) ::  attr_data  ! attribute data
  real (kind=pr), dimension (1:3) :: attr_data2
  integer, dimension (1:3) :: attr_data3

  integer     ::   error ! error flag
  integer(hsize_t), dimension(1) :: data_dims

  !   write (*,*) "fetching from ", filename," -> ", dsetname

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)

  ! Open an existing file.
  CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)
  ! Open an existing dataset.
  CALL h5dopen_f(file_id, dsetname, dset_id, error)


!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! open attribute (time)
!!!!!!!!!!!!!!!!!!!!!!!!!!
  aname = "time"
  CALL h5aopen_f(dset_id, aname, attr_id, error)

  ! Get dataspace and read
  CALL h5aget_space_f(attr_id, aspace_id, error)
  data_dims(1) = 1
  CALL h5aread_f( attr_id, H5T_NATIVE_DOUBLE, attr_data, data_dims, error)

  time = attr_data
  CALL h5aclose_f(attr_id, error) ! Close the attribute.
  CALL h5sclose_f(aspace_id, error) ! Terminate access to the data space.

!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! open attribute (domain_length)
!!!!!!!!!!!!!!!!!!!!!!!!!!
  aname2 = "domain_size"
  CALL h5aopen_f(dset_id, aname2, attr_id, error)

  ! Get dataspace and read
  CALL h5aget_space_f(attr_id, aspace_id, error)
  data_dims(1) = 3
  CALL h5aread_f( attr_id, H5T_NATIVE_DOUBLE, attr_data2, data_dims, error)

  xl = attr_data2(1)
  yl = attr_data2(2)
  zl = attr_data2(3)

  CALL h5aclose_f(attr_id, error) ! Close the attribute.
  CALL h5sclose_f(aspace_id, error) ! Terminate access to the data space.

!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! open attribute (sizes)
!!!!!!!!!!!!!!!!!!!!!!!!!!
  aname = "nxyz"
  CALL h5aopen_f(dset_id, aname, attr_id, error)

  ! Get dataspace and read
  CALL h5aget_space_f(attr_id, aspace_id, error)
  data_dims(1) = 3
  CALL h5aread_f( attr_id, H5T_NATIVE_INTEGER, attr_data3, data_dims, error)

  nx = attr_data3(1)
  ny = attr_data3(2)
  nz = attr_data3(3)

  CALL h5aclose_f(attr_id, error) ! Close the attribute.
  CALL h5sclose_f(aspace_id, error) ! Terminate access to the data space.

  CALL h5dclose_f(dset_id, error) ! End access to the dataset and release resources used by it.
  CALL h5fclose_f(file_id, error) ! Close the file.
  CALL h5close_f(error)  ! Close FORTRAN interface.

  !   write (*,'("time=",es12.4," domain=",3(es12.4,1x),2x,3(i3,1x))') time, xl,yl,zl, nx, ny, nz
end subroutine Fetch_attributes
