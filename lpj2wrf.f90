program postprocess

! program to postprocess LPJ-LMfire output for asynchronous coupling experiments (Raible/Velazquez UniBE).
! August 2017 - March 2018

! ifort -xHost -o postprocess postprocess.f90 -I/soge-home/users/cenv0668/local/netcdf/include -L/soge-home/users/cenv0668/local/netcdf/lib -lnetcdff -lnetcdf

use iso_fortran_env, only : int8,int16,int32,real32,real64
use netcdf

implicit none

integer, parameter :: sp = real32      ! 4 byte real (default)
integer, parameter :: dp = real64      ! 8 byte real
integer, parameter :: i1 = int8        ! 1 byte integer (signed, [-128,127])
integer, parameter :: i2 = int16       ! 2 byte integer
integer, parameter :: i4 = int32       ! 4 byte integer (default)

integer, parameter :: nc = 20  !number of cover classes

integer :: status
integer :: ifid
integer :: ofid
integer :: ffid
integer :: dimid
integer :: varid

character(100) :: jobfile
character(100) :: lpjfile
character(100) :: maskfile
character(100) :: climfile
character(100) :: outfile 

real(sp), allocatable, dimension(:,:,:) :: cover

real(sp), allocatable, dimension(:,:) :: icef
real(sp), allocatable, dimension(:,:) :: h2of
real(sp), allocatable, dimension(:,:) :: landf

real(sp), allocatable, dimension(:,:,:) :: covercatf
real(sp), allocatable, dimension(:,:,:) :: mLAI
real(sp), allocatable, dimension(:,:,:) :: fpar
real(sp), allocatable, dimension(:,:,:) :: almean

real(sp), allocatable, dimension(:,:,:,:) :: mLAIpft

real(sp), allocatable, dimension(:) :: maxLAI
real(sp), allocatable, dimension(:) :: phen

integer(i1), allocatable, dimension(:,:) :: covercat
real(sp), allocatable, dimension(:,:) :: snowmax

integer(i2), allocatable, dimension(:,:) :: ival
real(sp), allocatable, dimension(:,:) :: tmean

integer :: xlen
integer :: ylen
integer :: npft

integer :: i
integer :: x
integer :: y
integer :: m

logical, allocatable, dimension(:) :: present

real(sp) :: totalcover
real(sp) :: treecover
real(sp) :: treecover_decid
real(sp) :: treecover_evrgn
real(sp) :: grasscover
real(sp) :: vegcover

integer :: domtree
integer, dimension(1) :: pos

real(sp), dimension(9) :: PFTalbedo  !albedo at full leaf out
real(sp), dimension(9) :: alstem     !albedo of stems

namelist /albedo/ pftalbedo
namelist /infiles/ lpjfile,maskfile,climfile

real(sp), parameter :: alice   = 0.80  !from CLM
real(sp), parameter :: alsoil  = 0.25  !rough average over range of soil colors and hydration states
real(sp), parameter :: alwater = 0.08  !very rough average over range of solar zenith angles

real(sp), parameter :: alstem_woody = 0.16  !values from CLM
real(sp), parameter :: alstem_grass = 0.31

real(sp) :: alveg
real(sp) :: alunveg

logical, dimension(9) :: decid

!----------------------------------------

! variables required:
! 1. land cover fraction (category)
! 2. dominant land cover type (category)
! 3. FPAR (monthly)
! 4. LAI (monthly)
! 5. albedo (monthly)
! 6. soil temperature (annual mean)
! 7. snow albedo (annual max)

! notes
! in LPJ-LMfire, mean annual soil temperature is effectively air temperature, so we will copy that from the climate input
! albedo depends on snow cover, so need to archive this variable

! read in essential data

open(10,file='PFTalbedo.namelist')

read(10,nml=albedo)

close(10)

alstem(1:7) = alstem_woody
alstem(8:9) = alstem_grass

decid = .true.
decid(1) = .false.
decid(3) = .false.
decid(4) = .false.
decid(6) = .false.

!---

call getarg(1,jobfile)

open(10,file=jobfile)

read(10,nml=infiles)

close(10)

status = nf90_open(lpjfile,nf90_nowrite,ifid)
if (status /= nf90_noerr) call handle_err(status)

!------------

write(0,*)'reading'

status = nf90_inq_dimid(ifid,'lon',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ifid,dimid,len=xlen)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ifid,'lat',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ifid,dimid,len=ylen)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ifid,'pft',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ifid,dimid,len=npft)
if (status /= nf90_noerr) call handle_err(status)

allocate(cover(xlen,ylen,npft))
allocate(mLAI(xlen,ylen,12))
allocate(mLAIpft(xlen,ylen,npft,12))

allocate(fpar(xlen,ylen,12))
allocate(maxlai(npft))
allocate(phen(npft))
allocate(almean(xlen,ylen,12))
allocate(snowmax(xlen,ylen))

status = nf90_inq_varid(ifid,'cover',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ifid,varid,cover)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ifid,'mLAI',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ifid,varid,mLAIpft)
if (status /= nf90_noerr) call handle_err(status)

!----
!get mean temperature

allocate(ival(xlen,ylen))
allocate(tmean(xlen,ylen))

tmean = -9999.

status = nf90_open(climfile,nf90_nowrite,ffid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ffid,'tmp',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ffid,varid,ival)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(ffid)
if (status /= nf90_noerr) call handle_err(status)

where (ival /= -32768) tmean = real(ival) * 0.1

deallocate(ival)

!----
!get water and ice masks

allocate(icef(xlen,ylen))
allocate(h2of(xlen,ylen))
allocate(landf(xlen,ylen))

status = nf90_open(maskfile,nf90_nowrite,ffid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ffid,'landf',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ffid,varid,landf)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ffid,'icef',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ffid,varid,icef)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ffid,'waterf',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ffid,varid,h2of)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(ffid)
if (status /= nf90_noerr) call handle_err(status)

!----
! 1. calculate land cover categories

write(0,*)'calculating'

allocate(present(npft))

allocate(covercatf(xlen,ylen,nc))

covercatf = 0.
snowmax = 1.

where (cover < 0.) cover = 0.

mLAI = -9999.

do y = 1,ylen
  do x = 1,xlen
        
    present = .false.
  
    where (cover(x,y,:) > 0.) present = .true.
  
    treecover = sum(cover(x,y,1:7))
    treecover_decid = sum(cover(x,y,1:7),mask=decid(1:7))
    treecover_evrgn = sum(cover(x,y,1:7),mask=.not.decid(1:7))
    grasscover = sum(cover(x,y,8:9))
    totalcover = sum(cover(x,y,:))
    vegcover = max(totalcover - icef(x,y) - h2of(x,y),0.)
    
    pos = maxloc(cover(x,y,1:7),mask=cover(x,y,1:7) > 0.5)
    
    domtree = pos(1)
    
    if (treecover >= 0.75) then  !forests
    
      select case(domtree)

      case(3,6)  !cat 1: evergreen needleleaf forest
        
                covercatf(x,y,1) = vegcover
              
            case(1,4)  !cat 2: evergreen broadleaf forest  

                covercatf(x,y,2) = vegcover
                
            case(7)     !cat 3: deciduous needleleaf forest

                covercatf(x,y,3) = 0.5 * vegcover !assign half to needleleaf and half to broadleaf
                covercatf(x,y,4) = 0.5 * vegcover !assign half to needleleaf and half to broadleaf

            case (2,5) !cat 4: deciduous broadleaf forest

                covercatf(x,y,4) = vegcover !cover(x,y,2) + cover(x,y,5)
            
            case(0)        !cat 5: mixed forests (no dominant tree pft)

                covercatf(x,y,5) = vegcover

            end select
            
        else if (treecover > 0.25 .and. grasscover < 0.25) then  !shrublands
        
          if (treecover >= 0.5) then

        !cat 6: closed shrublands

                covercatf(x,y,6) = vegcover
                
            else

        !cat 7: open shrublands

                covercatf(x,y,7) = vegcover

      end if
        
    end if

    !cat 8: woody savannas
    
    if (cover(x,y,2) < 0.75 .and. cover(x,y,2) >= 0.5) then
      covercatf(x,y,8) = cover(x,y,2) + cover(x,y,9)
    end if
     
    !cat 9: savannas
    
    if (cover(x,y,2) < cover(x,y,9)) then
      covercatf(x,y,9) = vegcover
    end if
        
    !cat 10: grasslands
   
    if (any(present(1:5)) .and. grasscover > treecover) then
      covercatf(x,y,10) = grasscover
    end if
    
    !cat 11: permanent wetland
    !not present

    !cat 12: croplands
    !not present

    !cat 13: urban and built-up
    !not present

    !cat 14: mosaic of cropland and natural vegetation
    !not present

    !cat 15: snow and ice

    covercatf(x,y,15) = icef(x,y)
    
    !cat 16: barren or sparsely vegetated
    
    covercatf(x,y,16) = 1. - (vegcover + icef(x,y) + h2of(x,y))
    
!     write(0,*)covercatf(x,y,16),vegcover,icef(x,y),h2of(x,y)
!     read(*,*)
    
    !cat 17: water
        
    covercatf(x,y,17) = h2of(x,y)

    !tundra types
    
        if (all(.not.present(1:5))) then

            !cat 18: wooded tundra
    
            if (grasscover > 0.5 .and. treecover > 0.2) then
                covercatf(x,y,18) = vegcover
            end if

            !cat 19: mixed tundra
    
            if (grasscover > 0.5 .and. (treecover > 0. .and. treecover <= 0.2)) then
                covercatf(x,y,19) = vegcover
            end if
    
            !cat 20: bare ground tundra

            if (cover(x,y,8) <= 0.5 .and.all(.not.present(6:7))) then
              covercatf(x,y,20) = vegcover
          end if
            
        end if
        
        !----
        
        almean(x,y,:) = (1.-vegcover) * (icef(x,y) * alice + h2of(x,y) * alwater) + covercatf(x,y,16) * alsoil
        
        if (vegcover <= 0.) cycle

        maxlai = maxval(mLAIpft(x,y,:,:),dim=2)
                
        do m = 1,12
        
          mLAI(x,y,m) = sum(mLAIpft(x,y,:,m) * cover(x,y,:))
          
          where (maxlai /= 0.)
            phen = mLAIpft(x,y,:,m) / maxlai
          elsewhere
            phen = 0.
        end where
          
          fpar(x,y,m) = sum(cover(x,y,:) * phen)
          
         !mean albedo of the vegetated fraction  
          
          alveg = sum(cover(x,y,:) * phen * pftalbedo) + &
                  sum(cover(x,y,:) * 0.50 * (1. - phen) * alstem) +  &
                  sum(cover(x,y,:) * 0.50 * (1. - phen) * alsoil)
          
          !albedo of the unvegetated fraction
          
          almean(x,y,m) = vegcover * alveg + (1. - vegcover) * almean(x,y,m)

        end do

    !max snow albedo calculation
    !effectively anything can be completely covered by snow, but let's remove a bit to account for shading of snow under tree cover
    !1 - (fraction of evergreen trees + 25% fraction of deciduous trees(to account for stems/branches))
    !only on the vegetated fraction of the gridcell

        snowmax(x,y) = 1. - vegcover * (treecover_evrgn + 0.25 * treecover_decid)
        
  end do
end do


allocate(covercat(xlen,ylen))

covercat = -128

covercat = maxloc(covercatf,dim=3)

!----


!------------

status = nf90_close(ifid)
if (status /= nf90_noerr) call handle_err(status)

!------------
! mask the output variables for the unused buffer area

do i = 1,nc
  where (landf < 0.)
    covercatf(:,:,i) = -1.
  end where
end do

do i = 1,12
  where (landf < 0.)
    fpar(:,:,i)   = -9999.
    mLAI(:,:,i)   = -9999.
    almean(:,:,i) = -9999.
  end where
end do

where (landf < 0.)
  covercat = -128
  tmean    = -9999.
  snowmax  = -9999.
end where

!------------

write(0,*)'writing'

call getarg(2,outfile)

status = nf90_open(outfile,nf90_write,ofid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'lu_index',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,covercat)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'landusef',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,covercatf)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'LAI12m',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,mLAI)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'greenfrac',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,fpar)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'albedo12m',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,almean)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'snoalb',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,snowmax)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'soiltemp',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,tmean)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(ofid)
if (status /= nf90_noerr) call handle_err(status)

!------------

contains

!----------------------------

subroutine handle_err(status)

! Internal subroutine - checks error status after each netcdf call,
! prints out text message each time an error code is returned. 

  integer, intent (in) :: status
    
  if(status /= nf90_noerr) then 
    write(0,*)trim(nf90_strerror(status))
    stop
  end if

end subroutine handle_err

!----------------------------


end program postprocess
